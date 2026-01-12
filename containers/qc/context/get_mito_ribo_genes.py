#!/usr/bin/env python3
"""
get_mito_ribo_genes_mygene.py

MyGene.info replacement for BioMart using fetch_all + scroll_id to retrieve
all results (avoids the Elasticsearch from+size 10k limit).

Produces two CSVs optimized for spatial transcriptomics QC:
 - human_ribosomal_genes.csv (structural ribosomal proteins only: RPL*/RPS*)
 - human_mitochondrial_genes.csv (MT chromosome genes only, excluding pseudogenes)

Usage:
    python get_mito_ribo_genes_mygene.py [--cache-file mygene_cache.jsonl] [--max-records N] [--force-refresh]
"""
from pathlib import Path
import argparse
import json
import time
import sys
import urllib.request

import requests
import obonet
import networkx as nx
import pandas as pd

# ---------------- Configuration ----------------
OUTPUT_MITO_CSV = "human_mitochondrial_genes.csv"
OUTPUT_RIBO_CSV = "human_ribosomal_genes.csv"

GO_OBO_URL = "http://purl.obolibrary.org/obo/go.obo"
LOCAL_GO_OBO = "go.obo"

# GO root terms to expand - more focused for QC use case
RIBO_ROOTS = {
    "GO:0003735",  # structural constituent of ribosome (more specific than GO:0005840)
}
MITO_ROOTS = set()  # Not using GO terms for mito - using chromosome MT only for QC

# For QC purposes, we focus on:
# - Mitochondrial: Only genes encoded on chromosome MT (13 protein-coding genes + rRNAs/tRNAs)
# - Ribosomal: Only structural ribosomal proteins (RPL*/RPS* genes)
USE_CHROMOSOME_FOR_MITO = True
USE_GENE_NAME_PATTERN_FOR_RIBO = True

# MyGene.info query params defaults
MYGENE_BASE = "https://mygene.info/v3/query"
MYGENE_FIELDS = "ensembl.gene,symbol,genomic_pos,go,type_of_gene"
MYGENE_SPECIES = "human"

# polite delay between MyGene requests (seconds)
MYGENE_SLEEP = 0.2

# ---------------- Utilities ----------------
def download_file(url: str, dest: Path, show_size=True):
    dest = Path(dest)
    if dest.exists():
        return
    print(f"Downloading {url} -> {dest} ...")
    dest.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(url, str(dest))
    if show_size:
        try:
            print(f"Downloaded {dest} ({dest.stat().st_size:,} bytes)")
        except Exception:
            print(f"Downloaded {dest}")

def fetch_go_graph(obo_url=GO_OBO_URL, local_path=LOCAL_GO_OBO):
    local_path = Path(local_path)
    if not local_path.exists():
        download_file(obo_url, local_path)
    print(f"Parsing GO OBO from local file: {local_path}")
    t0 = time.time()
    graph = obonet.read_obo(local_path)
    dt = time.time() - t0
    print(f"Loaded GO OBO in {dt:.1f}s: {len(graph.nodes())} nodes, {len(graph.edges())} edges.")
    return graph

def expand_descendants(graph, root_go_ids):
    rev = graph.reverse()
    expanded = set()
    for gid in root_go_ids:
        if gid not in graph:
            print(f"Warning: GO id {gid} not found in ontology; skipping.")
            continue
        desc = nx.descendants(rev, gid)
        expanded.update(desc)
        expanded.add(gid)
    return set(x.upper() for x in expanded)

# ---------------- MyGene.info fetch & cache (scrolling) ----------------
def fetch_mygene_hits(cache_file: Path, batch_size:int=1000, max_records: int = None, force_refresh=False):
    """
    Fetch all MyGene.info hits for the species using fetch_all + scroll_id.
    - Initial request: q="*", species=MYGENE_SPECIES, fetch_all=true, size=batch_size.
    - API returns 'hits' and a '_scroll_id'. Subsequent requests: ?scroll_id=<id>
      return the next chunk until empty.
    - Caches results as JSONL to cache_file.
    """
    cache_file = Path(cache_file)

    # If cache exists and not forcing refresh, try to load it:
    if cache_file.exists() and not force_refresh:
        print(f"Loading cached MyGene hits from {cache_file}")
        hits = []
        with cache_file.open("r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    hits.append(json.loads(line))
                except Exception:
                    continue
        print(f"Loaded {len(hits)} cached hits.")
        if len(hits) > 0:
            return hits
        else:
            print("Cached file exists but contains 0 hits; will attempt to fetch fresh data using scroll API.")

    session = requests.Session()
    all_hits = []
    # Prepare initial request parameters (use fetch_all=true)
    initial_params = {
        "q": "*",
        "species": MYGENE_SPECIES,
        "fields": MYGENE_FIELDS,
        "size": batch_size,
        "ensemblonly": "true",
        "from": 0,
        "fetch_all": "true",  # triggers scroll behavior; returns _scroll_id
    }
    print(f"Debug: MyGene initial request params: {initial_params}")

    # initial request with retries
    tries = 0
    resp = None
    while tries < 4:
        try:
            resp = session.get(MYGENE_BASE, params=initial_params, timeout=60)
            resp.raise_for_status()
            break
        except requests.RequestException as e:
            tries += 1
            wait = 2 ** tries
            print(f"MyGene initial request failed (tries={tries}): {e}. Retrying in {wait}s...", file=sys.stderr)
            time.sleep(wait)
    if resp is None:
        raise RuntimeError("Failed initial MyGene request after retries.")
    j = resp.json()
    hits = j.get("hits", [])
    scroll_id = j.get("_scroll_id")
    if hits:
        all_hits.extend(hits)
    print(f"Fetched {len(all_hits)} hits (initial).", end="\r", flush=True)

    # If no scroll_id returned, then server either returned all results or zero; proceed accordingly.
    # Otherwise, iterate using scroll_id until no more hits or max_records reached.
    while scroll_id:
        if max_records is not None and len(all_hits) >= max_records:
            break
        # request next chunk with scroll_id
        tries = 0
        resp = None
        while tries < 4:
            try:
                resp = session.get(MYGENE_BASE, params={"scroll_id": scroll_id}, timeout=60)
                resp.raise_for_status()
                break
            except requests.RequestException as e:
                tries += 1
                wait = 2 ** tries
                print(f"MyGene scroll request failed (tries={tries}): {e}. Retrying in {wait}s...", file=sys.stderr)
                time.sleep(wait)
        if resp is None:
            print("Warning: failed to continue scroll; stopping early.", file=sys.stderr)
            break
        j = resp.json()
        hits = j.get("hits", [])
        scroll_id = j.get("_scroll_id") or None
        if not hits:
            break
        all_hits.extend(hits)
        print(f"Fetched {len(all_hits)} hits...", end="\r", flush=True)
        # polite sleep between scroll calls
        time.sleep(MYGENE_SLEEP)

    print(f"\nFinished fetching; total hits: {len(all_hits)}")
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    with cache_file.open("w", encoding="utf-8") as fh:
        for h in all_hits:
            fh.write(json.dumps(h, separators=(",", ":")) + "\n")
    print(f"Wrote cache to {cache_file}")

    if len(all_hits) == 0:
        print("\nWARNING: MyGene.info returned zero hits. Check network/access or parameters.")
    return all_hits

# ---------------- Convert MyGene hits to gene-GO rows ----------------
def is_ribosomal_structural_protein(gene_name):
    """
    Check if gene name matches canonical ribosomal protein pattern for QC.
    Returns True for RPL* (large subunit) and RPS* (small subunit) genes.
    Excludes pseudogenes and other ribosome-related but non-structural genes.
    """
    if not gene_name:
        return False
    gn = str(gene_name).strip().upper()
    # Match RPL or RPS followed by digits (optionally followed by A/B/C for paralogs)
    if gn.startswith("RPL") or gn.startswith("RPS"):
        remainder = gn[3:]  # Skip RPL/RPS prefix
        if not remainder:
            return False
        # Should start with digit
        if remainder[0].isdigit():
            # Check if it's not a pseudogene
            if "P" in remainder and remainder.index("P") > 0:
                # Could be a paralog (e.g., RPS27A) or pseudogene (e.g., RPL7P1)
                # Exclude if ends with P followed by digit (pseudogene pattern)
                if remainder.endswith("P") or (len(remainder) > 1 and remainder[-2] == "P" and remainder[-1].isdigit()):
                    return False
            return True
    return False

def mygene_hit_to_rows(hit):
    rows = []
    ensembl_id = None
    if "ensembl" in hit:
        en = hit["ensembl"]
        if isinstance(en, dict):
            en_gene = en.get("gene") or en.get("ensembl") or None
            if en_gene:
                ensembl_id = str(en_gene)
        elif isinstance(en, list) and en:
            first = en[0]
            if isinstance(first, dict):
                en_gene = first.get("gene") or first.get("ensembl")
                if en_gene:
                    ensembl_id = str(en_gene)
            else:
                ensembl_id = str(first)
    if not ensembl_id:
        _id = hit.get("_id")
        if isinstance(_id, str) and _id.startswith("ENSG"):
            ensembl_id = _id
    if not ensembl_id:
        return []

    gene_name = hit.get("symbol") or ""
    gene_type = hit.get("type_of_gene") or ""
    chromosome = ""
    gp = hit.get("genomic_pos")
    if gp:
        if isinstance(gp, list) and gp:
            gp0 = gp[0]
            chromosome = gp0.get("chr") or ""
        elif isinstance(gp, dict):
            chromosome = gp.get("chr") or ""
    go = hit.get("go")
    go_entries = []
    if go:
        if isinstance(go, dict):
            for cat in ("BP", "MF", "CC"):
                val = go.get(cat)
                if val:
                    if isinstance(val, list):
                        for elem in val:
                            gid = None
                            gname = None
                            if isinstance(elem, dict):
                                gid = elem.get("id") or elem.get("goid") or elem.get("GOID")
                                gname = elem.get("term") or elem.get("name")
                            else:
                                gid = str(elem)
                            if gid:
                                go_entries.append((str(gid).upper(), gname or ""))
                    elif isinstance(val, dict):
                        gid = val.get("id") or val.get("goid")
                        gname = val.get("term") or val.get("name")
                        if gid:
                            go_entries.append((str(gid).upper(), gname or ""))
                    else:
                        go_entries.append((str(val).upper(), ""))
        elif isinstance(go, list):
            for elem in go:
                if isinstance(elem, dict):
                    gid = elem.get("id") or elem.get("goid")
                    gname = elem.get("term") or elem.get("name")
                    if gid:
                        go_entries.append((str(gid).upper(), gname or ""))
                else:
                    go_entries.append((str(elem).upper(), ""))
        elif isinstance(go, str):
            go_entries.append((go.upper(), ""))
    normalized = []
    for gid, gname in go_entries:
        if not gid:
            continue
        g = gid.strip().upper()
        if g.isdigit():
            g = f"GO:{g.zfill(7)}"
        elif g.startswith("GO:"):
            parts = g.split("GO:")
            if parts and len(parts[-1].strip()) <= 7 and parts[-1].strip().isdigit():
                g = f"GO:{parts[-1].strip().zfill(7)}"
        normalized.append((g, gname or ""))
    if normalized:
        for gid, gname in normalized:
            rows.append({
                "ensembl_gene_id": ensembl_id,
                "gene_name": gene_name,
                "gene_type": gene_type,
                "chromosome": chromosome,
                "go_id": gid,
                "go_term": gname or ""
            })
    else:
        rows.append({
            "ensembl_gene_id": ensembl_id,
            "gene_name": gene_name,
            "gene_type": gene_type,
            "chromosome": chromosome,
            "go_id": None,
            "go_term": None
        })
    return rows

def hits_to_dataframe(hits):
    rows = []
    for h in hits:
        r = mygene_hit_to_rows(h)
        rows.extend(r)
    if not rows:
        return pd.DataFrame(columns=["ensembl_gene_id", "gene_name", "gene_type", "chromosome", "go_id", "go_term"])
    df = pd.DataFrame(rows)
    for c in ["ensembl_gene_id", "gene_name", "gene_type", "chromosome", "go_id", "go_term"]:
        if c not in df.columns:
            df[c] = None
    return df.astype(object)

# ---------------- Aggregation and output (robust) ----------------
def aggregate_go_annotations(df):
    agg = {}
    for _, row in df.iterrows():
        eid = str(row["ensembl_gene_id"]).strip()
        if not eid or eid.lower() == "nan":
            continue
        if eid not in agg:
            agg[eid] = {
                "gene_name": row["gene_name"] if pd.notna(row["gene_name"]) else "",
                "chromosome": row["chromosome"] if pd.notna(row["chromosome"]) else "",
                "gene_type": row.get("gene_type", "") if pd.notna(row.get("gene_type")) else "",
                "go_ids": set(),
                "go_terms": set(),
            }
        goid = row["go_id"]
        goterm = row["go_term"]
        if pd.notna(goid) and str(goid).strip():
            for part in str(goid).replace(",", ";").split(";"):
                p = part.strip()
                if not p:
                    continue
                p_up = p.upper()
                if p_up.startswith("GO:"):
                    agg[eid]["go_ids"].add(p_up)
                else:
                    if p.isdigit():
                        agg[eid]["go_ids"].add(f"GO:{p.zfill(7)}")
                    else:
                        agg[eid]["go_ids"].add(p_up)
        if pd.notna(goterm) and str(goterm).strip():
            for part in str(goterm).replace("|", ";").split(";"):
                t = part.strip()
                if t:
                    agg[eid]["go_terms"].add(t)
    return agg

def split_and_write_outputs(agg, ribo_expanded, mito_expanded, ensembl_meta, go_obo_url,
                            output_ribo=OUTPUT_RIBO_CSV, output_mito=OUTPUT_MITO_CSV):
    ribo_rows = []
    mito_rows = []

    for eid, info in agg.items():
        gene_name = info["gene_name"] or ""
        chrom = (info["chromosome"] or "").strip()
        gene_type = info.get("gene_type") or ""
        go_ids = sorted(info["go_ids"])
        go_terms = sorted(info["go_terms"])

        # Ribosomal: strict filtering for QC
        # Use gene name pattern (RPL*/RPS*) as primary filter
        # Optionally validate with GO term for additional confidence
        is_ribo_by_name = is_ribosomal_structural_protein(gene_name)
        matching_ribo = sorted(set(go_ids).intersection(ribo_expanded)) if ribo_expanded else []
        
        if USE_GENE_NAME_PATTERN_FOR_RIBO:
            # Primary: gene name pattern
            if is_ribo_by_name:
                evidence = ["gene_name_pattern:RPL*/RPS*"]
                if matching_ribo:
                    evidence.append(f"go_terms:{','.join(matching_ribo)}")
                ribo_rows.append({
                    "ensembl_gene_id": eid,
                    "gene_name": gene_name,
                    "chromosome": chrom,
                    "gene_type": gene_type,
                    "go_ids": ";".join(go_ids),
                    "go_terms": ";".join(go_terms),
                    "ribosomal_evidence": ";".join(evidence),
                    "ensembl_dataset": ensembl_meta,
                    "go_obo_source": go_obo_url,
                    "filter_note": "Structural ribosomal protein (RPL*/RPS* pattern) for QC",
                })
        else:
            # Fallback: GO term only
            if matching_ribo:
                ribo_rows.append({
                    "ensembl_gene_id": eid,
                    "gene_name": gene_name,
                    "chromosome": chrom,
                    "gene_type": gene_type,
                    "go_ids": ";".join(go_ids),
                    "go_terms": ";".join(go_terms),
                    "ribosomal_evidence": ";".join(matching_ribo),
                    "ensembl_dataset": ensembl_meta,
                    "go_obo_source": go_obo_url,
                    "filter_note": "GO:0003735 structural constituent of ribosome",
                })

        # Mitochondrial: strict filtering for QC
        # For spatial transcriptomics QC, we want ONLY genes encoded on MT chromosome
        # These are the 13 protein-coding genes + rRNA/tRNA genes
        is_mito_chr = (str(chrom).upper() == "MT")
        
        # Filter out pseudogenes even if on MT
        is_pseudogene = "pseudo" in gene_type.lower() if gene_type else False
        
        if is_mito_chr and not is_pseudogene:
            evidence = ["chromosome:MT"]
            if go_ids:
                evidence.append(f"go_annotation_present")
            mito_rows.append({
                "ensembl_gene_id": eid,
                "gene_name": gene_name,
                "chromosome": chrom,
                "gene_type": gene_type,
                "go_ids": ";".join(go_ids),
                "go_terms": ";".join(go_terms),
                "mitochondrial_evidence": ";".join(evidence),
                "ensembl_dataset": ensembl_meta,
                "go_obo_source": go_obo_url,
                "filter_note": "MT-encoded gene for QC (excludes nuclear-encoded mitochondrial proteins)",
            })

    ribo_columns = [
        "ensembl_gene_id", "gene_name", "chromosome", "gene_type", "go_ids", "go_terms",
        "ribosomal_evidence", "ensembl_dataset", "go_obo_source", "filter_note"
    ]
    mito_columns = [
        "ensembl_gene_id", "gene_name", "chromosome", "gene_type", "go_ids", "go_terms",
        "mitochondrial_evidence", "ensembl_dataset", "go_obo_source", "filter_note"
    ]

    if ribo_rows:
        ribo_df = pd.DataFrame(ribo_rows)[ribo_columns].sort_values("ensembl_gene_id").reset_index(drop=True)
    else:
        ribo_df = pd.DataFrame(columns=ribo_columns)

    if mito_rows:
        mito_df = pd.DataFrame(mito_rows)[mito_columns].sort_values("ensembl_gene_id").reset_index(drop=True)
    else:
        mito_df = pd.DataFrame(columns=mito_columns)

    print(f"Writing {len(ribo_df):,} ribosomal genes to {output_ribo} ...")
    print(f"  (Filtered for structural ribosomal proteins: RPL*/RPS* pattern)")
    ribo_df.to_csv(output_ribo, index=False)
    print(f"Writing {len(mito_df):,} mitochondrial genes to {output_mito} ...")
    print(f"  (Filtered for MT chromosome genes only, excluding nuclear-encoded and pseudogenes)")
    mito_df.to_csv(output_mito, index=False)

# ---------------- Main ----------------
def main():
    parser = argparse.ArgumentParser(
        description="Classify human genes into ribosomal/mitochondrial for spatial transcriptomics QC. "
                    "Optimized filters: RPL*/RPS* genes for ribosomal, MT chromosome genes for mitochondrial."
    )
    parser.add_argument("--cache-file", "-c", default="mygene_cache.jsonl", help="JSONL cache of raw MyGene hits")
    parser.add_argument("--batch-size", type=int, default=1000, help="MyGene.info page size (default 1000)")
    parser.add_argument("--max-records", type=int, default=None, help="Optional: stop after this many records (for testing)")
    parser.add_argument("--force-refresh", action="store_true", help="Force refresh of MyGene fetch (ignore cache)")
    parser.add_argument("--output-ribo", default=OUTPUT_RIBO_CSV)
    parser.add_argument("--output-mito", default=OUTPUT_MITO_CSV)
    parser.add_argument("--go-obo-url", default=GO_OBO_URL)
    parser.add_argument("--skip-go", action="store_true", 
                       help="Skip GO ontology fetch (faster; uses gene names and chromosome only)")
    args = parser.parse_args()

    print("=" * 70)
    print("QC Gene List Generator for Spatial Transcriptomics")
    print("=" * 70)
    print(f"Filter strategy:")
    print(f"  Ribosomal: RPL*/RPS* gene name pattern (structural proteins only)")
    print(f"  Mitochondrial: Chromosome MT only (excludes nuclear-encoded)")
    print("=" * 70)

    # Fetch GO ontology if needed for validation
    ribo_expanded = set()
    mito_expanded = set()
    if not args.skip_go and RIBO_ROOTS:
        go_graph = fetch_go_graph(args.go_obo_url, LOCAL_GO_OBO)
        if RIBO_ROOTS:
            print("Expanding ribosomal GO roots for validation ...")
            ribo_expanded = expand_descendants(go_graph, RIBO_ROOTS)
            print(f"  Ribosomal GO set size: {len(ribo_expanded):,}")
        if MITO_ROOTS:
            print("Expanding mitochondrial GO roots ...")
            mito_expanded = expand_descendants(go_graph, MITO_ROOTS)
            print(f"  Mitochondrial GO set size: {len(mito_expanded):,}")
    else:
        print("Skipping GO ontology (using gene name and chromosome filters only)")

    hits = fetch_mygene_hits(Path(args.cache_file), batch_size=args.batch_size,
                             max_records=args.max_records, force_refresh=args.force_refresh)

    df = hits_to_dataframe(hits)
    print(f"Converted MyGene hits to {len(df):,} gene-GO rows and {df['ensembl_gene_id'].nunique():,} unique Ensembl genes.")

    agg = aggregate_go_annotations(df)
    print(f"Aggregated {len(agg):,} unique Ensembl genes")

    ensembl_meta = f"mygene.info:v3;fetched={time.strftime('%Y-%m-%d')}"

    split_and_write_outputs(agg, ribo_expanded, mito_expanded, ensembl_meta, args.go_obo_url,
                            output_ribo=args.output_ribo, output_mito=args.output_mito)

    print("\n" + "=" * 70)
    print("QC Gene List Generation Complete")
    print("=" * 70)

if __name__ == "__main__":
    main()
