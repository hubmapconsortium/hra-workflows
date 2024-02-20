import re

_NON_WORD_REGEX = re.compile(r"\W+")
_NON_ALPHANUM_HYPHEN_REGEX = re.compile(r"[^a-z0-9-]+")
_ASCTB_ID_REGEX = re.compile(r"^ASCTB", re.IGNORECASE)
_ENSEMBL_ID_REGEX = re.compile(r"^ens", re.IGNORECASE)
_VALID_CELL_ID_REGEX = re.compile(r"^(CL|PCL):", re.IGNORECASE)


def create_temp_asctb_id(value: str) -> str:
    """Generate a temporary asctb id from a value.

    Args:
        value (str): The value to create an id for

    Returns:
        str: A temp id
    """
    if _ASCTB_ID_REGEX.match(value):
        return value

    value = value.strip().lower()
    value = re.sub(_NON_WORD_REGEX, "-", value)
    value = re.sub(_NON_ALPHANUM_HYPHEN_REGEX, "", value)
    return f"ASCTB-TEMP:{value}"


def create_cell_id(id: str) -> str:
    """Turn an id into a cell id.
    Cell ids start with CL: or PCL: otherwise it is turned into
    a temporary asctb id.

    Args:
        id (str): Original id

    Returns:
        str: Valid cell id
    """
    if _VALID_CELL_ID_REGEX.match(id):
        return id
    return create_temp_asctb_id(id)


def create_gene_id(id: str) -> str:
    """Turn an id into a gene id.

    Args:
        id (str): Original id

    Returns:
        str: Valid gene id
    """
    if _ENSEMBL_ID_REGEX.match(id):
        return f"ensembl:{id}"
    return create_temp_asctb_id(id)
