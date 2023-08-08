FROM hra-workflows/shared:main as shared
FROM hubmap/azimuth-annotate:1.3

COPY context/requirements-freeze.txt .
RUN pip install -r requirements-freeze.txt

COPY --from=shared /shared/ ./shared/
COPY context/* .

ENTRYPOINT ["python3", "/main.py"]