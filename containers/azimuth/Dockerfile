FROM hubmap/azimuth-annotate:1.3

COPY context/requirements-freeze.txt .
RUN pip install -r requirements-freeze.txt

COPY context/* .

ENTRYPOINT ["python3", "/main.py"]