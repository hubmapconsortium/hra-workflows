FROM hubmap/azimuth-annotate:1.3

COPY context/* .
RUN pip install -r requirements-freeze.txt

ENTRYPOINT ["python3", "/main.py"]