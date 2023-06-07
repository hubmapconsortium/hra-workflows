FROM hubmap/azimuth-annotate:1.3

COPY context/* .
RUN pip install -r requirements-freeze.txt

ENTRYPOINT ["python", "/main.py"]