FROM python:3

COPY context/requirements-freeze.txt .
RUN pip install -r requirements-freeze.txt

COPY context/* .

ENTRYPOINT ["python", "/main.py"]
