FROM python:3.10

COPY context/requirements-freeze.txt .
RUN pip install -r requirements-freeze.txt

COPY context/ .

CMD ["python", "/main.py"]
