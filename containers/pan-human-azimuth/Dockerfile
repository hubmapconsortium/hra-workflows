FROM python:3.10

COPY context/requirements-freeze.txt .
RUN pip install -r requirements-freeze.txt && pip cache purge

COPY context/ .

CMD ["python", "/main.py"]