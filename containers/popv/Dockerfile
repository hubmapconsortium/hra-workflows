FROM python:3

ENV ZENODO_MODELS_ID=7580707
COPY context/download-models.sh .
RUN ./download-models.sh

COPY context/download-ontology.sh .
RUN ./download-ontology.sh

COPY context/requirements-freeze.txt .
RUN pip install -r requirements-freeze.txt

COPY context/* .

ENTRYPOINT ["python", "/main.py"]
