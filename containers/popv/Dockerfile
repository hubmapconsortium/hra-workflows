FROM tensorflow/tensorflow:2.11.1

RUN apt-get -y update && apt-get -y install python3-dev g++ git wget

ENV MODELS_DIR=./models ZENODO_MODELS_ID=7580707
COPY context/download-models.sh .
RUN ./download-models.sh

COPY context/download-ontology.sh .
RUN ./download-ontology.sh

COPY context/requirements-freeze.txt .
RUN pip install -r requirements-freeze.txt

COPY context/* .

ENTRYPOINT ["python", "/main.py"]
