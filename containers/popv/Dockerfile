FROM tensorflow/tensorflow:2.11.1

RUN apt-get -y update && apt-get -y install python3-dev g++ git wget && apt-get clean

COPY context/download-ontology.sh .
RUN ./download-ontology.sh

COPY context/requirements-freeze.txt .
RUN pip install -r requirements-freeze.txt && pip cache purge

COPY context/ .

CMD ["python", "/main.py"]
