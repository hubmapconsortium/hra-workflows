FROM tensorflow/tensorflow:2.16.1

RUN apt-get -y update && apt-get -y install python3-dev g++ git wget && apt-get clean

COPY context/download-ontology.sh .
RUN ./download-ontology.sh

COPY context/requirements-freeze.txt .
RUN pip install -r requirements-freeze.txt && pip cache purge

COPY context/ .

# Popv unconditionally sets the number of threads to 10
# Ensure numba has the correct settings or it might throw an error
ENV NUMBA_NUM_THREADS=10

CMD ["python", "/main.py"]
