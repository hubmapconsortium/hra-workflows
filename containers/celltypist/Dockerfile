FROM python:3.10

RUN apt-get -y update && apt-get -y install cmake

COPY context/requirements-freeze.txt .
RUN pip install -r requirements-freeze.txt && pip cache purge

ENV CELLTYPIST_FOLDER="/models"
RUN python -c "import celltypist; celltypist.models.download_models()"

COPY context/ .

CMD ["python", "/main.py"]