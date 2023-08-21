FROM python:3

COPY context/requirements-freeze.txt .
RUN pip install -r requirements-freeze.txt

ENV CELLTYPIST_FOLDER="/models"
RUN python -c "import celltypist; celltypist.models.download_models()"

COPY context/* .

ENTRYPOINT ["python", "/main.py"]