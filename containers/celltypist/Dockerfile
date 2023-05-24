FROM python:3

COPY context/* .
RUN pip install -r requirements-freeze.txt

ENV CELLTYPIST_FOLDER="/models"
RUN python -c "import celltypist; celltypist.models.download_models()"

ENTRYPOINT ["python", "/main.py"]