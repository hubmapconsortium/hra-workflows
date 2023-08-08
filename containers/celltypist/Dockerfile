FROM hra-workflows/shared:main as shared
FROM python:3

COPY context/requirements-freeze.txt .
RUN pip install -r requirements-freeze.txt

ENV CELLTYPIST_FOLDER="/models"
RUN python -c "import celltypist; celltypist.models.download_models()"

COPY --from=shared /shared/ ./shared/
COPY context/* .

ENTRYPOINT ["python", "/main.py"]