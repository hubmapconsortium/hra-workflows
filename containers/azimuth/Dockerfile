FROM satijalab/azimuth:0.4.6

RUN R --no-echo -e "install.packages('rjson', repo='https://cloud.r-project.org')"

COPY context/requirements-freeze.txt .
RUN pip install -r requirements-freeze.txt

COPY context/ .

CMD ["python3", "/main.py"]