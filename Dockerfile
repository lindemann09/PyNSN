FROM python:3.6

WORKDIR /home

COPY requirements.txt ./
RUN pip install pip -U && \
	pip install --no-cache-dir -r requirements.txt

COPY pynsn/ pynsn/
COPY tests/ tests/
COPY README.md .
COPY setup.py .

CMD ls && \
    python setup.py install && \
    rm pynsn -rf && \
    pip list && \
    python -m unittest discover tests

