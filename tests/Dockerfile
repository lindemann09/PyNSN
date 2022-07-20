FROM python:3.8

WORKDIR /home

RUN pip install pip -U && \
	pip install numpy scipy pillow matplotlib svgwrite -U

COPY pynsn/ pynsn/
COPY tests/ tests/
COPY README.md .
COPY setup.py .

CMD python setup.py install && \
    rm pynsn -rf && \
    pip list && \
    python -m unittest discover tests

