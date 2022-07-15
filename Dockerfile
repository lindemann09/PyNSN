FROM jupyter/scipy-notebook

RUN pip install pip -U && \
	pip install pynsn
