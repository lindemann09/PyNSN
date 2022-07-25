.PHONY: install clean build

venv:
	python -m venv venv
	. venv/bin/activate && \
	pip install pip -U && \
	pip install setuptools wheel && \
	pip install -r requirements.txt && \
	pip list

build: venv
	. venv/bin/activate && \
	python setup.py sdist bdist_wheel

publish:
	twine check dist/*
	twine upload dist/*

docker_unittest: 
	docker build -t pynsn38 -f tests/Dockerfile-py38 . && \
	docker build -t pynsn310 -f tests/Dockerfile-py310 .
	docker run --rm pynsn38
	docker run --rm pynsn310

apiref: venv
	. venv/bin/activate && \
	cd documentation && \
	make html check_api

jupyter_examples: venv
	. venv/bin/activate && \
	pip install jupyter && \
	cd examples && \
	make html

tox: venv
	. venv/bin/activate && \
	python -m tox

clean:
	@rm -rf build \
		venv \
		dist \
		pynsn.egg-info \
		.tox \
		.pytest_cache
	cd documentation && make clean

