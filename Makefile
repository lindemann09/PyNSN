.PHONY: install clean build

venv:
	python3 -m venv venv
	. venv/bin/activate && \
	pip install pip -U && \
	pip install setuptools wheel && \
	pip install scipy numpy Pillow svgwrite && \
	pip install matplotlib && \
	pip list

build: venv
	. venv/bin/activate && \
	python3 setup.py sdist bdist_wheel

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
	pip install sphinx sphinx_rtd_theme sphinx_autodoc_typehints && \
	cd documentation && \
	make html

jupyter_examples: venv
	. venv/bin/activate && \
	pip install jupyter && \
	cd examples && \
	make html

tox: venv
	. venv/bin/activate && \
	python3 -m tox


clean:
	@rm -rf build \
		venv \
		dist \
		pynsn.egg-info \
		.tox \
		.pytest_cache

