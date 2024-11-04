.PHONY: install clean build


docker_unittest:
	docker build -t pynsn38 -f tests/Dockerfile-py38 . && \
	docker build -t pynsn310 -f tests/Dockerfile-py310 .
	docker run --rm pynsn38
	docker run --rm pynsn310

apiref:
	cd documentation && \
	make html check_api

jupyter_examples:
	cd examples && \
	make html

unittest:
	python -m unittest discover tests/

clean:
	@rm -rf build \
		dist \
		pynsn.egg-info \
		.tox \
		.pytest_cache \
		examples\pynsn
	cd documentation && make clean
	py3clean .
