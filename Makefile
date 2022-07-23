.PHONY: install clean build

build:
	python3 setup.py sdist bdist_wheel

install:
	python3 setup.py install

publish:
	twine check dist/*
	twine upload dist/*

docker_build: 
	docker build -t pynsn38 -f tests/Dockerfile-py38 .
	docker build -t pynsn310 -f tests/Dockerfile-py310 .

docker_unittest: 
	docker run --rm pynsn38
	docker run --rm pynsn310

clean:
		@rm -rf build \
			dist \
			pynsn.egg-info \
			.tox \
			.pytest_cache


