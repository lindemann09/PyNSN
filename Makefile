.PHONY: install clean build

build:
	python3 setup.py sdist bdist_wheel

install:
	python3 setup.py install

publish:
	twine check dist/*
	twine upload dist/*

unit_tests:
	python -m unittest discover tests

unit_test_docker:
	docker build -t pynsn_unit_test .
	docker run --rm pynsn_unit_test

clean:
		@rm -rf build \
			dist \
			pynsn.egg-info \
			.pytest_cache


