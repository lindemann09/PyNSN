.PHONY: install clean build

build:
	python3 setup.py sdist bdist_wheel

install:
	python3 setup.py install

publish:
	twine check dist/*
	twine upload dist/*

clean:
		@rm -rf build \
			dist \
			pynsn.egg-info 
