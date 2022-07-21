.PHONY: clean, venv, build

init:
	python3 -m venv venv
	. venv/bin/activate
	pip install matplotlib jupyter
	#pip install pynsn -U
	pip list
	rm *demo.html -f
	ln -s ../pynsn ./ -f


build_docs:
	make venv
	. venv/bin/activate
	ln -s ../pynsn ./ -f
	make pynsn_demo.html
	make make_object_arrays_demo.html
	make euro_flag_demo.html
	mv *.html ../docs
	rm -f pynsn

%.html:
	. venv/bin/activate
	jupyter nbconvert --execute --to html $*.ipynb

start_jupyter:
	. venv/bin/activate
	ln -s ../pynsn ./ -f
	jupyter notebook
	rm -f pynsn


clean:
	rm venv -rf
	rm -f pynsn

 
