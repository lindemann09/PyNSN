.PHONY: html start_jupyter

html:
	ln -s ../pynsn ./ -f
	pip install jupyter
	make pynsn_demo.html
	make make_object_arrays_demo.html
	make euro_flag_demo.html
	mv *.html ../docs
	rm -f pynsn

%.html:
	jupyter nbconvert --execute --to html $*.ipynb

start_jupyter:
	ln -s ../pynsn ./ -f
	jupyter notebook
	rm -f pynsn
 
