.PHONY: html start_jupyter

pynsn_link:
	ln -s ../pynsn ./ -f

html:
	make pynsn_link
	pip install jupyter
	make pynsn_demo.html
	make make_object_arrays_demo.html
	make euro_flag_demo.html
	cp *.html ../documentation/_static
	mv *.html ../docs
	rm -f pynsn

%.html:
	jupyter nbconvert --execute --to html $*.ipynb

start_jupyter:
	ln -s ../pynsn ./ -f
	jupyter notebook
	rm -f pynsn
 