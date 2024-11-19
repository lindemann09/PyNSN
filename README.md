PyNSN
=====

**Python Library for Creating Non-Symbolic Number Displays**

---

[![GitHub license](https://img.shields.io/github/license/lindemann09/PyNSN)](https://github.com/lindemann09/PyNSN/blob/master/LICENSE)
[![Python Version](https://img.shields.io/pypi/pyversions/pynsn?style=flat)](https://www.python.org)
[![PyPI](https://img.shields.io/pypi/v/pynsn?style=flat)](https://pypi.org/project/pynsn/)

(c) Oliver Lindemann (lindemann@cognitive-psychology.eu)

Project homepage: https://github.com/lindemann09/PyNSN


## Installing 

requires Python 3 (>=3.10)

```
python -m pip install pynsn
```

## Image formats

By default, PyNSN is able to write [SVG](https://en.wikipedia.org/wiki/Scalable_Vector_Graphics)
or [Pillow](https://pillow.readthedocs.io/en/stable/) images.
To generate [Pygame](https://www.pygame.org/news) or
[Matplotlib](https://matplotlib.org/stable/index.html) images or stimuli
for [Expyriment](http://expyriment.org), please install the respective
packages.

## Examples
* [making arrays](https://lindemann09.github.io/PyNSN/make_object_arrays_demo.html): manually creating nsn stimuluss and exporting picture file
  ([binder](https://mybinder.org/v2/gh/lindemann09/PyNSN/jupyter?labpath=examples%2Fmake_object_arrays_demo.ipynb))
* [random arrays](https://lindemann09.github.io/PyNSN/pynsn_demo.html): Creating random nsn stimuli
  ([binder](https://mybinder.org/v2/gh/lindemann09/PyNSN/jupyter?labpath=examples%2Fpynsn_demo.ipynb))
* matching visual features
* data base, sequences
* [Euro flag example](https://lindemann09.github.io/PyNSN/euro_flag_demo.html): using pictures as objects in array
  ([binder](https://mybinder.org/v2/gh/lindemann09/PyNSN/jupyter?labpath=examples%2Feuro_flag_demo.ipynb))
