PyNSN
=====

**Creating Non-Symbolic Number Displays**

---


*Released under the GNU General Public License v3* 

Oliver Lindemann (lindemann@cognitive-psychology.eu)

Project homepage: https://github.com/lindemann09/PyNSN


## Dependencies

* Python 3 (>=3.5)
* numpy (>=1.6)
* scipy (>=1.0)
* Pillow ()>=5.0)
* svgwrite (>=1.4)

### Optional requirements

Additional Python packages, which are optional and required only for 
some features:

* matplotlib (>=3.2)
* pygame (>=1.9)
* expyriment (>=0.9)
* PyQT5 (>=5.14)


## Installing via `pip`

```
python3 -m pip install pynsn
```

## Image formats

By default, PyNSN is able to write [SVG](https://en.wikipedia.org/wiki/Scalable_Vector_Graphics) 
or [Pillow](https://pillow.readthedocs.io/en/stable/) images. 
To generate [Pygame](https://www.pygame.org/news) or
[PyPlot](https://matplotlib.org/stable/index.html) images or stimuli 
for [Expyriment](http://expyriment.org), please install the respective
packages.


RUN GUI tool
-------------

*Note:* Using the PyNSN-GUI requires the installation of `PyQt5` 

```
python3 -m pynsn
```

or if installed correctly:

```
pynsn
```


