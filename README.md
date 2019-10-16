PyNSN
=====

**Creating Non-Symbolic Number Displays**

---


*Released under the GNU General Public License v3* 

Oliver Lindemann (lindemann@cognitive-psychology.eu)

Project homepage: https://github.com/lindemann09/PyNSN


Dependencies
------------

* Python 3 (>=3.5)
* numpy (>=1.6)
* scipy (>=1.0)
* pillow (>=5.0)

Additional Python packages, which are optional and required only for some 
features:

* PyQT4 (>=4.10)
* pygame (>=1.9)
* expyriment (>=0.8)


Installing via `pip`
--------------------

```
python3 -m pip install --index-url https://test.pypi.org/simple/ psnsn
```

To generate `pygame` or `expyriment` stimuli install the optional packages

```
python3 -m pip install --index-url https://test.pypi.org/simple/ psnsn[expyriment]
```



RUN GUI tool
-------------

Note: Using the PyNSN-GUI requires the installation `PyQt4 4.10` or larger

```
python3 -m pynsn.gui
```

or if installed correctly:

```
pynsn-gui
```


