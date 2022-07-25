===========================
Non-symbolic number stimuli
===========================

.. currentmodule:: pynsn

Object Array
=============
.. autosummary::
    :toctree: api/

    DotArray
    RectangleArray
    ArrayProperties


Dot Array
--------------

.. autoclass:: DotArray
   :members:
   :inherited-members:
   :undoc-members:

Rectangle Array
--------------

.. autoclass:: pynsn.RectangleArray
   :members:
   :inherited-members:
   :undoc-members:


Visual Properties
-----------------


.. autoclass:: ArrayProperties
   :members:
   :inherited-members:
   :undoc-members:


Object Shapes
=============

.. autosummary::
    :toctree: api/

    Dot
    Rectangle
    Point
    PictureFile



Dot
-----------------
.. autoclass:: pynsn.Dot
   :members:
   :inherited-members:
   :undoc-members:

Rectangle
-----------------
.. autoclass:: pynsn.Rectangle
   :members:
   :inherited-members:
   :undoc-members:

Point
-----------------
.. autoclass:: pynsn.Point
   :members:
   :inherited-members:
   :undoc-members:


PictureFile
-----------------

``PictureFile`` can be used as attribute of a `Rectangle`_ to use pictures in `Rectangle Array`_. The pictures will be
scaled to the size of the rectangle.

.. code-block:: py

   pict_object = pynsn.Rectangle(xy=(0,0), size=(80, 80),
                        attribute=nsn.PictureFile("mypict.png")



.. autoclass:: pynsn.PictureFile
   :members:
   :inherited-members:
   :undoc-members:

