===========================
Non-symbolic number stimuli
===========================

.. currentmodule:: pynsn

nsn stimulus
=============
.. autosummary::
    :toctree: api/

    DotArray
    RectangleArray
    ArrayProperties


nsn stimulus
--------------

.. autoclass:: DotArray
   :members:
   :inherited-members:
   :undoc-members:

nsn stimulus
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

``PictureFile`` can be used as attribute of a `Rectangle`_ to use pictures in `nsn stimulus`_. The pictures will be
scaled to the size of the rectangle.

.. code-block:: py

   pict_object = pynsn.Rectangle(xy=(0,0), size=(80, 80),
                        attribute=nsn.PictureFile("mypict.png")



.. autoclass:: pynsn.PictureFile
   :members:
   :inherited-members:
   :undoc-members:
