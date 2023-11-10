====================
Distributions
====================


Should be use to implement a random categorical variable (levels) that are not equally distributed.
To specify a equally distributes categorical variable for ``NSNFactory`` you can also simply use a tuple or list, e.g.:

.. code-block:: py

   example code = TODO
   example code = 23


.. currentmodule:: pynsn.distributions
.. autosummary::
    :toctree: api/

    Levels
    Uniform
    Normal
    Triangle
    Beta
    Normal2D


Discrete distribution
-----------------------


.. autoclass:: pynsn.distributions.Levels
   :members:
   :inherited-members:
   :undoc-members:

Continuous distributions
-------------------------


.. autoclass:: pynsn.distributions.Uniform
   :members:
   :inherited-members:
   :undoc-members:

.. autoclass:: pynsn.distributions.Normal
   :members:
   :inherited-members:
   :undoc-members:


.. autoclass:: pynsn.distributions.Triangle
   :members:
   :inherited-members:
   :undoc-members:

.. autoclass:: pynsn.distributions.Beta
   :members:
   :inherited-members:
   :undoc-members:

Multivariate distribution
--------------------------

.. autoclass:: pynsn.distributions.Normal2D
   :members:
   :inherited-members:
   :undoc-members:



