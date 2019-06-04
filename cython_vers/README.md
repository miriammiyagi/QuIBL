## QuIBL - Cython Version

As an alternative to the Python only version, this script is faster but requires Cython. To compile, type ``python setup.py build_ext --inplace``, then run ``QuIBL_cyth.py`` as normal.

You may receive benign Numpy version errors on compilation, these are normal since Cython uses a deprecated Numpy API.
