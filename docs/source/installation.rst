Installation instructions
=========================

Naming convention
-----------------

Please note that the package described here is going to be indexed on PyPI
as ``aerosolpython``.
However, throughout this manual it will be referenced by its
import name ``aeropy``


Dependencies
------------

``aeropy`` requires the use of `numpy <https://numpy.org/>`__,
`scipy <https://scipy.org/>`__, and `pandas <https://pandas.pydata.org/>`__.


Installation
------------

``aeropy`` is currently not yet available on PyPI, but it can be
installed by downloading the source code or cloning the  
`GitLab repository <https://gitlab.tuwien.ac.at/dominik.stolzenburg/aerosolpython>`__ 
and running the standard::

       python setup.py install

command or its usual variants (``python setup.py install --user``,
``python setup.py install --prefix=/PATH/TO/INSTALL/DIRECTORY``,
etc.).

The ``.tar.gz`` source file can also be installed using ``pip``.
In the directory of you current environment run the command::
        
      <absolute path to python.exe> -m pip install <path to tar.gz>

Note that ``<path to tar.gz>`` can be relative, absolute 
and even an online link.