Installation instructions
=========================

Naming convention
-----------------

Please note that the package described here is indexed on PyPI
as ``aerosolpy`` and can be pip installed and imported under this name. 

Earlier, unofficial versions were distributed as ``aeropy``. 
The usage of this is deprecated and does not work with this version of 
``aerosolpy``

Dependencies
------------

``aerosolpy`` requires the use of `numpy <https://numpy.org/>`__,
`scipy <https://scipy.org/>`__, and `pandas <https://pandas.pydata.org/>`__.


Installation
------------

``aerosolpy`` is available on `PyPI <https://pypi.org/project/aerosolpy/>`__, so it can be
installed using pip::

		pip install aerosolpy

Otherwise the source code can be downloaded or cloned from the  
`GitHub repository <https://github.com/DominikStolzenburg/aerosolpy>`__ 
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