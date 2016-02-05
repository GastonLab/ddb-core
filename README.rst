========
Overview
========

.. start-badges

.. list-table::
    :stub-columns: 1

    * - docs
      - |docs|
    * - tests
      - | |travis| |requires|
        | |coveralls|
        | |landscape|
    * - package
      - |version| |downloads| |wheel| |supported-versions| |supported-implementations|

.. |docs| image:: https://readthedocs.org/projects/ddb-core/badge/?style=flat
    :target: https://readthedocs.org/projects/ddb-core
    :alt: Documentation Status

.. |travis| image:: https://travis-ci.org/dgaston/ddb-core.svg?branch=master
    :alt: Travis-CI Build Status
    :target: https://travis-ci.org/dgaston/ddb-core

.. |requires| image:: https://requires.io/github/dgaston/ddb-core/requirements.svg?branch=master
    :alt: Requirements Status
    :target: https://requires.io/github/dgaston/ddb-core/requirements/?branch=master

.. |coveralls| image:: https://coveralls.io/repos/dgaston/ddb-core/badge.svg?branch=master&service=github
    :alt: Coverage Status
    :target: https://coveralls.io/r/dgaston/ddb-core

.. |landscape| image:: https://landscape.io/github/dgaston/ddb-core/master/landscape.svg?style=flat
    :target: https://landscape.io/github/dgaston/ddb-core/master
    :alt: Code Quality Status

.. |version| image:: https://img.shields.io/pypi/v/ddb-core.svg?style=flat
    :alt: PyPI Package latest release
    :target: https://pypi.python.org/pypi/ddb-core

.. |downloads| image:: https://img.shields.io/pypi/dm/ddb-core.svg?style=flat
    :alt: PyPI Package monthly downloads
    :target: https://pypi.python.org/pypi/ddb-core

.. |wheel| image:: https://img.shields.io/pypi/wheel/ddb-core.svg?style=flat
    :alt: PyPI Wheel
    :target: https://pypi.python.org/pypi/ddb-core

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/ddb-core.svg?style=flat
    :alt: Supported versions
    :target: https://pypi.python.org/pypi/ddb-core

.. |supported-implementations| image:: https://img.shields.io/pypi/implementation/ddb-core.svg?style=flat
    :alt: Supported implementations
    :target: https://pypi.python.org/pypi/ddb-core


.. end-badges

Dead Drift Bioinformatics Core Library.

* Free software: BSD license

Installation
============

::

    pip install ddb-core

Documentation
=============

https://ddb-core.readthedocs.org/

Development
===========

To run the all tests run::

    tox

Note, to combine the coverage data from all the tox environments run:

.. list-table::
    :widths: 10 90
    :stub-columns: 1

    - - Windows
      - ::

            set PYTEST_ADDOPTS=--cov-append
            tox

    - - Other
      - ::

            PYTEST_ADDOPTS=--cov-append tox
