DISC
====

|PyPI|

.. |PyPI| image:: https://img.shields.io/pypi/v/DISC.svg
    :target: https://pypi.org/project/disc
..
  |bioconda| |Docs| |Build Status| |Coverage| |Code Style| |Downloads|
.. |bioconda| image:: https://img.shields.io/badge/bioconda-blue.svg
    :target: http://bioconda.github.io/recipes/scvi/README.html
.. |Docs| image:: https://readthedocs.org/projects/scvi/badge/?version=latest
    :target: https://scvi.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
.. |Build Status| image:: https://travis-ci.org/YosefLab/scVI.svg?branch=master
    :target: https://travis-ci.org/YosefLab/scVI
.. |Coverage| image:: https://codecov.io/gh/YosefLab/scVI/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/YosefLab/scVI
.. |Code Style| image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/python/black
.. |Downloads| image:: https://pepy.tech/badge/scvi
   :target: https://pepy.tech/project/scvi

A scalable deep learning imputation method with semi-supervised learning

* Free software: Apache License 2.0

..
 * Documentation: https://scvi.readthedocs.io.


Requirements
----
- Python >= 3.6
- tensorflow-gpu >= 1.13.1
- pandas >= 0.23.0
- h5py >= 2.9.0

Quick Start
-----------

1. Install Python_. We typically use Python3.6 in Linux CentOS 7.

.. _Python: https://www.python.org/downloads/

2. Install TensorFlow_. You can install a gpu version of TensorFlow if you have an Nvidia GPU. We use tensorflow-gpu==1.13.1 here.

.. _TensorFlow: https://www.tensorflow.org/install/pip

3. Install DISC through pip:

    ``pip install DISC``

   Alternatively, you may clone this repository and run ``python setup.py install``.

4. Follow along with our Jupyter notebooks to quickly get familiar with DISC!

   a. Getting started:
       * `Running melanoma dataset`_

   b. Reproducing our results:
       * `Structure recovery comparing with FISH`_
       * `Downsampling recovery`_
       * `Clustering improvement`_
       * `Ultra-large dataset analysis`_

   c. Supplementary topics:
       * `Download datasets we used`_
       * `Other analysis scripts we used`_
..
   d. Advanced topics:

.. _`running melanoma dataset`: https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/data_loading.ipynb
.. _`Structure recovery comparing with FISH`: https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/data_loading.ipynb
.. _`Downsampling recovery`: https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/data_loading.ipynb
.. _`Clustering improvement`: https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/data_loading.ipynb
.. _`Ultra-large dataset analysis`: https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/data_loading.ipynb
.. _`Download datasets we used`: https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/data_loading.ipynb
.. _`Other analysis scripts we used`: https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/data_loading.ipynb



References
----------
..
 Romain Lopez, Jeffrey Regier, Michael Cole, Michael I. Jordan, Nir Yosef.
 **"Deep generative modeling for single-cell transcriptomics."**
 Nature Methods, 2018. `[pdf]`__
 
 .. __: https://rdcu.be/bdHYQ
 

