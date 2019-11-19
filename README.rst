DISC
====

|PyPI|

.. |PyPI| image:: https://img.shields.io/pypi/v/DISC.svg
    :target: https://pypi.org/project/disc

An accurate and scalable imputation algorithm based on semi-supervised deep learning for single-cell transcriptome

* Free software: Apache License 2.0

Requirements
------------

- Python_ >= 3.6
- tensorflow_ >= 1.13.1
- numpy_ >= 1.14.0
- pandas_ >= 0.21.0
- h5py_ >= 2.9.0

Installation
------------

**Installation with pip**

To install with ``pip``, run the following from a terminal::

  pip install DISC

**Installation from GitHub**

To clone the repository and install manually, run the following from a terminal::

  git clone git://github.com/iyhaoo/DISC.git

  cd DISC

  python setup.py install

Usage
-----

Quick Start
 1. Run DISC::

     disc \
     --dataset=matrix.loom \
     --out-dir=out_dir

    where ``matrix.loom`` is a `loom-formatted`_ raw count matrix with genes in rows and cells in columns and ``out_dir`` is the path of output directory.
 2. Results

    * ``log.tsv``: a tsv-formatted log file that records training states.
    * ``summary.pdf``: a pdf-formatted file that visualizes the fitting line and optimal point and it will be updated in real time when running.
    * ``summary.tsv``: a tsv-formatted file that shows the raw data of visualization.
    * ``result``: a directory for imputaion results as below:

      * ``imputation.loom``: a `loom-formatted`_ imputed matrix with genes in rows and cells in columns.
      * ``feature.loom``: a `loom-formatted`_ dimensionally reduced feature matrix provided by our method based on the imputed matrix above with feature in rows and cells in columns.
      * ``running_info.hdf5``: a `hdf5-formatted`_ saved some basic information about the input dataset such as library size, genes used for modelling and so on.

    * ``models``: a directory for trained models in every save interval

Tutorials
 1. Imputation

    * `Data pre-processing`_
    * `Imputation`_

 2. Reproducing our results:

    * `Structure recovery comparing with FISH`_
    * `Downsampling recovery`_
    * `Clustering improvement`_
    * `Ultra-large dataset analysis`_

 3. Supplementary topics:

    * `Download datasets we used`_
    * `Other analysis scripts we used`_



References
----------
..
 Romain Lopez, Jeffrey Regier, Michael Cole, Michael I. Jordan, Nir Yosef.
 **"Deep generative modeling for single-cell transcriptomics."**
 Nature Methods, 2018. `[pdf]`__
 
 .. __: https://rdcu.be/bdHYQ
 
History
-------

1.0.0 (2019-11-XX)
^^^^^^^^^^^^^^^^^^
* First release on PyPI_.


.. _Python: https://www.python.org/downloads/
.. _tensorflow: https://www.tensorflow.org/
.. _numpy: https://numpy.org/
.. _pandas: https://pandas.pydata.org/
.. _h5py: https://www.h5py.org/
.. _`hdf5-formatted`: https://www.hdfgroup.org/solutions/hdf5/
.. _`loom-formatted`: http://loompy.org/
.. _`Data pre-processing`: https://github.com/iyhaoo/DISC/blob/master/reproducibility/tutorials/data_preprocessing.ipynb
.. _`Imputation`: https://github.com/iyhaoo/DISC/blob/master/reproducibility/tutorials/imputation.md
.. _`Structure recovery comparing with FISH`: https://github.com/iyhaoo/DISC#
.. _`Downsampling recovery`: https://github.com/iyhaoo/DISC#
.. _`Clustering improvement`: https://github.com/iyhaoo/DISC#
.. _`Ultra-large dataset analysis`: https://github.com/iyhaoo/DISC#
.. _`Download datasets we used`: https://github.com/iyhaoo/DISC#
.. _`Other analysis scripts we used`: https://github.com/iyhaoo/DISC#
.. _PyPI: https://pypi.org/project/disc/
