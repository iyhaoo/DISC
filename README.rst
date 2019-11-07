DISC
=========

|PyPI|

.. |PyPI| image:: https://img.shields.io/pypi/v/DISC.svg
    :target: https://pypi.org/project/disc

A scalable deep learning imputation method with semi-supervised learning

* Free software: Apache License 2.0

..
 * Documentation: https://scvi.readthedocs.io.

Requirements
----
- Python >= 3.6
- tensorflow >= 1.13.1
- numpy >= 1.14.0
- pandas >= 0.21.0
- h5py >= 2.9.0

Quick Start
-----------

1. Install Python_. We typically use Python3.6 in Linux CentOS 7.

.. _Python: https://www.python.org/downloads/

2. Install DISC through pip:

``pip install DISC``

  Alternatively, you may clone this repository and run ``python setup.py install``.

3. Usage

  1.Quick Start

We use the melanoma dataset (GSE99330) here as an example.

The following code runs DISC:

``DISC --dataset=matrix.loom --out-dir={output_folder}``

where matrix.loom is a loom-formatted raw count matrix with genes in rows and cells in columns.

Results

output_folder contains the main output file (representing the mean parameter of ZINB distribution) as well as some additional matrices in TSV format:

mean.tsv is the main output of the method which represents the mean parameter of the ZINB distribution. This file has the same dimensions as the input file (except that the zero-expression genes or cells are excluded). It is formatted as a gene x cell matrix. Additionally, mean_norm.tsv file contains the library size-normalized expressions of each cell and gene. See normalize_total function from Scanpy for the details about the default library size normalization method used in DCA.

pi.tsv and dispersion.tsv files represent dropout probabilities and dispersion for each cell and gene. Matrix dimensions are same as mean.tsv and the input file.

reduced.tsv file contains the hidden representation of each cell (in a 32-dimensional space by default), which denotes the activations of bottleneck neurons.

Use -h option to see all available parameters and defaults.

4. Follow along with our Jupyter notebooks to quickly get familiar with DISC!

   1. Getting started:
       * `Running melanoma dataset`_

   2. Reproducing our results:
       * `Structure recovery comparing with FISH`_
       * `Downsampling recovery`_
       * `Clustering improvement`_
       * `Ultra-large dataset analysis`_

   3. Supplementary topics:
       * `Download datasets we used`_
       * `Other analysis scripts we used`_
..
   4. Advanced topics:

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
 
History
=========

1.0.0 (2019-11-XX)
------------------

* First release on PyPI.
