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
 2. Results:

    * ``log.tsv``: a tsv-formatted log file that records training states.
    * ``summary.pdf``: a pdf-formatted file that visualizes the fitting line and optimal point and it will be updated in real time when running.
    * ``summary.tsv``: a tsv-formatted file that shows the raw data of visualization.
    * ``result``: a directory for imputaion results as below:

      * ``imputation.loom``: a `loom-formatted`_ imputed matrix with genes in rows and cells in columns.
      * ``feature.loom``: a `loom-formatted`_ dimensionally reduced feature matrix provided by our method based on the imputed matrix above with feature in rows and cells in columns.
      * ``running_info.hdf5``: a `hdf5-formatted`_ saved some basic information about the input dataset such as library size, genes used for modelling and so on.

    * ``models``: a directory for trained models in every save interval

Data availability
  * `MELANOMA`_ :
      8,640 cells from the melanoma WM989 cell line were sequenced using Drop-seq, where 32,287 genes were detected.
      Its `loom-formatted`_ file after cell filtering is provide `here`__.

      .. __: https://raw.githubusercontent.com/iyhaoo/DISC/master/reproducibility/data/melanoma/raw.loom

  * `MELANOMA_FISH`_ :
      RNA FISH experiment of across 7,000-88,000 melanoma cells from the same cell line was conducted and 26 genes were detected.
      Its `loom-formatted`_ file is provide `here`__.

      .. __: https://raw.githubusercontent.com/iyhaoo/DISC/master/reproducibility/data/melanoma/fish.loom

  * `SSCORTEX`_ :
      Mouse somatosensory cortex of CD-1 mice at age of p28 and p29 were profiled by 10X where 7,477 cells were detected in total.
      Its `loom-formatted`_ file after cell filtering is provide `here`__.

      .. __: https://raw.githubusercontent.com/iyhaoo/DISC/master/reproducibility/data/sscortex/raw.loom

  * `SSCORTEX_FISH`_ :
      osmFISH experiment of 4,839 cells from somatosensory cortex, hippocampus and ventricle of a CD-1 mouse at age of p22 was conducted and 33 genes were detected.
      Its `loom-formatted`_ file of 4,388 cells from somatosensory cortex is provide `here`__.

      .. __: https://raw.githubusercontent.com/iyhaoo/DISC/master/reproducibility/data/sscortex/fish.loom

  * `PBMC`_ :
      2,700 freeze-thaw peripheral blood mononuclear cells (PBMC) from a healthy donor were profiled by 10X, where 32,738 genes were detect.
      Its `loom-formatted`_ file after cell filtering is provide `here`__.

      .. __: https://raw.githubusercontent.com/iyhaoo/DISC/master/reproducibility/data/pbmc3k/pbmc3k_filtered.loom

  * `CBMC`_ :
      Cord blood mononuclear cells were profiled by CITE-seq, where 8,005 human cells were detected in total.
      Its `loom-formatted`_ file after cell and with all human genes is provide `here`__.

      .. __: https://raw.githubusercontent.com/iyhaoo/DISC/master/reproducibility/data/cite_seq/raw.loom

  * `RETINA`_ :
      Retinas of mice at age of p14 were profiled in 7 different replicates on by Drop-seq, where 6,600, 9,000, 6,120, 7,650, 7,650, 8280, and 4000 (49,300 in total) STAMPs (single-cell transcriptomes attached to micro-particles) were collected with totally 24,658 genes detected.

  * `BRAIN_SPLIT`_ :
      156,049 mice nuclei from developing brain and spinal cord at age of p2 or p11 mice were profiled by SPLiT-seq, where 26,894 genes were detected

  * `BRAIN_1.3M`_ :
      1,306,127 cells from combined cortex, hippocampus, and subventricular zone of 2 E18 C57BL/6 mice were profiled by 10X, where 27998 genes were detected.


Tutorials
 1. Imputation

    * `Data pre-processing`_
    * `Run imputation`_

 2. Reproducing our results:

    * `Gene expression structures recovery validated by FISH`_
    * `Dropout event recovery`_
    * `Cell type identification improvement`_
    * `Ultra-large dataset analysis`_

..
 3. Supplementary topics:

References
----------
Yao He\ :sup:`#`, Hao Yuan\ :sup:`#`, Cheng Wu\ :sup:`#`, Zhi Xie\ :sup:`*`.
**"Reliable and efficient gene expression recovery in single-cell transcriptomes using DISC"**

History
-------

1.0.0 (2019-12-XX)
^^^^^^^^^^^^^^^^^^
* First release on PyPI_.


.. _Python: https://www.python.org/downloads/
.. _tensorflow: https://www.tensorflow.org/
.. _numpy: https://numpy.org/
.. _pandas: https://pandas.pydata.org/
.. _h5py: https://www.h5py.org/
.. _`hdf5-formatted`: https://www.hdfgroup.org/solutions/hdf5/
.. _`loom-formatted`: http://loompy.org/
.. _`Data pre-processing`: https://nbviewer.jupyter.org/github/iyhaoo/DISC/blob/master/reproducibility/tutorials/data_preprocessing.ipynb
.. _`Run imputation`: https://github.com/iyhaoo/DISC/blob/master/reproducibility/tutorials/run_imputation.md
.. _`Gene expression structures recovery validated by FISH`: https://rawcdn.githack.com/iyhaoo/DISC/d207c15306cd6aa73ce492953971d6e84c42fbc3/reproducibility/gene_expression/Gene_expression_structures_recovery_validated_by_FISH.nb.html
.. _`Dropout event recovery`: https://rawcdn.githack.com/iyhaoo/DISC/d207c15306cd6aa73ce492953971d6e84c42fbc3/reproducibility/gene_expression/Dropout_event_recovery.nb.html
.. _`Cell type identification improvement`: https://rawcdn.githack.com/iyhaoo/DISC/d207c15306cd6aa73ce492953971d6e84c42fbc3/reproducibility/cell_type_identification/Cell_type_identification_improvement.nb.html
.. _`Ultra-large dataset analysis`: https://github.com/iyhaoo/DISC#
.. _PyPI: https://pypi.org/project/disc/
.. _MELANOMA: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99330
.. _`the previous pipeline`: https://www.nature.com/articles/s41592-018-0033-z
.. _MELANOMA_FISH: https://www.dropbox.com/s/ia9x0iom6dwueix/fishSubset.txt?dl=0
.. _SSCORTEX: http://loom.linnarssonlab.org/dataset/cellmetadata/Mousebrain.org.level1/L1_Cortex2.loom
.. _SSCORTEX_FISH: http://linnarssonlab.org/osmFISH/availability/
.. _PBMC: https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html
.. _CBMC: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100866
.. _RETINA: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63472
.. _BRAIN_SPLIT: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110823
.. _BRAIN_1.3M: https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.3.0/1M_neurons
