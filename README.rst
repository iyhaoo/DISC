DISC
====

|PyPI|

.. |PyPI| image:: https://img.shields.io/pypi/v/DISC.svg
    :target: https://pypi.org/project/disc

A highly scalable and accurate inference of gene expression and structure for single-cell transcriptomes using semi-supervised deep learning.

* Free software: Apache License 2.0

Requirements
------------

- Python_ >=3.6
- TensorFlow_ >=1.13.1,<2.0.0
- numpy_ >=1.14.0
- pandas_ >=0.21.0
- h5py_ >=2.9.0

Installation
------------

- **Install TensorFlow**

  If you have an Nvidia GPU, be sure to install a version of TensorFlow that supports it first -- DISC runs much faster with GPU::

    pip install "tensorflow-gpu>= 1.13.1,<2.0.0"

  We typically tensorflow-gpu==1.13.1.

  Here are requirements for GPU version TensorFlow_::

    * Hardware
        * NVIDIA GPU card with CUDA Compute Capability 3.5 or higher.
    * Software
        * NVIDIA GPU drivers - CUDA 10.0 requires 410.x or higher.
        * CUDA Toolkit - TensorFlow_ supports CUDA 10.0 (TensorFlow >= 1.13.0)
        * CUPTI ships with the CUDA Toolkit.
        * cuDNN SDK (>= 7.4.1)

  See this__ for further information.

      .. __: https://www.tensorflow.org/install/gpu

- **Install DISC with pip**

  To install with ``pip``, run the following from a terminal::

    pip install disc

- **Install DISC from GitHub**

  To clone the repository and install manually, run the following from a terminal::

    git clone git://github.com/iyhaoo/DISC.git

    cd disc

    python setup.py install

Usage
-----

- **Quick Start**

  **(1). How to run DISC**::

      disc \
      --dataset=matrix.loom \
      --out-dir=out_dir

  where ``matrix.loom`` is a `loom-formatted`_ raw count matrix with genes in rows and cells in columns and ``out_dir`` is the target path for output folder.

  **(2). What DISC outputs**:

  * ``log.tsv``: records DISC training information.
  * ``summary.pdf``: shows the fitting line and optimal point and will be updated in real time when DISC is running.
  * ``summary.tsv``: records the raw data in ``summary.pdf``.
  * ``result``: imputaion result folder, which contains:

    * ``imputation.loom``: the imputed matrix with genes in rows and cells in columns.
    * ``feature.loom``: the feature matrix with feature in rows and cells in columns.
    * ``running_info.hdf5``: a `hdf5-formatted`_ file, contains some useful information of ``matrix.loom`` (e.g. library size, the expressed counts and cells for each genes, imputed genes, etc.).

  * ``models``: For every save interval, DISC freezes its parameters into this folder (in `pb`_ format).

- **Data availability**

  The sources of our data are listed here.

  * MELANOMA :
      8,640 cells from the melanoma WM989 cell line were sequenced
      using Drop-seq, where 32,287 genes were detected (`scRNA-seq`__).
      In addition, RNA FISH experiment of across 7,000-88,000 cells
      from the same cell line was conducted and 26 genes were detected (`FISH`__).

        .. __: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99330
        .. __: https://www.dropbox.com/s/ia9x0iom6dwueix/fishSubset.txt?dl=0

  * SSCORTEX :
      Mouse somatosensory cortex of CD-1 mice at age of p28 and p29
      were profiled by 10X where 7,477 cells were detected (`scRNA-seq`__).
      In addition, osmFISH experiment of 4,839 cells from somatosensory
      cortex, hippocampus and ventricle of a CD-1 mouse at age of p22 was
      conducted and 33 genes were detected (`FISH`__).

        .. __: http://loom.linnarssonlab.org/clone/Mousebrain.org.level1/L1_Cortex2.loom
        .. __: http://linnarssonlab.org/osmFISH/availability/

  * PBMC :
      2,700 freeze-thaw peripheral blood mononuclear cells (PBMC) from
      a healthy donor were profiled by 10X, where 32,738 genes
      were detect (`scRNA-seq`__).

        .. __: https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/frozen_pbmc_donor_a

  * CBMC :
      Cord blood mononuclear cells were profiled by CITE-seq, where
      8,005 human cells were detected in total (`scRNA-seq`__).

        .. __: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100866

  * JURKAT_293T :
      3258 jurkat cells (`scRNA-seq`__) and 2885 293T cells
      (`scRNA-seq`__) were profiled by 10X separately.
      This dataset has bulk RNA-seq data (`bulk RNA-seq`__).

        .. __: https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/jurkat
        .. __: https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/293t
        .. __: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129240


  * 10X_5CL :
      5,001 cells from 5 human lung adenocarcinoma cell lines H2228,
      H1975, A549, H838 and HCC827 were profiled by 10X (`scRNA-seq`__).
      This dataset has bulk RNA-seq data (`bulk RNA-seq`__).

        .. __: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126906
        .. __: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86337


  * BONE_MARROW :
      6,941 human bone marrow cells from sample MantonBM6 were profiled by 10X.
      The original single-cell RNA sequencing data provided by `HCA`__ was
      aligned to hg19, 6939 cells left after cell filtering (`scRNA-seq`__).
      This dataset has bulk RNA-seq data (`bulk RNA-seq`__).

        .. __: https://data.humancellatlas.org/explore/projects/cc95ff89-2e68-4a08-a234-480eca21ce79
        .. __: https://doc-04-6g-docs.googleusercontent.com/docs/securesc/rm132bl2k8nvnlftqa8a8d5p239lbngf/6o5dsruhjpmecgnkd0nn4b1ak3ss8ufd/1588554075000/07888005335114604629/01857410241295225190/1euh8YB8ThSLHJNQMTCuuKp_nRiME1KzN?e=download&authuser=0&nonce=7apqnnaq9bch8&user=01857410241295225190&hash=a60rd66gq56e0af1vc5ua60146t3gq7m
        .. __: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74246

  * RETINA :
      Retinas of mice at age of p14 were profiled in 7 different replicates
      on by Drop-seq, where 6,600, 9,000, 6,120, 7,650, 7,650, 8280, and
      4000 (49,300 in total) STAMPs (single-cell transcriptomes attached
      to micro-particles) were collected (`scRNA-seq`__). The dataset has
      `cell annotation`__.

        .. __: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63472
        .. __: http://mccarrolllab.org/wp-content/uploads/2015/05/retina_clusteridentities.txt

  * BRAIN_SPLiT :
      156,049 mice nuclei from developing brain and spinal cord at
      age of p2 or p11 mice were profiled by SPLiT-seq (`scRNA-seq`__).
      The cell annotation of this dataset is included in file
      GSM3017261_150000_CNS_nuclei.mat.gz at the same GEO page.

        .. __: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110823

  * BRAIN_1.3M :
      1,306,127 cells from combined cortex, hippocampus,
      and subventricular zone of 2 E18 C57BL/6 mice were
      profiled by 10X (`scRNA-seq`__).

        .. __: https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.3.0/1M_neurons

  We provide our pre-processed data here__.

    .. __: https://github.com/iyhaoo/DISC_data_availability

  +------------------+------------------+------------------+------------------+------------------+------------------+
  |Dataset           |Raw Data          |DS Data           |FISH Data         |Bulk Data         |Cell Annotation   |
  +==================+==================+==================+==================+==================+==================+
  |`MELANOMA`__      |`YES`__           |`0.5`__           |`YES`__           |NO                |NO                |
  +------------------+------------------+------------------+------------------+------------------+------------------+
  |`SSCORTEX`__      |`YES`__           |`0.5`__           |`YES`__           |NO                |NO                |
  +------------------+------------------+------------------+------------------+------------------+------------------+
  |`PBMC`__          |`YES`__           |`0.3`__, `0.5`__  |NO                |NO                |`YES`__           |
  +------------------+------------------+------------------+------------------+------------------+------------------+
  |`CBMC`__          |`YES`__           |`0.3`__, `0.5`__  |NO                |NO                |NO                |
  +------------------+------------------+------------------+------------------+------------------+------------------+
  |`JURKAT_293T`__   |`YES`__           |NO                |NO                |`YES`__           |NO                |
  +------------------+------------------+------------------+------------------+------------------+------------------+
  |`10X_5CL`__       |`YES`__           |NO                |NO                |`YES`__           |NO                |
  +------------------+------------------+------------------+------------------+------------------+------------------+
  |`BONE_MARROW`__   |`YES`__           |NO                |NO                |`YES`__           |`YES`__           |
  +------------------+------------------+------------------+------------------+------------------+------------------+
  |`RETINA`__        |`YES`__           |`0.3`__, `0.5`__  |NO                |NO                |`YES`__           |
  +------------------+------------------+------------------+------------------+------------------+------------------+
  |`BRAIN_SPLiT`__   |`YES`__           |`0.3`__, `0.5`__  |NO                |NO                |`YES`__           |
  +------------------+------------------+------------------+------------------+------------------+------------------+
  |BRAIN_1.3M        |NO (Too large)    |NO                |NO                |NO                |NO                |
  +------------------+------------------+------------------+------------------+------------------+------------------+

  .. __: https://github.com/iyhaoo/DISC_data_availability/tree/master/MELANOMA
  .. __: https://github.com/iyhaoo/DISC_data_availability/blob/master/MELANOMA/raw.loom
  .. __: https://github.com/iyhaoo/DISC_data_availability/tree/master/MELANOMA/ds_0.5
  .. __: https://github.com/iyhaoo/DISC_data_availability/blob/master/MELANOMA/fish.loom
  .. __: https://github.com/iyhaoo/DISC_data_availability/tree/master/SSCORTEX
  .. __: https://github.com/iyhaoo/DISC_data_availability/blob/master/SSCORTEX/raw.loom
  .. __: https://github.com/iyhaoo/DISC_data_availability/tree/master/SSCORTEX/ds_0.5
  .. __: https://github.com/iyhaoo/DISC_data_availability/blob/master/SSCORTEX/fish.loom
  .. __: https://github.com/iyhaoo/DISC_data_availability/tree/master/PBMC
  .. __: https://github.com/iyhaoo/DISC_data_availability/blob/master/PBMC/raw.loom
  .. __: https://github.com/iyhaoo/DISC_data_availability/tree/master/PBMC/ds_0.3
  .. __: https://github.com/iyhaoo/DISC_data_availability/tree/master/PBMC/ds_0.5
  .. __: https://github.com/iyhaoo/DISC_data_availability/blob/master/PBMC/cell_type.rds
  .. __: https://github.com/iyhaoo/DISC_data_availability/tree/master/CBMC
  .. __: https://github.com/iyhaoo/DISC_data_availability/blob/master/CBMC/raw.loom
  .. __: https://github.com/iyhaoo/DISC_data_availability/tree/master/CBMC/ds_0.3
  .. __: https://github.com/iyhaoo/DISC_data_availability/tree/master/CBMC/ds_0.5
  .. __: https://github.com/iyhaoo/DISC_data_availability/tree/master/JURKAT_293T
  .. __: https://github.com/iyhaoo/DISC_data_availability/blob/master/JURKAT_293T/raw.loom
  .. __: https://github.com/iyhaoo/DISC_data_availability/blob/master/JURKAT_293T/bulk.loom
  .. __: https://github.com/iyhaoo/DISC_data_availability/tree/master/10X_5CL
  .. __: https://github.com/iyhaoo/DISC_data_availability/blob/master/10X_5CL/raw.loom
  .. __: https://github.com/iyhaoo/DISC_data_availability/blob/master/10X_5CL/bulk.loom
  .. __: https://github.com/iyhaoo/DISC_data_availability/tree/master/BONE_MARROW
  .. __: https://github.com/iyhaoo/DISC_data_availability/blob/master/BONE_MARROW/raw.loom
  .. __: https://github.com/iyhaoo/DISC_data_availability/blob/master/BONE_MARROW/bulk.loom
  .. __: https://github.com/iyhaoo/DISC_data_availability/blob/master/BONE_MARROW/cell_type.rds
  .. __: https://github.com/iyhaoo/DISC_data_availability/tree/master/RETINA
  .. __: https://github.com/iyhaoo/DISC_data_availability/blob/master/RETINA/raw.loom.gz
  .. __: https://github.com/iyhaoo/DISC_data_availability/tree/master/RETINA/ds_0.3
  .. __: https://github.com/iyhaoo/DISC_data_availability/tree/master/RETINA/ds_0.5
  .. __: https://github.com/iyhaoo/DISC_data_availability/blob/master/RETINA/cell_type.rds
  .. __: https://github.com/iyhaoo/DISC_data_availability/tree/master/BRAIN_SPLiT
  .. __: https://github.com/iyhaoo/DISC_data_availability/blob/master/BRAIN_SPLiT
  .. __: https://github.com/iyhaoo/DISC_data_availability/tree/master/BRAIN_SPLiT/ds_0.3
  .. __: https://github.com/iyhaoo/DISC_data_availability/tree/master/BRAIN_SPLiT/ds_0.5
  .. __: https://github.com/iyhaoo/DISC_data_availability/blob/master/BRAIN_SPLiT/cell_type.rds

- Evaluations

  * Data Preparation, Imputation and Computational Resource Evaluation

    * Data Pre-processing

      +------------------+------------------+------------------+------------------+------------------+
      |`MELANOMA`__      |`SSCORTEX`__      |`PBMC`__          |`CBMC`__          |`JURKAT_293T`__   |
      +------------------+------------------+------------------+------------------+------------------+
      |`10X_5CL`__       |`BONE_MARROW`__   |`RETINA`__        |`BRAIN_SPLiT`__   |`BRAIN_1.3M`__    |
      +------------------+------------------+------------------+------------------+------------------+

      .. __: https://nbviewer.jupyter.org/github/iyhaoo/DISC/blob/master/reproducibility/Data%20Preparation%2C%20Imputation%20and%20Computational%20Resource%20Evaluation/Data%20Pre-processing/MELANOMA.ipynb
      .. __: https://nbviewer.jupyter.org/github/iyhaoo/DISC/blob/master/reproducibility/Data%20Preparation%2C%20Imputation%20and%20Computational%20Resource%20Evaluation/Data%20Pre-processing/SSCORTEX.ipynb
      .. __: https://raw.githack.com/iyhaoo/DISC/master/reproducibility/Data%20Preparation,%20Imputation%20and%20Computational%20Resource%20Evaluation/Data%20Pre-processing/PBMC.nb.html
      .. __: https://raw.githack.com/iyhaoo/DISC/master/reproducibility/Data%20Preparation,%20Imputation%20and%20Computational%20Resource%20Evaluation/Data%20Pre-processing/CBMC.nb.html
      .. __: https://raw.githack.com/iyhaoo/DISC/master/reproducibility/Data%20Preparation,%20Imputation%20and%20Computational%20Resource%20Evaluation/Data%20Pre-processing/JURKAT_293T.nb.html
      .. __: https://raw.githack.com/iyhaoo/DISC/master/reproducibility/Data%20Preparation,%20Imputation%20and%20Computational%20Resource%20Evaluation/Data%20Pre-processing/10X_5CL.nb.html
      .. __: https://raw.githack.com/iyhaoo/DISC/master/reproducibility/data_preparation_and_imputation/data_preprocessing_RETINA.nb.html
      .. __: https://nbviewer.jupyter.org/github/iyhaoo/DISC/blob/master/reproducibility/data_preparation_and_imputation/data_preprocessing_BRAIN_SPLiT.ipynb
      .. __: https://raw.githack.com/iyhaoo/DISC/master/reproducibility/data_preparation_and_imputation/data_preprocessing_10X_5CL.nb.html
      .. __: https://raw.githack.com/iyhaoo/DISC/master/reproducibility/data_preparation_and_imputation/data_preprocessing_BONE_MARROW.nb.html

    * Imputation

    * Computational Resource Evaluation

  * Data Structure Evaluation

    (1). Gene Expression Structures (FISH)
      Tutorial : `MELANOMA`__

       .. __: https://raw.githack.com/iyhaoo/DISC/master/reproducibility/gene_expression/Gene_expression_structures_recovery_validated_by_FISH_MELANOMA.nb.html

    (2). Gene and Cell Structures (Down-sampling)
      Tutorial : `MELANOMA`__

        .. __: https://raw.githack.com/iyhaoo/DISC/master/reproducibility/gene_expression/Dropout_event_recovery_MELANOMA.nb.html

    (S1). Spearman Correlation (Bulk)
      Tutorial : `JURKAT_293T`__

        .. __: https://raw.githack.com/iyhaoo/DISC/master/reproducibility/gene_expression/Dropout_event_recovery_MELANOMA.nb.html

    (S2). Identification of True Zeros (Down-sampling)
      Tutorial : `MELANOMA`__

        .. __: https://raw.githack.com/iyhaoo/DISC/master/reproducibility/gene_expression/Dropout_event_recovery_MELANOMA.nb.html

  * Down-stream Analysis:

    (1). Cell Type Identification (Down-sampling)
      Tutorial : `PBMC`__

        .. __: https://raw.githack.com/iyhaoo/DISC/master/reproducibility/gene_expression/Dropout_event_recovery_MELANOMA.nb.html

    (2). DEG Identification (Bulk)
      Tutorial : `JURKAT_293T`__

        .. __: https://raw.githack.com/iyhaoo/DISC/master/reproducibility/gene_expression/Dropout_event_recovery_MELANOMA.nb.html

    (3). Solution for Large Dataset Analysis
      Tutorial : `PBMC`__

        .. __: https://raw.githack.com/iyhaoo/DISC/master/reproducibility/cell_type_identification/Use_DISC_compressed_features_for_Seurat_clustering_PBMC.nb.html

    (S1). Trajectory Analysis
      Tutorial : `BONE_MARROW`__

        .. __: https://raw.githack.com/iyhaoo/DISC/master/reproducibility/cell_type_identification/Use_DISC_compressed_features_for_Seurat_clustering_PBMC.nb.html

  * Other Utility Scripts

    +------------------+------------------+------------------+
    |Script            |Output                               |
    +==================+==================+==================+
    |`Violin Plot`__   |`PBMC`__          |`RETINA`__        |
    +------------------+------------------+------------------+

        .. __: https://github.com/iyhaoo/DISC/blob/master/reproducibility/gene_expression/violin_plot.py
        .. __: https://github.com/iyhaoo/DISC/blob/master/reproducibility/results/PBMC/violin_plot.pdf
        .. __: https://github.com/iyhaoo/DISC/blob/master/reproducibility/results/RETINA/violin_plot.pdf

References
----------
Yao He\ :sup:`#`, Hao Yuan\ :sup:`#`, Cheng Wu\ :sup:`#`, Zhi Xie\ :sup:`*`.
**"DISC: a highly scalable and accurate inference of gene expression and structure for single-cell transcriptomes using semi-supervised deep learning"**

History
-------

1.0 (2019-12-16)
^^^^^^^^^^^^^^^^^^
* First release on PyPI_.


.. _Python: https://www.python.org/downloads/
.. _TensorFlow: https://www.tensorflow.org/
.. _numpy: https://numpy.org/
.. _pandas: https://pandas.pydata.org/
.. _h5py: https://www.h5py.org/
.. _`hdf5-formatted`: https://www.hdfgroup.org/solutions/hdf5/
.. _`Data availability`: https://github.com/iyhaoo/DISC_data_availability/
.. _`loom-formatted`: http://loompy.org/
.. _`pb`: https://www.tensorflow.org/guide/saved_model/
.. _`RDS-formatted`: https://stat.ethz.ch/R-manual/R-devel/library/base/html/readRDS.html
.. _`Run imputation`: https://github.com/iyhaoo/DISC/blob/master/reproducibility/data_preparation_and_imputation/run_imputation.md
.. _PyPI: https://pypi.org/project/disc/
