#1 Test if there is duplicated genes?
#python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/1_preprocess/loom_rename_duplicated.py \
#--input=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_10000_1k.loom




#2 Run

# 1k
# SAVER
/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/SAVER.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_50000_1k.loom 16 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/SAVER.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_100000_1k.loom 16 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/SAVER.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_500000_1k.loom 16 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/SAVER.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_all_1k.loom 16 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/SAVER.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_2.57M_1k.loom 16 1 -1

# MAGIC
/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/MAGIC.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_50000_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/MAGIC.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_100000_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/MAGIC.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_500000_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/MAGIC.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_all_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/MAGIC.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_2.57M_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

# DCA
/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/DCA.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_50000_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/DCA.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_100000_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/DCA.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_500000_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/DCA.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_all_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/DCA.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_2.57M_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

# scScope
/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scScope.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_50000_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scScope.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_100000_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scScope.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_500000_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scScope.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_all_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scScope.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_2.57M_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

# scVI
/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scVI.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_50000_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scVI.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_100000_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scVI.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_500000_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scVI.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_all_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scVI.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_2.57M_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

# DISC

/home/yuanhao/softwares/memusg python3 /home/yuanhao/DISC/disc/main.py \
--dataset=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_50000_1k.loom \
--out-dir=/home/yuanhao/DISC/performance/50k_1k \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1 \
--repeats=3 \
--library-size-factor=median \
--depth=16_8_1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/DISC/disc/main.py \
--dataset=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_100000_1k.loom \
--out-dir=/home/yuanhao/DISC/performance/100k_1k \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1 \
--repeats=3 \
--library-size-factor=median \
--depth=16_8_1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/DISC/disc/main.py \
--dataset=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_500000_1k.loom \
--out-dir=/home/yuanhao/DISC/performance/500k_1k \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1 \
--repeats=3 \
--library-size-factor=median \
--depth=16_8_1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/DISC/disc/main.py \
--dataset=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_all_1k.loom \
--out-dir=/home/yuanhao/DISC/performance/1.3m_1k \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1 \
--repeats=3 \
--library-size-factor=median \
--depth=16_8_1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/DISC/disc/main.py \
--dataset=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_2.57M_1k.loom \
--out-dir=/home/yuanhao/DISC/performance/2.6m_1k \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1 \
--repeats=3 \
--library-size-factor=median \
--depth=16_8_1



# 10k
# SAVER
/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/SAVER.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_50000_10k.loom 16 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/SAVER.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_100000_10k.loom 16 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/SAVER.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_500000_10k.loom 16 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/SAVER.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_all_10k.loom 16 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/SAVER.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_2.57M_10k.loom 16 1 -1

# MAGIC
/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/MAGIC.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_50000_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/MAGIC.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_100000_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/MAGIC.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_500000_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/MAGIC.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_all_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/MAGIC.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_2.57M_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

# DCA
/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/DCA.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_50000_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/DCA.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_100000_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/DCA.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_500000_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/DCA.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_all_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/DCA.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_2.57M_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

# scScope
/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scScope.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_50000_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scScope.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_100000_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scScope.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_500000_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scScope.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_all_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scScope.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_2.57M_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

# scVI
/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scVI.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_50000_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scVI.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_100000_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scVI.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_500000_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scVI.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_all_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scVI.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_2.57M_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

# DISC

/home/yuanhao/softwares/memusg python3 /home/yuanhao/DISC/disc/main.py \
--dataset=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_50000_10k.loom \
--out-dir=/home/yuanhao/DISC/performance/50k_10k \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1 \
--repeats=3 \
--library-size-factor=median \
--depth=16_8_1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/DISC/disc/main.py \
--dataset=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_100000_10k.loom \
--out-dir=/home/yuanhao/DISC/performance/100k_10k \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1 \
--repeats=3 \
--library-size-factor=median \
--depth=16_8_1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/DISC/disc/main.py \
--dataset=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_500000_10k.loom \
--out-dir=/home/yuanhao/DISC/performance/500k_10k \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1 \
--repeats=3 \
--library-size-factor=median \
--depth=16_8_1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/DISC/disc/main.py \
--dataset=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_all_10k.loom \
--out-dir=/home/yuanhao/DISC/performance/1.3m_10k \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1 \
--repeats=3 \
--library-size-factor=median \
--depth=16_8_1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/DISC/disc/main.py \
--dataset=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_2.57M_10k.loom \
--out-dir=/home/yuanhao/DISC/performance/2.6m_10k \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1 \
--repeats=3 \
--library-size-factor=median \
--depth=16_8_1





# 1k
# VIPER
/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_10000_1k.loom gene 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_50000_1k.loom gene 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_100000_1k.loom gene 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_500000_1k.loom gene 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_all_1k.loom gene 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_2.57M_1k.loom gene 1 -1

# deepImpute
/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_10000_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_50000_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_100000_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_500000_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_all_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_2.57M_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

#scImpute
/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_10000_1k.loom 16 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_50000_1k.loom 16 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_100000_1k.loom 16 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_500000_1k.loom 16 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_all_1k.loom 16 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_2.57M_1k.loom 16 1 -1




# 10k
# VIPER
/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_10000_10k.loom gene 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_50000_10k.loom gene 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_100000_10k.loom gene 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_500000_10k.loom gene 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_all_10k.loom gene 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/VIPER.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_2.57M_10k.loom gene 1 -1

# deepImpute
/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_10000_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_50000_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_100000_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_500000_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_all_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/deepImpute.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_2.57M_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

#scImpute
/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_10000_10k.loom 16 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_50000_10k.loom 16 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_100000_10k.loom 16 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_500000_10k.loom 16 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_all_10k.loom 16 1 -1

/home/yuanhao/softwares/memusg R --slave < /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/scImpute.r \
--args /home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_2.57M_10k.loom 16 1 -1



# MAGIC
/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/MAGIC.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_10000_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/MAGIC.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_50000_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/MAGIC.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_100000_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/MAGIC.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_500000_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/MAGIC.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_all_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/MAGIC.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_2.57M_1k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1


# MAGIC
/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/MAGIC.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_10000_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/MAGIC.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_50000_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/MAGIC.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_100000_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/MAGIC.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_500000_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/MAGIC.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_all_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1

/home/yuanhao/softwares/memusg python3 /home/yuanhao/single_cell/scripts/evaluation_pipeline/2_run_methods/single_run/MAGIC.py \
--loom=/home/yuanhao/data/fn/neuron1.3m/performance_test_set/1M_neurons_2.57M_10k.loom \
--min-expressed-cell=1 \
--min-expressed-cell-average-expression=-1


