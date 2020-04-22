# <ins>*M*</ins>icrobal <ins>*Ribo*</ins>some Profiling Analysis (<ins>*mRibo*</ins>)
Wrapper for running metagene analysis on microbial (Bacteria and Archaea) ribosome profiling data using modified scripts from Dr. Fuad Mohammad https://github.com/dgelsin/2018_Bacterial_Pipeline_riboseq

#### REQUIREMENTS:
```pip install DateTime
pip install multiprocess
pip install bcbio-gff
pip install biopython
conda install pandas
pip install numpy
pip install ipython
pip install collections-extended
pip install -U scikit-learn
```

**1. Must name reads.fastq file to DG1 and put it in the FASTQ directory**

**2. Must name alignment.sam file DG1_match.sam and put it in the alignment directory**


#### FOR GENERATING DENSITY FILES:

```python density_wrapper.py [FASTQ_directory] [GFF_file.gff] [GFF_DICTIONARY_file] [ALIGNMENT_directory]```

Example:

```python density_wrapper.py FASTQ/ annotations/Volcanii.gff annotations/Volcanii_dict_for_pause alignment/```



#### FOR METAGENE ANALYSIS AND PLOTTING:

```python main_wrapper.py [GFF_DICTIONARY_file]```

Example:

```python main_wrapper.py annotations/Volcanii_dict_for_pause```


#### Acknowledgements

Author of pipeline: Diego Rivera Gelsinger

Institution: Johns Hopkins, Department of Cell, Molecular, Developmental Biology, and Biophysics

For general feedback you can contact me at dgelsin1@jhu.edu.
