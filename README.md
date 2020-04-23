# <ins>*M*</ins>icrobal <ins>*Ribo*</ins>some Profiling Analysis (<ins>*mRibo*</ins>)
Wrapper for running metagene analysis on microbial (Bacteria and Archaea) ribosome profiling data using modified scripts from [Dr. Fuad Mohammad](https://github.com/dgelsin/2018_Bacterial_Pipeline_riboseq).

### REQUIREMENTS:
```pip install DateTime
pip install multiprocess
pip install bcbio-gff
pip install biopython
conda install pandas
pip install numpy
pip install ipython
pip install collections-extended
pip install -U scikit-learn
conda install bowtie
```

### GENERAL USAGE OF mRibo:

```
mRibo -U reads.fastq -R rRNA.fa -g genome.fa -a genome.gff -d genome_dict -C name_of_user -M name_of_microbe

mRibo core options:
	-U STR		ribosome profiling reads
	-R SRT		Ribosomal RNA and tRNA fasta file for filtering out corresponding reads
	-g STR		genome fasta file
	-a STR		genome annotation gtf/gff file
	-D INT		genome annotation dictionary file
	-C STR      	user name
	-M STR      	name of microbe being analyzed

Additional options:
	-I INT		min read length size for meta-gene analysis (default=10)
	-X INT		max read length size for meta-gene analysis (default=40)
	-p INT		number of threads to use for analysis (default=1)
	-c INT		number of cores to use for analysis (default=1)
	-lu INT		number of nucleotides to include upstream of the TSS (default=50)
	-ld INT		number of nucleotides to include downstream of the the TSS within the ORF, genes shorter than this are excluded (default=200)
	-m INT		style in which footprints are mapped for meta-gene analysis (5 prime or 3 prime) (default=5)
	-d STR		type of density to use for meta-gene plots (rpm or reads)  (default=rpm)
	-w STR		style in which to weigh density at positions within ORFs, either equal weight (yes) or not (no) (default=yes)
	-dg INT		distance between genes, genes closer than this are removed from start and stop, or can choose 'no' (default=50)
	-r INT		rpkm threshold to use for meta-gene analysis, genes below will be removed from average (default=10)
	-y INT		maximum value for y-axis of avggenes meta-gene plot (default=0)
	-s INT		A-site shift value for pause score meta-gene analysis (default=-11)
	-P STR		Style in which to conduct pause score meta-gene analysis, either amino acids (aa) or codons (codon) (default=aa)
	-Y INT		maximum value for y-axis of pause score meta-gene plot (default=0)

	--version | -v	show current mRibo version
```

Example:

```
mRibo -U /Users/DRG/Desktop/mRibo/mRibo/practice_files/FASTQ/DG1_subsample -R /Users/DRG/Desktop/mRibo/mRibo/practice_files/genomes/rRNA.fa -g /Users/DRG/Desktop/mRibo/mRibo/practice_files/genomes/genome.fa -a /Users/DRG/Desktop/mRibo/mRibo/annotations/Volcanii.gff -d /Users/DRG/Desktop/mRibo/mRibo/annotations/Volcanii_dict_for_pause -C DG -M Volcanii
```


### FOR ONLY GENERATING DENSITY FILES:

```
python density_wrapper.py [FASTQ_directory] [GFF_file.gff] [GFF_DICTIONARY_file] [ALIGNMENT_directory] [NAME_OF_USER] [NAME_OF_MICROBE] [MIN_LEN] [MAX_LEN] [NUM_THREADS] [NUM_CORES]
```

Example:

```
python density_wrapper.py FASTQ/ annotations/Volcanii.gff annotations/Volcanii_dict_for_pause alignment/ DG Volcanii 10 40 8 4
```


### FOR ONLY METAGENE ANALYSIS AND PLOTTING:

```
python main_wrapper.py [GFF_DICTIONARY_file] [NAME_OF_USER] [NAME_OF_MICROBE] [MIN_LEN] [MAX_LEN] [NUM_THREADS] [NUM_CORES] [LENGTH_UPSTREAM] [LENGTH_DOWNSTREAM] [ALIGNMENT_STYLE] [DENSITY_STYLE] [WEIGHT_STYLE] [DISTANCE_GENES] [RPKM_THRESHOLD] [YMAX_AVGGENES] [ASITE_SHIFT] [PAUSE_STYLE] [YMAX_PAUSE]
```

Example:

```
python main_wrapper.py annotations/Volcanii_dict_for_pause DG Volcanii 10 40 8 4 50 200 3 rpm yes 50 10 10 -11 aa 10
```


### Acknowledgements

Author of pipeline: [Diego Rivera Gelsinger](https://github.com/dgelsin)

Institution: Johns Hopkins, [Department of Cell, Molecular, Developmental Biology, and Biophysics](https://cmdb.jhu.edu/)

For general feedback you can contact me at [dgelsin1@jhu.edu](mailto:dgelsin1@jhu.edu). 
