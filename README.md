# <div align="center"><ins><strong><em>*M*</em></strong></ins>icrobal <ins><strong><em>*Ribo*</em></strong></ins>some Profiling Analysis (<ins><strong><em>*mRibo*</em></strong></ins>)</div>
Wrapper for running metagene analysis on microbial (Bacteria and Archaea) ribosome profiling data using modified scripts from [Dr. Fuad Mohammad](https://github.com/dgelsin/2018_Bacterial_Pipeline_riboseq).

### INSTALLATION:
To install locally:
```
git clone https://github.com/dgelsin/mRibo
```
then to make it executable from anywhere:
```
cp mRibo/mRibo/bin/* /usr/local/bin/
```

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
pip install plotly==3.10.0
pip install seaborn
conda install bowtie
```

### GENERAL USAGE OF mRibo:
```
mRibo -U reads.fastq -R rRNA.fa -g genome.fa -a genome.gff -d genome_dict -C name_of_user -M name_of_microbe -o output_directory_name

mRibo core options:
	-U STR		ribosome profiling reads
	-R SRT		Ribosomal RNA and tRNA fasta file for filtering out corresponding reads
	-g STR		genome fasta file
	-a STR		genome annotation gtf/gff file
	-D INT		genome annotation dictionary file
	-C STR      	user name
	-M STR      	name of microbe being analyzed
	-o STR          output directory

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
mRibo -U /Users/DRG/Desktop/mRibo/mRibo/practice_files/FASTQ/DG1_subsample -R /Users/DRG/Desktop/mRibo/mRibo/practice_files/genomes/rRNA.fa -g /Users/DRG/Desktop/mRibo/mRibo/practice_files/genomes/genome.fa -a /Users/DRG/Desktop/mRibo/mRibo/annotations/Volcanii.gff -d /Users/DRG/Desktop/mRibo/mRibo/annotations/Volcanii_dict_for_pause -C DG -M Volcanii -o mRibo_run1
```


### FOR ONLY GENERATING DENSITY FILES:
```
density_wrapper.py [FASTQ_directory] [GFF_file.gff] [GFF_DICTIONARY_file] [ALIGNMENT_directory] [NAME_OF_USER] [NAME_OF_MICROBE] [MIN_LEN] [MAX_LEN] [NUM_THREADS] [NUM_CORES]
```

Example:
```
density_wrapper.py input/FASTQ/ annotations/Volcanii.gff annotations/Volcanii_dict_for_pause input/alignment/ DG Volcanii 10 40 8 4
```


### FOR ONLY METAGENE ANALYSIS AND PLOTTING:
```
main_wrapper.py [GFF_DICTIONARY_file] [NAME_OF_USER] [NAME_OF_MICROBE] [MIN_LEN] [MAX_LEN] [NUM_THREADS] [NUM_CORES] [LENGTH_UPSTREAM] [LENGTH_DOWNSTREAM] [ALIGNMENT_STYLE] [DENSITY_STYLE] [WEIGHT_STYLE] [DISTANCE_GENES] [RPKM_THRESHOLD] [YMAX_AVGGENES] [ASITE_SHIFT] [PAUSE_STYLE] [YMAX_PAUSE]
```

Example:
```
main_wrapper.py annotations/Volcanii_dict_for_pause DG Volcanii 10 40 8 4 50 200 3 rpm yes 50 10 10 -11 aa 10
```

### OUTPUT:

All input and output files generated by mRibo will be placed within the output_directory that you named with the -o option.

In the input/ directory you will find a .sam bowtie alignment to the genome in the alignment/ directory, the reads after rRNA filtering used to map to the genome in the FASTQ/ directory, the rRNA bowtie index and .sam rRNA alignment in the bowtie_rRNA_alignment/ directory, and lastly the genome bowtie index in the bowtie_genome_alignment/ directory.

4 figures are generated as output from mRibo. They can be found in output_directory_name/output/reads/figures/user_name/

1. **aa_pausescore_.pdf** - a plot of pause scores for each amino acid at the A-site (large colored dot), P-site (gray dot), and E-site (black dot). Next are individual plots for each amino acid of read density along the length within ribosome footprints (top plot) with a heatmap underneath plotting footprint size (bottom plot). 
2. **asymmetry_score_.pdf** - boxplot of asymettry score for the library used in the analysis.
3. **avggene_.pdf** - meta-gene plots of ribosome density along ORFs (top plots). Underneath are heatmaps of ribosome footprint lengths and where they map along ORFs (bottom plots). The right plot corresponds to the start of ORFs where 0 is x nt upstream of the TSS (5' UTRs), set by the -lu option, and the TSS is marked by the number set by the -lu option. -ld option sets how far to plot downstream within ORFs. The left plot corresponds to the stop of ORFs where 0 is x nt upstream of the stop codon, set by the -ld option, and the stop codon is marked by the number set by the -ld option. -lu option sets how far to plot downstream of ORF stop codons (3' UTRs).
4. **frame_.pdf** - Various plots of reading frame distribution for all footprint lengths on average (left plot), each individual footprint (middle), and each gene (right plot). 0 corresponds to the 0t position in an in-frame ORF codon, 1 is +1 position, and 2 the +2 position.

Various data files are generated as output from mRibo for visualization and further analysis by the user.
1. **wig files** for viewing ribosome density per nt on a genome viewer are provided in output_directory_name/output/reads/density/density/wigfiles/
2. data used to generate avggene plots are provided in output_directory_name/output/reads/analysis/individual/user_name/avggenes/ for either the start of ORFs (**start_all.csv**) or stop of ORFs (**stop_all.csv**)
3. pause score values tp generate aa_pausescore figures are provided in output_directory_name/output/reads/analysis/individual/user_name/pause_score/ for either amino acids (**aa_scores.csv**) or individual codons (**codon_scores.csv**). Pause scores are provided for the -1, -2, A, P, E, +1, and +2 sites of the ribosome in both csv files. **aa_plot_values.csv** contains pause scores for all amino acids along the length of ORFs.
4. Expression values (rpkm) of genes are found in output_directory_name/output/reads/analysis/individual/user_name/***_genelist.csv**


### ACKNOWLEDGEMENTS:

Author of pipeline: [Diego Rivera Gelsinger](https://github.com/dgelsin)

Institution: Johns Hopkins, [Department of Cell, Molecular, Developmental Biology, and Biophysics](https://cmdb.jhu.edu/)

For general feedback you can contact me at [dgelsin1@jhu.edu](mailto:dgelsin1@jhu.edu). 
