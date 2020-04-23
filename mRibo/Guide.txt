######REQUIREMENTS######
pip install DateTime
pip install multiprocess
pip install bcbio-gff
pip install biopython
conda install pandas
pip install numpy
pip install ipython
pip install collections-extended
pip install -U scikit-learn

1. Must name reads.fastq file to DG1 and put it in the FASTQ directory
2. Must name alignment.sam file DG1_match.sam and put it in the alignment directory

######FOR GENERATING DENSITY FILES######

python density_wrapper.py [FASTQ_directory] [GFF_file.gff] [GFF_DICTIONARY file] [ALIGNMENT_directory]

example:
python density_wrapper.py FASTQ/ annotations/Volcanii.gff annotations/Volcanii_dict_for_pause alignment/



######FOR METAGENE ANALYSIS AND PLOTTING######

python main_wrapper.py [GFF_DICTIONARY]

Example:
python main_wrapper.py annotations/Volcanii_dict_for_pause