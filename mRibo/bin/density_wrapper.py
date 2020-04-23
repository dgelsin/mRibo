#!/usr/bin/env python2

import ribo_util
import ribo_main
import ribo_analysis
import ribo_plot
import sys

print   'density_wrapper.py [FASTQ_directory] [GFF_file.gff] [GFF_DICTIONARY_file] [ALIGNMENT_directory] [NAME_OF_USER] [NAME_OF_MICROBE] [MIN_LEN] [MAX_LEN] [NUM_THREADS] [NUM_CORES]'

creator=str(sys.argv[5])
microbe=str(sys.argv[6])
min_length=int(sys.argv[7])
max_length=int(sys.argv[8])
threads=int(sys.argv[9])
cores=int(sys.argv[10])

'''Settings and Inputs'''

library_creator = creator        #FM, KS, CW, Menkin, Li, etc... (initial of who made it)
organism        = microbe      #Coli, Subtilis, Tuberculosis etc...

inputs = {
    
    'files' : [library_creator + str(i) for i in range(1, 2)],       #Files to analyze
    # for data renaming: useful to rename files from seq_facility - can be ignored
    'order_name' : 'none',    # to rename/concat FASTQ if needed, else set to 'none' or ignore
    
    # select which functions to run: 'yes' and 'no' 
    'run_filtering': 'no',
    'run_bowtie'   : 'no',
    'run_density'  : 'yes',
    'run_readQC'   : 'no',
    
    # cuttoff for readsize and quality for filtering and density
    'minlength'    : min_length,
    'maxlength'    : max_length,
    'phred_cutoff' : 10,
                        
    # linker-1 for FM = CTGTAGGCACCATCAATAGATCGGAAGAGCACACGTCTGAACTCCAGTCA
    # for SM          = CTGTAGGCACCATCAATTCGTATGCCGTCTTCTGCTTG
    # Gross           = CTGTAGGCACCATCAATATCTCGTATGCCGTCTTCTGCTTG
    # Zoya            = ATCTCGTATGCCGTCTTCTGCTTG
    'linker'       :   'CTGTAGGCACCATCAATAGATCGGAAGAGCACACGTCTGAACTCCAGTCA',
                      
    
    # CPU information for multithreading applications, 
    'multiprocess' : 'yes',
    'threads'      : threads,
    'cores'        : cores, 
    
    }

'''Input directories'''
# Can be customized to your liking 

path_pc     = 'output/'
inpath      = path_pc + 'reads/'
path_script = 'scripts/'

paths_in = {
    #'fastq_download' : inpath  + 'FASTQ/downloaded/',
    'path_fastq'     : sys.argv[1],
    'path_gff'       : sys.argv[2],
    'path_gff_dict'  : sys.argv[3],   #will be made from GFF
    #'path_bowtie'    : path_script + 'bowtie/bowtie',
    #'btindex_ladder' : path_script + 'bowtie/indexes/ladder/ladder',
    #'btindex_trna'   : path_script + 'bowtie/indexes/'+organism+'/'+organism+'_tRNA',
    #'btindex_rrna'   : path_script + 'bowtie/indexes/'+organism+'/'+organism+'_rRNA',
    #'btindex_chr'    : path_script + 'bowtie/indexes/'+organism+'/'+organism+'_genome',
    }


### Output directories
paths_out = {
    'path_filter'       : inpath  + 'density/filtering_bowtie/filterdata/',
    #'path_ladder'       : inpath  + 'density/filtering_bowtie/alignments/ladder/',
    #'path_trna'         : inpath  + 'density/filtering_bowtie/alignments/tRNA/',
    'path_rrna'         : inpath  + 'density/filtering_bowtie/alignments/rRNA/',
    'path_chr'          : sys.argv[4],
    'path_temp'         : inpath  + 'density/filtering_bowtie/tmpds/',
    'path_density'      : inpath  + 'density/density/',
    'path_log'          : inpath  + 'density/logs/',
    'path_analysis_log' : inpath  + 'analysis/logs/',
    'path_analysis'     : inpath  + 'analysis/individual/',
    'path_figures'      : inpath  + 'figures/',
    }

gff_settings = {
    'path_out'         : 0,
    'feat_of_interest' : 'protein_coding',         #all, protein_coding, tRNA, rRNA
    'name_qual'        : 'Name',
    'name_qual_alt'    : 'gene_id',
    'biotype_qual'     : 'protein_coding',          #if biotype qualifier NA, biotype_qual = 'all'
    'aSD_seq'          : 'TCCTCC'
    }

# Modify FASTQ files downloaded from server, renaming and concatonating if necessary
#ribo_util.rename_FASTQ(inputs, library_creator, paths_in, paths_out)

# Check inputs, create output paths, and make gff dictionary if needed
step = 'density'
ribo_util.check_inputs(inputs, paths_in, step)
ribo_util.createpath(inputs, paths_out)

densityreads = ribo_main.run_density(inputs, paths_in, paths_out)