#!/usr/bin/env python2

### Settings and Inputs
# to exclude steps, change "y" to "n"
import ribo_util
import ribo_main
import ribo_analysis
import ribo_plot

import pandas as pd
import numpy as np
from math import log
from sklearn.linear_model import LinearRegression
import itertools

import plotly 
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import plot
import cPickle as pickle

import seaborn as sns
import matplotlib.pyplot as plt
from IPython.display import display

#from RNAstructure import RNAstructure

from datetime import datetime
from multiprocessing import Process, Pool
import os, time
import subprocess
import struct
import cPickle as pickle
import csv
from BCBio import GFF
from Bio import Seq
import itertools

import pandas as pd
import numpy as np
from IPython.display import display
import sys

creator=str(sys.argv[3])
microbe=str(sys.argv[4])
outpath=str(sys.argv[5])

library_creator = creator        #FM, KS, CW, Menkin, Li, etc...
organism        = microbe      #Coli, Subtilis, Tuberculosis etc...


### Input directories 
path_pc     = outpath
path_script = 'scripts/'
paths_in = {
    'path_gff'         : sys.argv[1],
    'path_gff_analysis': path_pc +'/'+organism+'_dictionary', 
    'path_badgenes'    : path_pc +'/'+organism+'bad_genes.csv',
    'SD_affinity'      : path_pc +'/'+organism+'bad_genes.csv',
    'path_gff_dict'    : path_pc +'/'+organism+'_dictionary',   #will be made from GFF
#    'SD_affinity'      : path_script + 'scripts/SD_CACCUCCU/sequences_octamers.csv',
#    'path_bowtie'      : path_script + 'scripts/ribo_seq/bowtie/bowtie',
    'path_new_seq'     : sys.argv[2],
    }


### Output directories
paths_out = {
    'path_log'       : path_pc  + 'density/logs/',
    'path_analysis'  : path_pc  + 'analysis/individual/',
    'path_figures'   : path_pc  + 'figures/',
    }



#def shine_dalgarno_affinity(aSD_seq, SD_seq):
    
#    '''Calculate diplex energy (Kcal/mol) using RNAstructure:
#    Dave Mathews Lab (http://rna.urmc.rochester.edu/RNAstructure.html)'''

#    aSD_seq = str(aSD_seq)
#    SD_seq  = str(SD_seq)
    
#    RNA_prog = RNAstructure(exe_path="/Users/fuad/anaconda/bin/RNAstructure/exe/")
#    energy = RNA_prog.DuplexFold(aSD_seq,SD_seq)
    
#    return energy


def GC_of_CDS(CDS_seq):
    
    '''Calculate nucleotide composition of each alias in GFF'''
    CDS_length = len(CDS_seq)
    
    G = 0 
    C = 0 
    A = 0 
    T = 0
    
    for nucleotide in CDS_seq:
        if nucleotide == 'G':
            G += 1
        elif nucleotide == 'C':
            C += 1
        elif nucleotide == 'A':
            A += 1
        elif nucleotide == 'T':
            T += 1
    
    G = float(G) / CDS_length * 100
    C = float(C) / CDS_length * 100
    A = float(A) / CDS_length * 100
    T = float(T) / CDS_length * 100
    
    return G, C, A, T

#def uorf_finder

    '''Find uORFs upstream of annotated genes
        - start from annotated start codons
            find upstream AUG within 100 nt from start in any frame
            find in frame stop and length 
            find SD strength          
    
    '''

    

def GFF_to_dict(paths_in, gff_settings):
    
    '''Parse gff into dict:
        - feat_of_interest = what to look for in gff (protein_coding, tRNA, rRNA, etc)
        - name_qual        = qualifier for alias/gene name (Name, gene_id)
        - name_qual_alt    = alternative qualifier, if none, set as 'none' 
        - biotype_qual     = qualifier for type of feature (biotype, etc)
        
        These values must correspont to values in the GFF'''
    
    '''Unload gff_settings'''
    
    path_out         = gff_settings['path_out']
    feat_of_interest = gff_settings['feat_of_interest']  #all, protein_coding, tRNA, rRNA
    name_qual        = gff_settings['name_qual']
    name_qual_alt    = gff_settings['name_qual_alt']
    remove_genes     = gff_settings['remove_genes']
#    aSD_seq          = gff_settings['aSD_seq']
    path_badgenes    = paths_in['path_badgenes']
    
    '''Output path can be defined, or use 0 to set as the annotation file for my main pipeline'''
    
    if path_out == 0:
        path_gff_dict = paths_in['path_gff_dict']
    else:
        path_gff_dict = path_out
    
    '''Parse GFF using BCBio'''
    
    GFFgen = GFF.parse(paths_in['path_gff'])
    feat_num = 0
    
    '''Define data arrays: will be used as columns for pandas DateFrame'''

    gff_dict   = {}
    chromosome = []
    aliaslist  = []
    startlist  = []
    stoplist   = []
    seqlist    = []
    typelist   = []
    strandlist = []
    startcodon = []
    stopcodon  = []
    SDaffinity = []
    G_content  = []
    C_content  = []
    A_content  = []
    T_content  = []
    
    aa_code, codon_code = ribo_util.get_genetic_code()
    aa_comp_dict = {}
   
    
    '''Make list of bad genes'''
            

# from Gene-Wei-Li 

#    bad_genes = pd.read_csv(path_badgenes)
#    bad_genes = bad_genes.to_dict(orient='list')
#    bad_genes = bad_genes['GeneName']
        
        

    '''Sift through GFF for relevant information'''
    
    for chromosome_number in range(1, 50):
        chr = next(GFFgen, None)
        
        if chr is None:
            break
        
        for feature in chr.features:
            chromosome_id = chr.id

            if feature.sub_features == []:
                    feat_num+=1
                    continue

            if remove_genes == 'yes':

                '''Skip over non-CDS annotations'''

                if not feature.sub_features[0].type == feat_of_interest:
                    feat_num+=1
                    continue        
                elif feature.qualifiers.has_key('pseudo') == True:
                    feat_num+=1
                    continue
                else:
                    feature_type = 'CDS'
            else: 

                '''Add feat type to GFF, noting pseudogenes'''

                if feature.qualifiers.has_key('pseudo') == True:
                    feature_type = 'pseudo'
                else: 
                    feature_type = feature.sub_features[0].type 


            '''Get feature name'''

            if name_qual in feature.qualifiers:
                feat_name = feature.qualifiers[name_qual][0]
            elif name_qual_alt in feature.qualifiers:
                feat_name = feature.qualifiers[name_qual_alt][0]
            else:
                feat_name = 'None'
                feat_num+=1
                continue

            '''Remove feature if bad'''

#            if remove_genes == 'yes':
#                if feat_name in bad_genes:
#                    feat_num+=1
#                    continue
#            else: 
#                if feat_name in bad_genes:
#                    feature_type = 'bad'


            '''Get start, end, and strand position'''

            start  = feature.location.start.position 
            end    = feature.location.end.position
            strand = feature.strand           

            '''Analyze features of interest for feat information'''   


            alias = feat_name

            '''Each strand is treated differently, + strand == 1'''

            if strand == 1:

                '''I save gene sequence + 50 bp from each end:
                makes it easier to analyze start and stop sequence 
                context without using whole genome sequence'''

                if start < 50:      # if gene is near the beginning of genome sequence:
                    sequence = 'N' * (50 - start)                  # TB GFF starts at 0, add N * 50
                    sequence = sequence + chr[0:end+50].seq        # gene sequence + 50nt at each end
                else: 
                    sequence = chr[start-50:end+50].seq            # gene sequence + 50nt at each end

                strand_val = '+'
                startcodon_pos = start
                stopcodon_pos  = end-1

                if start > 200: 
                    upstream_seq = chr[start-200:start+100].seq

            else:

                '''For minus strand, 'end' is start codon, 'start' is stop codon
                and sequence is reverse compliment of gene sequence.'''

                sequence_rc = chr[start-50:end+50].seq
                sequence = sequence_rc.reverse_complement()

                strand_val = '-'
                startcodon_pos = end-1
                stopcodon_pos  = start

                if end + 200 > len(chr.seq):
                    upstream_seq = 'none'
                else:
                    upstream_seq_rc = chr[end-100:end+200].seq
                    upstream_seq = upstream_seq_rc.reverse_complement()




            sequence = str(sequence)
            start_codon = sequence[50:53:1]
            stop_codon  = sequence[-53:-50]



            '''get sequence from start to stop for GC analysis'''

            CDS_seq = sequence[50:-50:1]


            G, C, A, T = GC_of_CDS(CDS_seq)

            '''Calculate SD affinity'''

#            SD_seq = sequence[30:50:1]    # analyze 20 nt upstream of start codons
#            SD_affinity = shine_dalgarno_affinity(aSD_seq, SD_seq)

            '''Append data to lists'''

            if alias == 'trmD':
                print sequence

            chromosome.append(chromosome_id)    
            typelist.append(feature_type)
            aliaslist.append(alias)
            seqlist.append(sequence)
            strandlist.append(strand_val)
            startlist.append(startcodon_pos)
            stoplist.append(stopcodon_pos)
            startcodon.append(start_codon)
            stopcodon.append(stop_codon)
#            SDaffinity.append(SD_affinity)
            G_content.append(G)
            C_content.append(C)
            A_content.append(A)
            T_content.append(T)

            feat_num+=1
        

        
    '''Append lists to gff_dict'''
    gff_dict['Chromosome'] = chromosome
    gff_dict['Alias']      = aliaslist
    gff_dict['Strand']     = strandlist
    gff_dict['Start']      = startlist
    gff_dict['Stop']       = stoplist
    gff_dict['Sequence']   = seqlist
    gff_dict['Start_Codon'] = startcodon
    gff_dict['Stop_Codon']  = stopcodon
    gff_dict['Type']        = typelist
#    gff_dict['SD_affinity'] = SDaffinity
    gff_dict['G_content'] = G_content
    gff_dict['C_content'] = C_content
    gff_dict['A_content'] = A_content
    gff_dict['T_content'] = T_content

           
    '''Pickle dict for use later'''
    ribo_util.makePickle(gff_dict,path_gff_dict)
    
    '''print dataframe, and save as .csv for use later'''
    ## Print GFF to check 
    gff_df = pd.DataFrame(gff_dict)
    display(gff_df)
    gff_df.to_csv(path_gff_dict + '.csv')
    
    return 

def SD_affinity_genome(paths_in):
    '''This function takes octamers of genomic sequence and calculates shine dalgarno affinity:
    output is a size separated dict that can be used like a density dict.'''
    # load octamers
    SD_affinity = paths_in['SD_affinity']
    
    affinity_list = pd.read_csv(SD_affinity)
    affinity_list = pd.Series(affinity_list.SD_affinity.values,index=affinity_list.Octamer).to_dict()
        
    length_range     = range(10, 46) 
    GFFgen = GFF.parse(paths_in['path_gff'])
    chr = GFFgen.next()
    feat_num = 0
    
    affinity_plus = []
    affinity_minus = []
    density_plus_sizesep = {}
    density_minus_sizesep = {}
    
    sequence    = chr.seq
    sequence_rc = sequence.reverse_complement()
    genome_size = len(sequence)
    
    position = 0
    for position in range (0, genome_size):
        if position < 8:
            motif    = 'AAAAAAAA'
            motif_rc = 'AAAAAAAA'
        elif genome_size - position < 8:
            motif_rc = 'AAAAAAAA'
            motif    = 'AAAAAAAA'
        else:
            
            motif    = sequence[position - 8: position].transcribe()
            motif_rc = sequence[position: position + 8].transcribe()
        
            motif_rc = motif_rc.reverse_complement()
        
        if len(motif) == 8 and len(motif_rc) == 8:
            SD_affinity_plus  = affinity_list[motif]  
            SD_affinity_minus = affinity_list[motif_rc]
        else:
            SD_affinity_plus  = 0.0  
            SD_affinity_minus = 0.0
        
        if position == 100000:
            print '100000'
        if position == 500000:
            print '500000'
        if position == 1000000:
            print '1000000'
        if position == 2000000:
            print '2000000'
        
        affinity_plus.append(SD_affinity_plus)
        affinity_minus.append(SD_affinity_minus)
    
    for length in length_range:
        density_plus_sizesep[length]  = affinity_plus
        density_minus_sizesep[length] = affinity_minus
        
    path_den   = inpath  + 'density/density/SD1/'
    ribo_util.makePickle(density_plus_sizesep,path_den+"plus_sizesep")
    ribo_util.makePickle(density_minus_sizesep,path_den+"minus_sizesep")

    return 

gff_settings = {
    'path_out'         : 0,
    'feat_of_interest' : 'CDS',         #all, protein_coding, tRNA, rRNA
    'name_qual'        : 'Name',
    'name_qual_alt'    : 'product',
    'remove_genes'     : 'yes',         #remove hard to align genes, if no, will add annotation 
#    'aSD_seq'          : 'TCCTCC'
    }


GFF_to_dict(paths_in, gff_settings)