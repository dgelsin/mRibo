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

gff_settings = {}
gff_settings['path_out']         = 0             # set to 0 to output in annotations folder
gff_settings['feat_of_interest'] = 'all'         # all, CDS, tRNA, rRNA: recommend using all
gff_settings['name_qual']        = 'Name'        # GFF gene qualifier
gff_settings['name_qual_alt']    = 'ID'          # Secondary GFF gene qualifier if no name is present
gff_settings['remove_genes']     = 'yes'         # remove hard to align genes listed in bad_genes.csv
gff_settings['gff_extra']         = 50           # additional sequence upstream and downstream of gene (set to 50)

GFF_conversion = ribo_util.GFF_to_dict(paths_in, gff_settings)