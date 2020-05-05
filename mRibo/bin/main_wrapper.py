#!/usr/bin/env python2

import ribo_util
import ribo_main
import ribo_analysis
import ribo_plot
import ribo_plot_with_legends
import sys

print   'main_wrapper.py [GFF_DICTIONARY_file] [NAME_OF_USER] [NAME_OF_MICROBE] [MIN_LEN] [MAX_LEN] [NUM_THREADS] [NUM_CORES] [LENGTH_UPSTREAM] [LENGTH_DOWNSTREAM] [ALIGNMENT_STYLE] [DENSITY_STYLE] [WEIGHT_STYLE] [DISTANCE_GENES] [RPKM_THRESHOLD] [YMAX_AVGGENES] [ASITE_SHIFT] [PAUSE_STYLE] [AMINO_ACID_NAMES] [CODON_NAMES] [YMAX_PAUSE]'

creator=str(sys.argv[2])
microbe=str(sys.argv[3])
min_length=int(sys.argv[4])
max_length=int(sys.argv[5])
threads=int(sys.argv[6])
cores=int(sys.argv[7])
length_upstream=int(sys.argv[8])
length_downstream=int(sys.argv[9])
alignment_style=int(sys.argv[10])
density_style=str(sys.argv[11])
weight_style=str(sys.argv[12])
distance_genes=int(sys.argv[13])
rpkm_threshold=int(sys.argv[14])
ymax_avggenes=int(sys.argv[15])
A_site_shift=int(sys.argv[16])
Pause_style=str(sys.argv[17])
Amino_acid_names=str(sys.argv[18])
Codon_names=str(sys.argv[19])
ymax_pause=int(sys.argv[20])

library_creator = creator        #FM, KS, CW, Menkin, Li, etc...
organism        = microbe      #Coli, Subtilis, Tuberculosis, Salmonella etc...

inputs = {}
#inputs['files']        = ['CW63']
inputs['files']        = [library_creator + str(i) for i in range(1, 2)]

# CPU information for multithreading applications
inputs['multiprocess'] = 'yes'
inputs['threads']      = threads   
inputs['cores']        = cores

# paths containing densities and annotations
path_pc     = 'output/'
inpath      = path_pc + 'reads/'
path_script = 'scripts/'

paths_in = {}
paths_in['path_gff_dict'] = sys.argv[1]

paths_out = {}
paths_out['path_density']      = inpath  + 'density/density/'
paths_out['path_log']          = inpath  + 'density/logs/'
paths_out['path_analysis_log'] = inpath  + 'analysis/logs/'
paths_out['path_analysis']     = inpath  + 'analysis/individual/'
paths_out['path_figures']      = inpath  + 'figures/'


# Check inputs, create output paths
step = 'analysis'
ribo_util.check_inputs(inputs, paths_out, step)
ribo_util.createpath(inputs, paths_out)

settings = {}
settings['minlength'] = min_length
settings['maxlength'] = max_length
settings['shift']     = 11
settings['gff_extra'] = length_upstream
settings['alignment'] = alignment_style

if not 'gff_dict' in globals(): 
    gff_dict, plus_dict, minus_dict = ribo_util.loadlargePickles(inputs, settings, paths_in, paths_out)

settings['length_out_ORF']  = length_upstream
settings['length_in_ORF']   = length_downstream        # genes shorter than this are excluded
settings['density_type']    = density_style    # 'reads' or 'rpm' 
settings['equal_weight']    = weight_style      # 'yes' or 'no', if yes, change density_type to reads -- faster
settings['next_gene']       = distance_genes         # genes closer than this are removed from start and stop, or 'no'
settings['threshold']       = rpkm_threshold       # RPKM, genes below will be removed from average
settings['subgroup'] = 'none'    #provide list of genes or 'none'

average_gene       = ribo_analysis.avggenes(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict)
asymmetry_analysis = ribo_analysis.asymmetry(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict)
frame_analysis     = ribo_analysis.frame(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict)

settings_plot = {}
settings_plot['HM_max'] = 1
settings_plot['ymax']   = ymax_avggenes # set to 0 for autoscaling
settings_plot['shift']  = 1

average_plot = ribo_plot.plot_avggene(inputs, paths_in, paths_out, settings, settings_plot)

frame_plot = ribo_plot.plot_frame(inputs, paths_in, paths_out, settings, settings_plot)

settings['A_site shift']    = A_site_shift
settings['plot_upstream']   = 40
settings['plot_downstream'] = 50
settings['start_trim']      = 50
settings['stop_trim']       = 20
settings['frameshift']      = 0
settings['motif_length']    = 9
settings['next_codon']      = 'no'

#motif_analysis = ribo_analysis.motif_pausescore(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict)
pausescore_analysis = ribo_analysis.pausescore(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict)

settings_plot['aa_or_codon'] = Pause_style 
settings_plot['amino_acid']  = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
settings_plot['codon']       = ['ATG','GTG','TTG']

settings_plot['ymax_dot']  = ymax_pause
settings_plot['ymax_line'] = 6
settings_plot['vmax_HM']   = 6

plot_pausescore = ribo_plot.plot_pausescore(inputs, paths_in, paths_out, settings, settings_plot)

genelists = ribo_analysis.genelist(inputs, paths_out, settings, gff_dict, plus_dict, minus_dict)

plot_pausescore = ribo_plot.plot_asymmetry_comp(inputs, paths_in, paths_out, settings)

average_plot = ribo_plot_with_legends.plot_avggene(inputs, paths_in, paths_out, settings, settings_plot)
#plot_pausescore = ribo_plot_with_legends.plot_pausescore(inputs, paths_in, paths_out, settings, settings_plot)