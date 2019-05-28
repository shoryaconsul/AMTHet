#%% This file runs MixClone using the synthetic data processed by datagen.py and
## feeds into the MixClone pipeline. It is similar to MixClone.py in the MixClone
## package in its function.

from os.path import normpath, dirname
from os.path import join as pjoin
from os import getcwd
#import argparse
from argparse import Namespace

from mixclone.model.run_model import run_model
from mixclone.postprocess.run_postprocess import *
from datagen import data_convert

fnum = 8
#cnvtype = "scnv"
indir = "ThetA_06_bootstrap/python/interval_files/scnv/"

#%% Preprocess synthetic data 
inname = data_convert(indir=indir,fnum=fnum)
inname_split = inname.split('\\')
cnvtype = inname_split[-2]
pathrel = '/'.join(inname_split[0:-3])
in_dir = inname_split[-3]

pathbase = dirname(getcwd())    # Move up one directory
#pathrel = "ThetA_06_bootstrap/python/"
#in_dir = "/interval_files/"
f_str = "/n3"
#inbase = pjoin(pathbase,normpath(pathrel+in_dir+cnvtype+f_str+str(fnum)))
inbase = normpath('/'.join(inname_split[0:-1])+f_str+str(fnum))
out_dir = "/output/"
outbase = normpath('/'.join(inname_split[0:-3])+out_dir+cnvtype+'/'+f_str+str(fnum))
outbase1 = pjoin(pathbase,normpath(pathrel+out_dir+cnvtype+f_str+str(fnum)))

#inpath = normpath(inbase+".MixClone.input.pkl")
#f_in = open(inpath,'r')
#f_in.close()

#%% Run MixClone Model
args_run_model = Namespace(input_filename_base = inbase, output_filename_base=outbase, 
                           max_copynumber=3, subclone_num=3, max_iters=30, stop_value=1e-6, 
                           baseline_thred=0.16)
run_model(args_run_model)

#%% Post-processing
import pickle as pkl


args_postprocess = Namespace(output_filename_base=outbase)

file_name = args_postprocess.output_filename_base + '.MixClone.output.pkl'    
infile = open(file_name, 'rb')

trainer = pkl.load(infile)
data = trainer.data
config_parameters = trainer.config_parameters
model_parameters = trainer.model_parameters
latent_variables = trainer.latent_variables
ll = trainer.ll

extract_paired_counts(data, args_postprocess.output_filename_base)
extract_segments(data, args_postprocess.output_filename_base)


extract_summary(model_parameters, config_parameters,
                ll, args_postprocess.output_filename_base)

infile.close()