# File to check MixClone files

import pickle as pkl
from os import getcwd
from os.path import normpath, dirname


pathbase = dirname(getcwd())
fname = "ThetA_06_bootstrap/python/interval_files/scnv/n21.MixClone.input.pkl"
fpath = normpath(pathbase+'/'+fname)
infile = open(fpath, 'rb')
data = pkl.load(infile)
print('DONE')
