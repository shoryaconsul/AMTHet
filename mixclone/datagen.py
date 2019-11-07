## This file synthesizes synthetic data in line with that produced for AltMin 
# and THetA. The output file from this can be used by the model part of MixClone.
# Written in Python 2.7

import pickle as pkl
import numpy as np
from numpy.random import exponential as exp, random, randint, randn

from mixclone import constants
from mixclone.preprocess.data import Data, Segment
from mixclone.preprocess.utils import get_BAF_counts

from os.path import normpath, dirname
from os.path import join as pjoin
from os import getcwd


def data_convert(indir, fnum, n, k):
    # Arguments
    ## indir: Input directory
    ## fnum: Iteration number
    ## n: Number of subclones
    ## k: Max copy number
    
    pathbase = dirname(getcwd())    # Move up one directory
#    pathrel = "ThetA_06_bootstrap/python/interval_files/"
#    cnvtype = "scnv" 
    pre_str = "/interval_count_n"+str(n)
    pre_gt = "/muC_n"+str(n)
    out_str = "/n"+str(n)
    
    #%% Constants
    snvrate = 1e-3/3  # SNVs every 1000 bases
    k = k # 3   # Max copy number
    cov_sd = 0.02   # Standard dev for variation in coverage
    readlen = 100
        
    #%% Reading and formatting data
    data = Data()
    fname = pjoin(pathbase,normpath(indir+pre_str+str(fnum)))
    gtname = pjoin(pathbase,normpath(indir+pre_gt+str(fnum)))
    f = open(fname,"r")
    next(f) # Skip header
    
    f_gt = open(gtname,"r")
    gt_str = f_gt.readlines()[1]
    f_gt.close()
    
    [mu_gt, cnv_gt] = gt_str.split()
    mu_gt = [float(x) for x in mu_gt.split(',')]   # Normal and tumor fractions
    cnv_gt = cnv_gt.split(':')[:-1]                 # Last element is empty
    cnv_gt_arr = [[int(y) for y in x.split(',')] for x in cnv_gt] # CNVs
    
    nt = len(cnv_gt_arr[0])      # Number of tumor subclones
    
    for s in f:
        sp = [int(x) for x in s.split()]
        chr_idx = sp[1]
        stpos = sp[2]   # Start pos
        endpos = sp[3]  # End pos
        nread = sp[4]   # Normal read num
        tread = sp[5]   # Tumor read num
        i = data.seg_num         # Current segment number
        
        chrom_name = str(chr_idx)
        seg_name = '_'.join([chrom_name, 'start', str(stpos), 'end', str(endpos)])
    
        # Identical to MixClone preprocess (data.py, load_segments)
        seg_i = Segment()
        seg_i.name = chrom_name
        seg_i.name = seg_name
        seg_i.chrom_idx = chr_idx
        seg_i.chrom_name = chrom_name
        seg_i.start = stpos
        seg_i.end = endpos
        seg_i.normal_reads_num = nread
        seg_i.tumor_reads_num = tread
        seg_i.log2_ratio = np.log2(1.0*tread/nread)
        
        #seglen = endpos-stpos
        currpos = stpos
        snv_pos =[]
        while currpos<endpos:
            incr = int(round(exp(1/snvrate)))
            if currpos+incr>=endpos:
                break
            else:
                currpos = currpos+incr
                snv_pos.append(currpos) 
        
        seg_baf = baf_gen(len(snv_pos))
        norm_ref = nread*seg_baf    # Ref allele count for normal
        norm_nref = nread-norm_ref  # Non ref allele count for normal
        
        
        seg_baf = baf_gen(len(snv_pos))
        cnv_seg = cnv_gt_arr[i]
        if sum(cnv_seg)== 2*nt:      # No CNVs in segment
            tumor_ref = tread*seg_baf # Ref allele count for tumor
            tumor_nref = tread-tumor_ref # Non ref allele count for tumor
        else:                       # CNVs in segment
            #
            baf_ref = np.array([mu_gt[0]*seg_baf])
            baf_nref = np.array([mu_gt[0]*(1-seg_baf)])
            for j in range(nt):
                cnv_j = cnv_seg[j] # Copy number for subclone
                if cnv_j == 0:
                    pass
                elif cnv_j == 1:
                    ref_flip = random>0.5 # Pick ref or nonref allele for CNV = 1
                    baf_ref = baf_ref+ref_flip*mu_gt[j+1]*seg_baf
                    baf_nref = baf_nref+(1-ref_flip)*mu_gt[j+1]*(1-seg_baf)
                elif cnv_j==2:
                    baf_ref = baf_ref+mu_gt[j+1]*seg_baf
                    baf_nref = baf_nref+mu_gt[j+1]*(1-seg_baf)
                else:
                    ref_num = randint(cnv_j) # Pick number of copies of ref allele
                    baf_ref = baf_ref+mu_gt[j+1]*ref_num*seg_baf
                    baf_nref = baf_nref+mu_gt[j+1]*(cnv_j-ref_num)*(1-seg_baf)
                
            seg_baf = baf_ref/(baf_ref+baf_nref) # Normalizing BAF
            tumor_ref = tread*baf_ref
            tumor_nref = tread*baf_nref
                    
        # Splitting reads between positions
        seglen = endpos-stpos+1
        read_mult = (1.0+cov_sd*randn(len(snv_pos)))*readlen/seglen
        read_mult_mat = np.matlib.repmat(read_mult,4,1)
        read_mult_mat = read_mult_mat.transpose()
        read_mult_mat = np.hstack((read_mult_mat,np.ones((len(snv_pos),2))))
        
        # Casting counts in same format as MixClone
        paired_counts = np.array([[], [], [], [], [], []], dtype=int).transpose()
        paired_counts = np.vstack((norm_ref,norm_nref,tumor_ref,tumor_nref,
                                   chr_idx*np.ones(len(snv_pos)),snv_pos)).transpose()
        paired_counts = np.round(np.multiply(read_mult_mat,paired_counts))
        BAF_counts = get_BAF_counts(paired_counts)
        
        seg_i.paired_counts = paired_counts
        seg_i.BAF_counts = BAF_counts
        
        data.segments.append(seg_i)
        data.seg_num += 1
            
    f.close()
    
    data.get_LOH_frac()
    
    ## Dumping output into pickle file
    outname = pjoin(pathbase,normpath(indir+out_str+str(fnum)+".MixClone.input.pkl"))
    outfile = open(outname, 'wb')
#    print outname
    pkl.dump(data, outfile, protocol=2)
    outfile.close()
    return outname

#%% Helper functions
def baf_gen(m):     # Generate BAFs for heterozygous sites
    baf_min = constants.BAF_N_MIN
    baf_max = constants.BAF_N_MAX
    
    return baf_min + (baf_max-baf_min)*random(m)

#%% Main
#indir = "ThetA_06_bootstrap/python/interval_files/scnv/"
#
#data=data_convert(indir=indir, fnum=0, n=2, k=3)

