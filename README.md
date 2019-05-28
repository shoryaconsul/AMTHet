# AMTHet
AMTHet - Inferring Intra-tumor hetoerogeneity

**Subdirectories:**
* mixclone: Contains files to run MixClone on our data
	1. datagen: Code to generate synthetic data from corresponding input file for THetA
	2. run_mixclone.py: Code to run MixClone pipeline from commnand line
* ThetA_06_bootstrap: Part of ThetA 0.6 required for our experiments
	1. interval_files: Contain input files for THetA. Generated input files for THetA and MixClone are stored here.
	2. output: Output files from MixClone and THetA. 
	3. PD4120a* files are used to run experiment on breast cancer sequencing data.

**Files:**
* alt_min_fn_upd.m: This function implements AMTHet.
* L_perm.m: Helper function for experiments
* synthdatagen_scnvsim.m: MATLAB script to generate input data for THetA with a process similar to SCNVSim.
* runMixClone_scnv.cmd: Windows batch script to run MixClone
* runTHeTA_MixClone_scnv.cmd: Windows batch script to run MixClone and THetA
* runTHeTA_n3k3_scnv.cmd: Windows batch script to run THetA
* runTHeTA_cancerseqdata.cmd: Windows batch script to run THetA on PD4120a 
* synthreadres.m: MATLAB script to read results from MixClone and THetA and run AMTHet on corresponding input files
* readres_PD4120a.m: MATLAB script to read results from THetA and run AMTHet on PD4120a 