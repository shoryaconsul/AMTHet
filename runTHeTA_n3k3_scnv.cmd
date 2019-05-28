@ECHO off

Rem This is a batch file to run THetA on the input files. The input directory is interval_Files and output directory is output.
Rem The tumors in these inputs have either only amplifications or only deletions in each interval.

set dirc=%~dp0ThetA_06_bootstrap\python
set dirin=%~dp0ThetA_06_bootstrap\python\interval_files\scnv
set prein=%~dp0ThetA_06_bootstrap\python\interval_files\scnv\interval_count_n2
set dirout=%~dp0ThetA_06_bootstrap\python\output\scnv
for %%f in (%prein%*) DO (
	python %dirc%\RunTHetA.py %%f -n 2 -k 3 -t 2 -d %dirout% --NO_INTERVAL_SELECTION
)
pause