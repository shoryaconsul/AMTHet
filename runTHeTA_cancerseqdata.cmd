@ECHO off

Rem This is a batch file to run THetA on the input files. The input directory is interval_Files and output directory is output.
Rem The tumors in these inputs have either only amplifications or only deletions in each interval.

set dirc=%~dp0ThetA_06_bootstrap\python
set dirin=%~dp0ThetA_06_bootstrap\python\interval_files
set prein=%~dp0ThetA_06_bootstrap\python\interval_files\PD4120a.n3
set dirout=%~dp0ThetA_06_bootstrap\python\output\
for %%f in (%prein%*) DO (
	python %dirc%\RunTHetA.py %dirout%\PD4120a.n2.withBounds -n 3 -k 3 -t 2 -d %dirout% --RESULTS %dirout%\PD4120a.n2.results
)
pause
