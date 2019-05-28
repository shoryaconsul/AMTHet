@ECHO off
setlocal enabledelayedexpansion

Rem This is a batch file to run THetA on the input files. The input directory is interval_Files and output directory is output.
Rem The tumors in these inputs have either only amplifications or only deletions in each interval.

set dirc=%~dp0ThetA_06_bootstrap\python
set dirm=%~dp0mixclone
set dirin=%~dp0ThetA_06_bootstrap\python\interval_files\scnv
set prein=%~dp0ThetA_06_bootstrap\python\interval_files\scnv\interval_count_n4
set dirout=%~dp0ThetA_06_bootstrap\python\output\scnv
for %%f in (%prein%*) DO (
	set ftmp=%%f
	set fn=!ftmp:~-1!
	python %dirm%\run_mixclone.py %dirin% 4 4 !fn!
)
pause