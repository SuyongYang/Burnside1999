The files for MATLAB are contained in the ZIP file distrib.zip.  You should 
unzip these to a directory of your choice, say D:\burnside, and maintain the
subdirectory structure given in the ZIP file.  That is, you should have 
subdirectories called D:\burnside\ber1993, D:\burnside\ce1992, 
D:\burnside\divisib, and D:\burnside\indivis.  

You should then add D:\burnside to your PATH in MATLAB.  This is accomplished
by opening MATLAB and using the File-Set Path-Path-Add to Path sequence to 
add the directory.  This command is permanent so you only have to go through 
it once.  It is a necessary step because the programs in the individual 
subdirectories all access the programs in the D:\burnside directory.  This 
step lets you avoid running the path command indicated in the notes every time
you run MATLAB.

To run one of the programs, say indivis.m, change directory in MATLAB, say by
executing 

cd d:\burnside\indivis

at the command prompt.  Then just execute

indivis

The output will appear on your screen.  

The other analagous programs are divisib.m, ber1993.m and ce1992.m.  They
control the other programs in their respective subdirectories.  See the 
notes for more details on how the programs work.  

The programs are used at your own risk.  I don't promise that they are 
error free, nor are they the responsibility of the World Bank.  If you find
them useful in your work, I only ask that you please cite the notes.  

Craig Burnside
January 2000