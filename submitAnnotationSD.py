import os
import glob
import subprocess
import numpy as npy
import pandas as pd
from myconfig import *
from subprocess import Popen, PIPE


def maxSubmitReached(max):
        p1 = Popen(["qstat", "-u", "fhormoz"], stdout=PIPE)
        p2 = Popen(["wc", "-l"], stdin=p1.stdout, stdout=PIPE)
        output = p2.communicate()[0]
        if(int(output) < max) :
                return False;
        else:
                return True;


TEMPLATE_SERIAL = """
#####################################
#!/bin/bash
#BSUB -n 1                #each  job run on 1 core
#BSUB -W 06:00            #job run 2 hour
#BSUB -J {name}
#BSUB -o {logfile}        #lsf output file
#BSUB -e {errfile}        #lsf error file
#BSUB -q short         	  #submit to "short" queue
#####################################
echo "------------------------------------------------------------------------"
echo "Job started on" `date`
echo "------------------------------------------------------------------------"
{script}
echo "------------------------------------------------------------------------"
echo "Job ended on" `date`
echo "------------------------------------------------------------------------"
"""

currentPath = os.getcwd();

for bedFiles in glob.glob(bedFolder + "/*.bed"):
	bedFile = bedFiles.replace(bedFolder,"").replace(".bed","");
	print bedFiles;
	annotationSdOUTName = annotationsFolder + "/" + annotationsName+ "/" + bedFile + "/" + bedFile + ".sd.out"; 
	if os.path.exists(annotationSdOUTName):
		continue; 
	print annotationSdOUTName;
	dataAnnotGW = npy.array([], dtype=npy.float64);
	for chrom in range(1,23,1):
		annotationsFile = annotationsFolder + "/" + annotationsName+ "/" + bedFile + "/" + bedFile + "." + str(chrom) + ".annot";
		freq1KGFile = "/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC." + str(chrom) + ".frq"; 

		annnotationsData = npy.genfromtxt(annotationsFile, dtype=None,skip_header=0, usecols=(2,4) )
		annnotationsDataFram = pd.DataFrame(data=annnotationsData[1:,:], columns=['SNP', 'ANNOT']);

		freq1KGData = npy.genfromtxt(freq1KGFile, dtype=None, skip_header=0, usecols=(1,4));
		freq1KGDataFrame = pd.DataFrame(data=freq1KGData[1:,:], columns=['SNP', 'MAF']);
	
		freq1KGannnotationsFrame = pd.merge(annnotationsDataFram, freq1KGDataFrame, on='SNP');

		freq1KGannnotationsFrame[['ANNOT','MAF']] = freq1KGannnotationsFrame[['ANNOT','MAF']].apply(pd.to_numeric)
		commonannnotations = freq1KGannnotationsFrame[freq1KGannnotationsFrame['MAF'] > 0.05];
		commonannnotationsVal = commonannnotations['ANNOT'];
		dataAnnotGW = npy.append(dataAnnotGW,commonannnotationsVal.as_matrix()); 
		print chrom;
	outfile = open(annotationSdOUTName, 'w');
	outfile.write(str(npy.std(dataAnnotGW))+"\n");
	outfile.close();
