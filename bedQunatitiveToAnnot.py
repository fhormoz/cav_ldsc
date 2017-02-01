import os
import glob
import pandas as pd
import numpy as npy
import subprocess
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
#BSUB -q short            #submit to "short" queue
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

	if not os.path.exists(annotationsFolder + "/" + annotationsName):
                os.makedirs(annotationsFolder + "/" + annotationsName);
        if not os.path.exists(annotationsFolder + "/" + annotationsName+ "/" + bedFile):
                os.makedirs(annotationsFolder + "/" + annotationsName+ "/" + bedFile);

	print bedFiles;
	print bedFile;
	print annotationsFolder + "/" + annotationsName+ "/" + bedFile + "/" + bedFile;
	
	for chrom in range(1,23,1):
		if os.path.exists(annotationsFolder + "/" + annotationsName+ "/" + bedFile + "/" + bedFile + '.' + str(chrom) + '.annot'):
			continue;
		bimDataM = npy.genfromtxt('/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC.'+str(chrom)+'.bim', dtype='str', skip_header=0)
		bimDataFram = pd.DataFrame(data=bimDataM, columns=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2']);

		bedData  = npy.genfromtxt(bedFiles, dtype='str',skip_header=0);
		bedDataFram = pd.DataFrame(data=bedData, columns=['CHR', 'BP', 'End', 'Annot']);

		merge = pd.merge(bimDataFram, bedDataFram, how='left', on=['CHR', 'BP']);
		merge = merge.fillna(0);

		merge.to_csv(annotationsFolder + "/" + annotationsName+ "/" + bedFile + "/" + bedFile + '.' + str(chrom) + '.annot', sep="\t", columns=['CHR', 'BP', 'SNP', 'CM', 'Annot'], index = False);
		
		scriptfile = currentPath + "/qsub/"+ bedFile + "_" + str(chrom);
		logfile    = currentPath + "/log/"+ bedFile + "_" + str(chrom) + ".log";
		errfile    = currentPath + "/err/"+ bedFile + "_" + str(chrom) + ".err";
		script = "python " + LDSCPath+ "ldsc.py" +\
			" --l2 --bfile /groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC." + str(chrom) +\
			" --ld-wind-cm 1 --annot " + annotationsFolder + "/" + annotationsName+ "/" + bedFile + "/" + bedFile + "." + str(chrom) +  ".annot" +\
			" --out " + annotationsFolder + "/" + annotationsName+ "/" + bedFile + "/" + bedFile + "." + str(chrom) +\
			" --print-snps " + currentPath + "/list.txt";

		scriptFILEHandler = open(scriptfile+'.qsub', 'wb');
		scriptFILEHandler.write(TEMPLATE_SERIAL.format(script=script, name="eCAVIAR", logfile=logfile, errfile=errfile, slots=1))
		scriptFILEHandler.close();
		subprocess.call('bsub < ' + scriptfile + '.qsub', shell=True)
		print chrom;
