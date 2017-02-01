import os
import optparse
import glob
import gzip
import subprocess
import numpy as np
from myconfigLD import *
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
#BSUB -W 03:00            #job run 5 hour
#BSUB -R "rusage[mem=10000]" 
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


def main(parser):
        (options, args) = parser.parse_args();
        annotLDPath    = options.annotLDPath;
        annotNewPath   = options.annotNewPath;
	LDSCResultFile = options.LDSCResultFile; 

	hg2 = 0;
	tau = np.loadtxt(open(LDSCResultFile,"r"), skiprows=1, usecols=(7,));
	for chrom in range(1,23):
		annotLDFileHandler  = gzip.open(annotLDPath + str(chrom) + ".annot.gz", "rb");
		annotNewFileHandler = open(annotNewPath + str(chrom) + ".annot");
		# Remove the header file
		annotLDFileHandler.readline();
		annotNewFileHandler.readline(); 
		print chrom;
		for dataBaseLineLDData, newAnnotData in zip(annotLDFileHandler, annotNewFileHandler):
			combineData = list(dataBaseLineLDData.split()[4:] + newAnnotData.split()[4:]);
			combineData = map(float, combineData);
			#combineDataSquare = [x**2 for x in combineData];
			hg2 = hg2 +  sum(combineData * tau);
			
	print sumFiles, hg2;

	outfile = open(LDSCResultFile.replace(".results", ".hg2"), 'w');
        outfile.write(str(hg2)+"\n");
        outfile.close();		


if __name__ == "__main__":
        parser = optparse.OptionParser("usage: %prog [options] ")
        parser.add_option("-n", "--NewAnnot", dest="annotNewPath",
                default="", type="string",
                help="specify the path of new annotation");
        parser.add_option("-l", "--LDBaseLinePath", dest="annotLDPath",
                default="", type="string",
                help="specify the Path to the BaseLine LD");
        parser.add_option("-r", "--LDSCResults", dest="LDSCResultFile",
                default="", type="string",
                help="specify the full path to LDSC result for each trait")

        main(parser);
