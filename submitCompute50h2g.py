import os
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

for bedFiles in glob.glob(bedFolder + "/M*.bed"):
	bedFile = bedFiles.replace(bedFolder,"").replace(".bed","");
	print bedFile;
	for sumFiles in glob.glob(summaryFolder+"/*.sumstats"):
		sumFile = sumFiles.replace(summaryFolder,"").replace(".sumstats","");
		resultFile = outFolder +  "/" + annotationsName + "/" + bedFile + "/" + sumFile + ".results";
		tau = np.loadtxt(open(resultFile,"r"), skiprows=1, usecols=(7,));
		annotationsLD  = "/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baselineLD/baselineLD.";
		annotationsNew = annotationsFolder + "/" + annotationsName + "/" + bedFile+"/"+bedFile+".";	
	

		scriptfile = currentPath + "/qsub/submitComputeh2g_"+ bedFile + "_"+ sumFile;
                logfile    = currentPath + "/log/submitComputeh2g"+ bedFile + "_" + sumFile + ".log";
                errfile    = currentPath + "/err/submitComputeh2g"+ bedFile + "_" + sumFile +   ".err";

                script =  "python " + currentPath + "compute50h2g.py" +\
                        " -n " + annotationsLD +\
                        " -l " + annotationsNew+\
                        " -r " + resultFile +"\n"

                print script;
                scriptFILEHandler = open(scriptfile+'.qsub', 'wb');
                scriptFILEHandler.write(TEMPLATE_SERIAL.format(script=script, name="Compute50h2g", logfile=logfile, errfile=errfile, slots=1))
                scriptFILEHandler.close();
                subprocess.call('bsub < ' + scriptfile + '.qsub', shell=True);
	exit(0);
