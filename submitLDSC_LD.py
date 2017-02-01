import os
import glob
import subprocess
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

for bedFiles in glob.glob(bedFolder + "/*.bed"):
	bedFile = bedFiles.replace(bedFolder,"").replace(".bed","");
	print bedFile;
	script = "";
	for sumFiles in glob.glob(summaryFolder+"/*.sumstats"):
		sumFile = sumFiles.replace(summaryFolder,"").replace(".sumstats","");
		print sumFile;
	
		if os.path.exists(outFolder +  "/" + annotationsName + "/" + bedFile + "/" + sumFile + ".results"):
			continue;
		else:
			print outFolder +  "/" + annotationsName + "/" + bedFile + "/" + sumFile + ".results";
			
		scriptfile = currentPath + "/qsub/submitLDSC_"+ bedFile + "_"+ sumFile;
        	logfile    = currentPath + "/log/submitLDSC"+ bedFile + "_" + sumFile + ".log";
        	errfile    = currentPath + "/err/submitLDSC"+ bedFile + "_" + sumFile +   ".err";
		if not os.path.exists(outFolder + "/" + annotationsName):
			os.makedirs(outFolder + "/" + annotationsName);
		if not os.path.exists(outFolder + "/" + annotationsName + "/" + bedFile):
                        os.makedirs(outFolder + "/" + annotationsName + "/" + bedFile);	
		script =  "python " + LDSCPath + "ldsc.py" +\
			" --h2 " + sumFiles +\
			" --ref-ld-chr /groups/price/ldsc/reference_files/1000G_EUR_Phase3/baselineLD/baselineLD.," + annotationsFolder + "/" + annotationsName + "/" + bedFile+"/"+bedFile+"."+\
			" --frqfile-chr /groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC." +\
			" --w-ld-chr /groups/price/ldsc/reference_files/1000G_EUR_Phase3/weights/weights.hm3_noMHC." +\
			" --overlap-annot --print-cov --print-coefficients --print-delete-vals " +\
			" --out " + outFolder +  "/" + annotationsName + "/" + bedFile + "/" + sumFile + "\n";
			 
		print script;
		scriptFILEHandler = open(scriptfile+'.qsub', 'wb');
		scriptFILEHandler.write(TEMPLATE_SERIAL.format(script=script, name="LDSC", logfile=logfile, errfile=errfile, slots=1))
		scriptFILEHandler.close();
		subprocess.call('bsub < ' + scriptfile + '.qsub', shell=True);
		
