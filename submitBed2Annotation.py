import os
import glob
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

	if not os.path.exists(annotationsFolder + "/" + annotationsName):
                os.makedirs(annotationsFolder + "/" + annotationsName);
	if not os.path.exists(annotationsFolder + "/" + annotationsName+ "/" + bedFile):
        	os.makedirs(annotationsFolder + "/" + annotationsName+ "/" + bedFile);

	for chrom in range(1,23,1):
		#l2.ldscore.gz
		if not os.path.exists(annotationsFolder + "/" + annotationsName+ "/" + bedFile + "/" + bedFile + "." + str(chrom) + ".l2.ldscore.gz"):
			print annotationsFolder + "/" + annotationsName+ "/" + bedFile + "/" + bedFile + "." + str(chrom) + ".l2.ldscore.gz";
		
			scriptfile = currentPath + "/qsub/"+ bedFile + "_" + str(chrom);
			logfile    = currentPath + "/log/"+ bedFile + "_" + str(chrom) + ".log";
			errfile	   = currentPath + "/err/"+ bedFile + "_" + str(chrom) + ".err";
			script = "python " + currentPath + "/bedToAnnot.py --bedfile-single " + bedFiles +\
				 " --bfile /groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC." + str(chrom) +\
				 " --annotfile " + annotationsFolder + "/" + annotationsName+ "/" + bedFile + "/tmp." + str(chrom) + ".annot\n" ;
	
			script = script + " cat " + annotationsFolder + "/" + annotationsName+ "/" + bedFile + "/tmp." + str(chrom)+".annot | uniq >" +\
					  annotationsFolder + "/" + annotationsName+ "/" + bedFile + "/" + bedFile + "." + str(chrom) +  ".annot\n" ;
	
			script = script + " rm " + annotationsFolder + "/" + annotationsName+ "/" + bedFile + "/tmp." + str(chrom)+".annot\n";
			script = script + "python " + LDSCPath+ "ldsc.py" +\
				" --l2 --bfile /groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC." + str(chrom) +\
				" --ld-wind-cm 1 --annot " + annotationsFolder + "/" + annotationsName+ "/" + bedFile + "/" + bedFile + "." + str(chrom) +  ".annot" +\
				" --out " + annotationsFolder + "/" + annotationsName+ "/" + bedFile + "/" + bedFile + "." + str(chrom) +\
			        " --print-snps " + currentPath + "/list.txt";
  
			scriptFILEHandler = open(scriptfile+'.qsub', 'wb');
       			scriptFILEHandler.write(TEMPLATE_SERIAL.format(script=script, name="eCAVIAR", logfile=logfile, errfile=errfile, slots=1))
			scriptFILEHandler.close();
       		 	subprocess.call('bsub < ' + scriptfile + '.qsub', shell=True)
		
