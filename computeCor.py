import os
import gzip
import glob
from myconfig import *
import numpy as np

line = gzip.open("/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baseline/baseline.1.annot.gz").readline();
numberCol = len(line.split());

dataALL = np.loadtxt("/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baseline/baseline.1.annot.gz", usecols=range(5,numberCol), skiprows=1);

for chrom in range(2,23):
	data = np.loadtxt("/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baseline/baseline." + str(chrom) + ".annot.gz", usecols=range(5,numberCol), skiprows=1);
	dataALL = np.concatenate((dataALL, data), axis = 0);
	print chrom;		

currentPath = os.getcwd();

print annotationsName;
print annotationsFolder;

for nameFile in glob.glob(annotationsFolder + annotationsName + "/*"):
	nameAnnot = nameFile.replace(annotationsFolder + annotationsName, "");
	print nameAnnot;	
	annotationCombine = nameFile  + "/" + nameAnnot + ".";
	
	dataSmallAnnot = np.atleast_2d(np.loadtxt(annotationCombine + "1.annot", usecols=(4,), skiprows=1)).T;
	for chrom in range(2,23):
		data = np.loadtxt(annotationCombine + str(chrom) + ".annot", usecols=(4,), skiprows=1);
		dataSmallAnnot = np.concatenate((dataSmallAnnot, np.atleast_2d(data).T), axis = 0);
		print annotationCombine + str(chrom) + ".annot";
	dataALL = np.concatenate((dataALL, dataSmallAnnot), axis=1);
	
np.savetxt("BaseLine_GTEx_cor.R",np.corrcoef(dataALL.transpose()));

