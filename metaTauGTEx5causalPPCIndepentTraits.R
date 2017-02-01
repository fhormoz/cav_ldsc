library("rmeta")

Mref       = 5961159

BedPath = "/groups/price/farhad/Data/Annotations/eQTL_CAVIAR_GTEx_v6_5causal/";
BedFiles = Sys.glob("//groups/price/farhad/Data/Annotations/eQTL_CAVIAR_GTEx_v6_5causal/*");
BedFiles = gsub("/groups/price/farhad/Data/Annotations/eQTL_CAVIAR_GTEx_v6_5causal/", "", BedFiles)
BedFiles = gsub("/","",BedFiles);

indep_traits = c('UKBiobank_Height3', 'UKBiobank_FVC', 'UKBiobank_Diastolic', 'UKBiobank_FEV1FV', 'PASS_Schizophrenia',
		 'UKBiobank_SmokingStatus', 'UKBiobank_Age_at_Menarche3', 'PASS_BMI1', 'PASS_Years_of_Education1',
		 'UKBiobank_Eczema', 'PASS_Anorexia', 'PASS_Crohns_Disease', 'PASS_LDL', 'PASS_Rheumatoid_Arthritis',
		 'UKBiobank_Age_at_Menopause2', 'PASS_Neuroticism', 'PASS_Coronary_Artery_Disease',
		 'PASS_Autism', 'PASS_Lupus', 'PASS_ENIGMA2_MeanPutamen');

cat(indep_traits);
cat("\n");
for (tissue in BedFiles) {
	sd_annot   = read.table(paste0(paste(paste(BedPath, tissue,sep="/"),tissue,sep="/"), ".sd.out"))$V1;
	#sd_annot   = read.table("/groups/price/farhad/Data/Annotations/eQTL_CAVIAR/LumpALL/LumpALL.sd.out")$V1	#you should put the sd of your annotations computed on the Mref SNPs with MAF>=5%
	list_trait = Sys.glob(paste0(paste("/groups/price/farhad/Results/LDSC_eQTL_GTEx_CAVIAR_5casaul_v6_7Nov2016/eQTL_CAVIAR_GTEx_v6_5causal/", tissue,sep="/"), "/*.log"));
	list_trait = gsub(".log", "", list_trait);
	meta_tau_star = NULL # col 1: tau*; col2: sd; col3: pval
	meta_tau = NULL #Meta-analysis of tau
	meta_eni = NULL
	for (trait in list_trait) {
	    traitName = gsub("/", "", gsub(paste("/groups/price/farhad/Results/LDSC_eQTL_GTEx_CAVIAR_5casaul_v6_7Nov2016/eQTL_CAVIAR_GTEx_v6_5causal/", tissue, sep="/"), "", trait));
	    if(traitName %in% indep_traits) {
	    log = read.table(paste(trait,".log",sep=""),h=F,fill=T)
	    h2g = as.numeric(as.character(log[which(log$V4=="h2:"),5]))
	    res = read.table(paste(trait,".results",sep=""),h=T) 
	    
	    meta_tau_star = rbind(meta_tau_star, c(
	      res[nrow(res),8]*Mref*sd_annot/h2g,
	      res[nrow(res),9]*Mref*sd_annot/h2g,
	      2*pnorm(-abs(res[nrow(res),10]))))
	   
	   meta_tau = rbind(meta_tau, c(
              res[nrow(res),8]*Mref/h2g,
              res[nrow(res),9]*Mref/h2g,
              2*pnorm(-abs(res[nrow(res),10]))))
	
	   meta_eni = rbind(meta_eni, c(
               res[nrow(res),5]*sd_annot/h2g,
               res[nrow(res),6]*sd_annot/h2g,
               2*pnorm(-abs(res[nrow(res),7]))))	
	 } else {
		#cat(traitName);
		#cat("\n");	
	 }
	}
	test_tau_star=meta.summaries(meta_tau_star[,1],meta_tau_star[,2],method="fix")
	cat(tissue);
	cat("\t");
	cat(c(test_tau_star$summary,test_tau_star$se.summary,2*pnorm(-abs(test_tau_star$summary/test_tau_star$se.summary))))
	cat("\t");
	
	test_tau=meta.summaries(meta_tau[,1],meta_tau[,2],method="fix")
	cat(c(test_tau$summary,test_tau$se.summary,2*pnorm(-abs(test_tau$summary/test_tau$se.summary))))
	cat("\t");

	test_eni=meta.summaries(meta_eni[,1],meta_eni[,2],method="fix")
        cat(c(test_eni$summary,test_eni$se.summary,2*pnorm(-abs(test_eni$summary/test_eni$se.summary))))
        cat("\t");

	cat(sum(meta_tau_star[,1]/meta_tau_star[,2]));
	cat("\n");
}
