library("rmeta")

Mref       = 5961159

BedPath = "/groups/price/farhad/Data/Annotations/eQTL_CAVIAR/";
BedFiles = Sys.glob("//groups/price/farhad/Data/Annotations/eQTL_CAVIAR/*");
BedFiles = gsub("/groups/price/farhad/Data/Annotations/eQTL_CAVIAR/", "", BedFiles)
BedFiles = gsub("/","",BedFiles);

tissue = "WrightBlood.binary_c_5"
baseLineAnnotCount = length(readLines(paste0(paste("/groups/price/farhad/Results/LDSC_eQTL_GTEx_CAVIAR/eQTL_CAVIAR/", tissue,sep="/"), "/PASS_Multiple_sclerosis.results")))-1;


for (annotIndex in 2:(baseLineAnnotCount-1)) {
	sd_annot   = as.matrix(read.table('/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baseline/baseline.sd.out'));
	list_trait = Sys.glob(paste0(paste("/groups/price/farhad/Results/LDSC_baseline/eQTL_CAVIAR/", tissue,sep="/"), "/*.log"));
	list_trait = gsub(".log", "", list_trait);
	meta_tau_star = NULL # col 1: tau*; col2: sd; col3: pval
	meta_enr_star = NULL #Meta-analysis of enrichment
	for (trait in list_trait) {
	    log = read.table(paste(trait,".log",sep=""),h=F,fill=T)
	    h2g = as.numeric(as.character(log[which(log$V4=="h2:"),5]))
	    res = read.table(paste(trait,".results",sep=""),h=T) 
	    
	    meta_tau_star = rbind(meta_tau_star, c(
	      res[annotIndex,8]*Mref*sd_annot[annotIndex-1]/h2g,
	      res[annotIndex,9]*Mref*sd_annot[annotIndex-1]/h2g,
	      2*pnorm(-abs(res[annotIndex,10]))))
	  
	   meta_enr_star = rbind(meta_enr_star, c(
              res[annotIndex,5],
              res[annotIndex,6],
              2*pnorm(-abs(res[annotIndex,7]))))	
	}
	test_tau=meta.summaries(abs(meta_tau_star[,1]),meta_tau_star[,2],method="random")
	cat(tissue);
	cat("\t");
	cat(c(test_tau$summary,test_tau$se.summary,2*pnorm(-abs(test_tau$summary/test_tau$se.summary))))
	cat("\t");
	test_enr=meta.summaries(meta_enr_star[,1],meta_enr_star[,2],method="random")
	cat(c(test_enr$summary,test_enr$se.summary,2*pnorm(-abs(test_enr$summary/test_enr$se.summary))))
	cat("\n");
}

tissue = "Meta_GTEx_Output_maxppc"
for (annotIndex in 2:(baseLineAnnotCount-1)) {
        sd_annot   = as.matrix(read.table('/groups/price/ldsc/reference_files/1000G_EUR_Phase3/baseline/baseline.sd.out'));
        list_trait = Sys.glob(paste0(paste("/groups/price/farhad/Results/LDSC_eQTL_GTEx_CAVIAR_5casaul_v6_7Nov2016/eQTL_CAVIAR_GTEx_v6_5causal/", tissue,sep="/"), "/*.log"));
        list_trait = gsub(".log", "", list_trait);
        meta_tau_star = NULL # col 1: tau*; col2: sd; col3: pval
        meta_enr_star = NULL #Meta-analysis of enrichment
        for (trait in list_trait) {
	    log = read.table(paste(trait,".log",sep=""),h=F,fill=T)
            h2g = as.numeric(as.character(log[which(log$V4=="h2:"),5]))
            res = read.table(paste(trait,".results",sep=""),h=T)

            meta_tau_star = rbind(meta_tau_star, c(
              res[annotIndex,8]*Mref*sd_annot[annotIndex-1]/h2g,
              res[annotIndex,9]*Mref*sd_annot[annotIndex-1]/h2g,
              2*pnorm(-abs(res[annotIndex,10]))))
           
           meta_enr_star = rbind(meta_enr_star, c(
              res[annotIndex,5],
              res[annotIndex,6],
              2*pnorm(-abs(res[annotIndex,7]))))
        }
        test_tau=meta.summaries(abs(meta_tau_star[,1]),meta_tau_star[,2],method="random")
        cat(tissue);
        cat("\t");
        cat(c(test_tau$summary,test_tau$se.summary,2*pnorm(-abs(test_tau$summary/test_tau$se.summary))))
        cat("\t");
        test_enr=meta.summaries(meta_enr_star[,1],meta_enr_star[,2],method="random")
        cat(c(test_enr$summary,test_enr$se.summary,2*pnorm(-abs(test_enr$summary/test_enr$se.summary))))
        cat("\n");
}
