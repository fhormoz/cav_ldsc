library("rmeta")

Mref       = 5961159

BedPath = "/groups/price/farhad/Data/Annotations/eQTL_CAVIAR/";
BedFiles = Sys.glob("//groups/price/farhad/Data/Annotations/eQTL_CAVIAR/*");
BedFiles = gsub("/groups/price/farhad/Data/Annotations/eQTL_CAVIAR/", "", BedFiles)
BedFiles = gsub("/","",BedFiles);

tissue = "LumpALL"
#for (tissue in BedFiles) {
	sd_annot   = read.table(paste0(paste(paste(BedPath, tissue,sep="/"),tissue,sep="/"), ".sd.out"))$V1;
	cat(sd_annot);
	cat("\n");
	#sd_annot   = read.table("/groups/price/farhad/Data/Annotations/eQTL_CAVIAR/LumpALL/LumpALL.sd.out")$V1	#you should put the sd of your annotations computed on the Mref SNPs with MAF>=5%
	list_trait = Sys.glob(paste0(paste("/groups/price/farhad/Results/LDSC_eQTL_GTEx_CAVIAR/eQTL_CAVIAR/", tissue,sep="/"), "/*.log"));
	list_trait = gsub(".log", "", list_trait);
	meta_tau_star = NULL # col 1: tau*; col2: sd; col3: pval
	for (trait in list_trait) {
	    if( (file.info(paste(trait,".log",sep=""))$size!=0)) {
	    log = read.table(paste(trait,".log",sep=""),h=F,fill=T)
	    h2g = as.numeric(as.character(log[which(log$V4=="h2:"),5]))
	    res = read.table(paste(trait,".results",sep=""),h=T) 
	    cat(res[nrow(res),8]*Mref*sd_annot/h2g);
	    cat("\t");
	    cat(res[nrow(res),9]*Mref*sd_annot/h2g);
	    cat("\n");   
	    meta_tau_star = rbind(meta_tau_star, c(
	      res[nrow(res),8]*Mref*sd_annot/h2g,
	      res[nrow(res),9]*Mref*sd_annot/h2g,
	      2*pnorm(-abs(res[nrow(res),10]))))
	   }
	}
        meta_tau_star
	test=meta.summaries(meta_tau_star[,1],meta_tau_star[,2],method="random")
	cat(tissue);
	cat("\t");
	cat(c(test$summary,test$se.summary,2*pnorm(-abs(test$summary/test$se.summary))))
	cat("\n");
#}
