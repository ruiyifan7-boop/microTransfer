## run_bootstrap_CI_management.R  (reviewer 2-item: bootstrap CI for Fig 4 BA/macro-F1)
## Class-stratified nonparametric bootstrap (2000 resamples) of leave-one-block-out held-out
## predictions; 95% percentile CI. set.seed(20260706).
set.seed(20260706)
root <- "global_harmonized"; outdir <- file.path(root,"enhancement","review"); dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
p <- read.delim(file.path(root,"enhancement","transferability","transferability_predictions.tsv"), check.names=FALSE)
p <- p[p$direction=="leave_one_block_out", c("kind","variant","sample_uid","truth","predicted","fold")]
bal <- function(truth,pred){lv<-levels(truth); pred<-factor(pred,levels=lv); mean(sapply(lv,function(z) mean(pred[truth==z]==z)), na.rm=TRUE)}
mf1 <- function(truth,pred){lv<-levels(truth); pred<-factor(pred,levels=lv); mean(sapply(lv,function(z){tp<-sum(truth==z&pred==z);fp<-sum(truth!=z&pred==z);fn<-sum(truth==z&pred!=z);pr<-ifelse(tp+fp==0,0,tp/(tp+fp));rc<-ifelse(tp+fn==0,0,tp/(tp+fn));ifelse(pr+rc==0,0,2*pr*rc/(pr+rc))}))}
boot <- function(df,Bn=2000){df$truth<-factor(df$truth);ob<-bal(df$truth,df$predicted);of<-mf1(df$truth,df$predicted)
  by<-split(seq_len(nrow(df)),df$truth);bb<-numeric(Bn);bf<-numeric(Bn)
  for(b in 1:Bn){ix<-unlist(lapply(by,function(ii) if(length(ii)) sample(ii,length(ii),replace=TRUE) else ii))
    tt<-factor(df$truth[ix],levels=levels(df$truth));pp<-df$predicted[ix];bb[b]<-bal(tt,pp);bf[b]<-mf1(tt,pp)}
  data.frame(n=nrow(df),balanced_accuracy=round(ob,3),BA_lo=round(quantile(bb,.025),3),BA_hi=round(quantile(bb,.975),3),
    macro_F1=round(of,3),F1_lo=round(quantile(bf,.025),3),F1_hi=round(quantile(bf,.975),3))}
res <- do.call(rbind, lapply(split(p, paste(p$kind,p$variant,sep="|")), boot))
res <- cbind(model=rownames(res),res); rownames(res)<-NULL; print(res,row.names=FALSE)
write.table(res, file.path(outdir,"S_bootstrap_CI_management.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
