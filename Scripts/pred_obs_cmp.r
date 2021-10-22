
a<-commandArgs(T)
data <- read.table(a[1])
dd <- unique(data)
outname <- a[2]

if(ncol(dd)==3){
    colnames(dd) <- c("loci","obs","fxr")
    pred_enc <- dd[,"fxr"]
    obs <- dd[,"obs"]
    pred_enc_floor <- floor(pred_enc)
    
    printdata_enc <- c()
    for( value_pred in seq(0,30)){
        this_value <- obs[which(pred_enc_floor == value_pred)]
        this_quantile <- c(quantile(this_value,probs=c(0.1,0.25,0.5,0.75,0.9)),mean(this_value),sd(this_value))
        printdata_enc <- rbind(printdata_enc,this_quantile)
    }
    rownames(printdata_enc) <- paste0("expcut",seq(0,30))
    colnames(printdata_enc) <- c("q10","q25","q50","q75","q90","mean","sd")
    
#    write.table(printdata_enc, file=paste0(outname,"_fxrQuantile.txt"),sep="\t",row.names=T,col.names=T,quote=F)
    outdata <-c(outname,
                        cor(log(obs+0.1),log(pred_enc+0.1),use="complete.obs"))
    write.table(t(outdata),file=paste0(outname,"_cor.txt"),sep="\t",quote=F,row.names=F,col.names=F)

}else{
    colnames(dd) <- c("loci","obs","raw","enc")
    pred_raw <- dd[,"raw"]
    pred_enc <- dd[,"enc"]
    obs <- dd[,"obs"]
    pred_raw_floor <- floor(pred_raw)
    pred_enc_floor <- floor(pred_enc)

    printdata_raw <- c()
    for( value_pred in seq(0,30)){
        this_value <- obs[which(pred_raw_floor == value_pred)]
        this_quantile <- c(quantile(this_value,probs=c(0.1,0.25,0.5,0.75,0.9)),mean(this_value),sd(this_value))
        printdata_raw <- rbind(printdata_raw,this_quantile)
    }
    rownames(printdata_raw) <- paste0("expcut",seq(0,30))
    colnames(printdata_raw) <- c("q10","q25","q50","q75","q90","mean","sd")
    write.table(printdata_raw, file=paste0(outname,"_rawQuantile.txt"),sep="\t",row.names=T,col.names=T,quote=F)
    
    printdata_enc <- c()
    for( value_pred in seq(0,30)){
        this_value <- obs[which(pred_enc_floor == value_pred)]
        this_quantile <- c(quantile(this_value,probs=c(0.1,0.25,0.5,0.75,0.9)),mean(this_value),sd(this_value))
        printdata_enc <- rbind(printdata_enc,this_quantile)
    }
    rownames(printdata_enc) <- paste0("expcut",seq(0,30))
    colnames(printdata_enc) <- c("q10","q25","q50","q75","q90","mean","sd")
    
#    write.table(printdata_enc, file=paste0(outname,"_encQuantile.txt"),sep="\t",row.names=T,col.names=T,quote=F)
    outdata <-c(outname,
                        cor(log(obs+0.1),log(pred_raw+0.1),use="complete.obs"),cor(log(obs+0.1),log(pred_enc+0.1),use="complete.obs"))
    write.table(t(outdata),file=paste0(outname,"_cor.txt"),sep="\t",quote=F,row.names=F,col.names=F)

}

