#' Predict a candidate gene for each region conditioning on the candidates selected in rest regions.
#'
#' @param
#' pickup Numeric. Region index.
#' region_chosen A vector of TRUE/FALSE indicates which regions are selected as conditioning regions
#' remaining A vector of genes chosen in the last iteration ((n-1)th round).
#' opt A list of parameter in Gibbs
#' @return Characters of a single gene 
#' @import
#' e1071
#' data.table
#' RLT
#' @export
#' @examples
#' predict_gene(pickup,region_chosen,remaining,opt,count))
#' 
rf2weight<-function(pickup,region_chosen,remaining,opt,count){
    cat("rf starts.. ",format(Sys.time()),"\n")
    ##cat(names(opt$evi$evi_list),"\n")
    rf_models <- lapply(names(opt$evi$evi_list),function(x){rf_training(pickup,region_chosen, remaining,opt,count,x)})
    rf_p <- Reduce(function(x,y){merge(x,y,by="official_name")},lapply(rf_models,function(x){x$rf_p}))
    pvaluecols <- names(rf_p)[grep("pvalue",names(rf_p),perl=T)]
    ##cat(pvaluecols,"\n")
    ##print(rf_p)
    rf_p[,risk:=prod(.SD),by="official_name",.SDcols=pvaluecols]
    rf_p[,extra_weight:=risk^opt$method_weight_power,]
    rf_p_o <- merge(opt$gene,rf_p,by="official_name")

    cat("rf ends... ",format(Sys.time()),"\n")    
    return(rf_p_o)
}


rf_training <- function(pickup,region_chosen,remaining,opt,count,evi_name){
    ## cat("processing",evi_name,format(Sys.time()),"\n")
    evi <-  opt$evi$evi_list[[evi_name]]
    exclude_samegene<-T
    index<-region_chosen
    index[pickup]<-F
    rows <- opt$nodes[opt$nodes %in% unlist(opt$region[pickup])]

    ## cat("format data",format(Sys.time()),"\n")
    positive<- data.table(official_name=remaining[!remaining %in% rows],label=1) %>% unique
    
    ## Filters: P<median, not in positive, not in current region
    negative <- data.table(official_name=sample(count[p<=median(p) & !(official_name %in% positive$official_name),][!official_name %in% rows,]$official_name,size=length(positive$official_name)),label="0") %>% unique

    tb <- data.table()
    tb <- rbind(positive,negative)
    tb$label <- as.factor(tb$label)
    tb <- tb[!duplicated(official_name),]

    training <- data.table()
    training <- merge(tb,evi,by="official_name")
    training <- training[!duplicated(official_name),]
    
    ## remove consistency columns
    ## cat("Column numbers of features: ",lapply(training,ncol) %>% unlist,"\n")
    ##training <- training[,lapply(training,function(y){ifelse(length(unique(y))==1,FALSE,TRUE)}) %>% unlist,with=F]
    ## ## cat("Column numbers of features after remove consistence:",lapply(training,ncol) %>% unlist,"\n")
    ## cat("training model",format(Sys.time()),"\n")
    rf_model <- RLT(training[,-1:-2],training$label,ntree=opt$ntrees,model="classification",track.obs=T, resample.prob = opt$resample_prob,use.cores = opt$threads)
    
    rownames(rf_model$ObsTrack) <- training$official_name
    rf_model$evi <- evi
    rf_model$feaname <- evi_name
    ## cat("predicting",format(Sys.time()),"\n")
    rf_model$rf_p <- rf_predict(rf_model,opt)
    ##rf_model$rf_p <- merge(rf_model$rf_p,tb,by="official_name",all.x=TRUE)

    ## cat("predicting done",format(Sys.time()),"\n")

    
    varimp_file <- paste(opt$res_path,"/",evi_name,".varimp",sep="")
    if(file.exists(varimp_file)){
        cat(rf_model$VarImp,file=varimp_file,sep="\t",append=T)
        cat("\n",file=varimp_file,sep="",append=T)
    }else{
        cat(colnames(rf_model$VarImp),file=varimp_file,sep="\t")
        cat("\n",file=varimp_file,sep="",append=T)
        cat(rf_model$VarImp,file=varimp_file,sep="\t",append=T)
        cat("\n",file=varimp_file,sep="",append=T)
    }
    
    return(rf_model)
}

rf_predict <- function(rf_model,opt){
    ## cat("predicting",format(Sys.time()),"\n")
    score <- predict.RLT(rf_model,rf_model$evi[,-1])
    ## cat("predicting done",format(Sys.time()),"\n")
    
    ## format ObsTrack
    ## cat("format obstrack",format(Sys.time()),"\n")
    restrowname <- opt$gene$official_name[! opt$gene$official_name %in% rownames(rf_model$ObsTrack)] %>% unique
    tree.rest <- matrix(rep(0,length(restrowname)*opt$ntrees),nrow = length(restrowname), ncol=opt$ntrees)
    rownames(tree.rest) <- restrowname
    if(length(restrowname)>0){
        mt <- rbind(rf_model$ObsTrack,tree.rest)
    }else{
        mt <- rf_model$ObsTrack
    }
    
    mt[mt>0]=1
    mt <- 1-mt
    mt <- mt[order(rownames(mt)),]

    ## format AllPrediction, use risk p, greater is better
    ## cat("format allprediction",format(Sys.time()),"\n")
    ## allprediction0 <- matrix(lapply(score$AllPrediction,function(x){x[,1]}) %>% unlist,nrow=dim(score$AllPrediction[[1]])[1])
    allprediction1 <- matrix(lapply(score$AllPrediction,function(x){x[,2]}) %>% unlist,nrow=dim(score$AllPrediction[[1]])[1])

    rownames(allprediction1) <- rf_model$evi$official_name

    ## validation
    ## apply(allprediction0,1,mean) %>% head
    ## apply(allprediction1,1,mean) %>% head

    ##score$ProbPredictio

    nodata.rowname <- rownames(mt)[! rownames(mt) %in%rownames(allprediction1)]
    nodata.mt <- matrix(rep(NA,length(nodata.rowname)*opt$ntrees),nrow=length(nodata.rowname),)
    rownames(nodata.mt) <- nodata.rowname

    if(length(nodata.rowname)>0){
        mt.predict <- rbind(allprediction1,nodata.mt)
    }else{
        mt.predict <- allprediction1
    }
    
    mt.predict <- mt.predict[match(rownames(mt),rownames(mt.predict)),]

    ## valitation
    ## sum(allprediction1[rownames(allprediction1)=="ADAMTS8",] == mt.predict[1,] )
    ## sum(rownames(mt.predict) ==rownames(mt))

    ## cat("multiplication",format(Sys.time()),"\n")
    mt.mean <- data.table(official_name = rownames(mt), sump=apply(mt*mt.predict,1,sum),n=apply(mt,1,sum))[,pvalue:=sump/n]

 
    ## validation
    ##sump1=apply(mt*mt.predict,1,sum)

    if(dim(mt.mean[n<50,])[1]>0){
        print(mt.mean[n<50,])
        stop("Three is at least gene with too few trees")
    }

    FUN <- match.fun(opt$fillNA)
    pvaluecolname <- names(mt.mean)[4]
    mt.mean[is.na(pvalue),pvalue:=FUN(mt.mean$pvalue[!is.na(mt.mean$pvalue)])]

    names(mt.mean) <- c("official_name",paste(names(mt.mean)[-1],rf_model$feaname,sep="_"))
    ## cat("Done",format(Sys.time()),"\n")
    return(mt.mean)
}

