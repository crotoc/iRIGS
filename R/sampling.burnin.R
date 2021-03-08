#' The main sampling function.
#'
#' @param opt A list of parameters in Gibbs script
#' mode String. burnin/post_burnin
#' with_region Logical. Whether the region risk is evaluated
#' @return A list.
#' count  A data.table including all posterior probability and region risk p.
#' chosen A data.table ncluding the genes chosen in each round.
#' region_risk A data.table including the TRUE/FALSE in each round.
#' @import data.table doParallel foreach
#' @export
#' @examples
#' burnin <- sampling(opt,mode="burnin",with_region=FALSE)
#' 
sampling_burnin <- function(opt){
    cat("STEP: burnin","\n")
    cat("WITHREGION: ",with_region,"\n")

####################################################################
#######		initiate parameter
#####################################################################
    param <- list()
    param$n <- opt$burnin_round    
    cat("ROUND TOTAL: ",n,"\n")
    param$thres <- opt$dif_threshold;
    param$dif <- param$thres+1;
    
    if(opt$burnin_start=="random"){
        cat("RANDOM GENES AS START: ",n,"\n")
        param$remaining<-unlist(lapply(opt$region,function(x) sample(x,1)))
    }else if(opt$burnin_start=="nearest"){
        ## remaining<-opt$remaining
    }

    cat("\nStart genes: \n ",remaining,"\n\n")

    param$count <- opt$gene[,c(1,3)]
    param$count[,num0:=0]
    param$count[,num1:=0]
    param$count[,round:=0]

    param$region_chosen<-rep(T,length(opt$region)); ##region_risk<-NULL

#####    
    out <- sampling_burnin_circle(param,opt)
    
    out <- list()
    out$count <- out$count[,.SD[order(-p)],by=SNP]
    out$count <- out$count[,rank:=1:.N,by=SNP]
    out$count[rank==1,grp:="positive"]
    threshold = summary(out$count$p)[2]
    out$count[p<=threshold,grp:="negative"]
    if(with_region==TRUE){
        out$region_risk <- region_risk %>% data.table
        names(out$region_risk) <- names(opt$region)
        out$region_risk_p <- data.table(SNP=names(opt$region),region_p=out$region_risk[,colSums(.SD)]/opt$after_burnin_round)
        out$count <- merge(out$count, out$region_risk_p,by="SNP",all.x=TRUE)
    } else if(with_region==FALSE){
        out$count[,region_p:=1]
    }
    print(out)    
    out
}



    
sampling_burnin_circle <- function(param,opt){    
    registerDoParallel(opt$threads)
    method_weight_all <- opt$gene[,c("SNP","official_name")]

    if(opt$progress_bar){pb <- txtProgressBar(min = 1,max = n, style = 3)}
    
    for(circle in 1:param$n){
        
        ## Progress monitoring
        if(opt$progress_bar){ setTxtProgressBar(pb, circle)}
        ## stop conditions
        if(param$dif<param$thres) break

        chosen <- NULL
        region_risk <- NULL
        
        ## If choosing region is true, extract the region logical vector
        if(opt$with_region){
            region_chosen<-unlist(foreach(pickup=1:length(opt$region)) %dopar% region_chosen_update(pickup,param$region_chosen,param$remaining,opt))
            region_risk<-rbind(region_risk,region_chosen)
        }

        if(circle == 1){cat("\nBegin to predicting genes...\n")}

        ## Predict genes
        remaining<-foreach(pickup=1:length(opt$region),.combine=c) %dopar% predict_gene(pickup,param$region_chosen,param$remaining,opt,circle)
        
        param$count <- count_processing(param$count,remaining,opt)

        if( circle>1){param$dif <- param$count$dif %>% sum %>% sqrt}
    }

    stopImplicitCluster()
    return(list(count=param$count,remining=remaining,region_risk=region_risk))
}
