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
sampling <- function(opt_tmp,run_step,with_region){


    cat("\nMODE: ",run_step,"\n")
    cat("STEP: ",opt_tmp$run_mode,"\n")
    cat("WITHREGION: ",with_region,"\n")
    n=0
    
    if(run_step=="burnin"){
        n <- opt_tmp$burnin_round
    }else{
        n <- opt_tmp$after_burnin_round
    }

    cat("ROUND TOTAL: ",n,"\n")
    
    thres<-opt_tmp$dif_threshold; dif<-thres+1;

    if(length(opt_tmp$remaining)==0){
        cat("RANDOM GENES AS START: ",n,"\n")
        remaining<-unlist(lapply(opt_tmp$region,function(x) sample(x,1)))
    }else{
        cat("USE lastet BURNIN GENES AS START: ",n,"\n")
        remaining<-opt_tmp$remaining
    }

    cat("\nStart genes: \n ",remaining,"\n\n")

    count <- opt_tmp$gene[,c(1,3)]
    count[,num0:=0]
    count[,num1:=0]
    count[,round:=0]

    if(opt_tmp$run_mode=="identification" && opt_tmp$method != "NULL"){
        cat("\nuse count table from training step: \n ","\n\n")
        count1 <- opt_tmp$training_count
        print(count1)
    }
    
    region_chosen<-rep(T,length(opt_tmp$region)); ##region_risk<-NULL
    chosen <- NULL
    region_risk <- NULL

    nblock <- ceiling(length(opt_tmp$region)/opt_tmp$threads)
    blocks <- data.table(start=(seq(1,opt_tmp$threads)-1)*nblock+1)
    blocks[,end:=start+nblock-1]
    blocks <- blocks[start <= length(opt_tmp$region),]
    blocks <-blocks[end>length(opt_tmp$region),end:=length(opt_tmp$region)]
    
    registerDoParallel(opt_tmp$threads)
    
    if(opt_tmp$progress_bar){
        pb <- txtProgressBar(min = 1,max = n, style = 3)
    }

    method_weight_all <- opt_tmp$gene[,c("SNP","official_name")]

    for(circle in 1:n)
    {

        ## Progress monitoring
        if(opt_tmp$progress_bar){
            setTxtProgressBar(pb, circle)
        }

        
        if(dif<thres && run_step =="burnin") break


        ## If choosing region is true, extract the region logical vector
        if(with_region){
            ## region_chosen<-unlist(foreach(i=1:nrow(blocks)) %dopar% unlist(lapply(blocks[i,]$start:blocks[i,]$end,function(x){region_chosen_update(x,region_chosen,remaining,opt_tmp)})))
            region_chosen<-unlist(foreach(pickup=1:length(opt_tmp$region)) %dopar% region_chosen_update(pickup,region_chosen,remaining,opt_tmp))
            region_risk<-rbind(region_risk,region_chosen)
        }


        ##Calculate the weights from different methods every certerns rounds.
        if(((opt_tmp$method_interval==1 && circle%%opt_tmp$method_interval==0) || (opt_tmp$method_interval!=1 && circle%%opt_tmp$method_interval==1)) && opt_tmp$run_mode=="identification" && opt_tmp$method =="svm"){
            ##SVM
            cat("\nRound ",circle," use SVM","\n")
            
            opt_tmp$method_weight <- foreach(pickup=1:length(opt_tmp$region),.combine = rbind) %dopar%  svm2weight(pickup,region_chosen,remaining,opt_tmp,count1)
            method_weight_all=merge(method_weight_all,opt_tmp$method_weight[,c("SNP","official_name","extra_weight")],by=c("SNP","official_name"))
            names(method_weight_all)[2+ceiling(circle/opt_tmp$method_interval)]=paste("svm",ceiling(circle/opt_tmp$method_interval),sep="")
            
        }else if(((opt_tmp$method_interval==1 && circle%%opt_tmp$method_interval==0) || (opt_tmp$method_interval!=1 && circle%%opt_tmp$method_interval==1)) && opt_tmp$run_mode=="identification" && opt_tmp$method =="rf"){
            ##RAMDOM FOREST
            cat("\nRound ",circle," use RF","\n")
            pickup <- sample(1:length(opt_tmp$region),1)
            opt_tmp$method_weight <- rf2weight(pickup,region_chosen,remaining,opt_tmp,count1)
            ## cat("merge begin")
            method_weight_all=merge(method_weight_all,opt_tmp$method_weight[,c("official_name","SNP","extra_weight")],by=c("official_name","SNP"))
            ## cat("merge end")
            names(method_weight_all)[2+ceiling(circle/opt_tmp$method_interval)]=paste("method",ceiling(circle/opt_tmp$method_interval),sep="")
        }
       
        if(circle == 1){
            cat("\nBegin to predicting genes...\n")
        }

        ## cat("\nBegin to predicting genes...\n")
        remaining<-foreach(pickup=1:length(opt_tmp$region),.combine=c) %dopar% predict_gene(pickup,region_chosen,remaining,opt_tmp,circle)
        ## cat("\nEnd to predicting genes...\n")
        
        ## remaining<-unlist(foreach(i=1:nrow(blocks)) %dopar% unlist(lapply(blocks[i,]$start:blocks[i,]$end,function(x){predict_gene(x,region_chosen,remaining,opt_tmp)})))
  
        ## region_risk<-rbind(region_risk,region_chosen)
        if(run_step == "after_burnin"){
            chosen<-rbind(chosen,remaining)
        }


        count <- count_processing(count,remaining,opt_tmp)

        if(run_step =="burnin" && circle>1){
            dif <- count$dif %>% sum %>% sqrt
        }

        
        if(opt_tmp$run_mode=="identification" && opt_tmp$method !="NULL"){
            count1 <- count_processing(count1,remaining,opt_tmp)
        }
    }

    stopImplicitCluster()
    
    out <- list()
    
    count <- count[,.SD[order(-p)],by=SNP]
    count <- count[,rank:=1:.N,by=SNP]

    count[rank==1,grp:="positive"]
    threshold = summary(count$p)[2]
    count[p<=threshold,grp:="negative"]
    
    out$count <- count
    print(out)

    if(opt_tmp$run_mode=="identification" && opt_tmp$method !="NULL"){
        out$count1 <- count1
    }
    
    out$remaining <- remaining
    
    if(run_step == "after_burnin"){
        out$chosen <- chosen %>% data.table
        names(out$chosen) <- names(opt_tmp$region)
    }

    if(opt_tmp$method!="NA"){
        out$method_weight=method_weight_all
    }

    if(with_region==TRUE){
        out$region_risk <- region_risk %>% data.table
        names(out$region_risk) <- names(opt_tmp$region)
        out$region_risk_p <- data.table(SNP=names(opt_tmp$region),region_p=out$region_risk[,colSums(.SD)]/opt_tmp$after_burnin_round)
        out$count <- merge(out$count, out$region_risk_p,by="SNP",all.x=TRUE)
    } else if(with_region==FALSE){
        out$count[,region_p:=1]
    }
    out
}



count_processing <- function(count,remaining,opt_tmp){
    count[,num1:=0]
    for(j in 1:length(opt_tmp$region)){
        count[SNP==names(opt_tmp$region[j]) & official_name==remaining[j],num1:=1]
    }
    
    count[,num1:=num0+num1]
    count[,round:=round+1]
    if(count$round[1]>1){
        count[,p0:=num0/(count$round-1)]
        count[,p:=num1/count$round]
        count[,dif:=(p-p0)^2]
    }
    count[,num0:=num1]
    return(count)
}
