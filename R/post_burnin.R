#' Two-step sampling bundle consist of a burnin step and a post burnin step.
#'
#' @param opt A list of parameters in Gibbs script
#' mode burnin/post_burnin mode
#' with_region Logical. Whether estimate region risk. 
#' @return A list containing results of sampling
#' @import data.table
#' @export
#' @examples
#' opt$training_gene<-post_burnin(opt,evi_mode="generic",with_region = FALSE,step="training")
#' 
post_burnin <- function(opt,evi_mode,with_region,run_mode="training"){

    opt_tmp <- opt
    opt_tmp$run_mode <- run_mode
    opt_tmp$evi_mode <- evi_mode


    if(opt_tmp$evi_mode == "generic"){
        cat("Use generic weight Matrix: \n")
        opt_tmp$extra_weight <- opt$generic$generic_weight
        opt_tmp$extra_processed <- opt$generic$generic_processed
        cat("Weight Matrix: \n")
        print(opt_tmp$extra_weight)
    }else if(opt_tmp$evi_mode == "evi"){
        if(opt_tmp$method == "svm"){
            cat("SVM use evi weight Matrix: \n")
            opt_tmp$extra_weight <- opt_tmp$evi$evi_weight
            opt_tmp$extra_processed <- opt_tmp$evi$evi_processed
            opt_tmp$training_count <- copy(opt_tmp$training_gene$count)
        }
        if(opt_tmp$method == "rf"){
            opt_tmp$training_count <- copy(opt_tmp$training_gene$count)
            opt_tmp$extra_evi_list <- opt_tmp$evi$evi_list
        }
        if(opt_tmp$method == "NULL"){
            opt_tmp$extra_weight <- opt_tmp$evi$evi_weight
            opt_tmp$extra_processed <- opt$evi$evi_processed
            cat("Weight Matrix: \n")
            print(opt_tmp$extra_weight)
        }
    }

    ##TEST
    #run_mode = "training"
    #evi_mode = "generic"
    
    
    cat("MODE: ",opt_tmp$evi_mode,"\n")
    cat("RUN_MODE: ",opt_tmp$run_mode,"\n")



    
    ##Burnin
    if(run_mode=="training"){
        burnin <- sampling(opt_tmp,run_step="burnin",with_region=FALSE)
        opt_tmp$remaining <- burnin$remaining
    }else{
        cat("Skip burnin step for identification step and use traning remaining and count...")
        opt_tmp$remaining <- opt$training_gene$remaining
    }
    
    ##after_burnin
    after_burnin <- sampling(opt_tmp,run_step="after_burnin",with_region=with_region)
    
    after_burnin
}
