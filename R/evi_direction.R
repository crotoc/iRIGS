
#' Conver the direction of evidence according to their column names.
#'
#' @param ori A data.table including the evidence data and column names containing the directons greater/less.
#' @return A data.table with correct directions (Greater is better)
#' @import data.table
#' @export
#' @examples
#' evi <- evi_direction(evi) 
#' 

evi_direction <- function(ori){
    if(names(ori)[1]!="official_name"){
        stop("The first column of evidence data has to be offical_name...!\n")
    }else{
        cat("check the first column of generic evi file ...PASS\n")
    }
    
    if(length(grep("binary|greater|less",names(ori)[-1],perl=T,ignore.case=T))!=(dim(ori)[2]-1)){
        stop("Every column name must contain binary|greater|less!!! ")
    } else {
        cat("Checking the generic evi file header...PASS\n")
    }

    myless <- grepl("less",colnames(ori)[-1],perl=T,ignore.case=T)
    mygreater <- grepl("greater",colnames(ori)[-1],perl=T,ignore.case=T)
    cat("Reverse ", sum(myless), paste(colnames(ori)[-1][myless],sep=";")," features\n")
    cat("Leave", sum(mygreater), paste(colnames(ori)[-1][mygreater],sep=";")," features\n")
    
    for(i in 2:ncol(ori)){
        ## if(grepl("binary",colnames(ori)[i],perl=T,ignore.case=T)){
        ##     N<- ori[,i,with=F] %>% sum
        ##     ori[[i]] <- as.numeric(ori[[i]])
        ##     ori[get(names(ori)[i])==1,(names(ori)[i]):=qnorm((N/nrow(ori)),lower.tail=F)]
        ##     ori[get(names(ori)[i])==0,(names(ori)[i]):=qnorm((1-N/nrow(ori)),lower.tail=F)]
        ## }
        if(grepl("less",colnames(ori)[i],perl=T,ignore.case=T)){
            ori[[i]] = -ori[[i]]
        }
    }
    
    ori <- ori[order(official_name),]
    ori <- ori[!duplicated(official_name),]
    
}



