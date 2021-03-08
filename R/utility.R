#' The main sampling function.
#'
#' @param
#' opt A list of parameters in Gibbs script
#' remaining A vector from predict_gene
#' count A table include the sampling count of each round of iteration
#' 
#'
#' @return A talbe
#' count  A data.table including all posterior probability and region risk p.
#' 
#'
#' @import
#' data.table
#' pacman
#' optparse
#' magrittr
#' stringi
#' devtools
#' roxygen2
#'
#' @export
#'
#' @examples
#' 
#' 

count_processing <- function(count,remaining,opt){
    count[,num1:=0]
    for(j in 1:length(opt$region)){
        count[SNP==names(opt$region[j]) & official_name==remaining[j],num1:=1]
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

