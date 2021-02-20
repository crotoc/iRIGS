#' Returning the network string.
#'
#' @param opt A list of parameters in Gibbs script
#' @return A network string based on the parameter inputed in Gibbs.
#' @export
#' @examples
#' opt$transp <- load_network(opt)

load_network <- function(opt){
####===================== loading data ====================####
    network <- opt$network
    rp <- opt$res_prob
    cand_network<-c("hi","go","coexpr","pina","string","gpc")
    if(!is.element(network,cand_network)) stop("Wrong network option: must be one of ",paste(cand_network,collapse=", "),"...!\n")  
    if(!any(rp==c(0.1,0.3,0.5,0.7)))  stop("Wrong restart probability option: must be one of 0.1, 0.3, 0.5, 0.7...!\n")  

###======================= specify different network =============================###
    go_path<-"/scratch/cgg/chenr6/Database/Gibbs_data_from_quan/"
    pina_path<-"/fs0/wangq6/Breast_Cancer_causal_gene/RWR_gibbs/PINA_propagation_probabilities/"
    coexpr_path<-"/fs0/wangq6/Temp/temp_for_QiangWei/co-expression/propagation_probability/"
    string_path<-"/home/wangq6/cgg/Breast_Cancer_causal_gene/STRING/"
    gpc_path<-"/fs0/jiy1/gene_extract/network/"
    hi_path<-"/home/wangq6/cgg/Feixiong_CVD/propogation_matrix/"

    if(network=="go")      transp<-paste(go_path,"go_propogation_probality_rp_",rp,".RData",sep="")
    if(network=="pina")    transp<-paste(pina_path,"propogation_probality_rp_",rp,".RData",sep="")
    if(network=="coexpr")  transp<-paste(coexpr_path,"coexpression_propogation_probality_rp_",rp,".RData",sep="")
    if(network=="string")  transp<-paste(string_path,"STRING_propogation_probality_rp_",rp,".RData",sep="")
    if(network=="gpc")     transp<-paste(gpc_path,"gpc_propogation_probality_rp_",rp,".RData",sep="")
    if(network=="hi")      transp<-paste(hi_path,"propogation_probality_rp_",rp,".RData",sep="")
    transp
}
