#' Efficient function to calculate Matthew's Correlation Co-efficient (MCC)
#'
#' @param predicted predicted continuous values as vector
#' @param actual actual values (0,1)
#' @param cutoff cut-off for calculation, made for p-values so backwards to 
#' usual cut-offs. i.e. if prediciton<cut-off, it is a true prediction
#' @return Matthew's Correlation Co-efficient (MCC)
#' @keywords internal
calc_mcc<-function(predicted,actual,cutoff){
    predicted<-ifelse(predicted<cutoff,1,0)
    FP <- as.double(sum( predicted & !actual,na.rm = TRUE))
    FN <- as.double(sum(!predicted &  actual,na.rm = TRUE))
    act_pos <- as.double(sum(actual))
    act_neg <- as.double(length(actual)-act_pos)
    TN <- act_neg - FP
    TP <- act_pos - FN
    message("TP: ",TP," FP: ",FP," TN: ",TN," FN: ",FN)
    if(is.na(TN)) TN <- 0
    if(is.na(TP)) TP <- 0
    if(is.na(FN)) FN <- 0
    if(is.na(FP)) FP <- 0
    MCC <- (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    if(is.na(MCC)) MCC <- 0
    #return Type 1 error as test
    return(list("MCC"=MCC,"type_1_err"=(FP/(FP+TN))))
}