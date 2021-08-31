#' Create FC groups with counts based on the number of genes wanted
#'
#' @inheritParams error_hierarchicell
#' @return named vector containing the FC value and the number of genes to be 
#' derived for it
#' @keywords internal
get_deg_groups<-
    function(n_genes,fc_range_min,fc_range_max,seed,deg_fc_interval) {
    #if user wanted set the seed
    if(!isFALSE(seed))
        set.seed(seed)
    all_fcs <- seq(from=fc_range_min,to=fc_range_max,by=deg_fc_interval)
    counts <- as.vector(table(sample(all_fcs, size = n_genes, replace = T)))
    names(counts) <- all_fcs
    return(counts)
}