#' @title Computing Type 1 Error for Single Cell Expression Case-Control Differential
#'   Expression Analysis
#'
#' @name compute_error
#'
#' @description Computes type 1 error for single cell data that is cell-type specifc,
#'   hierarchical, and compositonal. This function computes type 1 error with the
#'   single-cell differential expression analysis tool 'MAST', using random
#'   effects to account for the correlation structure that exists among measures
#'   from cells within an individual. The type 1 error calculations will borrow
#'   information from the input data (or the package default data) to simulate
#'   data under a variety of pre-determined conditions. These conditions include
#'   foldchange, number of genes, number of samples (i.e., independent
#'   experimental units), and the mean number of cells per individual.
#'
#' @details Prior to running the \code{\link{error_hierarchicell}} function, it
#'   is important to run the \code{\link{filter_counts}} function followed by
#'   the \code{\link{compute_data_summaries}} function to build an R object that
#'   is in the right format for the following simulation function to properly
#'   work.
#'
#' @note Data should be \strong{only for cells of the specific cell-type} you
#'   are interested in simulating or computing type 1 error for. Data should also
#'   contain as many unique sample identifiers as possible. If you are inputing
#'   data that has less than 5 unique values for sample identifier (i.e.,
#'   independent experimental units), then the empirical estimation of the
#'   inter-individual heterogeneity is going to be very unstable. Finding such a
#'   dataset will be difficult at this time, but, over time (as experiments grow
#'   in sample size and the numbers of publically available single-cell RNAseq
#'   datasets increase), this should improve dramatically.
#'

NULL

#'@title Compute Type 1 Error for Single Cell Expression Case-Control Analysis
#'
#'@rdname error_hierarchicell
#'
#'@description Computes type 1 error for single cell data that is cell-type
#'  specifc, hierarchical, and compositonal. This function computes type 1 error
#'  with the single-cell differential expression analysis tool 'MAST', using
#'  random effects to account for the correlation structure that exists among
#'  measures from cells within an individual. The type 1 error calculations will
#'  borrow information from the input data (or the package default data) to
#'  simulate data under a variety of pre-determined conditions. These conditions
#'  include foldchange, number of genes, number of samples (i.e., independent
#'  experimental units), and the mean number of cells per individual.
#'
#'@details Prior to running the \code{\link{error_hierarchicell}} function, it
#'  is important to run the \code{\link{filter_counts}} function followed by the
#'  \code{\link{compute_data_summaries}} function to build an R object that is
#'  in the right format for the following simulation function to properly work.
#'
#'@note Data should be \strong{only for cells of the specific cell-type} you are
#'  interested in simulating or computing type 1 error for. Data should also
#'  contain as many unique sample identifiers as possible. If you are inputing
#'  data that has less than 5 unique values for sample identifier (i.e.,
#'  independent experimental units), then the empirical estimation of the
#'  inter-individual heterogeneity is going to be very unstable. Finding such a
#'  dataset will be difficult at this time, but, over time (as experiments grow
#'  in sample size and the numbers of publically available single-cell RNAseq
#'  datasets increase), this should improve dramatically.
#'
#'@param data_summaries an R object that has been output by the package's
#'  compute_data_summaries function. No default
#'
#'@param method a name. The method for differential expression to be used for
#'  the computation of error. Possible methods include: MAST with random effects
#'  ("MAST_RE"), MAST ("MAST"), MAST with batch effect correction
#'  ("MAST_Combat"), GLM assuming a tweedie distribution ("GLM_tweedie"), GLMM
#'  assuming a tweedie distribution ("GLMM_tweedie"), generalized estimating
#'  equations ("GEE1"), ROTS ("ROTS"), or Monocle ("Monocle"). Defaults to
#'  "MAST_RE" which is the currently recommended analysis pipeline for
#'  single-cell data. See \code{\link{de_methods}} for more details on each of
#'  the methods.
#'
#'@param n_genes an integer. The number of genes you would like to simulate for
#'  your dataset for both differentially and non-differentially expressed genes.
#'  Too large of a number may cause memory failure and may slow the simulation 
#'  down tremendously. We recommend an integer less than 40,000.
#'  Defaults to 10,000.
#'
#'@param n_per_group an integer. The number of independent samples per
#'  case/control group for simulation. Creates a balanced design, for unbalanced
#'  designs, specify n_cases and n_controls separately. If not specifying a
#'  foldchange, the number of cases and controls does not matter. Defaults to 3.
#'
#'@param n_cases an integer. The number of independent control samples for
#'  simulation. Defaults to n_per_group.
#'
#'@param n_controls an integer. The number of independent case samples for
#'  simulation. Defaults to n_per_group.
#'
#'@param cells_per_control an integer. The mean number of cells per
#'  control you would like to simulate. Too large of a number may cause
#'  memory failure and may slow the simulation down tremendously. We recommend
#'  an integer less than 300, but more is possible. We note that anything
#'  greater than 100, brings marginal improvements in type 1 error. Defaults to
#'  100.
#'
#'@param cells_per_case an integer. The mean number of cells per
#'  case you would like to simulate. Too large of a number may cause
#'  memory failure and may slow the simulation down tremendously. We recommend
#'  an integer less than 300, but more is possible. We note that anything
#'  greater than 100, brings marginal improvements in type 1 error. Defaults to
#'  100.
#'
#'
#'@param ncells_variation_type either "Poisson", "NB", or "Fixed". Allows the
#'  number of cells per individual to be fixed at exactly the specified number
#'  of cells per individual, vary slightly with a poisson distribution with a
#'  lambda equal to the specified number of cells per individual, or a negative
#'  binomial with a mean equal to the specified number of cells and dispersion
#'  size equal to one.Defaults to "Poisson".
#'
#'
#'@param pval a number. The significance threshold (alpha) to use for
#'  significance. Defaults to 0.05. Can also be a vector of pvalue - up to a
#'  length of 5.
#'
#'@param alter_dropout_cases a numeric proportion between 0 and 1. The
#'  proportion by which you would like to simulate decreasing the amount of
#'  dropout between case control groups. For example, if you would like to
#'  simulate a decrease in the amount of dropout in your cases by twenty
#'  percent, then 0.2 would be appropriate. This component of the simulation
#'  allows the user to adjust the proportion of dropout if they believe the
#'  stochastic expression of a gene will differ between cases and controls. For
#'  a two-part hurdle model, like MAST implements, this will increase your
#'  ability to detect differences. Defaults to 0.
#'
#'@param decrease_dropout a numeric proportion between 0 and 1. The proportion
#'  by which you would like to simulate decreasing the amount of dropout in your
#'  data. For example, if you would like to simulate a decrease in the amount of
#'  dropout in your data by twenty percent, then 0.2 would be appropriate. This
#'  component of the simulation allows the user to adjust the proportion of
#'  dropout if they believe future experiments or runs will have improved
#'  calling rates (due to improved methods or improved cell viability) and
#'  thereby lower dropout rates. Defaults to 0.
#'
#'@param fc_range_min a number between 1 and 10. The minimum amount of fold 
#'  change to simulate a difference in expression between case and control 
#'  groups. This is used to derive the differentially expressed genes. The
#'  fold change changes genes in either direction, so a fold change of 2 would
#'  cause the mean expression in cases to be either twice the amount or half the
#'  amount for any particular gene. Defaults 1.1 which arguably would be hard to
#'  by a model since the FC difference between case and controls is so low..
#'  
#'@param fc_range_max a number between 1 and 10. The maximum amount of fold 
#'  change to simulate a difference in expression between case and control 
#'  groups. This is used to derive the differentially expressed genes. The
#'  fold change changes genes in either direction, so a fold change of 2 would
#'  cause the mean expression in cases to be either twice the amount or half the
#'  amount for any particular gene. Defaults 3. n_genes will be randomly 
#'  sampled between fc_range_min and fc_range_max at an interval of 0.1.
#'  
#'@param deg_fc_interval the incremental value between fc_range_min and 
#'  fc_range_max used to derive the number of fold change groups. This uniform 
#'  sampling of fold change values is used to create DEGs from a distribution 
#'  with the specified fold change. Default is 0.5
#'  
#'@param seed Use this to ensure the same count datasets are generated in 
#'  subsequent runs, this is helpful when comparing models so we know the exact
#'  same dataset was used. Default False, so no seed is set. If integer value
#'  passed, seed is set using it.
#'@param checkpoint add location to save results at checkpoints in the code so 
#'  if the run is cancelled it can be pisked up where it left off. Default of 
#'  NULL will not save checkpoints. This becomes very useful as the number of 
#'  cells and number of samples increases down to mixture models' long runtime.
#'  
#'@param force_new if checkpoint set to a folder, force_new set to TRUE will 
#'  rerun all the code even if a checkpoint has been saved. Default is FALSE. 
#'  Even if FALSE used and checkpoint saved, if new input parameters are used 
#'  affected parts of the function will be rerun.  
#'
#'@return The estimated error under the specified conditions when using 'MAST'
#'  with random effects to account for the correlation structure that exists
#'  among measures from cells within an individual.
#'
#'
#'@export

error_hierarchicell <- function(data_summaries,
                                   method = "MAST_RE",
                                   n_genes = 10000,
                                   n_per_group = 3,
                                   n_cases = n_per_group,
                                   n_controls = n_per_group,
                                   cells_per_case = 100,
                                   cells_per_control = 100,
                                   ncells_variation_type = "Poisson",
                                   pval = 0.05,
                                   decrease_dropout = 0,
                                   alter_dropout_cases = 0,
                                   fc_range_min = 1.1,
                                   fc_range_max = 10,
                                   deg_fc_interval = 0.5,
                                   seed = FALSE,
                                   checkpoint = NULL,
                                   force_new = FALSE){
  if(!is.null(checkpoint)){
    #### Do a bit of QC to get the full path ####
    ## Expand "~" into full path bc isn't recognized in certain envs(eg Python)
    checkpoint <- path.expand(checkpoint)
    ## Expand relative path "./" into absolute path bc it's less ambiguous
    checkpoint <- gsub("^[.]/", paste0(getwd(), "/"), checkpoint)
    ## Remove trailing / if present
    checkpoint <- gsub("/$", "", checkpoint)
  }
  
  if (method == "MAST_RE") {


    if (!requireNamespace(c("MAST","SummarizedExperiment","lme4"),quietly = TRUE)){
      stop("The packages 'MAST', 'lme4', 'fitdistrplus', and \n
           'SummarizedExperiment' are required. Please install them.\n
           It may be a problem with dependencies for these packages,\n
           type '!requireNamespace(\"MAST\")' to see if this is the issue.",
           call. = FALSE)
    } else {
      if(!is.numeric(fc_range_min)||fc_range_min<1)
        stop("fc_range_min should be between 1-10")
      if(!is.numeric(fc_range_max)||fc_range_max<1.1)
        stop("fc_range_max should be between 1.1-10")
      if(fc_range_min>=fc_range_max)
        stop("fc_range_max should be greater than fc_range_min")
      if(!is.numeric(deg_fc_interval)||deg_fc_interval<=0||deg_fc_interval>(fc_range_max-fc_range_min))
        stop("deg_fc_interval should be numeric and smaller than the difference between fc_range_max and fc_range_min")
      
      if (cells_per_control < 50 | cells_per_case < 50){
        message("----------------------------------------------")
        message("Mean number of cells per individual is less than 50.\n The probability of complete separation will start to increase.")
        message("----------------------------------------------")
      }
      #simulate the non differentially expressed genes
      #if checkpoint saved and not forced, don't rerun
      checkpoint_params <- 
        paste(method,n_genes,n_per_group,n_cases,n_controls,cells_per_case,
                cells_per_control,ncells_variation_type,pval,decrease_dropout,
                alter_dropout_cases,seed,sep="_")
      if(!force_new && !(!is.null(checkpoint) && 
          file.exists(paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds")))){
        #if user wanted set the seed
        if(!isFALSE(seed))
          set.seed(seed)
        non_deg <- suppressMessages(simulate_hierarchicell(data_summaries,
                                                             n_genes = n_genes,
                                                             n_per_group = n_per_group,
                                                             n_cases = n_cases,
                                                             n_controls = n_controls,
                                                             cells_per_case = cells_per_case,
                                                             cells_per_control = cells_per_control,
                                                             ncells_variation_type = ncells_variation_type,
                                                             foldchange = 1, #fc=1 means not DEG
                                                             decrease_dropout = decrease_dropout,
                                                             alter_dropout_cases = alter_dropout_cases))
        
        #non degs
        genecounts_ndeg <- as.matrix(t(non_deg[,c(-1,-2,-3)]))
        coldata_ndeg <- non_deg[,1:3]
        coldata_ndeg$Status <- as.factor(coldata_ndeg$Status)
        genecounts_ndeg <- genecounts_ndeg[which(apply(genecounts_ndeg, 1, mean) > 5), ]
        genecounts_ndeg <- genecounts_ndeg[,rownames(coldata_ndeg)]
        log2counts_ndeg <- log2(genecounts_ndeg + 1)
        
        fData_ndeg <- data.frame(primerid=rownames(genecounts_ndeg))
        sca_ndeg <- suppressMessages(MAST::FromMatrix(exprsArray=log2counts_ndeg, cData=coldata_ndeg, fData=fData_ndeg))
        cdr2_ndeg <- colSums(SummarizedExperiment::assay(sca_ndeg)>0)
        SummarizedExperiment::colData(sca_ndeg)$ngeneson <- scale(cdr2_ndeg)
        SummarizedExperiment::colData(sca_ndeg)$Status <-
          factor(SummarizedExperiment::colData(sca_ndeg)$Status)
        SummarizedExperiment::colData(sca_ndeg)$DonorID <-
          factor(SummarizedExperiment::colData(sca_ndeg)$DonorID)
        zlmCond_ndeg <- suppressMessages(MAST::zlm(~ ngeneson + Status + (1 | DonorID),
                                                   sca_ndeg, method='glmer',ebayes = F,
                                                   strictConvergence = FALSE))
        
        summaryCond_ndeg <- suppressMessages(MAST::summary(zlmCond_ndeg,
                                                           doLRT='StatusControl'))
        summaryDt_ndeg <- summaryCond_ndeg$datatable
        fcHurdle_ndeg <- merge(summaryDt_ndeg[summaryDt_ndeg$contrast=='StatusControl'
                                              & summaryDt_ndeg$component=='logFC', c(1,7,5,6,8)],
                               summaryDt_ndeg[summaryDt_ndeg$contrast=='StatusControl'
                                              & summaryDt_ndeg$component=='H', c(1,4)],
                               by = 'primerid')
        
        fcHurdle_ndeg <- stats::na.omit(as.data.frame(fcHurdle_ndeg))
        
        if(!is.null(checkpoint))
          saveRDS(fcHurdle_ndeg,
                    file=paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds"))
      }else{
        #load checkpoint data
        fcHurdle_ndeg <- 
          readRDS(paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds"))
      }
      
      
      #simulate the differentially expressed genes
      deg_groups <-
        get_deg_groups(n_genes,fc_range_min,fc_range_max,seed,deg_fc_interval)
      #separate call for each
      degs_sep_fcs <- vector(mode="list",length=length(deg_groups))
      names(degs_sep_fcs) <- names(deg_groups)
      #if user wanted set the seed
      if(!isFALSE(seed))
        set.seed(seed)
      #cycle through fc choices
      for (fc_i in as.numeric(names(deg_groups))){
        if(!force_new && !(!is.null(checkpoint) && 
                           file.exists(paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds")))){
          degs <- 
            suppressMessages(simulate_hierarchicell(data_summaries,
                                                    n_genes = deg_groups[[as.character(fc_i)]],
                                                    n_per_group = n_per_group,
                                                    n_cases = n_cases,
                                                    n_controls = n_controls,
                                                    cells_per_case = 
                                                      cells_per_case,
                                                    cells_per_control = 
                                                      cells_per_control,
                                                    ncells_variation_type = 
                                                      ncells_variation_type,
                                                    foldchange = fc_i,
                                                    decrease_dropout = 
                                                      decrease_dropout,
                                                    alter_dropout_cases = 
                                                      alter_dropout_cases))
          
          #degs
          genecounts_degs <- as.matrix(t(degs[,c(-1,-2,-3)]))
          coldata_degs <- degs[,1:3]
          coldata_degs$Status <- as.factor(coldata_degs$Status)
          genecounts_degs <- genecounts_degs[which(apply(genecounts_degs, 1, mean) > 5), ]
          genecounts_degs <- genecounts_degs[,rownames(coldata_degs)]
          log2counts_degs <- log2(genecounts_degs + 1)
          
          fData_degs <- data.frame(primerid=rownames(genecounts_degs))
          sca_degs <- suppressMessages(MAST::FromMatrix(exprsArray=log2counts_degs, cData=coldata_degs, fData=fData_degs))
          cdr2_degs <- colSums(SummarizedExperiment::assay(sca_degs)>0)
          SummarizedExperiment::colData(sca_degs)$ngeneson <- scale(cdr2_degs)
          SummarizedExperiment::colData(sca_degs)$Status <-
            factor(SummarizedExperiment::colData(sca_degs)$Status)
          SummarizedExperiment::colData(sca_degs)$DonorID <-
            factor(SummarizedExperiment::colData(sca_degs)$DonorID)
          zlmCond_degs <- suppressMessages(MAST::zlm(~ ngeneson + Status + (1 | DonorID),
                                                     sca_degs, method='glmer',ebayes = F,
                                                     strictConvergence = FALSE))
          
          summaryCond_degs <- suppressMessages(MAST::summary(zlmCond_degs,
                                                             doLRT='StatusControl'))
          summaryDt_degs <- summaryCond_degs$datatable
          fcHurdle_degs <- merge(summaryDt_degs[summaryDt_degs$contrast=='StatusControl'
                                                & summaryDt_degs$component=='logFC', c(1,7,5,6,8)],
                                 summaryDt_degs[summaryDt_degs$contrast=='StatusControl'
                                                & summaryDt_degs$component=='H', c(1,4)],
                                 by = 'primerid')
          if(!is.null(checkpoint))
            saveRDS(fcHurdle_degs,
                    file=paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds"))
        }else{
          #load checkpoint data
          fcHurdle_degs <- 
            readRDS(paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds"))
        }  
        
        degs_sep_fcs[[as.character(fc_i)]] <- 
          stats::na.omit(as.data.frame(fcHurdle_degs))
        
      }
      #just combine the results
      fcHurdle_degs <- rbindlist(degs_sep_fcs,idcol="fc")
      
      if (length(pval) == 1){
        #create actuals
        act_neg <- rep(0,length(fcHurdle_ndeg[,6]))
        #dt
        act_pos <- rep(1,length(fcHurdle_degs[[6]]))
        #combine pos neg
        pred <- c(fcHurdle_ndeg[,6],fcHurdle_degs[[6]])
        actuals <- c(act_neg,act_pos)
        
        signif <- ifelse(fcHurdle_ndeg[,6] < pval, 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval," is: ", rate))
        return(calc_mcc(predicted=pred,actual=actuals,cutoff=pval))

      } else if (length(pval) == 2) {

        signif <- ifelse(fcHurdle[,6] < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

      } else if (length(pval) == 3) {

        signif <- ifelse(fcHurdle[,6] < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))


      } else if (length(pval) == 4) {

        signif <- ifelse(fcHurdle[,6] < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[4], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[4]," is: ", rate))

      } else if (length(pval) == 5) {

        signif <- ifelse(fcHurdle[,6] < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[4], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[4]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[5], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[5]," is: ", rate))

      } else {

        message("Too many pvalues, shorten vector of pvalues to 5 or less")

      }

    }

  } else if (method == "MAST") {
    if (!requireNamespace(c("MAST","SummarizedExperiment"),quietly = TRUE)){
      stop("The packages 'MAST', 'fitdistrplus', and \n
           'SummarizedExperiment' are required. Please install them.\n
           It may be a problem with dependencies for these packages,\n
           type '!requireNamespace(\"MAST\")' to see if this is the issue.",
           call. = FALSE)
    } else {

      if(!is.numeric(fc_range_min)||fc_range_min<1)
        stop("fc_range_min should be between 1-10")
      if(!is.numeric(fc_range_max)||fc_range_max<1.1)
        stop("fc_range_max should be between 1.1-10")
      if(fc_range_min>=fc_range_max)
        stop("fc_range_max should be greater than fc_range_min")
      if(!is.numeric(deg_fc_interval)||deg_fc_interval<=0||deg_fc_interval>(fc_range_max-fc_range_min))
        stop("deg_fc_interval should be numeric and smaller than the difference between fc_range_max and fc_range_min")
      


      if (cells_per_control < 50 | cells_per_case < 50){
        message("----------------------------------------------")
        message("Mean number of cells per individual is less than 50.\n The probability of complete separation will start to increase.")
        message("----------------------------------------------")
      }

      #simulate the non differentially expressed genes
      #if checkpoint saved and not forced, don't rerun
      checkpoint_params <- 
        paste(method,n_genes,n_per_group,n_cases,n_controls,cells_per_case,
              cells_per_control,ncells_variation_type,pval,decrease_dropout,
              alter_dropout_cases,seed,sep="_")
      if(!force_new && !(!is.null(checkpoint) && 
                         file.exists(paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds")))){
        #if user wanted set the seed
        if(!isFALSE(seed))
          set.seed(seed)
        non_deg <- suppressMessages(simulate_hierarchicell(data_summaries,
                                                           n_genes = n_genes,
                                                           n_per_group = n_per_group,
                                                           n_cases = n_cases,
                                                           n_controls = n_controls,
                                                           cells_per_case = cells_per_case,
                                                           cells_per_control = cells_per_control,
                                                           ncells_variation_type = ncells_variation_type,
                                                           foldchange = 1, #fc=1 means not DEG
                                                           decrease_dropout = decrease_dropout,
                                                           alter_dropout_cases = alter_dropout_cases))
        
        #non degs
        genecounts_ndeg <- as.matrix(t(non_deg[,c(-1,-2,-3)]))
        coldata_ndeg <- non_deg[,1:3]
        coldata_ndeg$Status <- as.factor(coldata_ndeg$Status)
        genecounts_ndeg <- genecounts_ndeg[which(apply(genecounts_ndeg, 1, mean) > 5), ]
        genecounts_ndeg <- genecounts_ndeg[,rownames(coldata_ndeg)]
        log2counts_ndeg <- log2(genecounts_ndeg + 1)
        
        fData_ndeg <- data.frame(primerid=rownames(genecounts_ndeg))
        sca_ndeg <- suppressMessages(MAST::FromMatrix(exprsArray=log2counts_ndeg, cData=coldata_ndeg, fData=fData_ndeg))
        cdr2_ndeg <- colSums(SummarizedExperiment::assay(sca_ndeg)>0)
        SummarizedExperiment::colData(sca_ndeg)$ngeneson <- scale(cdr2_ndeg)
        SummarizedExperiment::colData(sca_ndeg)$Status <-
          factor(SummarizedExperiment::colData(sca_ndeg)$Status)
        SummarizedExperiment::colData(sca_ndeg)$DonorID <-
          factor(SummarizedExperiment::colData(sca_ndeg)$DonorID)
        zlmCond_ndeg <- suppressMessages(MAST::zlm(~ ngeneson + Status,
                                                   sca_ndeg, method='glm',ebayes = F))
        
        summaryCond_ndeg <- suppressMessages(MAST::summary(zlmCond_ndeg,
                                                           doLRT='StatusControl'))
        summaryDt_ndeg <- summaryCond_ndeg$datatable
        fcHurdle_ndeg <- merge(summaryDt_ndeg[summaryDt_ndeg$contrast=='StatusControl'
                                              & summaryDt_ndeg$component=='logFC', c(1,7,5,6,8)],
                               summaryDt_ndeg[summaryDt_ndeg$contrast=='StatusControl'
                                              & summaryDt_ndeg$component=='H', c(1,4)],
                               by = 'primerid')
        
        fcHurdle_ndeg <- stats::na.omit(as.data.frame(fcHurdle_ndeg))
        
        if(!is.null(checkpoint))
          saveRDS(fcHurdle_ndeg,
                  file=paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds"))
      }else{
        #load checkpoint data
        fcHurdle_ndeg <- 
          readRDS(paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds"))
      }
      
      #simulate the differentially expressed genes
      deg_groups <-
        get_deg_groups(n_genes,fc_range_min,fc_range_max,seed,deg_fc_interval)
      #separate call for each
      degs_sep_fcs <- vector(mode="list",length=length(deg_groups))
      names(degs_sep_fcs) <- names(deg_groups)
      #if user wanted set the seed
      if(!isFALSE(seed))
        set.seed(seed)
      #cycle through fc choices
      for (fc_i in as.numeric(names(deg_groups))){
        if(!force_new && !(!is.null(checkpoint) && 
                           file.exists(paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds")))){
          degs <- 
            suppressMessages(simulate_hierarchicell(data_summaries,
                                                    n_genes = deg_groups[[as.character(fc_i)]],
                                                    n_per_group = n_per_group,
                                                    n_cases = n_cases,
                                                    n_controls = n_controls,
                                                    cells_per_case = 
                                                      cells_per_case,
                                                    cells_per_control = 
                                                      cells_per_control,
                                                    ncells_variation_type = 
                                                      ncells_variation_type,
                                                    foldchange = fc_i,
                                                    decrease_dropout = 
                                                      decrease_dropout,
                                                    alter_dropout_cases = 
                                                      alter_dropout_cases))
          
          #degs
          genecounts_degs <- as.matrix(t(degs[,c(-1,-2,-3)]))
          coldata_degs <- degs[,1:3]
          coldata_degs$Status <- as.factor(coldata_degs$Status)
          genecounts_degs <- genecounts_degs[which(apply(genecounts_degs, 1, mean) > 5), ]
          genecounts_degs <- genecounts_degs[,rownames(coldata_degs)]
          log2counts_degs <- log2(genecounts_degs + 1)
          
          fData_degs <- data.frame(primerid=rownames(genecounts_degs))
          sca_degs <- suppressMessages(MAST::FromMatrix(exprsArray=log2counts_degs, cData=coldata_degs, fData=fData_degs))
          cdr2_degs <- colSums(SummarizedExperiment::assay(sca_degs)>0)
          SummarizedExperiment::colData(sca_degs)$ngeneson <- scale(cdr2_degs)
          SummarizedExperiment::colData(sca_degs)$Status <-
            factor(SummarizedExperiment::colData(sca_degs)$Status)
          SummarizedExperiment::colData(sca_degs)$DonorID <-
            factor(SummarizedExperiment::colData(sca_degs)$DonorID)
          zlmCond_degs <- suppressMessages(MAST::zlm(~ ngeneson + Status,
                                                     sca_degs, method='glm',ebayes = F))
          
          summaryCond_degs <- suppressMessages(MAST::summary(zlmCond_degs,
                                                             doLRT='StatusControl'))
          summaryDt_degs <- summaryCond_degs$datatable
          fcHurdle_degs <- merge(summaryDt_degs[summaryDt_degs$contrast=='StatusControl'
                                                & summaryDt_degs$component=='logFC', c(1,7,5,6,8)],
                                 summaryDt_degs[summaryDt_degs$contrast=='StatusControl'
                                                & summaryDt_degs$component=='H', c(1,4)],
                                 by = 'primerid')
          
          if(!is.null(checkpoint))
            saveRDS(fcHurdle_degs,
                    file=paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds"))
        }else{
          #load checkpoint data
          fcHurdle_degs <- 
            readRDS(paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds"))
        }  
        
        degs_sep_fcs[[as.character(fc_i)]] <- 
          stats::na.omit(as.data.frame(fcHurdle_degs))
        
      }
      #just combine the results
      fcHurdle_degs <- rbindlist(degs_sep_fcs,idcol="fc")
      
      
      if (length(pval) == 1){
        #create actuals
        act_neg <- rep(0,length(fcHurdle_ndeg[,6]))
        #dt
        act_pos <- rep(1,length(fcHurdle_degs[[6]]))
        #combine pos neg
        pred <- c(fcHurdle_ndeg[,6],fcHurdle_degs[[6]])
        actuals <- c(act_neg,act_pos)
        
        signif <- ifelse(fcHurdle_ndeg[,6] < pval, 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval," is: ", rate))
        return(calc_mcc(predicted=pred,actual=actuals,cutoff=pval))

      } else if (length(pval) == 2) {

        signif <- ifelse(fcHurdle[,6] < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

      } else if (length(pval) == 3) {

        signif <- ifelse(fcHurdle[,6] < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))


      } else if (length(pval) == 4) {

        signif <- ifelse(fcHurdle[,6] < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[4], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[4]," is: ", rate))

      } else if (length(pval) == 5) {

        signif <- ifelse(fcHurdle[,6] < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[4], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[4]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[5], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[5]," is: ", rate))

      } else {

        message("Too many pvalues, shorten vector of pvalues to 5 or less")

      }

    }


  } else if (method == "MAST_Combat") {
    if (!requireNamespace(c("MAST","SummarizedExperiment","sva"),quietly = TRUE)){
      stop("The packages 'MAST', 'fitdistrplus',  \n
           'SummarizedExperiment', and 'sva' are required. Please install them.\n
           It may be a problem with dependencies for these packages,\n
           type '!requireNamespace(\"MAST\")' to see if this is the issue.",
           call. = FALSE)
    } else {

      if(!is.numeric(fc_range_min)||fc_range_min<1)
        stop("fc_range_min should be between 1-10")
      if(!is.numeric(fc_range_max)||fc_range_max<1.1)
        stop("fc_range_max should be between 1.1-10")
      if(fc_range_min>=fc_range_max)
        stop("fc_range_max should be greater than fc_range_min")
      if(!is.numeric(deg_fc_interval)||deg_fc_interval<=0||deg_fc_interval>(fc_range_max-fc_range_min))
        stop("deg_fc_interval should be numeric and smaller than the difference between fc_range_max and fc_range_min")
      


      if (cells_per_control < 50 | cells_per_case < 50){
        message("----------------------------------------------")
        message("Mean number of cells per individual is less than 50.\n The probability of complete separation will start to increase.")
        message("----------------------------------------------")
      }

      #simulate the non differentially expressed genes
      
      #if checkpoint saved and not forced, don't rerun
      checkpoint_params <- 
        paste(method,n_genes,n_per_group,n_cases,n_controls,cells_per_case,
              cells_per_control,ncells_variation_type,pval,decrease_dropout,
              alter_dropout_cases,seed,sep="_")
      if(!force_new && !(!is.null(checkpoint) && 
                         file.exists(paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds")))){
        #if user wanted set the seed
        if(!isFALSE(seed))
          set.seed(seed)
        non_deg <- suppressMessages(simulate_hierarchicell(data_summaries,
                                                           n_genes = n_genes,
                                                           n_per_group = n_per_group,
                                                           n_cases = n_cases,
                                                           n_controls = n_controls,
                                                           cells_per_case = cells_per_case,
                                                           cells_per_control = cells_per_control,
                                                           ncells_variation_type = ncells_variation_type,
                                                           foldchange = 1, #fc=1 means not DEG
                                                           decrease_dropout = decrease_dropout,
                                                           alter_dropout_cases = alter_dropout_cases))
        
        #non degs
        genecounts_ndeg <- as.matrix(t(non_deg[,c(-1,-2,-3)]))
        coldata_ndeg <- non_deg[,1:3]
        coldata_ndeg$Status <- as.factor(coldata_ndeg$Status)
        genecounts_ndeg <- genecounts_ndeg[which(apply(genecounts_ndeg, 1, mean) > 5), ]
        genecounts_ndeg <- genecounts_ndeg[,rownames(coldata_ndeg)]
        log2counts_ndeg <- log2(genecounts_ndeg + 1)
        log2counts_ndeg <- sva::ComBat(log2counts_ndeg,coldata_ndeg$DonorID)
        
        fData_ndeg <- data.frame(primerid=rownames(genecounts_ndeg))
        sca_ndeg <- suppressMessages(MAST::FromMatrix(exprsArray=log2counts_ndeg, cData=coldata_ndeg, fData=fData_ndeg))
        cdr2_ndeg <- colSums(SummarizedExperiment::assay(sca_ndeg)>0)
        SummarizedExperiment::colData(sca_ndeg)$ngeneson <- scale(cdr2_ndeg)
        SummarizedExperiment::colData(sca_ndeg)$Status <-
          factor(SummarizedExperiment::colData(sca_ndeg)$Status)
        SummarizedExperiment::colData(sca_ndeg)$DonorID <-
          factor(SummarizedExperiment::colData(sca_ndeg)$DonorID)
        zlmCond_ndeg <- suppressMessages(MAST::zlm(~ ngeneson + Status,
                                              sca_ndeg, method='glm',ebayes = F))
        
        summaryCond_ndeg <- suppressMessages(MAST::summary(zlmCond_ndeg,
                                                           doLRT='StatusControl'))
        summaryDt_ndeg <- summaryCond_ndeg$datatable
        fcHurdle_ndeg <- merge(summaryDt_ndeg[summaryDt_ndeg$contrast=='StatusControl'
                                              & summaryDt_ndeg$component=='logFC', c(1,7,5,6,8)],
                               summaryDt_ndeg[summaryDt_ndeg$contrast=='StatusControl'
                                              & summaryDt_ndeg$component=='H', c(1,4)],
                               by = 'primerid')
        
        fcHurdle_ndeg <- stats::na.omit(as.data.frame(fcHurdle_ndeg))
        
        if(!is.null(checkpoint))
          saveRDS(fcHurdle_ndeg,
                  file=paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds"))
      }else{
        #load checkpoint data
        fcHurdle_ndeg <- 
          readRDS(paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds"))
      }
      
      #simulate the differentially expressed genes
      deg_groups <-
        get_deg_groups(n_genes,fc_range_min,fc_range_max,seed,deg_fc_interval)
      #separate call for each
      degs_sep_fcs <- vector(mode="list",length=length(deg_groups))
      names(degs_sep_fcs) <- names(deg_groups)
      #if user wanted set the seed
      if(!isFALSE(seed))
        set.seed(seed)
      #cycle through fc choices
      for (fc_i in as.numeric(names(deg_groups))){
        if(!force_new && !(!is.null(checkpoint) && 
                           file.exists(paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds")))){
          degs <- 
            suppressMessages(simulate_hierarchicell(data_summaries,
                                                    n_genes = deg_groups[[as.character(fc_i)]],
                                                    n_per_group = n_per_group,
                                                    n_cases = n_cases,
                                                    n_controls = n_controls,
                                                    cells_per_case = 
                                                      cells_per_case,
                                                    cells_per_control = 
                                                      cells_per_control,
                                                    ncells_variation_type = 
                                                      ncells_variation_type,
                                                    foldchange = fc_i,
                                                    decrease_dropout = 
                                                      decrease_dropout,
                                                    alter_dropout_cases = 
                                                      alter_dropout_cases))
          
          #degs
          genecounts_degs <- as.matrix(t(degs[,c(-1,-2,-3)]))
          coldata_degs <- degs[,1:3]
          coldata_degs$Status <- as.factor(coldata_degs$Status)
          genecounts_degs <- genecounts_degs[which(apply(genecounts_degs, 1, mean) > 5), ]
          genecounts_degs <- genecounts_degs[,rownames(coldata_degs)]
          log2counts_degs <- log2(genecounts_degs + 1)
          log2counts_degs <- sva::ComBat(log2counts_degs,coldata_degs$DonorID)
          
          fData_degs <- data.frame(primerid=rownames(genecounts_degs))
          sca_degs <- suppressMessages(MAST::FromMatrix(exprsArray=log2counts_degs, cData=coldata_degs, fData=fData_degs))
          cdr2_degs <- colSums(SummarizedExperiment::assay(sca_degs)>0)
          SummarizedExperiment::colData(sca_degs)$ngeneson <- scale(cdr2_degs)
          SummarizedExperiment::colData(sca_degs)$Status <-
            factor(SummarizedExperiment::colData(sca_degs)$Status)
          SummarizedExperiment::colData(sca_degs)$DonorID <-
            factor(SummarizedExperiment::colData(sca_degs)$DonorID)
          zlmCond_degs <- suppressMessages(MAST::zlm(~ ngeneson + Status,
                                                     sca_degs, method='glm',ebayes = F))
          
          summaryCond_degs <- suppressMessages(MAST::summary(zlmCond_degs,
                                                             doLRT='StatusControl'))
          summaryDt_degs <- summaryCond_degs$datatable
          fcHurdle_degs <- merge(summaryDt_degs[summaryDt_degs$contrast=='StatusControl'
                                                & summaryDt_degs$component=='logFC', c(1,7,5,6,8)],
                                 summaryDt_degs[summaryDt_degs$contrast=='StatusControl'
                                                & summaryDt_degs$component=='H', c(1,4)],
                                 by = 'primerid')
          if(!is.null(checkpoint))
            saveRDS(fcHurdle_degs,
                    file=paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds"))
        }else{
          #load checkpoint data
          fcHurdle_degs <- 
            readRDS(paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds"))
        }  
        
        degs_sep_fcs[[as.character(fc_i)]] <- 
          stats::na.omit(as.data.frame(fcHurdle_degs))
        
      }
      #just combine the results
      fcHurdle_degs <- rbindlist(degs_sep_fcs,idcol="fc")

      if (length(pval) == 1){
        #create actuals
        act_neg <- rep(0,length(fcHurdle_ndeg[,6]))
        #dt
        act_pos <- rep(1,length(fcHurdle_degs[[6]]))
        #combine pos neg
        pred <- c(fcHurdle_ndeg[,6],fcHurdle_degs[[6]])
        actuals <- c(act_neg,act_pos)
        
        signif <- ifelse(fcHurdle_ndeg[,6] < pval, 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval," is: ", rate))
        return(calc_mcc(predicted=pred,actual=actuals,cutoff=pval))

      } else if (length(pval) == 2) {

        signif <- ifelse(fcHurdle[,6] < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

      } else if (length(pval) == 3) {

        signif <- ifelse(fcHurdle[,6] < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))


      } else if (length(pval) == 4) {

        signif <- ifelse(fcHurdle[,6] < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[4], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[4]," is: ", rate))

      } else if (length(pval) == 5) {

        signif <- ifelse(fcHurdle[,6] < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[4], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[4]," is: ", rate))

        signif <- ifelse(fcHurdle[,6] < pval[5], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[5]," is: ", rate))

      } else {

        message("Too many pvalues, shorten vector of pvalues to 5 or less")

      }

    }


  } else if (method == "GLM_tweedie") {
    if (!requireNamespace(c("glmmTMB"),quietly = TRUE)){
      stop("The 'glmmTMB package is required. Please install it.\n
           It may be a problem with dependencies for these packages,\n
           type '!requireNamespace(\"glmmTMB\")' to see if this is the issue.",
           call. = FALSE)
    } else {

      message("This function is slow and requires a lot of memory.")

      if(!is.numeric(fc_range_min)||fc_range_min<1)
        stop("fc_range_min should be between 1-10")
      if(!is.numeric(fc_range_max)||fc_range_max<1.1)
        stop("fc_range_max should be between 1.1-10")
      if(fc_range_min>=fc_range_max)
        stop("fc_range_max should be greater than fc_range_min")
      if(!is.numeric(deg_fc_interval)||deg_fc_interval<=0||deg_fc_interval>(fc_range_max-fc_range_min))
        stop("deg_fc_interval should be numeric and smaller than the difference between fc_range_max and fc_range_min")
      

      if (cells_per_control < 50 | cells_per_case < 50){
        message("----------------------------------------------")
        message("Mean number of cells per individual is less than 50.\n The probability of complete separation will start to increase.")
        message("----------------------------------------------")
      }

      #simulate the non differentially expressed genes
      
      #if checkpoint saved and not forced, don't rerun
      checkpoint_params <- 
        paste(method,n_genes,n_per_group,n_cases,n_controls,cells_per_case,
              cells_per_control,ncells_variation_type,pval,decrease_dropout,
              alter_dropout_cases,seed,sep="_")
      if(!force_new && !(!is.null(checkpoint) && 
                         file.exists(paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds")))){
        #if user wanted set the seed
        if(!isFALSE(seed))
          set.seed(seed)
        non_deg <- suppressMessages(simulate_hierarchicell(data_summaries,
                                                           n_genes = n_genes,
                                                           n_per_group = n_per_group,
                                                           n_cases = n_cases,
                                                           n_controls = n_controls,
                                                           cells_per_case = cells_per_case,
                                                           cells_per_control = cells_per_control,
                                                           ncells_variation_type = ncells_variation_type,
                                                           foldchange = 1, #fc=1 means not DEG
                                                           decrease_dropout = decrease_dropout,
                                                           alter_dropout_cases = alter_dropout_cases))
  
        genecounts_ndeg <- as.matrix(t(non_deg[,c(-1,-2,-3)]))
        coldata_ndeg <- non_deg[,1:3]
        coldata_ndeg$Status <- as.factor(coldata_ndeg$Status)
        coldata_ndeg$DonorID <- as.factor(coldata_ndeg$DonorID)
        genecounts_ndeg <- genecounts_ndeg[which(apply(genecounts_ndeg, 1, mean) > 5), ] + 1
        genecounts_ndeg <- log(genecounts_ndeg)
        genecounts_ndeg <- t(genecounts_ndeg[,rownames(coldata_ndeg)])
        allcells_ndeg <- cbind(coldata_ndeg,genecounts_ndeg)
  
        fitfixed_ndeg <- lapply(4:ncol(allcells_ndeg),
                           function(x){glmmTMB::glmmTMB(allcells_ndeg[,x] ~ Status,
                                                        data=allcells_ndeg,
                                                        family=glmmTMB::tweedie,
                                                        ziformula= ~0)})
        summaries_ndeg <- lapply(fitfixed_ndeg, summary)
        pvalues_ndeg <- as.numeric(unlist(lapply(summaries_ndeg, function(x){stats::coef(x)$cond[2,4]})))
        
        if(!is.null(checkpoint))
          saveRDS(pvalues_ndeg,
                  file=paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds"))
      }else{
        #load checkpoint data
        pvalues_ndeg <- 
          readRDS(paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds"))
      }
      
      #simulate the differentially expressed genes
      deg_groups <-
        get_deg_groups(n_genes,fc_range_min,fc_range_max,seed,deg_fc_interval)
      #separate call for each
      degs_sep_fcs <- vector(mode="list",length=length(deg_groups))
      names(degs_sep_fcs) <- names(deg_groups)
      #if user wanted set the seed
      if(!isFALSE(seed))
        set.seed(seed)
      #cycle through fc choices
      for (fc_i in as.numeric(names(deg_groups))){
        if(!force_new && !(!is.null(checkpoint) && 
                           file.exists(paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds")))){
          degs <- 
            suppressMessages(simulate_hierarchicell(data_summaries,
                                                    n_genes = deg_groups[[as.character(fc_i)]],
                                                    n_per_group = n_per_group,
                                                    n_cases = n_cases,
                                                    n_controls = n_controls,
                                                    cells_per_case = 
                                                      cells_per_case,
                                                    cells_per_control = 
                                                      cells_per_control,
                                                    ncells_variation_type = 
                                                      ncells_variation_type,
                                                    foldchange = fc_i,
                                                    decrease_dropout = 
                                                      decrease_dropout,
                                                    alter_dropout_cases = 
                                                      alter_dropout_cases))
          genecounts_degs <- as.matrix(t(degs[,c(-1,-2,-3)]))
          coldata_degs <- degs[,1:3]
          coldata_degs$Status <- as.factor(coldata_degs$Status)
          coldata_degs$DonorID <- as.factor(coldata_degs$DonorID)
          genecounts_degs <- genecounts_degs[which(apply(genecounts_degs, 1, mean) > 5), ] + 1
          genecounts_degs <- log(genecounts_degs)
          genecounts_degs <- t(genecounts_degs[,rownames(coldata_degs)])
          allcells_degs <- cbind(coldata_degs,genecounts_degs)
          
          fitfixed_degs <- lapply(4:ncol(allcells_degs),
                                  function(x){glmmTMB::glmmTMB(allcells_degs[,x] ~ Status,
                                                               data=allcells_degs,
                                                               family=glmmTMB::tweedie,
                                                               ziformula= ~0)})
          summaries_degs <- lapply(fitfixed_degs, summary)
          pvalues_degs <- as.numeric(unlist(lapply(summaries_degs, function(x){stats::coef(x)$cond[2,4]})))
          if(!is.null(checkpoint))
            saveRDS(pvalues_degs,
                    file=paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds"))
        }else{
          #load checkpoint data
          pvalues_degs <- 
            readRDS(paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds"))
        }    
        degs_sep_fcs[[as.character(fc_i)]] <- pvalues_degs
      }
      #just combine the results
      pvalues_degs <- unlist(degs_sep_fcs)
      
      
      
      if (length(pval) == 1){
        #create actuals
        act_neg <- rep(0,length(pvalues_ndeg))
        act_pos <- rep(1,length(pvalues_degs))
        #combine pos neg
        pred <- c(pvalues_ndeg,pvalues_degs)
        actuals <- c(act_neg,act_pos)
        
        signif <- ifelse(pvalues_ndeg < pval, 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval," is: ", rate))
        return(calc_mcc(predicted=pred,actual=actuals,cutoff=pval))

      } else if (length(pval) == 2) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

      } else if (length(pval) == 3) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(pvalues < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))


      } else if (length(pval) == 4) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(pvalues < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))

        signif <- ifelse(pvalues < pval[4], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[4]," is: ", rate))

      } else if (length(pval) == 5) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(pvalues < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))

        signif <- ifelse(pvalues < pval[4], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[4]," is: ", rate))

        signif <- ifelse(pvalues < pval[5], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[5]," is: ", rate))

      } else {

        message("Too many pvalues, shorten vector of pvalues to 5 or less")

      }


    }


  } else if (method == "GLMM_tweedie") {
    if (!requireNamespace(c("glmmTMB"),quietly = TRUE)){
      stop("The 'glmmTMB package is required. Please install it.\n
           It may be a problem with dependencies for these packages,\n
           type '!requireNamespace(\"glmmTMB\")' to see if this is the issue.",
           call. = FALSE)
    } else {

      message("This function is slow and requires a lot of memory.")

      if(!is.numeric(fc_range_min)||fc_range_min<1)
        stop("fc_range_min should be between 1-10")
      if(!is.numeric(fc_range_max)||fc_range_max<1.1)
        stop("fc_range_max should be between 1.1-10")
      if(fc_range_min>=fc_range_max)
        stop("fc_range_max should be greater than fc_range_min")
      if(!is.numeric(deg_fc_interval)||deg_fc_interval<=0||deg_fc_interval>(fc_range_max-fc_range_min))
        stop("deg_fc_interval should be numeric and smaller than the difference between fc_range_max and fc_range_min")
      


      if (cells_per_control < 50 | cells_per_case < 50){
        message("----------------------------------------------")
        message("Mean number of cells per individual is less than 50.\n The probability of complete separation will start to increase.")
        message("----------------------------------------------")
      }

      #simulate the non differentially expressed genes
      
      #if checkpoint saved and not forced, don't rerun
      checkpoint_params <- 
        paste(method,n_genes,n_per_group,n_cases,n_controls,cells_per_case,
              cells_per_control,ncells_variation_type,pval,decrease_dropout,
              alter_dropout_cases,seed,sep="_")
      if(!force_new && !(!is.null(checkpoint) && 
                         file.exists(paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds")))){
        
        #if user wanted set the seed
        if(!isFALSE(seed))
          set.seed(seed)
        non_deg <- suppressMessages(simulate_hierarchicell(data_summaries,
                                                           n_genes = n_genes,
                                                           n_per_group = n_per_group,
                                                           n_cases = n_cases,
                                                           n_controls = n_controls,
                                                           cells_per_case = cells_per_case,
                                                           cells_per_control = cells_per_control,
                                                           ncells_variation_type = ncells_variation_type,
                                                           foldchange = 1, #fc=1 means not DEG
                                                           decrease_dropout = decrease_dropout,
                                                           alter_dropout_cases = alter_dropout_cases))
        
        genecounts_ndeg <- as.matrix(t(non_deg[,c(-1,-2,-3)]))
        coldata_ndeg <- non_deg[,1:3]
        coldata_ndeg$Status <- as.factor(coldata_ndeg$Status)
        coldata_ndeg$DonorID <- as.factor(coldata_ndeg$DonorID)
        genecounts_ndeg <- genecounts_ndeg[which(apply(genecounts_ndeg, 1, mean) > 5), ] + 1
        genecounts_ndeg <- log(genecounts_ndeg)
        genecounts_ndeg <- t(genecounts_ndeg[,rownames(coldata_ndeg)])
        allcells_ndeg <- cbind(coldata_ndeg,genecounts_ndeg)
        
        fitmixed_ndeg <- lapply(4:ncol(allcells_ndeg),
                                function(x){glmmTMB::glmmTMB(allcells_ndeg[,x] ~ Status+ (1 | DonorID),
                                                             data=allcells_ndeg,
                                                             family=glmmTMB::tweedie,
                                                             ziformula= ~0)})
        summaries_ndeg <- lapply(fitmixed_ndeg, summary)
        pvalues_ndeg <- as.numeric(unlist(lapply(summaries_ndeg, function(x){stats::coef(x)$cond[2,4]})))
        
        if(!is.null(checkpoint))
          saveRDS(pvalues_ndeg,
                  file=paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds"))
      }else{
        #load checkpoint data
        pvalues_ndeg <- 
          readRDS(paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds"))
      }
      
      #simulate the differentially expressed genes
      deg_groups <-
        get_deg_groups(n_genes,fc_range_min,fc_range_max,seed,deg_fc_interval)
      #separate call for each
      degs_sep_fcs <- vector(mode="list",length=length(deg_groups))
      names(degs_sep_fcs) <- names(deg_groups)
      #if user wanted set the seed
      if(!isFALSE(seed))
        set.seed(seed)
      #cycle through fc choices
      for (fc_i in as.numeric(names(deg_groups))){
        if(!force_new && !(!is.null(checkpoint) && 
                           file.exists(paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds")))){
          degs <- 
            suppressMessages(simulate_hierarchicell(data_summaries,
                                                    n_genes = deg_groups[[as.character(fc_i)]],
                                                    n_per_group = n_per_group,
                                                    n_cases = n_cases,
                                                    n_controls = n_controls,
                                                    cells_per_case = 
                                                      cells_per_case,
                                                    cells_per_control = 
                                                      cells_per_control,
                                                    ncells_variation_type = 
                                                      ncells_variation_type,
                                                    foldchange = fc_i,
                                                    decrease_dropout = 
                                                      decrease_dropout,
                                                    alter_dropout_cases = 
                                                      alter_dropout_cases))
          genecounts_degs <- as.matrix(t(degs[,c(-1,-2,-3)]))
          coldata_degs <- degs[,1:3]
          coldata_degs$Status <- as.factor(coldata_degs$Status)
          coldata_degs$DonorID <- as.factor(coldata_degs$DonorID)
          genecounts_degs <- genecounts_degs[which(apply(genecounts_degs, 1, mean) > 5), ] + 1
          genecounts_degs <- log(genecounts_degs)
          genecounts_degs <- t(genecounts_degs[,rownames(coldata_degs)])
          allcells_degs <- cbind(coldata_degs,genecounts_degs)
          
          fitmixed_degs <- lapply(4:ncol(allcells_degs),
                                  function(x){glmmTMB::glmmTMB(allcells_degs[,x] ~ Status+ (1 | DonorID),
                                                               data=allcells_degs,
                                                               family=glmmTMB::tweedie,
                                                               ziformula= ~0)})
          summaries_degs <- lapply(fitmixed_degs, summary)
          pvalues_degs <- as.numeric(unlist(lapply(summaries_degs, function(x){stats::coef(x)$cond[2,4]})))
          
          if(!is.null(checkpoint))
            saveRDS(pvalues_degs,
                    file=paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds"))
        }else{
          #load checkpoint data
          pvalues_degs <- 
            readRDS(paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds"))
        }  
        degs_sep_fcs[[as.character(fc_i)]] <- pvalues_degs
      }
      #just combine the results
      pvalues_degs <- unlist(degs_sep_fcs)

      if (length(pval) == 1){
        #create actuals
        act_neg <- rep(0,length(pvalues_ndeg))
        act_pos <- rep(1,length(pvalues_degs))
        #combine pos neg
        pred <- c(pvalues_ndeg,pvalues_degs)
        actuals <- c(act_neg,act_pos)
        
        signif <- ifelse(pvalues_ndeg < pval, 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval," is: ", rate))
        return(calc_mcc(predicted=pred,actual=actuals,cutoff=pval))

      } else if (length(pval) == 2) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

      } else if (length(pval) == 3) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(pvalues < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))


      } else if (length(pval) == 4) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(pvalues < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))

        signif <- ifelse(pvalues < pval[4], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[4]," is: ", rate))

      } else if (length(pval) == 5) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(pvalues < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))

        signif <- ifelse(pvalues < pval[4], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[4]," is: ", rate))

        signif <- ifelse(pvalues < pval[5], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[5]," is: ", rate))

      } else {

        message("Too many pvalues, shorten vector of pvalues to 5 or less")

      }

    }


  } else if (method == "GEE1") {
    if (!requireNamespace(c("geepack"),quietly = TRUE)){
      stop("The 'glmmTMB package is required. Please install it.\n
           It may be a problem with dependencies for these packages,\n
           type '!requireNamespace(\"geepack\")' to see if this is the issue.",
           call. = FALSE)
    } else {

      message("This function is slow and requires a lot of memory.")

      if(!is.numeric(fc_range_min)||fc_range_min<1)
        stop("fc_range_min should be between 1-10")
      if(!is.numeric(fc_range_max)||fc_range_max<1.1)
        stop("fc_range_max should be between 1.1-10")
      if(fc_range_min>=fc_range_max)
        stop("fc_range_max should be greater than fc_range_min")
      if(!is.numeric(deg_fc_interval)||deg_fc_interval<=0||deg_fc_interval>(fc_range_max-fc_range_min))
        stop("deg_fc_interval should be numeric and smaller than the difference between fc_range_max and fc_range_min")
      


      if (cells_per_control < 50 | cells_per_case < 50){
        message("----------------------------------------------")
        message("Mean number of cells per individual is less than 50.\n The probability of complete separation will start to increase.")
        message("----------------------------------------------")
      }

      #simulate the non differentially expressed genes
      
      #if checkpoint saved and not forced, don't rerun
      checkpoint_params <- 
        paste(method,n_genes,n_per_group,n_cases,n_controls,cells_per_case,
              cells_per_control,ncells_variation_type,pval,decrease_dropout,
              alter_dropout_cases,seed,sep="_")
      if(!force_new && !(!is.null(checkpoint) && 
                         file.exists(paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds")))){
        #if user wanted set the seed
        if(!isFALSE(seed))
          set.seed(seed)
        non_deg <- suppressMessages(simulate_hierarchicell(data_summaries,
                                                           n_genes = n_genes,
                                                           n_per_group = n_per_group,
                                                           n_cases = n_cases,
                                                           n_controls = n_controls,
                                                           cells_per_case = cells_per_case,
                                                           cells_per_control = cells_per_control,
                                                           ncells_variation_type = ncells_variation_type,
                                                           foldchange = 1, #fc=1 means not DEG
                                                           decrease_dropout = decrease_dropout,
                                                           alter_dropout_cases = alter_dropout_cases))
        
        genecounts_ndeg <- as.matrix(t(non_deg[,c(-1,-2,-3)]))
        coldata_ndeg <- non_deg[,1:3]
        coldata_ndeg$Status <- as.factor(coldata_ndeg$Status)
        coldata_ndeg$DonorID <- as.factor(coldata_ndeg$DonorID)
        genecounts_ndeg <- genecounts_ndeg[which(apply(genecounts_ndeg, 1, mean) > 5), ] + 1
        genecounts_ndeg <- log(genecounts_ndeg)
        genecounts_ndeg <- t(genecounts_ndeg[,rownames(coldata_ndeg)])
        allcells_ndeg <- cbind(coldata_ndeg,genecounts_ndeg)
  
        fitgee_ndeg <- lapply(4:ncol(allcells_ndeg),
                           function(x){geepack::geeglm(allcells_ndeg[,x] ~ Status,
                                                       data=allcells_ndeg,
                                                       family=stats::gaussian(link="identity"),
                                                       id = DonorID,
                                                       corstr="exchangeable")})
        summaries_ndeg <- lapply(fitgee_ndeg, summary)
        pvalues_ndeg <- as.numeric(unlist(lapply(summaries_ndeg, function(x){stats::coef(x)[2,4]})))
        
        if(!is.null(checkpoint))
          saveRDS(pvalues_ndeg,
                  file=paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds"))
      }else{
        #load checkpoint data
        pvalues_ndeg <- 
          readRDS(paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds"))
      }
      
      #simulate the differentially expressed genes
      deg_groups <-
        get_deg_groups(n_genes,fc_range_min,fc_range_max,seed,deg_fc_interval)
      #separate call for each
      degs_sep_fcs <- vector(mode="list",length=length(deg_groups))
      names(degs_sep_fcs) <- names(deg_groups)
      #if user wanted set the seed
      if(!isFALSE(seed))
        set.seed(seed)
      #cycle through fc choices
      for (fc_i in as.numeric(names(deg_groups))){
        if(!force_new && !(!is.null(checkpoint) && 
                           file.exists(paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds")))){
          degs <- 
            suppressMessages(simulate_hierarchicell(data_summaries,
                                                    n_genes = deg_groups[[as.character(fc_i)]],
                                                    n_per_group = n_per_group,
                                                    n_cases = n_cases,
                                                    n_controls = n_controls,
                                                    cells_per_case = 
                                                      cells_per_case,
                                                    cells_per_control = 
                                                      cells_per_control,
                                                    ncells_variation_type = 
                                                      ncells_variation_type,
                                                    foldchange = fc_i,
                                                    decrease_dropout = 
                                                      decrease_dropout,
                                                    alter_dropout_cases = 
                                                      alter_dropout_cases))
          genecounts_degs <- as.matrix(t(degs[,c(-1,-2,-3)]))
          coldata_degs <- degs[,1:3]
          coldata_degs$Status <- as.factor(coldata_degs$Status)
          coldata_degs$DonorID <- as.factor(coldata_degs$DonorID)
          genecounts_degs <- genecounts_degs[which(apply(genecounts_degs, 1, mean) > 5), ] + 1
          genecounts_degs <- log(genecounts_degs)
          genecounts_degs <- t(genecounts_degs[,rownames(coldata_degs)])
          allcells_degs <- cbind(coldata_degs,genecounts_degs)
          
          fitgee_degs <- lapply(4:ncol(allcells_degs),
                                function(x){geepack::geeglm(allcells_degs[,x] ~ Status,
                                                            data=allcells_degs,
                                                            family=stats::gaussian(link="identity"),
                                                            id = DonorID,
                                                            corstr="exchangeable")})
          summaries_degs <- lapply(fitgee_degs, summary)
          pvalues_degs <- as.numeric(unlist(lapply(summaries_degs, function(x){stats::coef(x)[2,4]})))
          
          if(!is.null(checkpoint))
            saveRDS(pvalues_degs,
                    file=paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds"))
        }else{
          #load checkpoint data
          pvalues_degs <- 
            readRDS(paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds"))
        }  
        degs_sep_fcs[[as.character(fc_i)]] <- pvalues_degs
      }
      #just combine the results
      pvalues_degs <- unlist(degs_sep_fcs)
      
      
      if (length(pval) == 1){
        #create actuals
        act_neg <- rep(0,length(pvalues_ndeg))
        act_pos <- rep(1,length(pvalues_degs))
        #combine pos neg
        pred <- c(pvalues_ndeg,pvalues_degs)
        actuals <- c(act_neg,act_pos)
        
        signif <- ifelse(pvalues_ndeg < pval, 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval," is: ", rate))
        return(calc_mcc(predicted=pred,actual=actuals,cutoff=pval))

      } else if (length(pval) == 2) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

      } else if (length(pval) == 3) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(pvalues < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))


      } else if (length(pval) == 4) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(pvalues < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))

        signif <- ifelse(pvalues < pval[4], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[4]," is: ", rate))

      } else if (length(pval) == 5) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(pvalues < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))

        signif <- ifelse(pvalues < pval[4], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[4]," is: ", rate))

        signif <- ifelse(pvalues < pval[5], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[5]," is: ", rate))

      } else {

        message("Too many pvalues, shorten vector of pvalues to 5 or less")

      }


    }
  } else if (method == "ROTS") {

    if (!requireNamespace(c("ROTS","tidyr"),quietly = TRUE)){
      stop("The packages 'ROTS', 'fitdistrplus', and \n
           'tidyr' are required. Please install them.\n
           It may be a problem with dependencies for these packages,\n
           type '!requireNamespace(\"ROTS\")' to see if this is the issue.",
           call. = FALSE)
    } else {

      if(!is.numeric(fc_range_min)||fc_range_min<1)
        stop("fc_range_min should be between 1-10")
      if(!is.numeric(fc_range_max)||fc_range_max<1.1)
        stop("fc_range_max should be between 1.1-10")
      if(fc_range_min>=fc_range_max)
        stop("fc_range_max should be greater than fc_range_min")
      if(!is.numeric(deg_fc_interval)||deg_fc_interval<=0||deg_fc_interval>(fc_range_max-fc_range_min))
        stop("deg_fc_interval should be numeric and smaller than the difference between fc_range_max and fc_range_min")
      


      if (cells_per_control < 50 | cells_per_case < 50){
        message("----------------------------------------------")
        message("Mean number of cells per individual is less than 50.\n The probability of complete separation will start to increase.")
        message("----------------------------------------------")
      }

      #simulate the non differentially expressed genes
      #if checkpoint saved and not forced, don't rerun
      checkpoint_params <- 
        paste(method,n_genes,n_per_group,n_cases,n_controls,cells_per_case,
              cells_per_control,ncells_variation_type,pval,decrease_dropout,
              alter_dropout_cases,seed,sep="_")
      if(!force_new && !(!is.null(checkpoint) && 
                         file.exists(paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds")))){
        #if user wanted set the seed
        if(!isFALSE(seed))
          set.seed(seed)
        non_deg <- suppressMessages(simulate_hierarchicell(data_summaries,
                                                           n_genes = n_genes,
                                                           n_per_group = n_per_group,
                                                           n_cases = n_cases,
                                                           n_controls = n_controls,
                                                           cells_per_case = cells_per_case,
                                                           cells_per_control = cells_per_control,
                                                           ncells_variation_type = ncells_variation_type,
                                                           foldchange = 1, #fc=1 means not DEG
                                                           decrease_dropout = decrease_dropout,
                                                           alter_dropout_cases = alter_dropout_cases))
        
        genecounts_ndeg <- as.matrix(t(non_deg[,c(-1,-2,-3)]))
        coldata_ndeg <- non_deg[,1:3]
        coldata_ndeg$Status <- as.factor(coldata_ndeg$Status)
        coldata_ndeg$DonorID <- as.factor(coldata_ndeg$DonorID)
        genecounts_ndeg <- genecounts_ndeg[which(apply(genecounts_ndeg, 1, mean) > 5), ]
        genecounts_ndeg <- genecounts_ndeg[,rownames(coldata_ndeg)]
        coldata_ndeg$Status <- ifelse(coldata_ndeg$Status == "Control",0,1)
        coldata_ndeg$DonorID <- as.factor(coldata_ndeg$DonorID)
        genecounts_ndeg <- genecounts_ndeg[,rownames(coldata_ndeg)]
        results_ndeg <- suppressMessages(ROTS::ROTS(data = genecounts_ndeg, groups = coldata_ndeg$Status, B = 1000, K = 500))
        results_ndeg <- suppressMessages(ROTS::summary.ROTS(results_ndeg, num.genes=nrow(genecounts_ndeg)))
        results_ndeg <- as.data.frame(results_ndeg)
        results_ndeg <- stats::na.omit(results_ndeg)
        pvalues_ndeg <- as.numeric(results_ndeg$pvalue)
        if(!is.null(checkpoint))
          saveRDS(pvalues_ndeg,
                  file=paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds"))
      }else{
        #load checkpoint data
        pvalues_ndeg <- 
          readRDS(paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds"))
      }
      
      #simulate the differentially expressed genes
      deg_groups <-
        get_deg_groups(n_genes,fc_range_min,fc_range_max,seed,deg_fc_interval)
      #separate call for each
      degs_sep_fcs <- vector(mode="list",length=length(deg_groups))
      names(degs_sep_fcs) <- names(deg_groups)
      #if user wanted set the seed
      if(!isFALSE(seed))
        set.seed(seed)
      #cycle through fc choices
      for (fc_i in as.numeric(names(deg_groups))){
        if(!force_new && !(!is.null(checkpoint) && 
                           file.exists(paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds")))){
          degs <- 
            suppressMessages(simulate_hierarchicell(data_summaries,
                                                    n_genes = deg_groups[[as.character(fc_i)]],
                                                    n_per_group = n_per_group,
                                                    n_cases = n_cases,
                                                    n_controls = n_controls,
                                                    cells_per_case = 
                                                      cells_per_case,
                                                    cells_per_control = 
                                                      cells_per_control,
                                                    ncells_variation_type = 
                                                      ncells_variation_type,
                                                    foldchange = fc_i,
                                                    decrease_dropout = 
                                                      decrease_dropout,
                                                    alter_dropout_cases = 
                                                      alter_dropout_cases))
          genecounts_degs <- as.matrix(t(degs[,c(-1,-2,-3)]))
          coldata_degs <- degs[,1:3]
          coldata_degs$Status <- as.factor(coldata_degs$Status)
          coldata_degs$DonorID <- as.factor(coldata_degs$DonorID)
          genecounts_degs <- genecounts_degs[which(apply(genecounts_degs, 1, mean) > 5), ]
          genecounts_degs <- genecounts_degs[,rownames(coldata_degs)]
          coldata_degs$Status <- ifelse(coldata_degs$Status == "Control",0,1)
          coldata_degs$DonorID <- as.factor(coldata_degs$DonorID)
          genecounts_degs <- genecounts_degs[,rownames(coldata_degs)]
          results_degs <- suppressMessages(ROTS::ROTS(data = genecounts_degs, groups = coldata_degs$Status, B = 1000, K = 500))
          results_degs <- suppressMessages(ROTS::summary.ROTS(results_degs, num.genes=nrow(genecounts_degs)))
          results_degs <- as.data.frame(results_degs)
          results_degs<- stats::na.omit(results_degs)
          pvalues_degs <- as.numeric(results_degs$pvalue)
          if(!is.null(checkpoint))
            saveRDS(pvalues_degs,
                    file=paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds"))
        }else{
          #load checkpoint data
          pvalues_degs <- 
            readRDS(paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds"))
        }  
        degs_sep_fcs[[as.character(fc_i)]] <- pvalues_degs
      }
      #just combine the results
      pvalues_degs <- unlist(degs_sep_fcs)

      if (length(pval) == 1){
        #create actuals
        act_neg <- rep(0,length(pvalues_ndeg))
        act_pos <- rep(1,length(pvalues_degs))
        #combine pos neg
        pred <- c(pvalues_ndeg,pvalues_degs)
        actuals <- c(act_neg,act_pos)
        
        signif <- ifelse(pvalues_ndeg < pval, 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval," is: ", rate))
        return(calc_mcc(predicted=pred,actual=actuals,cutoff=pval))

      } else if (length(pval) == 2) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

      } else if (length(pval) == 3) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(pvalues < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))


      } else if (length(pval) == 4) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(pvalues < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))

        signif <- ifelse(pvalues < pval[4], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[4]," is: ", rate))

      } else if (length(pval) == 5) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(pvalues < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))

        signif <- ifelse(pvalues < pval[4], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[4]," is: ", rate))

        signif <- ifelse(pvalues < pval[5], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[5]," is: ", rate))

      } else {

        message("Too many pvalues, shorten vector of pvalues to 5 or less")

      }

      }


  } else if (method == "Monocle") {
    if (!requireNamespace(c("monocle","tidyr"),quietly = TRUE)){
      stop("The packages 'monocle', 'fitdistrplus', and \n
           'tidyr' are required. Please install them.\n
           It may be a problem with dependencies for these packages,\n
           type '!requireNamespace(\"monocle\")' to see if this is the issue.",
           call. = FALSE)
    } else {

      if(!is.numeric(fc_range_min)||fc_range_min<1)
        stop("fc_range_min should be between 1-10")
      if(!is.numeric(fc_range_max)||fc_range_max<1.1)
        stop("fc_range_max should be between 1.1-10")
      if(fc_range_min>=fc_range_max)
        stop("fc_range_max should be greater than fc_range_min")
      if(!is.numeric(deg_fc_interval)||deg_fc_interval<=0||deg_fc_interval>(fc_range_max-fc_range_min))
        stop("deg_fc_interval should be numeric and smaller than the difference between fc_range_max and fc_range_min")
      


      if (cells_per_control < 50 | cells_per_case < 50){
        message("----------------------------------------------")
        message("Mean number of cells per individual is less than 50.\n The probability of complete separation will start to increase.")
        message("----------------------------------------------")
      }

      #simulate the non differentially expressed genes
      
      #if checkpoint saved and not forced, don't rerun
      checkpoint_params <- 
        paste(method,n_genes,n_per_group,n_cases,n_controls,cells_per_case,
              cells_per_control,ncells_variation_type,pval,decrease_dropout,
              alter_dropout_cases,seed,sep="_")
      if(!force_new && !(!is.null(checkpoint) && 
                         file.exists(paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds")))){
        #if user wanted set the seed
        if(!isFALSE(seed))
          set.seed(seed)
        non_deg <- suppressMessages(simulate_hierarchicell(data_summaries,
                                                           n_genes = n_genes,
                                                           n_per_group = n_per_group,
                                                           n_cases = n_cases,
                                                           n_controls = n_controls,
                                                           cells_per_case = cells_per_case,
                                                           cells_per_control = cells_per_control,
                                                           ncells_variation_type = ncells_variation_type,
                                                           foldchange = 1, #fc=1 means not DEG
                                                           decrease_dropout = decrease_dropout,
                                                           alter_dropout_cases = alter_dropout_cases))
        
        genecounts_ndeg <- as.matrix(t(non_deg[,c(-1,-2,-3)]))
        coldata_ndeg <- non_deg[,1:3]
        coldata_ndeg$Status <- as.factor(coldata_ndeg$Status)
        genecounts_ndeg <- genecounts_ndeg[which(apply(genecounts_ndeg, 1, mean) > 5), ]
        genecounts_ndeg <- genecounts_ndeg[,rownames(coldata_ndeg)]
        coldata_ndeg$Status <- ifelse(coldata_ndeg$Status == "Control",0,1)
        coldata_ndeg$DonorID <- as.factor(coldata_ndeg$DonorID)
        genecounts_ndeg <- genecounts_ndeg[,rownames(coldata_ndeg)]
        features_ndeg <- data.frame(rownames(genecounts_ndeg),"Function")
        colnames(features_ndeg) <- c("gene_short_name","Function")
        rownames(features_ndeg) <- features_ndeg$gene_short_name
        features_ndeg <- methods::new("AnnotatedDataFrame", data = features_ndeg)
        pheno_ndeg <- methods::new("AnnotatedDataFrame", data = coldata_ndeg)
        celldat_ndeg <- monocle::newCellDataSet(genecounts_ndeg,phenoData = pheno_ndeg, featureData = features_ndeg, expressionFamily = VGAM::negbinomial.size())
        celldat_ndeg <- monocle::detectGenes(celldat_ndeg, min_expr = 0.1)
        celldat_ndeg <- BiocGenerics::estimateSizeFactors(celldat_ndeg)
        celldat_ndeg <- BiocGenerics::estimateDispersions(celldat_ndeg)
        results_ndeg <- monocle::differentialGeneTest(celldat_ndeg, fullModelFormulaStr = "~Status")
        results_ndeg <- stats::na.omit(results_ndeg)
        pvalues_ndeg <- as.numeric(results_ndeg$pval)
        
        if(!is.null(checkpoint))
          saveRDS(pvalues_ndeg,
                  file=paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds"))
      }else{
        #load checkpoint data
        pvalues_ndeg <- 
          readRDS(paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds"))
      }
      
      
      #simulate the differentially expressed genes
      deg_groups <-
        get_deg_groups(n_genes,fc_range_min,fc_range_max,seed,deg_fc_interval)
      #separate call for each
      degs_sep_fcs <- vector(mode="list",length=length(deg_groups))
      names(degs_sep_fcs) <- names(deg_groups)
      #if user wanted set the seed
      if(!isFALSE(seed))
        set.seed(seed)
      #cycle through fc choices
      for (fc_i in as.numeric(names(deg_groups))){
      #Monocole somtimes fails at the following call:
      #celldat_degs <- BiocGenerics::estimateDispersions(celldat_degs)
      #There is no fix for this, see this thread:
      #https://github.com/cole-trapnell-lab/monocle-release/issues/38
      #or here:
      #https://github.com/cole-trapnell-lab/monocle-release/issues/5
      #having cells (or genes) with all 0 counts could cause it but that isn't
      #happening here. Workaround is just to abondon run and keep going with 
      #others
        if(!force_new && !(!is.null(checkpoint) && 
                           file.exists(paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds")))){  
          degs <- 
            suppressMessages(simulate_hierarchicell(data_summaries,
                                                    n_genes = deg_groups[[as.character(fc_i)]],
                                                    n_per_group = n_per_group,
                                                    n_cases = n_cases,
                                                    n_controls = n_controls,
                                                    cells_per_case = 
                                                      cells_per_case,
                                                    cells_per_control = 
                                                      cells_per_control,
                                                    ncells_variation_type = 
                                                      ncells_variation_type,
                                                    foldchange = fc_i,
                                                    decrease_dropout = 
                                                      decrease_dropout,
                                                    alter_dropout_cases = 
                                                      alter_dropout_cases))
          genecounts_degs <- as.matrix(t(degs[,c(-1,-2,-3)]))
          coldata_degs <- degs[,1:3]
          coldata_degs$Status <- as.factor(coldata_degs$Status)
          genecounts_degs <- genecounts_degs[which(apply(genecounts_degs, 1, mean) > 5), ]
          genecounts_degs <- genecounts_degs[,rownames(coldata_degs)]
          coldata_degs$Status <- ifelse(coldata_degs$Status == "Control",0,1)
          coldata_degs$DonorID <- as.factor(coldata_degs$DonorID)
          genecounts_degs <- genecounts_degs[,rownames(coldata_degs)]
          features_degs <- data.frame(rownames(genecounts_degs),"Function")
          colnames(features_degs) <- c("gene_short_name","Function")
          rownames(features_degs) <- features_degs$gene_short_name
          features_degs <- methods::new("AnnotatedDataFrame", data = features_degs)
          pheno_degs <- methods::new("AnnotatedDataFrame", data = coldata_degs)
          celldat_degs <- monocle::newCellDataSet(genecounts_degs,phenoData = pheno_degs, featureData = features_degs, expressionFamily = VGAM::negbinomial.size())
          celldat_degs <- monocle::detectGenes(celldat_degs, min_expr = 0.1)
          celldat_degs <- BiocGenerics::estimateSizeFactors(celldat_degs)
          estimateDisp_err<-
            tryCatch(
              celldat_degs <- BiocGenerics::estimateDispersions(celldat_degs),
              error = function(e) e,
              warning = function(w) w)
          if(is(estimateDisp_err, "error"))
            break
          results_degs <- monocle::differentialGeneTest(celldat_degs, fullModelFormulaStr = "~Status")
          results_degs <- stats::na.omit(results_degs)
          pvalues_degs <- as.numeric(results_degs$pval)
          if(!is.null(checkpoint))
            saveRDS(pvalues_degs,
                    file=paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds"))
        }else{
          #load checkpoint data
          pvalues_degs <- 
            readRDS(paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds"))
        }  
        degs_sep_fcs[[as.character(fc_i)]] <- pvalues_degs
      }
      #just combine the results
      pvalues_degs <- unlist(degs_sep_fcs)

      if (length(pval) == 1){
        #create actuals
        act_neg <- rep(0,length(pvalues_ndeg))
        act_pos <- rep(1,length(pvalues_degs))
        #combine pos neg
        pred <- c(pvalues_ndeg,pvalues_degs)
        actuals <- c(act_neg,act_pos)
        
        signif <- ifelse(pvalues_ndeg < pval, 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval," is: ", rate))
        return(calc_mcc(predicted=pred,actual=actuals,cutoff=pval))

      } else if (length(pval) == 2) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

      } else if (length(pval) == 3) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(pvalues < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))


      } else if (length(pval) == 4) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(pvalues < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))

        signif <- ifelse(pvalues < pval[4], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[4]," is: ", rate))

      } else if (length(pval) == 5) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(pvalues < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))

        signif <- ifelse(pvalues < pval[4], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[4]," is: ", rate))

        signif <- ifelse(pvalues < pval[5], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[5]," is: ", rate))

      } else {

        message("Too many pvalues, shorten vector of pvalues to 5 or less")

      }

      }

  }  else if (method == "Pseudobulk_mean") {
    if (!requireNamespace(c("DESeq2","tidyr"),quietly = TRUE)){
      stop("The packages 'DESeq2', 'fitdistrplus', and \n
           'tidyr' are required. Please install them.\n
           It may be a problem with dependencies for these packages,\n
           type '!requireNamespace(\"DESeq2\")' to see if this is the issue.",
           call. = FALSE)
    } else {

      if(!is.numeric(fc_range_min)||fc_range_min<1)
        stop("fc_range_min should be between 1-10")
      if(!is.numeric(fc_range_max)||fc_range_max<1.1)
        stop("fc_range_max should be between 1.1-10")
      if(fc_range_min>=fc_range_max)
        stop("fc_range_max should be greater than fc_range_min")
      if(!is.numeric(deg_fc_interval)||deg_fc_interval<=0||deg_fc_interval>(fc_range_max-fc_range_min))
        stop("deg_fc_interval should be numeric and smaller than the difference between fc_range_max and fc_range_min")
      

      #simulate the non differentially expressed genes
      
      #if checkpoint saved and not forced, don't rerun
      checkpoint_params <- 
        paste(method,n_genes,n_per_group,n_cases,n_controls,cells_per_case,
              cells_per_control,ncells_variation_type,pval,decrease_dropout,
              alter_dropout_cases,seed,sep="_")
      if(!force_new && !(!is.null(checkpoint) && 
                         file.exists(paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds")))){
        #if user wanted set the seed
        if(!isFALSE(seed))
          set.seed(seed)
        non_deg <- suppressMessages(simulate_hierarchicell(data_summaries,
                                                           n_genes = n_genes,
                                                           n_per_group = n_per_group,
                                                           n_cases = n_cases,
                                                           n_controls = n_controls,
                                                           cells_per_case = cells_per_case,
                                                           cells_per_control = cells_per_control,
                                                           ncells_variation_type = ncells_variation_type,
                                                           foldchange = 1, #fc=1 means not DEG
                                                           decrease_dropout = decrease_dropout,
                                                           alter_dropout_cases = alter_dropout_cases))
        
        genecounts_ndeg <- as.matrix(non_deg[,c(-1,-2,-3)])
        genecounts_ndeg <- genecounts_ndeg[ ,which(apply(genecounts_ndeg, 2, mean) > 5)]
        genecounts_ndeg <- cbind(non_deg[,1:2],genecounts_ndeg)
        computemeans_ndeg <- function(x){tapply(x,genecounts_ndeg[,2],mean)}
        cellmeans_ndeg <- sapply(genecounts_ndeg[,c(-1,-2)],computemeans_ndeg)
        rownames(cellmeans_ndeg) <- paste0(rownames(cellmeans_ndeg),"_Mean")
        coldata_ndeg <- as.data.frame(cbind(rownames(cellmeans_ndeg),rownames(cellmeans_ndeg)))
        colnames(coldata_ndeg) <- c("SampleID","ToSep")
        coldata_ndeg <- tidyr::separate(coldata_ndeg,ToSep,c("Status", "Donor_Number", "Mean"), sep="_")
        rownames(coldata_ndeg) <- coldata_ndeg$SampleID
        coldata_ndeg$Status <- as.factor(coldata_ndeg$Status)
        coldata_ndeg$Status <- stats::relevel(coldata_ndeg$Status, "Control")
        cellmeans_ndeg <- round(t(cellmeans_ndeg),0)
        cellmeans_ndeg <- cellmeans_ndeg[, rownames(coldata_ndeg)]
        
        dsd_ndeg <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(countData = cellmeans_ndeg, colData = coldata_ndeg, design = ~ Status))
        normFactors_ndeg <- matrix(rep(1,(nrow(dsd_ndeg)*ncol(dsd_ndeg))),ncol=ncol(dsd_ndeg),nrow=nrow(dsd_ndeg),dimnames=list(1:nrow(dsd_ndeg),1:ncol(dsd_ndeg)))
        normFactors_ndeg <- normFactors_ndeg / exp(rowMeans(log(normFactors_ndeg)))
        DESeq2::normalizationFactors(dsd_ndeg) <- normFactors_ndeg
        dsd_ndeg <- suppressMessages(DESeq2::DESeq(dsd_ndeg))
        res_ndeg <- as.data.frame(DESeq2::results(dsd_ndeg))
        pvalues_ndeg <- as.numeric(res_ndeg$pvalue)
        
        if(!is.null(checkpoint))
          saveRDS(pvalues_ndeg,
                  file=paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds"))
      }else{
        #load checkpoint data
        pvalues_ndeg <- 
          readRDS(paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds"))
      }
      
      #simulate the differentially expressed genes
      deg_groups <-
        get_deg_groups(n_genes,fc_range_min,fc_range_max,seed,deg_fc_interval)
      #separate call for each
      degs_sep_fcs <- vector(mode="list",length=length(deg_groups))
      names(degs_sep_fcs) <- names(deg_groups)
      #if user wanted set the seed
      if(!isFALSE(seed))
        set.seed(seed)
      #cycle through fc choices
      for (fc_i in as.numeric(names(deg_groups))){
        if(!force_new && !(!is.null(checkpoint) && 
                           file.exists(paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds")))){
          degs <- 
            suppressMessages(simulate_hierarchicell(data_summaries,
                                                    n_genes = deg_groups[[as.character(fc_i)]],
                                                    n_per_group = n_per_group,
                                                    n_cases = n_cases,
                                                    n_controls = n_controls,
                                                    cells_per_case = 
                                                      cells_per_case,
                                                    cells_per_control = 
                                                      cells_per_control,
                                                    ncells_variation_type = 
                                                      ncells_variation_type,
                                                    foldchange = fc_i,
                                                    decrease_dropout = 
                                                      decrease_dropout,
                                                    alter_dropout_cases = 
                                                      alter_dropout_cases))
          genecounts_degs <- as.matrix(degs[,c(-1,-2,-3)])
          genecounts_degs <- genecounts_degs[ ,which(apply(genecounts_degs, 2, mean) > 5)]
          genecounts_degs <- cbind(degs[,1:2],genecounts_degs)
          computemeans_degs <- function(x){tapply(x,genecounts_degs[,2],mean)}
          cellmeans_degs <- sapply(genecounts_degs[,c(-1,-2)],computemeans_degs)
          rownames(cellmeans_degs) <- paste0(rownames(cellmeans_degs),"_Mean")
          coldata_degs <- as.data.frame(cbind(rownames(cellmeans_degs),rownames(cellmeans_degs)))
          colnames(coldata_degs) <- c("SampleID","ToSep")
          coldata_degs <- tidyr::separate(coldata_degs,ToSep,c("Status", "Donor_Number", "Mean"), sep="_")
          rownames(coldata_degs) <- coldata_degs$SampleID
          coldata_degs$Status <- as.factor(coldata_degs$Status)
          coldata_degs$Status <- stats::relevel(coldata_ndeg$Status, "Control")
          cellmeans_degs <- round(t(cellmeans_degs),0)
          cellmeans_degs <- cellmeans_degs[, rownames(coldata_degs)]
          
          dsd_degs <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(countData = cellmeans_degs, colData = coldata_degs, design = ~ Status))
          normFactors_degs <- matrix(rep(1,(nrow(dsd_degs)*ncol(dsd_degs))),ncol=ncol(dsd_degs),nrow=nrow(dsd_degs),dimnames=list(1:nrow(dsd_degs),1:ncol(dsd_degs)))
          normFactors_degs <- normFactors_degs / exp(rowMeans(log(normFactors_degs)))
          DESeq2::normalizationFactors(dsd_degs) <- normFactors_degs
          dsd_degs <- suppressMessages(DESeq2::DESeq(dsd_degs))
          res_degs <- as.data.frame(DESeq2::results(dsd_degs))
          pvalues_degs <- as.numeric(res_degs$pvalue)
          if(!is.null(checkpoint))
            saveRDS(pvalues_degs,
                    file=paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds"))
        }else{
          #load checkpoint data
          pvalues_degs <- 
            readRDS(paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds"))
        }  
        degs_sep_fcs[[as.character(fc_i)]] <- pvalues_degs
      }
      #just combine the results
      pvalues_degs <- unlist(degs_sep_fcs)
      
      
      if (length(pval) == 1){
        #create actuals
        act_neg <- rep(0,length(pvalues_ndeg))
        act_pos <- rep(1,length(pvalues_degs))
        #combine pos neg
        pred <- c(pvalues_ndeg,pvalues_degs)
        actuals <- c(act_neg,act_pos)
        
        signif <- ifelse(pvalues_ndeg < pval, 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval," is: ", rate))
        return(calc_mcc(predicted=pred,actual=actuals,cutoff=pval))

      } else if (length(pval) == 2) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

      } else if (length(pval) == 3) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(pvalues < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))


      } else if (length(pval) == 4) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(pvalues < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))

        signif <- ifelse(pvalues < pval[4], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[4]," is: ", rate))

      } else if (length(pval) == 5) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(pvalues < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))

        signif <- ifelse(pvalues < pval[4], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[4]," is: ", rate))

        signif <- ifelse(pvalues < pval[5], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[5]," is: ", rate))

      } else {

        message("Too many pvalues, shorten vector of pvalues to 5 or less")

      }

    }

  }else if (method == "Pseudobulk_sum") {
    if (!requireNamespace(c("DESeq2","tidyr"),quietly = TRUE)){
      stop("The packages 'DESeq2', 'fitdistrplus', and \n
           'tidyr' are required. Please install them.\n
           It may be a problem with dependencies for these packages,\n
           type '!requireNamespace(\"DESeq2\")' to see if this is the issue.",
           call. = FALSE)
    } else {

      if(!is.numeric(fc_range_min)||fc_range_min<1)
        stop("fc_range_min should be between 1-10")
      if(!is.numeric(fc_range_max)||fc_range_max<1.1)
        stop("fc_range_max should be between 1.1-10")
      if(fc_range_min>=fc_range_max)
        stop("fc_range_max should be greater than fc_range_min")
      if(!is.numeric(deg_fc_interval)||deg_fc_interval<=0||deg_fc_interval>(fc_range_max-fc_range_min))
        stop("deg_fc_interval should be numeric and smaller than the difference between fc_range_max and fc_range_min")
      


      #simulate the non differentially expressed genes
      #if checkpoint saved and not forced, don't rerun
      checkpoint_params <- 
        paste(method,n_genes,n_per_group,n_cases,n_controls,cells_per_case,
              cells_per_control,ncells_variation_type,pval,decrease_dropout,
              alter_dropout_cases,seed,sep="_")
      if(!force_new && !(!is.null(checkpoint) && 
                         file.exists(paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds")))){
        #if user wanted set the seed
        if(!isFALSE(seed))
          set.seed(seed)
        non_deg <- suppressMessages(simulate_hierarchicell(data_summaries,
                                                           n_genes = n_genes,
                                                           n_per_group = n_per_group,
                                                           n_cases = n_cases,
                                                           n_controls = n_controls,
                                                           cells_per_case = cells_per_case,
                                                           cells_per_control = cells_per_control,
                                                           ncells_variation_type = ncells_variation_type,
                                                           foldchange = 1, #fc=1 means not DEG
                                                           decrease_dropout = decrease_dropout,
                                                           alter_dropout_cases = alter_dropout_cases))
        
        genecounts_ndeg <- as.matrix(non_deg[,c(-1,-2,-3)])
        genecounts_ndeg <- genecounts_ndeg[ ,which(apply(genecounts_ndeg, 2, mean) > 5)]
        genecounts_ndeg <- cbind(non_deg[,1:2],genecounts_ndeg)
        computesums_ndeg <- function(x){tapply(x,genecounts_ndeg[,2],sum)}
        cellsums_ndeg <- sapply(genecounts_ndeg[,c(-1,-2)],computesums_ndeg)
        rownames(cellsums_ndeg) <- paste0(rownames(cellsums_ndeg),"_sum")
        coldata_ndeg <- as.data.frame(cbind(rownames(cellsums_ndeg),rownames(cellsums_ndeg)))
        colnames(coldata_ndeg) <- c("SampleID","ToSep")
        coldata_ndeg <- tidyr::separate(coldata_ndeg,ToSep,c("Status", "Donor_Number", "sum"), sep="_")
        rownames(coldata_ndeg) <- coldata_ndeg$SampleID
        coldata_ndeg$Status <- as.factor(coldata_ndeg$Status)
        coldata_ndeg$Status <- stats::relevel(coldata_ndeg$Status, "Control")
        cellsums_ndeg <- round(t(cellsums_ndeg),0)
        cellsums_ndeg <- cellsums_ndeg[, rownames(coldata_ndeg)]
  
        
        dsd_ndeg <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(countData = cellsums_ndeg, colData = coldata_ndeg, design = ~ Status))
        normFactors_ndeg <- matrix(rep(1,(nrow(dsd_ndeg)*ncol(dsd_ndeg))),ncol=ncol(dsd_ndeg),nrow=nrow(dsd_ndeg),dimnames=list(1:nrow(dsd_ndeg),1:ncol(dsd_ndeg)))
        normFactors_ndeg <- normFactors_ndeg / exp(rowMeans(log(normFactors_ndeg)))
        DESeq2::normalizationFactors(dsd_ndeg) <- normFactors_ndeg
        dsd_ndeg <- suppressMessages(DESeq2::DESeq(dsd_ndeg))
        res_ndeg <- as.data.frame(DESeq2::results(dsd_ndeg))
        pvalues_ndeg <- as.numeric(res_ndeg$pvalue)
        if(!is.null(checkpoint))
          saveRDS(pvalues_ndeg,
                  file=paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds"))
      }else{
        #load checkpoint data
        pvalues_ndeg <- 
          readRDS(paste0(checkpoint,"/","nDEG_",checkpoint_params,".rds"))
      }
      
      
      #simulate the differentially expressed genes
      deg_groups <-
        get_deg_groups(n_genes,fc_range_min,fc_range_max,seed,deg_fc_interval)
      #separate call for each
      degs_sep_fcs <- vector(mode="list",length=length(deg_groups))
      names(degs_sep_fcs) <- names(deg_groups)
      #if user wanted set the seed
      if(!isFALSE(seed))
        set.seed(seed)
      #cycle through fc choices
      for (fc_i in as.numeric(names(deg_groups))){
        if(!force_new && !(!is.null(checkpoint) && 
                           file.exists(paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds")))){
          degs <- 
            suppressMessages(simulate_hierarchicell(data_summaries,
                                                    n_genes = deg_groups[[as.character(fc_i)]],
                                                    n_per_group = n_per_group,
                                                    n_cases = n_cases,
                                                    n_controls = n_controls,
                                                    cells_per_case = 
                                                      cells_per_case,
                                                    cells_per_control = 
                                                      cells_per_control,
                                                    ncells_variation_type = 
                                                      ncells_variation_type,
                                                    foldchange = fc_i,
                                                    decrease_dropout = 
                                                      decrease_dropout,
                                                    alter_dropout_cases = 
                                                      alter_dropout_cases))
          genecounts_degs <- as.matrix(degs[,c(-1,-2,-3)])
          genecounts_degs <- genecounts_degs[ ,which(apply(genecounts_degs, 2, mean) > 5)]
          genecounts_degs <- cbind(degs[,1:2],genecounts_degs)
          computesums_degs <- function(x){tapply(x,genecounts_degs[,2],sum)}
          cellsums_degs <- sapply(genecounts_degs[,c(-1,-2)],computesums_degs)
          rownames(cellsums_degs) <- paste0(rownames(cellsums_degs),"_sum")
          coldata_degs <- as.data.frame(cbind(rownames(cellsums_degs),rownames(cellsums_degs)))
          colnames(coldata_degs) <- c("SampleID","ToSep")
          coldata_degs <- tidyr::separate(coldata_degs,ToSep,c("Status", "Donor_Number", "sum"), sep="_")
          rownames(coldata_degs) <- coldata_degs$SampleID
          coldata_degs$Status <- as.factor(coldata_degs$Status)
          coldata_degs$Status <- stats::relevel(coldata_ndeg$Status, "Control")
          cellsums_degs <- round(t(cellsums_degs),0)
          cellsums_degs <- cellsums_degs[, rownames(coldata_degs)]
          
          dsd_degs <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(countData = cellsums_degs, colData = coldata_degs, design = ~ Status))
          normFactors_degs <- matrix(rep(1,(nrow(dsd_degs)*ncol(dsd_degs))),ncol=ncol(dsd_degs),nrow=nrow(dsd_degs),dimnames=list(1:nrow(dsd_degs),1:ncol(dsd_degs)))
          normFactors_degs <- normFactors_degs / exp(rowMeans(log(normFactors_degs)))
          DESeq2::normalizationFactors(dsd_degs) <- normFactors_degs
          dsd_degs <- suppressMessages(DESeq2::DESeq(dsd_degs))
          res_degs <- as.data.frame(DESeq2::results(dsd_degs))
          pvalues_degs <- as.numeric(res_degs$pvalue)
          if(!is.null(checkpoint))
            saveRDS(pvalues_degs,
                    file=paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds"))
        }else{
          #load checkpoint data
          pvalues_degs <- 
            readRDS(paste0(checkpoint,"/","DEG_",fc_i,"_",checkpoint_params,".rds"))
        }  
        degs_sep_fcs[[as.character(fc_i)]] <- pvalues_degs
      }
      #just combine the results
      pvalues_degs <- unlist(degs_sep_fcs)
      

      if (length(pval) == 1){
        #create actuals
        act_neg <- rep(0,length(pvalues_ndeg))
        act_pos <- rep(1,length(pvalues_degs))
        #combine pos neg
        pred <- c(pvalues_ndeg,pvalues_degs)
        actuals <- c(act_neg,act_pos)
        
        signif <- ifelse(pvalues_ndeg < pval, 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval," is: ", rate))
        return(calc_mcc(predicted=pred,actual=actuals,cutoff=pval))

      } else if (length(pval) == 2) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

      } else if (length(pval) == 3) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(pvalues < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))


      } else if (length(pval) == 4) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(pvalues < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))

        signif <- ifelse(pvalues < pval[4], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[4]," is: ", rate))

      } else if (length(pval) == 5) {

        signif <- ifelse(pvalues < pval[1], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[1]," is: ", rate))

        signif <- ifelse(pvalues < pval[2], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[2]," is: ", rate))

        signif <- ifelse(pvalues < pval[3], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[3]," is: ", rate))

        signif <- ifelse(pvalues < pval[4], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[4]," is: ", rate))

        signif <- ifelse(pvalues < pval[5], 1, 0)
        rate <- mean(signif)
        message(paste0("Type 1 error for ",pval[5]," is: ", rate))

      } else {

        message("Too many pvalues, shorten vector of pvalues to 5 or less")

      }

    }

  }

}
