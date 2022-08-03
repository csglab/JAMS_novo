suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggseqlogo))
suppressPackageStartupMessages(library(ggforce))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(reticulate))
options( error = traceback, nwarnings = 10000 )
set.seed(5)

## https://rstudio.github.io/reticulate/articles/python_packages.html
## https://rstudio.github.io/reticulate/reference/use_python.html
## https://rstudio.github.io/reticulate/articles/package.html
## https://rstudio.github.io/reticulate/

## install python packages
# py_install("numba")
# virtualenv_create("r-reticulate")
# virtualenv_install("r-reticulate", "numba", version = "0.55.2" )
# use_virtualenv("r-reticulate")
# Used for testing with Rstudio, normally commented
# setwd("/home/ahcorcha/repos/tools/JAMS_novo") # ahcorcha
# ./data/CTCF_demo/02_formatted_data/smallest_demo
# /home/ahcorcha/repos/tools/JAMS_novo/data/CTCF_demo/05_motif_discovery/IN/NRF1_HUMAN_GM12878_ENCFF910PTH
########################################################   IN and load data ####
option_list = list(
  make_option(c("-e", "--experiment"), type="character",
              default="test",
              help="Experiment ID", metavar="character"),

  make_option(c("-f", "--flanking"), type="integer", default=20,
              help="length of flanking sequence around the motif", 
              metavar="character"),
  
  make_option(c("-l", "--pfm_length"), type="integer", default=8,
              help="", metavar="character"),
  
  make_option(c("-d", "--input_dir"), type="character", metavar="character",
              default="/home/ahcorcha/repos/tools/JAMS_novo/data/CTCF_demo/05_motif_discovery/IN/NRF1_HUMAN_GM12878_ENCFF910PTH",
              help="Input directory with PFM, methylation counts etc ..."),

  make_option(c("-i", "--max_iterations"), type="character", metavar="integer",
              default=100,
              help="Input directory with PFM, methylation counts etc ..."),  
    
  make_option(c("-o", "--output_dir"), type="character",
              default="./data/CTCF_demo/05_motif_discovery/runs",
              help="", metavar="character"),
  
  make_option(c("-x", "--shifting_pos"), type="integer",
              default=2,
              help="", metavar="character"),  
  
  make_option(c("-z", "--inf_pct"), type="double",
              default=0.4,
              help="threshold percentage of seq coeffs for shifting", 
              metavar="character"),  
  
  make_option(c("-m", "--exclude_meth"), type="logical",
              action = "store", default = "FALSE",
              help="", metavar="character"),
  
  make_option(c("-s", "--script_path"), type="character",
              default="/home/ahcorcha/repos/tools/JAMS_novo/src",
              help="", metavar="character")
  );


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser); rm(option_list, opt_parser)

opt$exclude_meth <- as.logical( toupper(opt$exclude_meth) )

## ahcorcha
# opt$exclude_meth <- TRUE

source( paste0( opt$script_path, "/Methyl_ChIP_ftns.R") )
## Source_python is required for source de_novo_discovery_ftns.R
source_python(paste0( opt$script_path, "/motif_discovery.py" ) )
source( paste0( opt$script_path, "/de_novo_discovery_ftns.R" ) )

experiment <- paste0( "de_novo_motif_finding_", opt$experiment, 
                      "_flanking_", opt$flanking )

opt$output_dir <- paste0( opt$output_dir, "/de_novo_motif_", experiment )
dir.create( opt$output_dir, showWarnings = FALSE )
prefix <- paste0( opt$output_dir, "/", experiment )

sink( paste0( prefix, "_log.txt" ) )
cat(paste0( experiment, "\n"))
cat(paste0( "Start wall time: ", Sys.time(), "\n"))

if ( opt$exclude_meth ) {
  cat( paste0( "Training TF models without intra-motif methylation, ",
               "these coefficients will be added to *coefficients_with_FDR.txt ", 
               "for conformity with other functions\n" ) )
}


cat( paste0( "\n\nLoading data (motif length = ", opt$pfm_length, ") ...\n" ) )
dat_all <- load_dat( opt$input_dir, pfm_length = opt$pfm_length )
num_peaks <- nrow( dat_all$x.A.all )

## Random (but constant from iteration to iteration) peaks for motif heatmap 
rnd_num <- sort( sample.int( num_peaks, 7500 ) )

pos_limits <- c( (opt$flanking+1), 
               ( ncol(dat_all$x.C.all)-opt$flanking-opt$pfm_length ) )


###########################################################   Pre iteration ####
## The starting position here makes the center of the motif is at the center of the peaks

## Start at the peak's center, intuitive
# start_pos <- rep_len(x = 101 , # Start at peak's middle
#                      length.out = nrow(dat_all$x.Met.all)) # Number of ChIP-seq peaks

## Start at random positions, for testing
start_pos <- floor( runif( num_peaks,
                           min = pos_limits[1],
                           max = pos_limits[2] ) )


start_pos_list <-list( start_pos )
pfm_length_change_flag <- FALSE
prev_mean_abs_pos_change <- 1000
shift_pos <- 10
iteration <- 1

## The furthest right we can go, position at the start of the left flank
possible_position <- ncol( dat_all$x.C.all ) - 2 * opt$flanking - opt$pfm_length
## Pre compute position specific predictors 
### Create list of length upper_limit_pos, each entry is the dat_all for the 
### correspondent position.
cat( paste0( "Pre-calc. predictors per position (motif length = ", 
             opt$pfm_length, ") ... \n" ) )

predictors_list <- pre_calc_by_pos_dat( this_dat_all = dat_all, 
                                        possible_position = possible_position, 
                                        flanking = opt$flanking, 
                                        pfm_length = opt$pfm_length )

cat(paste0( "Start iterations wall time: ", Sys.time(), "\n"))

while( shift_pos != 0 ){
  
    if ( as.integer(iteration) > as.integer(opt$max_iterations) ){ 
      cat( "Max number of iterations reached\n")
    break }
  
  ## If the motif length had to be expanded 
  ##  (can't be expanded on the first iteration)
  if ( pfm_length_change_flag ){
    
    
    if ( shift_pos == "expand" ){
      cat( paste0( "\n\nExpanding motif at both ends by: ", 
                   opt$shifting_pos, "\n" ) )
      start_pos <- start_pos - opt$shifting_pos
      opt$pfm_length <- opt$pfm_length + opt$shifting_pos * 2
      
    } else{
      cat( paste0( "Expand motif to: ", shift_pos, "\n" ) )
      
      if( shift_pos < 0 ){
        start_pos <- start_pos - opt$shifting_pos
        opt$pfm_length <- opt$pfm_length + opt$shifting_pos }
      
      if( shift_pos > 0 ){
        opt$pfm_length <- opt$pfm_length + opt$shifting_pos }
    }
    
    
    cat(paste0( "Loading data (motif length = ", opt$pfm_length, ") ...\n" ) )
    
    rm( dat_all, predictors_list ); gc()
    dat_all <- load_dat( opt$input_dir, pfm_length = opt$pfm_length )
    
    pos_limits <- c( (opt$flanking+1), 
                   ( ncol(dat_all$x.C.all)-opt$flanking-opt$pfm_length ) )
    
    ## The furthest right we can go, position at the start of the left flank
    possible_position <- ncol( dat_all$x.C.all ) - 2 * opt$flanking - opt$pfm_length
    
    
    start_pos[ start_pos < pos_limits[1] ] <- pos_limits[1]
    start_pos[ start_pos > pos_limits[2] ] <- pos_limits[2]
    
    ## Pre compute position specific predictors 
    ### Create list of length upper_limit_pos, each entry is the dat_all for the 
    ### correspondent position.
    cat( paste0( "Pre-calc. predictors per position (motif length = ", 
                 opt$pfm_length, ") ... \n" ) )
  
    predictors_list <- pre_calc_by_pos_dat( this_dat_all = dat_all, 
                                            possible_position = possible_position, 
                                            flanking = opt$flanking, 
                                            pfm_length = opt$pfm_length )
    pfm_length_change_flag <- FALSE
  }
  
  
  ###############################################################   Iteration ####
  # opt$max_iterations <- 20
  for (i in iteration:as.integer(opt$max_iterations)) {
    iteration <- iteration + 1
    ### Starting iteration
    # i <- 1
    cat( paste0( "Iteration: ", i, "\n" ) )
    prefix_iteration <- paste0( prefix, "_iteration_", format_iteration(i) )
    
    #############################################################   Train GLM ####
    cat( "Training GLM ...\n" )
    this_glm <- train_GLM_at_shifted_pos( flanking = opt$flanking,
                                          pfm_length = opt$pfm_length,
                                          dat_all = dat_all,
                                          start_pos = start_pos,
                                          exclude_meth = opt$exclude_meth )

    ############### Evaluate every position within +/- 200 bps of peak center ####
    cat("Evaluate every position within +/- 200 bps of peak center ...\n")
    pdwn_coeffs <- as.data.frame( coefficients( summary( this_glm ) ) )
    
    ## Get the pulldown coefficients
    pdwn_coeffs <- pdwn_coeffs[ grepl(pattern = ":t$",
                                x = rownames(pdwn_coeffs)),]
    
    ## Get the name of the variables included in th model
    X_var_names <- gsub( ":t", "", rownames( pdwn_coeffs ) )
    
    ## Make sure those variables are in the correct order for the matrix*vector
    ## Mult. the PULLDOWN coeffs vector and the variable matrix for each position
    c_pdwn_predicted <- sapply( X = predictors_list, 
                                FUN = eval_coeffs, 
                                pdn_coeff = pdwn_coeffs$Estimate, 
                                X_names = X_var_names )
      
    rownames(c_pdwn_predicted) <- rownames(predictors_list[[1]])
    c_pdwn_predicted <- as.data.frame( c_pdwn_predicted )
    
    ## Per row, get column name of max value (max TF signal).
    new_start_pos <- colnames(c_pdwn_predicted)[max.col(c_pdwn_predicted, 
                                                        ties.method = "first")]
  
    new_start_pos <- as.numeric( gsub( "V", "", new_start_pos ) )
    new_start_pos <- new_start_pos + opt$flanking ## ahcorcha
    
    ## Change start_pos to the ones with max TF signal
    len <- length(start_pos_list)
    start_pos_list[[len+1]] <- new_start_pos
    
    # compare new_start_pos with start_pos
    mean_abs_pos_change <- mean( abs( start_pos - new_start_pos ) )
    
    cat( paste0("   Mean absolute position change: ", round(mean_abs_pos_change, 4), "\n") )
    
    median_abs_pos_change <- median( abs( start_pos - new_start_pos ) )
    cat( paste0(" Median absolute position change: ", round(median_abs_pos_change, 4), "\n") )
    
    abs_change_pos <- ( start_pos - new_start_pos )
    abs_change_pos_df <- as.data.frame( abs_change_pos )
    
    phist <- ggplot( abs_change_pos_df, aes( x = abs_change_pos ) ) + 
                     geom_histogram( aes(y=..count..), binwidth = 2, 
                                     colour="black", fill="white") +
                     xlim( -( ncol(dat_all$x.C.all)-2*opt$flanking), 
                            ( ncol(dat_all$x.C.all)-2*opt$flanking) ) +
                     labs(x = "Position change", y = "Density") +
                     theme_light()
    
    start_pos <- new_start_pos 
    
    #############################   Visualize motif start pos over iterations ####
    ###### Save run's information: write coefficients / draw logo and DNA coeffs
    cat( "Visualization ...\n" )
    dna_acc_plot_name <- paste0( prefix_iteration, "_dna_coefficients.pdf" )
  
    p_dna_coeffs <- plot_dna_acc_coefficients( this_glm )

    motif_coefs <- write.sequence.model.av.met( seq_fit = this_glm, opt$exclude_meth, opt$pfm_length )
  
    write.table( x =  motif_coefs[[1]], quote = FALSE, sep = "\t",
                 col.names = TRUE, row.names = TRUE,
                 file = paste0( prefix_iteration, "_coefficients_with_FDR.txt") )
  
    p_motif_coefs <- ( motif_coefs[[2]] / p_dna_coeffs )
    
    p_motif_coefs <- p_motif_coefs +
               plot_annotation( title = paste0(experiment, ", iteration: ", i ) ) +
               plot_layout( heights = c(2, 0.75) )
    
    sample_start_pos <- start_pos[rnd_num]
    
    motif_ht <- motif_pos_heatmap( sample_start_pos, n_cols = ncol(dat_all$x.C.all), opt$pfm_length, i )
    
    p_motif_ht <- grid.grabExpr( draw(motif_ht, heatmap_legend_side = "bottom") )
    
    layout <- "AADD
               AADD
               BBDD
               CCDD"
    
    this.tittle <- paste0( experiment, "\n",
                           "flanking = ", opt$flanking,
                           ", motif length = ", opt$pfm_length,
                           ", iteration = ", i,
                           ", n peaks = ", nrow( dat_all$x.C.all ) )
    
    if(opt$exclude_meth) { 
      this.tittle <- paste0(this.tittle, ", TF model without intra-motif methylation") }
    
    p1 <- p_motif_coefs + phist +  p_motif_ht  +
          plot_layout( design = layout ) +
          plot_annotation(title = this.tittle )
  
    ggsave( filename = paste0( prefix_iteration, "_logo_acc_coeffs_motif_ht.pdf" ),
            p1, height = 10, width = 14 )
    
    ggsave( filename = paste0( prefix_iteration, "_logo_acc_coeffs_motif_ht.png" ),
            p1, height = 10, width = 14 )
    
    ############################################ Conditions to end iterations ####
    delta_mean_pos_change <- abs( mean_abs_pos_change - prev_mean_abs_pos_change )
    prev_mean_abs_pos_change <- mean_abs_pos_change
    
    if ( delta_mean_pos_change < 0.1 ){
      
      cat( paste0( "Expand threshold percentage: ", opt$inf_pct, "\n" ))
      
      ## returns positions to be shifted +/- 1,2,3,4,5 or expand
      shift_pos <- pos_to_wiggle( pdwn_coeffs, opt$shifting_pos, opt$inf_pct )
      
      if( shift_pos == 0 ){ 
        cat("Low information content in each end\n")
        break }
      
      if ( shift_pos != 0 ){
        cat( "Expand to one side ...\n" )
        pfm_length_change_flag <- TRUE
        break }
      
    }
  }
}

######################## Format and write start positions across iterations ####
cat( "Format and write start positions across iterations ...\n" )
start_pos_list_df <- as.data.frame( start_pos_list )
colnames(start_pos_list_df) <- paste0( "iteration_", 0:(ncol(start_pos_list_df)-1) )
rownames(start_pos_list_df) <- rownames(predictors_list[[1]])
start_pos_list_df$peak_id <-rownames(start_pos_list_df)
start_pos_list_df$peak_chr <- gsub( ":.*", "", start_pos_list_df$peak_id )

last_cols <- c( ( ncol(start_pos_list_df) - 1 ), ncol(start_pos_list_df) )

start_pos_list_df <- cbind( start_pos_list_df[,last_cols], start_pos_list_df[,-last_cols] )

write.csv( x = start_pos_list_df, row.names = FALSE, 
           file = paste0( prefix, "_start_position_across_iterations.csv" ) )

######################################################## Write out TF motif ####
cat( "Write out TF motif ...\n" )
write.table( x = get_motif(pdwn_coeffs, magnitud = "Estimate" ),
           file = paste0( prefix, "_motif_from_scaled_estimate.tab" ), 
           sep = "\t", quote = FALSE, col.names = FALSE )

write.table( x = get_motif(pdwn_coeffs, magnitud = "z value" ),
           file = paste0( prefix, "_motif_from_scaled_z_value.tab" ), 
           sep = "\t", quote = FALSE, col.names = FALSE )


############################################### Write out complete motif 
cat( "Write out complete motif  ...\n" )
motif <- as.matrix(get_motif(pdwn_coeffs, magnitud = "z value" ))
p_logo <- ggseqlogo( data = motif, method = "custom", seq_type = "dna" ) +
                     ylab( "z value" )

ggsave(filename = paste0( prefix, "_complete_motif_logo.pdf" ), plot = p_logo)

write.pfm.cisbp( filename = paste0( prefix, "_complete_motif.pfm.txt" ), 
                 motif = motif, motif_name = opt$experiment )


############################################### Write out trimmed motif
cat( "Write out trimmed motif ...\n" )
motif <- trim_motif( motif, base_th = 0.2 )
p_logo <- ggseqlogo( data = motif, method = "custom", seq_type = "dna" ) +
                     ylab( "z value" )

ggsave(filename = paste0( prefix, "_trimmed_motif_logo.pdf" ), plot = p_logo)

write.pfm.cisbp( filename = paste0( prefix, "_trimmed_motif.pfm.txt" ), 
                 motif = motif, motif_name = opt$experiment )


####################################################################### End ####
warning()
cat(paste0( "End wall time: ", Sys.time(), "\n"))
sink()





################################################## Use only 90% for testing ####
# ninety <- sample(1:nrow(dat_all$acc), floor( nrow(dat_all$acc)*0.9) )
# dat_all$acc <- dat_all$acc[ninety, ]
# dat_all$x.Met.all <- dat_all$x.Met.all[ninety, ]
# dat_all$x.A.all <- dat_all$x.A.all[ninety, ]
# dat_all$x.C.all <- dat_all$x.C.all[ninety, ]
# dat_all$x.G.all <- dat_all$x.G.all[ninety, ]
# dat_all$x.T.all <- dat_all$x.T.all[ninety, ]
# dat_all$x.CpG.all <- dat_all$x.CpG.all[ninety, ]
# dat_all$x.CG.all <- dat_all$x.CG.all[ninety, ]
# dat_all$target <- dat_all$target[ninety, ]
# dat_all$x.M.all <- dat_all$x.M.all[ninety, ]
# dat_all$x.W.all <- dat_all$x.W.all[ninety, ]

# get reverse complement predictors
# predictors_list <- vapply( X = predictors_list,
#                            FUN = rev_complement_predictor,
#                            FUN.VALUE = list(possible_position),
#                            pfm_length = opt$pfm_length )
# 
# predictors_list <- vapply( X = predictors_list,
#                            FUN = rev_complement_predictor,
#                            FUN.VALUE = list(possible_position),
#                            pfm_length = opt$pfm_length )
# predictors_list <- rev_compl_predictors_list
# rev_compl_predictors_list <- predictors_list
# ht_path <- paste0( opt$output_dir, "/ht.pdf")
# vis_sequence(predictors_list[[101]], opt$pfm_length, ht_path )
# ht_path <- paste0( opt$output_dir, "/ht_rev.pdf")
# vis_sequence(rev_compl_predictors_list[[101]], opt$pfm_length, ht_path )