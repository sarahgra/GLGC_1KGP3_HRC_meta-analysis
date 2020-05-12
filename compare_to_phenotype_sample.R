#!/usr/bin/Rscript

library(data.table)
library(dplyr)
library(argparse)
library(GenABEL)
library(ggplot2)

parser <- ArgumentParser()
parser$add_argument("-r", "--results_file", default="all_MGI_results.txt", help="This file should contain the file names of the results files, formatted as columns with 1) individual 2) risk_score")
args <- parser$parse_args()


risk_score_files <- fread(args$results_file, header=F)

#Add necessary steps to generate a data table with LDL values and covariates
#Adjust LDL for meds
ldl_transformed$LDL <- ifelse(ldl_transformed$statin_med == "Yes", ldl_transformed$LDL/0.7, ldl_transformed$LDL)


#Table to store results
ALL_R2_results <- data.table(NULL)

for (ancestry in c("Caucasian", "African American")){
    ldl_ancestry <- ldl_transformed[ldl_transformed$Race == ancestry,]

    n_rows <- length(risk_score_files$V1) + 1
    R2_results <- data.table(Risk_score=c("Covariates_only", risk_score_files$V1), adj_R2=numeric(n_rows), Lower_95_CI=numeric(n_rows), Upper_95_CI=numeric(n_rows))
    
    covar_only <- summary(lm(as.formula(paste0('LDL_mean~gender+BATCH+Age+PC1+PC2+PC3+PC4')), data=ldl_ancestry))
    
    #This will shuffle the data up, and select the same number of people but there will be some duplicates (so not the exact same data table), calculate linear model again, check adj R2.  Called percentile bootstrapping
    #Textbook example with percentile bootstrapping: http://pages.stat.wisc.edu/~larget/stat302/chap3.pdf
    ci <- c()
    for(j in 1:1000){
        tmp <- summary(lm(as.formula(paste0('LDL_mean~gender+BATCH+Age+PC1+PC2+PC3+PC4')), data=ldl_ancestry[sample(nrow(ldl_ancestry), replace=T),]))
        ci <- c(ci,tmp$adj.r.squared)
      if(j%%100==0) print(j)
    }
    ci <- ifelse(ci<0, 0, ci)
    
    R2_results[R2_results$Risk_score == "Covariates_only", "adj_R2"] <- covar_only$adj.r.squared[1]
    R2_results[R2_results$Risk_score == "Covariates_only", "Lower_95_CI"] <- quantile(ci, c(0.025))
    R2_results[R2_results$Risk_score == "Covariates_only", "Upper_95_CI"] <- quantile(ci, c(0.975))
    
    
    for (risk_score_file_n in risk_score_files$V1){
        print(paste("Now analyzing:", risk_score_file_n))
    
        risk_score <- fread(risk_score_file_n)
        if ("risk_score_sum" %in% names(risk_score)){
            setnames(risk_score, "risk_score_sum", "score")
        }
        risk_score <- left_join(risk_score, ldl_ancestry, by=c("individual"="Encrypted_PatientID"))
        
        #Normalize the risk score
        risk_score$normalized <- rntransform(risk_score$score)
        risk_score <- risk_score[!(is.na(risk_score$score)|is.na(risk_score$LDL_mean)),]
        
        predicted_LDL_with_covars <- summary(lm(as.formula(paste0('LDL_mean~gender+BATCH+Age+PC1+PC2+PC3+PC4+normalized')), data=risk_score))
        ci <- c()
        for(j in 1:1000){
            tmp <- summary(lm(as.formula(paste0('LDL_mean~gender+BATCH+Age+PC1+PC2+PC3+PC4+normalized')), data=risk_score[sample(nrow(risk_score), replace=T),]))
         ci <- c(ci,tmp$adj.r.squared)
          if(j%%100==0){ print(j)}
        }
        ci <- ifelse(ci<0, 0, ci)
        
        R2_results[R2_results$Risk_score == risk_score_file_n, "adj_R2"] <- predicted_LDL_with_covars$adj.r.squared
        R2_results[R2_results$Risk_score == risk_score_file_n, "Lower_95_CI"] <- quantile(ci, c(0.025))
        R2_results[R2_results$Risk_score == risk_score_file_n, "Upper_95_CI"] <- quantile(ci, c(0.975))
    }
    write.table(R2_results, paste0("all_adjR2_results_different_models_", sub(" ", "", ancestry), ".txt"), col.names=T, row.names=F, sep="\t", quote=F)
}




