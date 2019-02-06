#' @importFrom data.table dcast copy melt
#' @importFrom outliers chisq.out.test
#' @importFrom MASS fitdistr


#' @export
keep_only_proteotypic <- function(input_dt) {

return(input_dt[which(grepl("^1/", input_dt$ProteinName)),])

}


#' @export
long2wide <- function(input_dt, level = "PeptideIon") {
  
if (!level %in% c("PeptideIon", "Protein")) {
stop("Please select a valid type. Options:  \"PeptideIon(default)\", \"Protein\"")
}
 
  if (level == "PeptideIon") {
    output_dt <-
      dcast(input_dt, ProteinName + PeptideIon ~ filename, value.var = "Intensity")
  } else if (level == "Protein") {
    output_dt <-
      dcast(input_dt, ProteinName ~ run_id, value.var = "Quant_sum")
  }
  
  return(output_dt)
  
}


#' @export
wide2long <- function(input_dt) {
  
  output_dt <-
    melt(
      input_dt,
      id.vars = c("PeptideIon", "ProteinName", "FullName", "NTT", "NMC"),
      measure.vars = names(input_dt)[which(grepl("mean_", names(input_dt)))],
      variable.name = "run_id",
      value.name = "Intensity",
      na.rm = T
    )
  
  return(output_dt)
  
}


#' @export
keep_only_semi <- function(input_dt) {
  
  semi_cons_peptides_long <- copy(input_dt)
  semi_cons_peptides_long[which(semi_cons_peptides_long$NTT == 2),]$Intensity <- 0
  
  return(semi_cons_peptides_long)
  
}


#' @export
pept2prot <- function(input_pept_long, intensity_normalization = "sqrt") {
  
    #for CPTAC itraq breast cancer and benchmarking data, this seems that sqrt works better.
    if (intensity_normalization == "sqrt") {
      output_prot_long <-
        input_pept_long[, .(Quant_sum = sum_na(sqrt(Intensity + 1)),
                            NumPeptides = length(which(Intensity > 0))), by = .(ProteinName, run_id)]
      
    } else if (intensity_normalization == "log") {
      output_prot_long <-
        input_pept_long[, .(Quant_sum = sum_na(log2(Intensity + 1)),
                            NumPeptides = length(which(Intensity > 0))), by = .(ProteinName, run_id)]
      
    } else if (intensity_normalization == "raw") {
      output_prot_long <-
        input_pept_long[, .(Quant_sum = sum_na(Intensity),
                            NumPeptides = length(which(Intensity > 0))), by = .(ProteinName, run_id)]
      
    }
    
    
    return(output_prot_long)
    
  }


#' @export
normalize_data <- function(input_dt, input_index, normalization = "mediancenter") {
    
    output_dt <- copy(input_dt)
    
    output_dt[, (input_index) := log2(output_dt[, input_index, with = F])]
    
    run_median <-
      sapply(output_dt[, input_index, with = F], median_na)
    
    for (i in 1:length(input_index)) {
      output_dt[, names(output_dt)[input_index[i]]] <-
        output_dt[, input_index[i], with = F] - run_median[i] + median(run_median)
      output_dt[, names(output_dt)[input_index[i]]] <-
        2 ^ output_dt[, input_index[i], with = F]
    }
    
    return(output_dt)
    
  }


#' @export
calculate_stats <- function(input_dt, input_index) {
  
  index_semi <- which(input_dt$NTT < 2)
  index_missed <- which(input_dt$NMC > 0)
  index_normal <- which(input_dt$NTT == 2 & input_dt$NMC == 0)
  
  numPept <- apply(input_dt[, input_index, with = F], 2, function(x)
      length(which(x > 0)))
  
  numSemi <-
    apply(input_dt[index_semi, input_index, with = F], 2, function(x)
      length(which(x > 0)))
  numMissed <-
    apply(input_dt[index_missed, input_index, with = F], 2, function(x)
      length(which(x > 0)))
  numNormal <-
    apply(input_dt[index_normal, input_index, with = F], 2, function(x)
      length(which(x > 0)))
  
  TIC <- apply(input_dt[, input_index, with = F], 2, sum_na)
  
  proteomic_measure <- data.frame(cbind(names(input_dt)[input_index], numPept, numSemi, numMissed, numNormal, TIC))
  names(proteomic_measure) <-
    c("filename",
      "numPept",
      "numSemi",
      "numMissed",
      "numNormal",
      "TIC")
  
  return(proteomic_measure)
  
}


#' @export
merge_replicates <- function(input_dt, input_anno) {
  
  output_dt <- copy(input_dt)
  
  list_samples <- unique(input_anno$sample_id)
  
  message("It starts to merge replicates...")
  
  for (i in 1:length(list_samples)) {
    cat("Processing ", i, " sample: ", list_samples[i], "\n")
    output_dt[, paste0("mean_", list_samples[i]) := apply(.SD, 1, mean_na), .SDcols = input_anno[input_anno$sample_id %in% list_samples[i],]$filename]
  }
  
  message("Done with merging replicates...")
  
  return(output_dt)
  
}


#' @export
generate_iPIS_matrix <- function(input_all_dt, input_semi_dt) {
  
  index_temp <- which(grepl("^mean", names(input_all_dt)))

  m_iPIS <- copy(input_all_dt)
  
  #for(i in 1:length(index_temp)) {
  #m_iPIS[, index_temp[i], with=F] <- 1 - input_semi_dt[, index_temp[i], with=F] / input_all_dt[, index_temp[i], with=F]
  #}
  
  m_iPIS[, index_temp] <- 1 - input_semi_dt[, index_temp, with = F] / input_all_dt[, index_temp, with = F]
  
  names(m_iPIS)[index_temp] <- sapply(strsplit(names(m_iPIS)[index_temp], "mean_"), "[[", 2)

  index_order <- order(apply(m_iPIS[, index_temp, with = F], 1, function(x) length(which(is.na(x)))))

  m_iPIS <- m_iPIS[index_order, ]

#m_iPIS <- m_iPIS[order(apply(m_iPIS[, index_temp, with = F], 1, function(x) length(which(is.na(x))))), ]

  return(m_iPIS)
  
}


#' @export
old_calculate_PIN <- function(input_iPIS) {
  
  index_temp <- seq(2, dim(input_iPIS)[2], 1)
  
  m_PIN <- data.frame(
    sample_id = names(input_iPIS)[index_temp],
    PIN = apply(input_iPIS[, index_temp, with =
                             F], 2, mean_na),
    pval = 1.0
  )
  
  m_PIN <- m_PIN[order(m_PIN$PIN), ]
  
  #wenguang: re-scale PIN... not sure whether it is a good idea...
  #m_PIN$PIN <- -log2(0.98 - m_PIN$PIN)
  
  for (i in 1:(dim(m_PIN)[1] - 1)) {
    
    if (grepl("lowest",
              chisq.out.test(m_PIN$PIN[i:dim(m_PIN)[1]], opposite = F)$alternative)) {
      m_PIN[i,]$pval <-
        chisq.out.test(m_PIN$PIN[i:dim(m_PIN)[1]], opposite = F)$p.value
    } else {
      m_PIN[i,]$pval <-
        chisq.out.test(m_PIN$PIN[i:dim(m_PIN)[1]], opposite = T)$p.value
    }
    
    if (i > 1) {
      if (m_PIN[i,]$pval < m_PIN[i - 1,]$pval) {
        m_PIN[i,]$pval <- m_PIN[i - 1,]$pval + m_PIN[i,]$pval * 0.05
        #m_PIN[i, ]$pval <- m_PIN[i, ]$pval*0.6 + m_PIN[i-1, ]$pval * 0.4
      }
    }
    
  }
  
  
  #wenguang: need to be tested, and will be released later...
  
  if (1 == 1) {
    #index_sig <- which(m_PIN$pval < 0.01)
    #index_insig <- which(m_PIN$pval >= 0.01)
    
    index_sig <- which(m_PIN$pval < 0.05)
    index_insig <- which(m_PIN$pval >= 0.05)
    
    for (i in 1:length(index_sig)) {
      if (grepl("lowest",
                chisq.out.test(m_PIN$PIN[c(index_sig[i], index_insig)], opposite = F)$alternative)) {
        m_PIN[index_sig[i],]$pval <-
          chisq.out.test(m_PIN$PIN[c(index_sig[i], index_insig)], opposite = F)$p.value
      } else {
        m_PIN[index_sig[i],]$pval <-
          chisq.out.test(m_PIN$PIN[c(index_sig[i], index_insig)], opposite = T)$p.value
      }
      
      
    }
    
    #p.adjust(raw$pval, method="hommel")
    
  }
  
  return(m_PIN)
  
}





#' @export
calculate_PIN <- function(input_iPIS, input_remove_zero_rows) {
  
  index_temp <- seq(2, dim(input_iPIS)[2], 1)

#wenguang: CAUTION!!! This will calculate PINs without proteins detected consistently with only semi-tryptic peptides across all measurements (i.e. the row is full with zero).
#wenguang: by default it is off. It was only turned on, that yansheng plasma samples were analyzed, as those rows took up a large fraction of total proteins identified.
if(input_remove_zero_rows == TRUE) {
  input_iPIS <- input_iPIS[ - which( apply(input_iPIS[, index_temp, with=F], 1, function(x) length(which(x == 0.0))) == length(index_temp) ), ]
}
  
  m_PIN <- data.frame(
    sample_id = names(input_iPIS)[index_temp],
    PIN = apply(input_iPIS[, index_temp, with=F], 2, mean_na),
    pval = 1.0
  )
  
  m_PIN <- m_PIN[order(m_PIN$PIN), ]
  
  fit_weibull <- fitdistr(m_PIN$PIN, densfun = "weibull")
  m_PIN$pval <-   pweibull(m_PIN$PIN,
                           shape = fit_weibull$estimate[1],
                           scale = fit_weibull$estimate[2])
  
  length_insig_last <- dim(m_PIN)[1]
  count <- 1
  index_insig <- which(m_PIN$pval >= 0.05)
  
  while (length(index_insig) < length_insig_last) {
    
    cat("Iteration for fitting null distribution: ", count, "\n")
    

#wenguang: note that for error like "Error in fitdistr(m_PIN[index_insig, ]$PIN, densfun = "weibull" optimization failed", you may want to try 
#          fit_weibull <- fitdistr(m_PIN[index_insig,]$PIN, densfun = "weibull", lower=c(10, 0.6))

if( !is(try(fitdistr(m_PIN[index_insig, ]$PIN, densfun = "weibull"), silent=T), "try-error") ) {

    fit_weibull <- fitdistr(m_PIN[index_insig,]$PIN, densfun = "weibull")

} else if( !is(try(fitdistr(m_PIN[index_insig, ]$PIN, densfun = "weibull", lower=c(10, 0.6)), silent=T), "try-error") ) {

   fit_weibull <- fitdistr(m_PIN[index_insig, ]$PIN, densfun = "weibull", lower=c(10, 0.6))

} else if( !is(try(fitdistr(m_PIN[index_insig, ]$PIN, densfun = "weibull", lower=c(1, 0.1)), silent=T), "try-error") ) {

   fit_weibull <- fitdistr(m_PIN[index_insig, ]$PIN, densfun = "weibull", lower=c(1, 0.1))

}

    m_PIN$pval  <- pweibull(m_PIN$PIN,
                           shape = fit_weibull$estimate[1],
                           scale = fit_weibull$estimate[2])
    
    length_insig_last <- length(index_insig)
    count <- count + 1
    
index_insig <- which(m_PIN$pval >= 0.05)
#index_insig <- which(m_PIN$pval >= 0.02)
    
  }
  
  #wenguang: need to be tested, and will be released later...
  
  if (1 == 1) {
    #index_sig <- which(m_PIN$pval < 0.02)
    #index_insig <- which(m_PIN$pval >= 0.02)
    
    index_sig <- which(m_PIN$pval < 0.05)
    index_insig <- which(m_PIN$pval >= 0.05)



    for (i in 1:length(index_sig)) {
      
      if (grepl("lowest",
                chisq.out.test(
                  m_PIN$PIN[c(index_sig[i], index_insig)],
                  variance = var(m_PIN[index_insig,]$PIN),
                  opposite = F
                )$alternative)) {
        m_PIN[index_sig[i],]$pval <-
          chisq.out.test(m_PIN$PIN[c(index_sig[i], index_insig)],
                         variance = var(m_PIN[index_insig,]$PIN),
                         opposite = F)$p.value
      } else {
        m_PIN[index_sig[i],]$pval <-
          chisq.out.test(m_PIN$PIN[c(index_sig[i], index_insig)],
                         variance = var(m_PIN[index_insig,]$PIN),
                         opposite = T)$p.value
      }
      
      
    }
    
    #p.adjust(raw$pval, method="hommel")
    
  }
  
  return(m_PIN)
  
}







#' @export
calculate_PIN_new <- function(input_iPIS, input_remove_zero_rows) {

  insig_pval <- 0.05

  index_temp <- seq(2, dim(input_iPIS)[2], 1)

#wenguang: CAUTION!!! This will calculate PINs without proteins detected consistently with only semi-tryptic peptides across all measurements (i.e. the row is full with zero).
#wenguang: by default it is off. It was only turned on, that yansheng plasma samples were analyzed, as those rows took up a large fraction of total proteins identified.
  if(input_remove_zero_rows == TRUE) {
    input_iPIS <- input_iPIS[ - which( apply(input_iPIS[, index_temp, with=F], 1, function(x) length(which(x == 0.0))) == length(index_temp) ), ]
  }
  
  m_PIN <- data.frame(
    sample_id = names(input_iPIS)[index_temp],
    PIN = apply(input_iPIS[, index_temp, with=F], 2, mean_na),
    pval = 1.0
  )

  m_PIN <- m_PIN[order(m_PIN$PIN), ]
  
  fit_weibull <- fitdistr(m_PIN$PIN, densfun = "weibull")
  m_PIN$pval <-   pweibull(m_PIN$PIN,
                           shape = fit_weibull$estimate[1],
                           scale = fit_weibull$estimate[2])
  
  length_insig_last <- dim(m_PIN)[1]
  count <- 1
  index_insig <- which(m_PIN$pval >= insig_pval)

best_fit_count <- 1
best_fit_D <- 0
best_fit_bool <- FALSE
  
  while (length(index_insig) < length_insig_last) {
    
    cat("Iteration for fitting null distribution: ", count, "\n")
    

#wenguang: note that for error like "Error in fitdistr(m_PIN[index_insig, ]$PIN, densfun = "weibull" optimization failed", you may want to try 
#          fit_weibull <- fitdistr(m_PIN[index_insig,]$PIN, densfun = "weibull", lower=c(10, 0.6))

if( !is(try(fitdistr(m_PIN[index_insig, ]$PIN, densfun = "weibull"), silent=T), "try-error") ) {

    fit_weibull <- fitdistr(m_PIN[index_insig,]$PIN, densfun = "weibull")

} else if( !is(try(fitdistr(m_PIN[index_insig, ]$PIN, densfun = "weibull", lower=c(10, 0.6)), silent=T), "try-error") ) {

   fit_weibull <- fitdistr(m_PIN[index_insig, ]$PIN, densfun = "weibull", lower=c(10, 0.6))

} else if( !is(try(fitdistr(m_PIN[index_insig, ]$PIN, densfun = "weibull", lower=c(1, 0.1)), silent=T), "try-error") ) {

   fit_weibull <- fitdistr(m_PIN[index_insig, ]$PIN, densfun = "weibull", lower=c(1, 0.1))

}

    m_PIN$pval  <- pweibull(m_PIN$PIN,
                           shape = fit_weibull$estimate[1],
                           scale = fit_weibull$estimate[2])
    
    length_insig_last <- length(index_insig)
    
index_insig <- which(m_PIN$pval >= insig_pval)
    
print(ks.test(m_PIN[index_insig, ]$PIN, rweibull(n=10000, shape=fit_weibull$estimate[1], scale=fit_weibull$estimate[2])))

fit_test <- ks.test(m_PIN[index_insig, ]$PIN, rweibull(n=10000, shape=fit_weibull$estimate[1], scale=fit_weibull$estimate[2]))

if( fit_test$statistic > best_fit_D) {
best_fit_D <- fit_test$statistic
best_fit_count <- count

if( fit_test$p.value < 0.01 ) best_fit_bool <- TRUE

}


    count <- count + 1

  }


if(best_fit_bool == "TRUE") {

cat("It fits best when iteration", best_fit_count, "times.\n" )

  m_PIN <- data.frame(
    sample_id = names(input_iPIS)[index_temp],
    PIN = apply(input_iPIS[, index_temp, with=F], 2, mean_na),
    pval = 1.0
  )



 m_PIN <- m_PIN[order(m_PIN$PIN), ]
  
  fit_weibull <- fitdistr(m_PIN$PIN, densfun = "weibull")
  m_PIN$pval <-   pweibull(m_PIN$PIN,
                           shape = fit_weibull$estimate[1],
                           scale = fit_weibull$estimate[2])
  
  length_insig_last <- dim(m_PIN)[1]
  index_insig <- which(m_PIN$pval >= insig_pval)

reloop_count <- 1

  while (reloop_count <= best_fit_count) {
    
    cat("Iteration for fitting null distribution: ", reloop_count, "\n")
    

#wenguang: note that for error like "Error in fitdistr(m_PIN[index_insig, ]$PIN, densfun = "weibull" optimization failed", you may want to try 
#          fit_weibull <- fitdistr(m_PIN[index_insig,]$PIN, densfun = "weibull", lower=c(10, 0.6))

if( !is(try(fitdistr(m_PIN[index_insig, ]$PIN, densfun = "weibull"), silent=T), "try-error") ) {

    fit_weibull <- fitdistr(m_PIN[index_insig,]$PIN, densfun = "weibull")

} else if( !is(try(fitdistr(m_PIN[index_insig, ]$PIN, densfun = "weibull", lower=c(10, 0.6)), silent=T), "try-error") ) {

   fit_weibull <- fitdistr(m_PIN[index_insig, ]$PIN, densfun = "weibull", lower=c(10, 0.6))

} else if( !is(try(fitdistr(m_PIN[index_insig, ]$PIN, densfun = "weibull", lower=c(1, 0.1)), silent=T), "try-error") ) {

   fit_weibull <- fitdistr(m_PIN[index_insig, ]$PIN, densfun = "weibull", lower=c(1, 0.1))

}

    m_PIN$pval  <- pweibull(m_PIN$PIN,
                           shape = fit_weibull$estimate[1],
                           scale = fit_weibull$estimate[2])
    
    length_insig_last <- length(index_insig)
    
index_insig <- which(m_PIN$pval >= insig_pval)
#index_insig <- which(m_PIN$pval >= 0.02)
    
print(ks.test(m_PIN[index_insig, ]$PIN, rweibull(n=10000, shape=fit_weibull$estimate[1], scale=fit_weibull$estimate[2])))

fit_test <- ks.test(m_PIN[index_insig, ]$PIN, rweibull(n=10000, shape=fit_weibull$estimate[1], scale=fit_weibull$estimate[2]))

    reloop_count <- reloop_count + 1

  }


  if (1 == 0) {
    #index_sig <- which(m_PIN$pval < 0.02)
    #index_insig <- which(m_PIN$pval >= 0.02)
    
    index_sig <- which(m_PIN$pval < 0.01)
    index_insig <- which(m_PIN$pval >= 0.01)
    
    for (i in 1:length(index_sig)) {
      
      if (grepl("lowest",
                chisq.out.test(
                  m_PIN$PIN[c(index_sig[i], index_insig)],
                  variance = var(m_PIN[index_insig,]$PIN),
                  opposite = F
                )$alternative)) {
        m_PIN[index_sig[i],]$pval <-
          chisq.out.test(m_PIN$PIN[c(index_sig[i], index_insig)],
                         variance = var(m_PIN[index_insig,]$PIN),
                         opposite = F)$p.value
      } else {
        m_PIN[index_sig[i],]$pval <-
          chisq.out.test(m_PIN$PIN[c(index_sig[i], index_insig)],
                         variance = var(m_PIN[index_insig,]$PIN),
                         opposite = T)$p.value
      }
      
      
    }
    
    #p.adjust(raw$pval, method="hommel")
    
  }


} else {


  m_PIN <- data.frame(
    sample_id = names(input_iPIS)[index_temp],
    PIN = apply(input_iPIS[, index_temp, with=F], 2, mean_na),
    pval = 1.0
  )
  
  m_PIN <- m_PIN[order(m_PIN$PIN), ]
  
  fit_weibull <- fitdistr(m_PIN$PIN, densfun = "weibull")
  m_PIN$pval <-   pweibull(m_PIN$PIN,
                           shape = fit_weibull$estimate[1],
                           scale = fit_weibull$estimate[2])
  
  length_insig_last <- dim(m_PIN)[1]
  count <- 1
  index_insig <- which(m_PIN$pval >= 0.05)
  
  while (length(index_insig) < length_insig_last) {
    
    cat("Iteration for fitting null distribution: ", count, "\n")

#wenguang: note that for error like "Error in fitdistr(m_PIN[index_insig, ]$PIN, densfun = "weibull" optimization failed", you may want to try 
#          fit_weibull <- fitdistr(m_PIN[index_insig,]$PIN, densfun = "weibull", lower=c(10, 0.6))

if( !is(try(fitdistr(m_PIN[index_insig, ]$PIN, densfun = "weibull"), silent=T), "try-error") ) {

    fit_weibull <- fitdistr(m_PIN[index_insig,]$PIN, densfun = "weibull")

} else if( !is(try(fitdistr(m_PIN[index_insig, ]$PIN, densfun = "weibull", lower=c(10, 0.6)), silent=T), "try-error") ) {

   fit_weibull <- fitdistr(m_PIN[index_insig, ]$PIN, densfun = "weibull", lower=c(10, 0.6))

} else if( !is(try(fitdistr(m_PIN[index_insig, ]$PIN, densfun = "weibull", lower=c(1, 0.1)), silent=T), "try-error") ) {

   fit_weibull <- fitdistr(m_PIN[index_insig, ]$PIN, densfun = "weibull", lower=c(1, 0.1))

}

    m_PIN$pval  <- pweibull(m_PIN$PIN,
                           shape = fit_weibull$estimate[1],
                           scale = fit_weibull$estimate[2])
    
    length_insig_last <- length(index_insig)
    count <- count + 1
    
index_insig <- which(m_PIN$pval >= 0.05)
#index_insig <- which(m_PIN$pval >= 0.02)
    
  }
  
  #wenguang: need to be tested, and will be released later...
  
  if (1 == 1) {
    #index_sig <- which(m_PIN$pval < 0.02)
    #index_insig <- which(m_PIN$pval >= 0.02)
    
if( dim(m_PIN)[1] < 100 ) {
    index_sig <- which(m_PIN$pval < 0.04)
    index_insig <- which(m_PIN$pval >= 0.04)
} else {
    index_sig <- which(m_PIN$pval < 0.02)
    index_insig <- which(m_PIN$pval >= 0.02)
}


if( length(index_sig) > 0 ) {

    for (i in 1:length(index_sig)) {
      
      if (grepl("lowest",
                chisq.out.test(
                  m_PIN$PIN[c(index_sig[i], index_insig)],
                  variance = var(m_PIN[index_insig,]$PIN),
                  opposite = F
                )$alternative)) {
        m_PIN[index_sig[i],]$pval <-
          chisq.out.test(m_PIN$PIN[c(index_sig[i], index_insig)],
                         variance = var(m_PIN[index_insig,]$PIN),
                         opposite = F)$p.value
      } else {
        m_PIN[index_sig[i],]$pval <-
          chisq.out.test(m_PIN$PIN[c(index_sig[i], index_insig)],
                         variance = var(m_PIN[index_insig,]$PIN),
                         opposite = T)$p.value
      }
      
      
    }

}
    
    #p.adjust(raw$pval, method="hommel")
    
  }



}

  return(m_PIN)

}
