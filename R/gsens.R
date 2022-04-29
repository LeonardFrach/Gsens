load.lib = c('lavaan', 'dplyr', 'parallel')
install.lib <- load.lib[!load.lib %in% installed.packages()] # Install missing libraries
sapply(load.lib, require, character = TRUE) # Load libraries

#' @title Adjusting for genetic confounding using PGS for the outcome

#' Adjusting for genetic confounding in exposure-outcome associations using the polygenic score for the outcome.
#' This is the recommended function for most scenarios, and the only function that has been extended to the multiple exposure case.
 
#' @param data Either a data frame of raw data or a covariance/correlation matrix. 
#' Must contain the variables G (PGS for Y), Y (outcome) and X (exposure) or instead X1, X2, ..., Xi if multiple exposures are used.
#' If data is a cov/cor matrix, the additional argument sample.nobs is required.
#' @param sample.nobs Optional. Number of observations used for the data.
#' Only necessary if input is a covariance or correlation matrix.
#' @param h2 Heritability estimate of the outcome (Y).
#' Can be chosen to be any external value, e.g. SNP- or twin-heritability estimates.
#' @param estimator Optional. Type of estimator used.
#' Can be any estimator for continuous data implemented in lavaan, e.g. "ML" or "GLS" (default = "GLS").
#' @param se Optional. Method to compute standard errors, can be any method implemented in lavaan,
#' e.g. "robust.sem" or "bootstrap". If se = "bootstrap", bootstrapping *cannot* be used for estimation of confidence intervals,
#' which is however strongly recommended (default = "standard").
#' @param boot.ci.type Bootstrapping method used to compute the confidence intervals,
#' can be any method implemented in lavaan. Default is "perc".
#' @param fmi Optional. If full-information maximum likelihood estimation is desired (estimator = "ML", missing "ML"),
#' and fmi = TRUE, missingness is reported (default = FALSE).
#' @param print Optional. Can be one of c("exposure", "mediation", "confounding", "overlap", "all" or "summary").
#' If print = "all", all relevant parameters will be printed.
#' If print = "summary", the lavaan output of the summary() function will be printed (default = "summary").
#' @param ... Additional arguments passed from lavaan, including 'bootstrap' (number of bootstraps for CIs, default = 1000)
#' and 'constraints' (constraints on the model).
 
#' @return Estimates for the adjusted exposure-outcome associations, exposure-mediated genetic effects,
#' genetic confounding and genetic overlap.

#' @examples
#' df <- data.frame(G, X1, X2, Y) 
#' gsensY(df, h2 = 0.5);

#' @author Leonard Frach & Jean-Baptiste Pingault
#' @export
#' @import lavaan
#' @import parallel
#' @import dplyr

#' @references Frach, L., Rijsdijk, F., Dudbridge, F. & Pingault, J. B. (in preparation).
#' Adjusting for genetic confounding using polygenic scores within structural equation models
#' @references Pingault, J. B., Rijsdijk, F., Schoeler, T., Choi, S. W., Selzam, S., Krapohl, E., ... & Dudbridge, F. (2021).
#' Genetic sensitivity analysis: Adjusting for genetic confounding in epidemiological associations. PLoS genetics, 17(6), e1009590. https://doi.org/10.1371/journal.pgen.1009590

gsensY = function(data,
                  sample.nobs = NULL,
                  h2,
                  estimator = "GLS",
                  se = "standard",
                  boot.ci.type = "perc",
                  fmi = F, 
                  print = "summary", ...) {
    
    
    if (dim(data)[1] != dim(data)[2]) {
        message("Using raw data as input.")
        
        data <- as.data.frame(data)
        
        # no correlation matrix
        cor <- NULL
        
        # standardize all variables
        data <- data %>% mutate_all(scale)
        message("Your variables have been standardized.")
        
    } else if (dim(data)[1] == dim(data)[2] & all(diag(data) == 1)) {
        message("Using correlation matrix as input.")
        
        cor <- data
        data <- as.data.frame(data)
        
    } else if (dim(data)[1] == dim(data)[2] & !all(diag(data) == 1)) {
        message("Converting covariance matrix to correlation matrix.")
        
        cor <- cov2cor(data)
        data <- as.data.frame(cor)
    }
    
    
    # select exposures
    namesX <- colnames(select(data, dplyr::starts_with("X")))
    NX <- length(namesX) 
    
    
    # check if the data frame contains more variables than needed
    if (dim(data)[2] > (NX + 2)) {
        stop("Your data contain unused variables.")
    }
    
    
    # create list of required column names (alternatively use lapply) 
    for (i in 1:NX) {
        namesX[i] <- paste("X", i, sep = "")
        names <- c(namesX, "Y", "G", "X")
    }
    
    
    # check if col names are in data frame. X is optional, X1 would also be acceptable when using only one exposure.
    if (all((names %in% colnames(data))[1:length(names) - 1]) == T |
        all((names %in% colnames(data))[2:length(names)]) == T) {
    } else {
        stop("Your data do not contain all the required variable names: G, Y and X/X1")
    }
    
    
    # create covariance structure and label for the model
    
    covstruc <- outer(namesX, namesX, function(x, y) paste(x, "~~", y))
    labelsa = paste0("a", 1:NX) 
    labelsb = paste0("b", 1:NX)
    labelsm = paste0("m", 1:NX)
    labels_gc = paste0("gc_", labelsb)
    labels_go = paste0("go_", labelsb)
    
    # if correlation matrix is used
    if (!is.null(cor)) {
        data <- NULL
    }
    
    
    ### main lavaan model ###
    
    gC <- numeric() 
    
    
    ## model specification
    
    model = c(
        paste("Y ~", paste0(c(labelsb,"c"),"*", c(namesX,"GG"), collapse = " + ")),    # Y depends on X and true polygenic score
        paste(namesX,"~", paste0(labelsa,"*","GG")),                                   # X1-Xi depend on true polygenic score
        covstruc[lower.tri(covstruc, diag = TRUE)],                                    # covariance structure 
        "GG =~ l*G" , 
        "GG ~~ 1*GG",
        "G ~~ me*G",
        "Y ~~ Y",
        # total mediation effect
        paste("m :=", paste0(labelsa,"*",labelsb, collapse = " + ")),
        paste0(labelsm, " := ", labelsa,"*",labelsb),                                  # specific mediation pathways for each G->Xi->Y
        # heritability constraints
        paste0("h := ","t(matrix(","c(", 
               paste0(labelsa, collapse = ", "),")))"," %*% ","matrix(c(", 
               paste0(labelsb, collapse = ", "),"))"," + c"),
        paste('h == sqrt(',h2,')'),
        for (i in 1:NX) {                                                              # genetic confounding for each Xi->Y association 
            gC <- c(gC, paste0(labels_gc[i], " := ", paste0(labelsa[i],"*", labelsb[-i],"*", labelsa[-i],
                                                            collapse = " + "), " + ", labelsa[i], "*c")) 
        },
        gC,
        paste0(labels_go, " := ", labelsa,"*", labelsa,"*", labelsb, " + ", labels_gc) # genetic overlap for each Xi->Y association 
    )
    
    
    ## model estimation
    
    # number of cores for parallelization, if more than one core available, number of cores minus 1
    nCores <- ifelse(parallel::detectCores() == 1, 1, parallel::detectCores() - 1)
    
    if(se == "bootstrap") {
        message("Bootstrapping will be performed for the standard errors. This might take a while.")
    }
    
    fit <- lavaan(model, data = data, sample.nobs = sample.nobs, sample.cov = cor,
                  estimator = estimator, se = se, ...)
    
    
    ## parameter estimates
    
    # if FIML is desired, proportion of missingness can be reported
    if (fmi) {
        pe <- parameterEstimates(fit, standardized = T, boot.ci.type = boot.ci.type, fmi = T)
        
    } else if (se != "bootstrap") {
        pe <- parameterEstimates(fit, boot.ci.type = boot.ci.type, standardized = T)
        
    } else {
        pe <- parameterEstimates(fit, standardized = T)
    }
    
    if(se == "bootstrap") {
        warning("Bootstrapping was NOT used for estimation of confidence intervals, which is however, strongly recommended. Consider changing your se = 'bootstrap' argument.")
    }
    
    results <- data.frame(rbind(
        pe[pe$label %in% labelsb,],
        pe[pe$label %in% labelsm,],
        pe[pe$label %in% "m",],
        pe[pe$label %in% labels_gc,],
        pe[pe$label %in% labels_go,]
    ))[, c(5:10)]
    
    results <- results %>%
        mutate_if(is.numeric, round, 3) # round all numeric variables
    
    
    # name the effects of the exposures
    for (i in 1:NX) {
        rownames(results)[i] <- c(paste("Adjusted Bx", i, "y", sep = ""))
    }
    
    # name the genetic effects mediated by the exposures
    for (i in 1:NX) {
        rownames(results)[NX + i] <- c(paste("Mediation m", i, sep = ""))
    }
    
    # total mediation
    rownames(results)[2*NX + 1] <- "Total mediation"
    
    # name the genetic confounding for each exposure-outcome association
    for (i in 2:(NX + 1)) {
        rownames(results)[(2*NX + i)] <- c(paste("Genetic confounding Bx", (i - 1), "y", sep = ""))
    }  
    
    # name the genetic overlap for each exposure-outcome association
    for (i in 2:(NX + 1)) {
        rownames(results)[3*NX + i] <- c(paste("Genetic overlap x", (i - 1), "y", sep = ""))
    }  
    
    # rename estimate column to est.std as in standardizedSolution()
    colnames(results)[1] <- "est.std"
    
    
    results$pvalue = formatC(2*pnorm(-abs(results$z)), digits = 3)
    
    
    ## Print results
    
    if (print == "summary") {summary(fit, standardized = T)}
    if (print == "exposure") {print(results[1:NX, ])} 
    if (print == "mediation") {print(results[(NX + 1):(3*NX - 1), ])} 
    if (print == "confounding") {print(results[(3*NX):(4*NX - 1), ])} 
    if (print == "overlap") {print(results[(4*NX):(5*NX - 1), ])} 
    if (print == "all") {print(results)} 
    
}  



#' @title Adjusting for genetic confounding using PGS for the exposure

#' Adjusting for genetic confounding in exposure-outcome associations using the polygenic score for the exposure.

#' @param rxy The observed phenotypic correlation between exposure X and outcome Y.
#' @param rgx The observed correlation between the polygenic score for X and exposure X.
#' @param rgy The observed correlation between the polygenic score for X and outcome Y.
#' @param n Sample size
#' @param h2 The additive genetic variance explained in exposure X under the scenario of interest.
#' @param print Optional. Enables the examination of model parameters (default = FALSE).
#' @param constrain Optional. Argument to constrain model parameters (default = NULL).

#' @return Estimates for the adjusted exposure-outcome associations, exposure-mediated genetic effects,
#' genetic confounding and genetic overlap.

#' @author Jean-Baptiste Pingault, Tabea Schoeler & Frank Dudbridge
#' @export
#' @import lavaan
#' @import dplyr

#' @references Pingault, J.-B., O’Reilly, P. F., Schoeler, T., Ploubidis, G. B., Rijsdijk, F., & Dudbridge, F. (2018).
#' Using genetic data to strengthen causal inference in observational research. Nature Reviews Genetics, 19(9), 566–580.
#' https://doi.org/10.1038/s41576-018-0020-3  
#' @references Pingault, J. B., Rijsdijk, F., Schoeler, T., Choi, S. W., Selzam, S., Krapohl, E., ... & Dudbridge, F. (2021).
#' Genetic sensitivity analysis: Adjusting for genetic confounding in epidemiological associations. PLoS genetics, 17(6), e1009590.
#' https://doi.org/10.1371/journal.pgen.1009590

gsensX = function(rxy,
                  rgx,
                  rgy,
                  n,
                  h2,
                  constrain = NULL,
                  print = FALSE) {
    mat = matrix(c(1, rgx, rgy, rgx, 1, rxy, rgy, rxy, 1), ncol = 3, nrow = 3)
    colnames(mat) = c("G","X","Y"); rownames(mat) = c("G","X","Y")
    
    model1 <- paste('
                  Y ~ bxy*X + bgy*GG    # Y depends on X and true polygenic score
                  X ~ bgx*GG            # X depends on true polygenic score
                  GG =~ l*G             # true polygenic score is estimated by G
                  GG ~~ 1*GG            # true polygenic score will be standardised
                  Y ~~ Y                # residual error of Y
                  X ~~ X                # residual error of X
                  G ~~ vg*G             # measurement error in G due to SNP selection & sampling error
                  bgx == sqrt(',h2,')
                  conf := bgx*bgy       # Confounding effect in standardized model
                  total := conf + bxy
                  ', constrain,'        # optional constraints
                  ')
    
    
    fit1 = lavaan(model1, sample.cov = mat, sample.nobs = n, estimator = "GLS")
    if (print) {  summary(fit1)}
    pe = parameterestimates(fit1)
    pe$pvalue = formatC(2*pnorm(-abs(pe$z)), digits = 5)
    results = data.frame(rbind(
        pe[pe$label == "bxy",],
        pe[pe$label == "conf",],
        pe[pe$label == "total",]
    ))[,5:10]
    results = results %>%
        mutate_if(is.numeric, round, 3) # round all numeric variables
    
    rownames(results) = c("Adjusted Bxy","Genetic confounding","Total effect")
    results
}


## Two polygenic scores ##

#' @title Adjusting for genetic confounding using PGS for the exposure and the outcome

#' Adjusting for genetic confounding in exposure-outcome associations using the polygenic score for the exposure and the outcome

#' @param rxy The observed phenotypic correlation between exposure X and outcome Y.
#' @param rg1x The observed correlation between the polygenic score for X and exposure X.
#' @param rg2y The observed correlation between the polygenic score for X and outcome Y.
#' @param rg1y The correlation between the outcome Y and the observed polygenic score for X.
#' @param rg2x The correlation between the exposure X and the observed polygenic score for Y.
#' @param rg1g2 The correlation between the two observed polygenic scores.
#' @param n Sample size
#' @param h2.x The additive genetic variance explained in exposure X under the scenario of interest.
#' @param h2.y The additive genetic variance explained in exposure Y under the scenario of interest.
#' @param print Optional. Enables the examination of model parameters (default = FALSE).
#' @param constrain Optional. Argument to constrain model parameters (default = NULL).

#' @return Estimates for the adjusted exposure-outcome associations, exposure-mediated genetic effects,
#' genetic confounding and genetic overlap.

#' @author Jean-Baptiste Pingault, Tabea Schoeler & Frank Dudbridge
#' @export
#' @import lavaan
#' @import dplyr

#' @references Pingault, J.-B., O’Reilly, P. F., Schoeler, T., Ploubidis, G. B., Rijsdijk, F., & Dudbridge, F. (2018).
#' Using genetic data to strengthen causal inference in observational research. Nature Reviews Genetics, 19(9), 566–580.
#' https://doi.org/10.1038/s41576-018-0020-3  
#' @references Pingault, J. B., Rijsdijk, F., Schoeler, T., Choi, S. W., Selzam, S., Krapohl, E., ... & Dudbridge, F. (2021).
#' Genetic sensitivity analysis: Adjusting for genetic confounding in epidemiological associations. PLoS genetics, 17(6), e1009590.
#' https://doi.org/10.1371/journal.pgen.1009590


gsensXY = function(rxy,
                   rg1x,
                   rg2y,
                   rg1y,
                   rg2x,
                   rg1g2,
                   n,
                   h2.x,
                   h2.y,
                   constrain = NULL, 
                   print = FALSE) {
    
    lower = c(
        1,
        rg1g2,1,
        rg1x,rg2x,1,
        rg1y,rg2y,rxy,1)
    
    mat = getCov(lower, names = c("G1","G2","X","Y"))
    
    model1 <- paste('
                  Y ~ bxy*X + bg1y*GG1 + bg2y*GG2   # Y depends on X and true polygenic scores
                  X ~ bg1x*GG1 + bg2x*GG2           # X depends on true polygenic scores
                  GG1 =~ lg1*G1                     # true polygenic score is estimated by G1
                  GG2 =~ lg2*G2                     # true polygenic score is estimated by G2
                  GG1 ~~ 1*GG1                      # true polygenic score will be standardised
                  GG2 ~~ 1*GG2                      # true polygenic score will be standardised
                  Y ~~ vy*Y                         # residual error of Y
                  X ~~ vx*X                         # residual error of X
                  GG1 ~~ bg1g2*GG2                  # correlation of true polygenic scores
                  G1 ~~ vg1*G1                      # measurement error in G1 due to SNP selection & sampling error
                  G2 ~~ vg2*G2                      # measurement error in G2 due to SNP selection & sampling error
                  
                  # heritability constraints
                  bg1x + bg1g2*bg2x == sqrt(',h2.x,')
                  bg2y + (bg2x + bg1g2*bg1x)*bxy + bg1g2*bg1y == sqrt(',h2.y,') # SNP heritability constraints
                  
                  # Confounding and total effects
                  conf := bg1x*bg1y + bg2x*bg2y + bg1x*bg2y*bg1g2 + bg2x*bg1y*bg1g2
                  total := conf + bxy
                  ', constrain,' #optional constraints
                  
                  ')
    
    fit1 = lavaan(model1, sample.cov = mat, sample.nobs = n, estimator = "GLS")
    if (print) {  summary(fit1)}
    pe = parameterestimates(fit1)
    pe$pvalue = formatC(2*pnorm(-abs(pe$z)), digits = 5)
    results = data.frame(rbind(
        pe[pe$label == "bxy",],
        pe[pe$label == "conf",],
        pe[pe$label == "total",]
    ))[, 5:10]
    
    results = results %>%
        mutate_if(is.numeric, round, 3) # round all numeric variables
    
    rownames(results) = c("Adjusted Bxy","Genetic confounding","Total effect")
    results
    
}



