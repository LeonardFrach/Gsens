#' Adjusting for genetic confounding using PGS for the outcome
#'
#' Adjusting for genetic confounding in exposure--outcome associations using the polygenic score for the outcome.
#' This is the recommended function for most scenarios, and the only function that has been extended to the multiple exposure case.

#' @param data Optional. Used when input contains raw data (data frame).
#' @param sample.cov Optional. If input is a covariance or correlation matrix, the lavaan arguments `sample.cov = df_cov` and `sample.nobs = n` are required, where `df_cov` should be the name of your input matrix and `n` the sample size.
#' @param h2 Heritability estimate of the outcome (Y).
#' Can be chosen to be any external value, e.g. SNP- or twin-heritability estimates.
#' @param exposures Vector of variable name(s) of the exposure(s).
#' Example: `exposures = c("x1", "x2")`
#' @param outcome Name of the outcome variable.
#' @param pgs Name of the polygenic score variable (pgs corresponding to the outcome).
#' @param model Optional. If additional variables beyond exposures, PGS and outcome are used, e.g., measured confounders,
#' these additional relations can be specified and will be added to the GsensY model. Example: `model = c("X1 ~ Z", "X2 ~ Z")`. 
#' If `model` is not specified, the standard GsensY model will be estimated. 
#' @param ... Additional arguments passed from lavaan, including `se` (estimation method for the standard errors),
#' `estimator` (estimator used for model, default is ML), `bootstrap` (number of bootstraps for CIs, default = 1000),
#' `sample.nobs` (Number of observations for estimation using summary data, not recommended), and more.
#' See the lavaan documentation for details, e.g., [lavaan::lavaan()] or [lavaan::lavOptions()]

#' @return The Gsens model output will be returned as a lavaan object.
#' For example, the `summary()` or [lavaan::parameterEstimates()] functions can be used for more detailed outputs, e.g. for standardized estimates.

#' @examples
#' \dontrun{
#' df <- data.frame(X1, X2, X3, Y, PGS_outcome)
#' df_cov <- cov(df)
#' gsens_out <- gsensY(data = df, h2 = 0.5, exposures = c("X1", "X2", "X3"), outcome = "Y", pgs = "PGS_outcome")
#' gsens_out_cov <- gsensY(sample.cov = df_cov, sample.nobs = 5000, h2 = 0.5, exposures = c("X1", "X2", "X3"), outcome = "Y", pgs = "PGS_outcome")
#' 
#' ## print GsensY results
#' gsens_out@external$gsensY
#' ## any other standard lavaan output also available
#' summary(gsens_out)
#' }
#'

#' @author Leonard Frach & Jean-Baptiste Pingault
#' @export
#' @import lavaan
#' @importFrom dplyr mutate across

#' @references
#' Frach, L., Rijsdijk, F., Hannigan, L. J., Dudbridge, F. & Pingault, J.-B. (in preparation).
#' Genetic sensitivity analysis: estimating genetic confounding and environmentally mediated genetic effects using multiple exposures.
#'
#' Pingault, J. B., Rijsdijk, F., Schoeler, T., Choi, S. W., Selzam, S., Krapohl, E., ... & Dudbridge, F. (2021).
#' Genetic sensitivity analysis: Adjusting for genetic confounding in epidemiological associations.
#' *PLoS genetics, 17*(6), e1009590. \doi{10.1371/journal.pgen.1009590}

gsensY = function(data = NULL,
                  sample.cov = NULL,
                  h2,
                  exposures,
                  outcome,
                  pgs,
                  model = NULL,
                  loading = FALSE,
                  lg2 = NULL,
                  ...) {
    
    #check for correlation flag
    dotdotdot <- list(...)
    correlation.flag <- FALSE
    if (!is.null(dotdotdot$correlation) && dotdotdot$correlation) {
        correlation.flag <- TRUE
    }
    
    
    # create covariance structure and label for the model
    
    covstruc <- outer(exposures,
                      exposures,
                      function(x, y) paste(x, "~~", y))
    
    # Residual (co)variances
    covstruc <- covstruc[lower.tri(covstruc, diag = !correlation.flag)]
    
    NX <- length(exposures)
    labelsa   <- paste0("a", 1:NX)
    labelsb   <- paste0("b", 1:NX)
    labelsm   <- paste0("m", 1:NX)
    labels_gc <- paste0("gc_", labelsb)
    labels_go <- paste0("go_", labelsb)
    labels_eX <- paste0("eX", 1:NX)
    
    labels_d <- c()
    
    
    ### main lavaan model ###
    
    # empty variable for genetic confounding
    gC <- c()
    
    # check if additional model arguments were passed on
    if (!is.null(model)) {
        model <- model
    } else {
        model <- NULL
    }
    
    ## model specification
    
    # check for additional causal effects
    add <- do.call(paste, c(expand.grid(exposures, exposures), sep = " ~ "))
    
    if (any(add %in% model)) {
        delta <- TRUE
    } else {
        delta <- FALSE
    }
    
    # add labels for causal paths
    if (delta) {
        add_delta <- add[add %in% model]
        interm_m <- strsplit(add_delta, "~") %>% unlist()
        model <- model[model != add_delta]
    }
    
    m_delta <- c()
    
    #if (loading == TRUE) {
        
    # Correlation of PGS and outcome
    if (is.null(data)) {
        df_cor <- cov2cor(sample.cov)
        rPGS <- df_cor[rownames(df_cor) == outcome, colnames(df_cor) == pgs]
        
    } else {
        rPGS <- as.numeric(cor(data[, outcome], data[, pgs]))
    }
    
    lg2 <- rPGS^2 / h2
    
    model <- c(
        paste(outcome, "~", paste0(c(labelsb,"c"),"*", c(exposures,"GG"), collapse = " + ")),    # Y depends on X and true polygenic score
        paste(exposures,"~", paste0(labelsa,"*","GG")),                                   # X1-Xi depend on true polygenic score
        
        # loading constraint
        ifelse(loading, paste("GG =~ ", (sqrt(lg2)),"*", pgs), paste("GG =~ lg*", pgs)), # if loading constraint should be used, add first term, otherwise add second
        #ifelse(loading, paste("GG =~ ", (sqrt(rPGS^2 / h2)),"*", pgs), ),
        
        # heritability constraint
        ifelse(!loading, paste(h2, " == (h^2)/VarY"), ""), 
        
        # variance of latent genetic factor
        "GG ~~ 1*GG",
        
        # Residual variance of PGS
        ifelse(loading, paste(pgs, " ~~ ", (1 - lg2), "*", pgs), paste(pgs, " ~~ vp*", pgs)), # Using the loading constraint, the residual variance of the pgs can be freely estimated which would produce slightly larger SEs
        #ifelse(loading, paste(pgs, " ~~ ", (1 - (rPGS^2 / h2)), "*", pgs), ), # this can be freely estimated which would produce slightly larger SEs
        
        # Residual variance of the outcome
        ifelse(!correlation.flag, paste(outcome, " ~~ ResVarY*", outcome), ""),
        
        # total mediation effect (add causal effects if specified)
        if (delta) {
            for (i in 1:length(add_delta)) {
                m_delta[i] <- gsub("~", paste0("~ d", i, " * "), add_delta[i])
            }
            
            # split string to get the labels
            x <- strsplit(m_delta, " * ") %>% unlist()
            
            for (i in 1:length(add_delta)) {
                labels_d <- c(labels_d,
                              paste0(labelsa[exposures %in% interm_m[2*i]], "*",
                                     x[5*i - 2], "*", labelsb[exposures %in% interm_m[2*i - 1]]))
            }
            
            covstruc <- covstruc[covstruc != gsub("~", "~~", add_delta)]
            
            paste("m :=", paste(paste0(labelsa,"*", labelsb, collapse = " + "), "+", paste0(labels_d, collapse = " + ")))
            
        } else {
            paste("m :=", paste0(labelsa,"*",labelsb, collapse = " + "))
            
        },
        
        # residual (co)variances 
        covstruc,
        # causal effect(s)
        m_delta,
        
        paste0(labelsm, " := ", labelsa,"*",labelsb),    # specific mediation pathways for each G->Xi->Y
        # heritability constraints
        
        # add external model specifications
        model,
        
        ## MODEL CONSTRAINTS:
        paste0("h := m + c"), # heritability definition
        ifelse(correlation.flag | !loading, paste(h2, " == (h^2)/VarY"), ""), # constraint might still be needed for correlation feature
        
        # Observed variance of Y
        if (is.null(data)) {
            paste("VarY := ", sample.cov[rownames(sample.cov) == outcome, colnames(sample.cov) == outcome]) # scaling to get the population variance not needed for now -> * (length(!is.na(data[, outcome])) - 1) / length(!is.na(data[, outcome]))),
            
        } else {
            
            paste("VarY := ", as.numeric(var(data[, outcome], na.rm = TRUE))) # scaling to get the population variance not needed for now -> * (length(!is.na(data[, outcome])) - 1) / length(!is.na(data[, outcome]))),
            
        },
        
        # genetic confounding for each Xi->Y association
        if (NX > 1) {
            for (i in 1:NX) {     
                
                gC <- c(gC, paste0(labels_gc[i], " := ",
                                   paste0(labelsa[i],"*", labelsb[-i],"*", labelsa[-i],
                                          collapse = " + "), " + ", labelsa[i], "*c"))
            }
            
            if (delta) {
                for (i in length(add_delta)) {
                    
                    gC[exposures %in% interm_m[2*i - 1]] <- paste0(gC[exposures %in% interm_m[2*i - 1]], " + ", 
                                                                   paste0(labelsa[exposures %in% interm_m[2*i]], "*",
                                                                          x[5*i - 2], "*c"))
                }
            }
            
            gC
        } else {
            
            (gC <- paste0(labels_gc, " := ",
                          paste0(labelsa, "*c")))
        },
        
        paste0(labels_go, " := ", labelsa,"*", labelsa,"*", labelsb, " + ", labels_gc) # genetic overlap for each Xi->Y association
    )
        

    ## model estimation options
    
    # run model
    fit_mod <- lavaan(model = model, data = data, sample.cov = sample.cov, ...)
    
    
    ## parameter estimates options
    pe <- parameterEstimates(fit_mod)
    
    
    results <- data.frame(rbind(
        pe[pe$label %in% labelsb,],
        pe[pe$label %in% labelsm,],
        pe[pe$label %in% "m",],
        pe[pe$label %in% labels_gc,],
        pe[pe$label %in% labels_go,]
    ))[, c(5:dim(pe)[2])]
    
    
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
    for (i in 1:NX) {
        rownames(results)[(2*NX + 1 + i)] <- c(paste("Genetic confounding Bx", i, "y", sep = ""))
    }
    
    # name the genetic overlap for each exposure-outcome association
    for (i in 1:NX) {
        rownames(results)[3*NX + 1 + i] <- c(paste("Genetic overlap x", i, "y", sep = ""))
    }
    
    results$pvalue <- format(2*pnorm(-abs(results$z)), digits = 3, scientific = TRUE)
    results <- results |> dplyr::mutate(dplyr::across(!"pvalue", \(x) round(x, 3)))

    
    ## store results in the lavaan object
    fit_mod@external$gsensY <- results
    
    return(fit_mod) # model output can be used, e.g. using summary() function
    
}


#' Adjusting for genetic confounding using PGS for the exposure
#'
#' Adjusting for genetic confounding in exposure-outcome associations using the polygenic score for the exposure.

#' @param rxy The observed phenotypic correlation between exposure X and outcome Y.
#' @param rgx The observed correlation between the polygenic score for X and exposure X.
#' @param rgy The observed correlation between the polygenic score for X and outcome Y.
#' @param n Sample size
#' @param h2 The additive genetic variance explained in exposure X under the scenario of interest.
#' @param print Optional. Enables the examination of model parameters (default = FALSE).
#' @param constrain Optional. Argument to constrain model parameters (default = NULL).

#' @return Estimates for the adjusted exposure-outcome association,
#' genetic confounding and the total effect.

#' @author Jean-Baptiste Pingault, Tabea Schoeler & Frank Dudbridge
#' @export
#' @import lavaan
#' @importFrom dplyr mutate_all

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
  mat <- matrix(c(1, rgx, rgy, rgx, 1, rxy, rgy, rxy, 1), ncol = 3, nrow = 3)
  colnames(mat) <- c("G","X","Y")
  rownames(mat) <- c("G","X","Y")

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


  fit1 <- lavaan(model1, sample.cov = mat, sample.nobs = n, estimator = "GLS")
  if (print) {  summary(fit1)}
  
  pe <- parameterestimates(fit1)
  pe$pvalue <- formatC(2*pnorm(-abs(pe$z)), digits = 5)
  results <- data.frame(rbind(
    pe[pe$label == "bxy",],
    pe[pe$label == "conf",],
    pe[pe$label == "total",]
  ))[,5:10]
  results <- dplyr::mutate_all(results, round, 3) # round all numeric variables

  rownames(results) <- c("Adjusted Bxy","Genetic confounding","Total effect")
  results
}


## Two polygenic scores ##

#' Adjusting for genetic confounding using PGS for the exposure and the outcome
#'
#' Adjusting for genetic confounding in exposure-outcome associations using the polygenic score for the exposure and the outcome.

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

#' @return Estimates for the adjusted exposure-outcome association,
#' genetic confounding and the total effect.

#' @author Jean-Baptiste Pingault, Tabea Schoeler & Frank Dudbridge
#' @export
#' @import lavaan
#' @importFrom dplyr mutate_all

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

  lower <- c(
    1,
    rg1g2,1,
    rg1x,rg2x,1,
    rg1y,rg2y,rxy,1)

  mat <- getCov(lower, names = c("G1","G2","X","Y"))

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

  fit1 <- lavaan(model1, sample.cov = mat, sample.nobs = n, estimator = "GLS")
  if (print) {  summary(fit1)}
  
  pe <- parameterestimates(fit1)
  pe$pvalue <- formatC(2*pnorm(-abs(pe$z)), digits = 5)
  results <- data.frame(rbind(
    pe[pe$label == "bxy",],
    pe[pe$label == "conf",],
    pe[pe$label == "total",]
  ))[, 5:10]

  results <- dplyr::mutate_all(results, round, 3) # round all numeric variables

  rownames(results) <- c("Adjusted Bxy","Genetic confounding","Total effect")
  results

}



