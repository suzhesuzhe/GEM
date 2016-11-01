#' Effect Size Calculation
#' @name effectSize
#' 
#' @description Calculates the effect size of a moderator when there are only two treatment groups. More details please see:
#' \cite{Kraemer, H. C. (2013). Discovering, comparing, and combining moderators of treatment on outcome after randomized clinical trials: a parametric approach. Statistics in medicine, 32(11), 1964-1973.}
#'
#' @param response A vector giving the outcome for all subjects
#' @param treatment A vector giving the treatment group index for all subjects
#' @param moderator A vector giving the moderator
#'
#' @return \code{eff_size} the calculated effect size for the moderator 
#' 
#' 
#' @examples
#' #constructing the covariance matrix
#' co <- matrix(0.2, 10, 10)
#' diag(co) <- 1
#' dataEx <- data_generator1(d = 0.3, R2 = 0.5, v2 = 1, n = 3000,
#'                     co = co, beta1 = rep(1,10),inter = c(0,0))
#' #fit the GEM
#' dat <- dataEx[[1]]
#' model_nu <- gem_fit(dat = dat, method = "nu")
#' augmentData <- model_nu[[4]]
#' es <- effectSize(augmentData$y, augmentData$trt, augmentData$Z) 
#' #this should be the same with effect size calculated by the gem_fit function
#' @export


effectSize <- function(response, treatment, moderator)
{
    mod_scale <- scale(moderator, center=TRUE, scale=TRUE)
    treat_recode <- as.numeric(as.character(factor(treatment, labels=c(-0.5, 0.5))))
    mm <- lm(response~treat_recode + mod_scale + treat_recode * mod_scale)
    
    eff_size <- as.numeric(mm[[1]][4]/2/sqrt(summary(mm)[[6]]^2+mm[[1]][3]^2+mm[[1]][4]^2/4))
    
    return(eff_size)
    
}