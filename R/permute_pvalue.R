
#' Calculation of permutation p-value 
#' @name permute_pvalue
#' 
#' @description Calculates the permutation p-value for a fitted GEM. See more detail in \cite{E Petkova, T Tarpey, Z Su, and RT Ogden. Generated effect modifiers (GEMs) in randomized clinical trials. Biostatistics, (First published online: July 27, 2016). doi: 10.1093/biostatistics/kxw035.}
#'
#'
#' @param dat Data frame with first column as the treatment index, second column as the outcome, and the 
#' remaining columns as the covariates design matrix. The elements of the treatment index take \eqn{K} distinct values, where \eqn{K} is the number of treatment groups. 
#' The outcome has to be a continuous variable.
#' @param permuteN Number of permutation 
#' @param method Choice of the criterion that the generated effect modifier optimizes. This is a string in 
#' \code{c("nu","de","F")}, which corresponde to the numerator, denominator and F-statistics
#' criteria respectively. The default method is the F-statistics method.
#' 
#' 
#' @return \code{perm_p} Permutation p-value for the data and choosen criterior 
#' @return \code{p} A vector of calculated p-value for the original and permuted data set under the choosen criterior 
#' 
#' 
#' @examples
#' #constructing the covariance matrix
#' co <- matrix(0.2, 10, 10)
#' diag(co) <- 1
#' #simulate a data set
#' dataEx <- data_generator1(d = 0.3, R2 = 0.5, v2 = 1, n = 300,
#'                         co = co, beta1 = rep(1,10),inter = c(0,0))
#' #calculate the permuted p value
#' dat <- dataEx[[1]]
#' permute_pvalue(dat, permuteN = 200, method = "nu")
#' @export


permute_pvalue <- function(dat, permuteN, method = "F")
{
    
    colnames(dat)[1] <- "trt"
	K <- length(unique(dat$trt))
    dat[,1] <- factor(dat[,1],labels=c(1:K))
    N <- ddply(dat, .(dat$trt), nrow)[,2]
    mm <- gem_fit(dat,method)
    p <- mm[[3]]
  
    for(i in 1:permuteN)
    {
        datt <- dat
    	datt[,1] <- sample(dat$trt)
        mmm <- gem_fit(datt,method)
        p <- c(p, mmm[[3]])        
    }
    p_value <- sum(p[2:(1+permuteN)] <= p[1])/permuteN
    return(list("perm_p"=p_value,
    	         "p"=p))
    
}
