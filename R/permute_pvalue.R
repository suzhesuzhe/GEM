
#' Permuted P Value Calculation
#' @name permute_pvalue
#' 
#' @description Calculation of permuted p value for a fitted GEM.
#'
#' @param dat Dataframe with first column as the treatment index, second column as response, and the 
#' remaining columns as covariates design matrix. 
#' @param method Choice of criterior on which the GEM is based. This can be a string in 
#' \code{c("nu","de","F")}, which corresponde to the numerator, denominator and F-statistics
#' methods respectively.
#' 
#' 
#' @return \code{eff_size} the calculated permuted p value for the data and choosen criterior 
#' 
#' 
#' @examples
#' #constructing the covariance matrix
#' co <- matrix(0.2, 10, 10)
#' diag(co) <- 1
#' dataEx <- data_generator1(d = 0.3, R2 = 0.5, v2 = 1, n = 3000, co = co, beta1 = rep(1,10),inter = c(0,0))
#' #fit the GEM
#' dat <- dataEx[[1]]
#' pvalue_permuteCpp(dat,"nu")
#' @export


permute_pvalue <- function(dat,method)
{
    
    colnames(dat)[1] <- "trt"
	K <- length(unique(dat$trt))
    dat[,1] <- factor(dat[,1],labels=c(1:K))
    N <- ddply(dat, .(trt), nrow)[,2]
    mm <- gem_fit(dat,type)
    p <- mm[[3]]
  
    for(i in 1:1000)
    {
        datt <- dat
    	datt[,1] <- sample(dat$trt)
        mmm <- gem_fit(datt,type)
        p <- c(p, mmm[[3]])        
    }
    p_value <- sum(p[2:1001] <= p[1])/1000
    return(list("pvalue"=p_value,
    	         "p"=p))
    
}
