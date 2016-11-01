
#' GEM Fit
#' @name gem_fit
#' @import plyr
#' @import ggplot2
#' @import RcppArmadillo
#' @useDynLib GEM
#' @importFrom Rcpp sourceCpp
#' 
#' @description The main algorithm in GEM package for calculating the coefficients of the linear combination of the covariates
#' to generate a GEM. See more detail in \cite{E Petkova, T Tarpey, Z Su, and RT Ogden. Generated effect modifiers (GEMs) in randomized clinical trials. Biostatistics, (First published online: July 27, 2016). doi: 10.1093/biostatistics/kxw035.}
#'
#' @param dat Data frame with first column as the treatment index, second column as outcome, and the 
#' remaining columns as covariates design matrix. 
#' @param method Choice of criterior on which the GEM is based. This can be a string in 
#' \code{c("nu","de","F")}, which corresponde to the numerator, denominator and F-statistics
#' methods respectively. The default method is the F-statistics method.

#' @return 
#' \enumerate{
#' 		\item \code{method} An indicator for which method should be used for fitting the GEM 
#' 		\item \code{gemObject} The calculated weights for combining the covariates
#' 		\item \code{p_value} The p-value for the interaction term in model \eqn{Y  = a + trt + Z + trt*Z + \epsilon}
#' 		\item \code{Augmented_Data} The input data argumented with the combined moderator as the last column 
#' 		\item \code{effect.size} The effect size of the GEM if there are only two treatment groups
#' 		\item \code{plot} A scatter plot of Y versus the GEM Z with fitted lines and grouped by treatment
#' 		
#' }
#' 
#' 
#' @examples
#' #constructing the covariance matrix
#' co <- matrix(0.2, 10, 10)
#' diag(co) <- 1
#' dataEx <- data_generator1(d = 0.3, R2 = 0.5, v2 = 1, n = 3000,
#'                         co = co, beta1 = rep(1,10),inter = c(0,0))
#' #fit the GEM
#' dat <- dataEx[[1]]
#' model_nu <- gem_fit(dat = dat, method = "nu")
#' model_de <- gem_fit(dat = dat, method = "de")
#' model_F <- gem_fit(dat = dat, method = "F")
#' @export
gem_fit <- function(dat, method = "F") 
{
       colnames(dat)[1:2] <- c("trt","y")
       k <- length(unique(dat$trt))
       dat_list <- dlply(dat,.(dat$trt),function(x){as.matrix(x[,-1])})

       gemObject <- gemCpp(dat_list, method)
       dat$Z <- as.matrix(dat[,c(-1,-2)]) %*% as.matrix(gemObject[[1]])
       
       mod <- fastLm(y~factor(trt)*Z,dat = dat)
       p_value <- summary(mod)[[2]][4,4]
       p <- ggplot(dat, aes_string(x = "Z", y = "y", group = factor("trt"), color = factor("trt")))+
		       geom_point() + geom_smooth(method="lm",fullrange = T) 
        if ( k==2)
        {
            es <- effectSize(dat$y, dat$trt, dat$Z)
        }else
        {
            es <- NA
        }
	   results <- list("method"= method,
                       "gem"=gemObject,
                       "p_value" = p_value,
                       "Augmented_Data" = dat,
                       "effectSize" = es,
                       "plot" = p)
	return(results)	
	
}
