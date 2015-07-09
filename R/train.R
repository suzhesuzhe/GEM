#' fit gem on certain data set
#' @name train
#' 
#' @description The main algorithm for calculating a weight for a list of moderators
#' to generate a optimal combined moderator based on three different criterias.
#'
#' @param dat A dataframe with first column as the treatment index, second column as response, and the 
#' remaining columns as predictors design matrix. 
#' @param type a description of the method to be used in the model. This can be a character string in 
#' \code{c("nu","de","F")}, which are corresponding to the numerator, denominator and F-statistics
#' methods respectively.
#' 
#' 
#' @details These are two functions for simulating data set under nonGem scenario, for both 
#' in-sample and out-sample types respectively. The data is simulated as a linear combination of 
#' certain predictors, which are multivariate Guassian distributed and have preset coefficients as the 
#' functions' argument. Since they are designed for simulating data to test with calculting gem, 
#' those data generators take the simplest setting for two treatment groups with same number 
#' of observations and same proportion of explained variance \eqn{R^2}. 
#' 
#' The \code{nongemData_inSample} will generate only one reponse for each observation under certain treatment 
#' group, while \code{nongemData_outSample} will generate two responses under both of the treatments. In general 
#' \code{nongemData_inSample} is used to simulate a sample data to train the gem algorithm and \code{nongemData_outSample}
#' is used for simulating the population.
#' 
#' As the gem algorithm is based on the linear model coefficients and the number of subjects
#' in each treatment group, see\url{http://en.wikipedia.org/wiki/Integer_overflow}, 
#' the function \code{gemData_inSample} calculates the theoretical value of the alpha under three gem methods. 
# 
#' @return Output of the functions is a list of length three and each of them contain some key aspect of 
#' information about the gem fit
#' 
#' For the first list component \code{gem}
#' \enumerate{
#' 		\item \code{type} A indicator showing which method is used for fitting the GEM 
#' 		\item \code{alpha} The calculated weight for the sequnce of moderators
#' 		\item \code{p_value} The p value for the interaction term when fitting the reponse
#' 		versus treatment interacting the combined moderator
#' 		\item \code{effect.size} The effect size of the combined moderator
#' }
#' 
#' For the second list component, \code{permuted_pvalue} is a adjust p value for the interaction term
#' see also \url{www.rblogger.com}
#' 
#' For the third component, they are other versions of the calculated alpha weight, to be more specific 
#' \enumerate{
#' 		\item \code{gamma_gem} A vector for the intercept term and the combined moderator. The linear
#' 		combination of these two terms by the vector will be a decision rule for treatment assignment
#' 		\item \code{eta_gem} A vector for the intercept term and the original sequence of moderators. The linear
#' 		combination of those terms by the vector will be a decision rule for treatment assignment
#' 		\item \code{intersectGem} A number indicting the decision rule of the combined moderator, larger or
#' 		smaller than this number will result in selecting certain treatment assignment
#' }
#' 
#' @examples
#' bet <- vector("list",2)
#' bet[[1]] <- c(1,1,1,1,1)
#' bet[[2]] <- c(1,2,3,4,5)
#' co <- matrix(0.2,5,5)
#' diag(co) <- 1
#' popu1 <- nongemData_inSample(200,co,0.1,bet,c(0,1))
#' mod_nu <- train(popu1[[1]],"nu")
#' @export
train <- function(dat,type)
{
	gemObject <- alphaGem(dat, type)
	
	colnames(dat)[1] <- "trt"
	K = length(unique(dat[,1]))
    dat[,1] <- factor(dat[,1],labels=0:(K-1))
    p <- gemObject[[3]]
  
    for(i in 1:1000)
    {
        datt <- dat
    	datt[,1] <- factor(sample(dat[,1],length(dat[,1]),replace=F))
        mmm <- alphaGem(datt, type)
        p <- c(p, mmm[[3]])
        
    }
    p_value <- sum(p[2:1001] <= p[1])/1000

	
	addingZObject <- addingZ(dat, gemObject)
	etaAndGamma <- gemTransform(addingZObject)

	return(list(
		"gem" = gemObject,
		"permuted_pvalue"=p_value,
		"etaAndGamma" = etaAndGamma))	
}


