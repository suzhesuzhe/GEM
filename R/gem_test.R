
#' Implement Fitted GEM criterior on a Data Set
#' @name gem_test
#' 
#' @description Calculates the value of the treatment decision based on the information from a fitted GEM model.The information
#' is contained in the \code{gemObject}, which is obtained from the \code{gem_fit} function. With \code{gem_test_sample} the value 
#' of decision based on a GEM model is obtained for a test sample that for each subject has observed outcome under only one treatment 
#' condition (this would be the situation when the test sample is a "real" data set). With \code{gem_test_simsample} the value 
#' of decision is calculated when the test sample has the outcome under all treatment conditions for all subjects
#' (this would be the situation when simulated data is used).
#'
#' @param dat Data frame with first column as the treatment index, second column as outcome, and the 
#' remaining columns as covariates design matrix. The treatment index could only have two values ana the outcome should be of continuous type.
#' @param gemObject A list containing the fitted GEM information, which could be the second element \code{gemObject} 
#' of the ouput from the \code{gem_fit} function or a list with the same structure.
#' @param XFrame Design matrix of the predictors in the simulated sample with known outcomes under both treatment conditions
#' @param y0 Outcome vector for all subjects under the first treatment assignment
#' @param y1 Outcome vector for all subjects under the second treatment assignment

#' 
#' @details The treatment decision rule estimated by the \code{gem_fit} function can be applied to a new (real) data set to estimate its value.  
#' It can also be applied to a simulated data set, where the outcome is known under both conditions, to study its performance of such treatment decision rule.  
#' These two functions correspond to those two situations and compute the population average benefit (or the value of the decision rule based on a GEM model).
#'  
#' @return \code{PAB_gem} Population average benefit of a the treatment regime based on a GEM model
#' @return \code{PAB_unres} Population average benefit of a treatment regime based on an unrestricted linear model
#' @return \code{opt_gem} The optimal treatment assignment for each subject

#' @examples
#' #constructing the covariance matrix
#' co <- matrix(0.2, 10, 10)
#' diag(co) <- 1
#' dataEx <- data_generator1(d = 0.3, R2 = 0.5, v2 = 1, n = 300, 
#'                          co = co, beta1 = rep(1,10),inter = c(0,0))
#' #fit the GEM
#' dat <- dataEx[[1]]
#' model_nu <- gem_fit(dat = dat, method = "nu")
#' 
#' #calculate the population average benefit in the data sample
#' gem_test_sample(dat,model_nu[[2]])
#' 
#' #calculate the population average benefit when outcome under both treatment conditions 
#' #is known, usually in a simulated sample
#' bigData <- data_generator3(n = 1000,co = co,bet =dataEx[[2]], inter = c(0,0)) 
#' gem_test_simsample(bigData[[1]],bigData[[2]],bigData[[3]],model_nu[[2]])
#' @export
gem_test_sample<- function(dat, gemObject)
{

	dat[,1] <- factor(dat[,1],labels=c(0,1))
	XFrame <- dat[,c(-1,-2)]
	
	alpha <- gemObject$alpha
	gamma0 <- gemObject$GammaByTrt[1,]
	gamma1 <- gemObject$GammaByTrt[2,]
	 
	#restricted model
	#here we use XFrame design Matrix without sweeping the mean due to the meanning of those coefficient
	optTrt_gem <- as.numeric(gamma1[1] + as.matrix(XFrame) %*% as.matrix(gamma1[2]*alpha) >=
	        gamma0[1] + as.matrix(XFrame) %*% as.matrix(gamma0[2]*alpha)) 
	
	PAB_gem <- sum(dat[dat[,1]==optTrt_gem,2])/sum(dat[,1]==optTrt_gem)

	#unrestricted model
	unresCoeff0 <- gemObject$betaByTrt[1,]
	unresCoeff1 <- gemObject$betaByTrt[2,]

	#here we use XFrame design Matrix without sweeping the mean due to the meanning of those coefficient
	optTrt_unres <-  as.numeric(as.matrix(cbind(1,XFrame)) %*% as.matrix(as.numeric(unresCoeff1)) >= 
		            as.matrix(cbind(1,XFrame)) %*% as.matrix(as.numeric(unresCoeff0)))
	
	PAB_unres <-  sum(dat[dat[,1]==optTrt_unres,2])/sum(dat[,1]==optTrt_unres)

	result <- list(
		"PAB_gem" = PAB_gem,
		"PAB_unres" = PAB_unres,
		"opt_gem" = optTrt_gem)
	return(result)
	
}

#' @name gem_test
#' @export
gem_test_simsample <- function(y0, y1,XFrame, gemObject)
{

	alpha <- gemObject$alpha
	gamma0 <- gemObject$GammaByTrt[1,]
	gamma1 <- gemObject$GammaByTrt[2,]

	#here we use XFrame design Matrix without sweeping the mean due to the meanning of those coefficient
	optTrt_gem <- as.numeric(gamma1[1] + as.matrix(XFrame) %*% as.matrix(gamma1[2]*alpha) >=
	        gamma0[1] + as.matrix(XFrame) %*% as.matrix(gamma0[2]*alpha))
	
	PAB_gem <- mean(y0*(1-optTrt_gem) + y1 * optTrt_gem)

	########################
	unresCoeff0 <- gemObject$betaByTrt[1,]
	unresCoeff1 <- gemObject$betaByTrt[2,]
	
	#here we use XFrame design Matrix without sweeping the mean due to the meanning of those coefficient
	optTrt_unres <-  as.numeric(as.matrix(cbind(1,XFrame)) %*% as.matrix(as.numeric(unresCoeff1)) >= 
		            as.matrix(cbind(1,XFrame)) %*% as.matrix(as.numeric(unresCoeff0))) 
	
	PAB_unres <- mean(y0*(1-optTrt_unres) + y1 * optTrt_unres)

	result <- list(
		"PAB_gem" = PAB_gem,
		"PAB_unres" = PAB_unres,
		"opt_gem" = optTrt_gem)
	return(result)
	
}


