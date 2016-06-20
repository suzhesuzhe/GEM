
#' Implement Fitted GEM criterior on a Data Set
#' @name gem_test
#' 
#' @description Implement a fitted GEM on a testing dataset, calculate the population average benefit
#'
#' @param dat Dataframe with first column as the treatment index, second column as response, and the 
#' remaining columns as covariates design matrix. 
#' @param gemObject Object as the second element of the ouput list from gem_fit function
#' @param XFrame Design matrix for the predictors of the out-of-sample data set
#' @param y0 Response for each observation under the first treatment assignment
#' @param y1 Response for each observation under the second treatment assignment

#' 
#' @details A GEM fit provides a treatment assignment strategy for two types of data: each subject has only 
#' one observation of response under certain treatment or has observations of response under both treatments. 
#' These two gem_test functions are both used to calculate the population average benefit, however, they are designed
#' for those two types of data respectively.
#'  
#' @return \code{PAB_gem} Population average benefit based on the treatment regime from the GEM fit
#' @return \code{PAB_unres} Population average benefit based on on restricted linear model
#' 
#' @examples
#' #constructing the covariance matrix
#' co <- matrix(0.2, 10, 10)
#' diag(co) <- 1
#' dataEx <- data_generator1(d = 0.3, R2 = 0.5, v2 = 1, n = 3000, co = co, beta1 = rep(1,10),inter = c(0,0))
#' #fit the GEM
#' dat <- dataEx[[1]]
#' model_nu <- gem_fit(dat = dat, method = "nu")
#' 
#' #calculate the population average benefit for dat itself
#' gem_test_insample(dat,model_nu[[2]])
#' 
#' #calculate the population average benefit for a second type data
#' bigData <- data_generator3(n = 10000,co = co,bet =dataEx[[2]], inter = c(0,0)) 
#' gem_test_outsample(bigData[[1]],bigData[[2]],bigData[[3]],model_nu[[2]])
#' @export

#' @rdname gem_test_insample
#' @export
gem_test_insample<- function(dat, gemObject)
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
		"PAB_unres" = PAB_unres)
	return(result)
	
}

#' @rdname gem_test_outsample
#' @export
gem_test_outsample <- function(y0, y1,XFrame, gemObject)
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
		"PAB_unres" = PAB_unres)
	return(result)
	
}


