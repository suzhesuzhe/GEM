
#' Implement gem criteria on a data set
#' @name test
#' @import ggplot2
#' 
#' @description Implement a fitted GEM on a testing dataset, calculte the combined moderator, the population 
#' average benefit, the effect size and plot the combined moderator interacting with the treatment.
#'
#' @param dat A dataframe with first column as the treatment index, second column as response, and the 
#' remaining columns as predictors design matrix. Here you should know the number of predictors in the data
#' set that used to fit the \code{trainingObject} and the \code{dat} should have same number of predictors.
#' @param trainingObject A list that return from \code{train}
#' @param XFrame Design matrix for the predictors of the out-of-sample data set
#' @param y0 Response for each observation under the first treatment assignment
#' @param y1 Response for each observation under the second treatment assignment

#' 
#' @details These are two functions for validating a gem criteria on a data set, for both 
#' in-sample and out-sample types respectively. 
#' 
#' The \code{test_inSample} take the data set which
#' only one reponse is available for each observation and the population average benefit is calculated
#' as the mean of responses whose pre-assigned treatment coincides with the optimal treatment.
#' 
#'  \code{test_outSample} validate on a data set that two responses for both treatments are available, 
#' and for this function the only output is the population average benefit, which is the mean of the 
#' response under optimal treatment assignmnet
#' 
#' The effect size is calulated in terms of the combined moderator.
#' 
#' @return Output of these two functions for the data generator under Gem case is not the same,
#' and you can find explanation for each of them in the following:
#' #' \enumerate{
#' 		\item \code{PAB_gem} Population average benefit, see details
#' 		\item \code{es_gem} Effect size of the combined moderator, see details
#' 		\item \code{argumentedData} The input dataframe plus one more column as the combined moderator
#' 		\item \code{combinedModeratorPlot} A ggplot2 object that plots the combined moderator versus response 
#' 		with two treatment group marked with different colors
#' }

#' @examples
#' numb_p <-5
#' bet <- vector("list",2)
#' bet[[1]] <- c(1,1,1,1,1)
#' bet[[2]] <- c(1,2,3,4,5)
#' co <- matrix(0.2,numb_p,numb_p)
#'       diag(co) <- 1
#'
#' popu1 <- nongemData_inSample(200,co,0.1,bet,c(0,1))
#' popu2 <- nongemData_outSample(20000,co,0.2,bet,c(0,1))
#'
#' mod1 <- train(popu1[[1]],"nu")
#'
#' test1 <- test_inSample(popu1[[1]],mod1)
#' test2 <- test_outSample(popu2[[3]],popu2[[1]],popu2[[2]],mod1)
#' @export
test_inSample <- function(dat, trainingObject) #the cri is centered or not can be checked by the third element of criApObject
{

	constructorObject <- dataConstructor(dat)
	K=constructorObject[[1]]
	p=constructorObject[[2]]
	N=constructorObject[[3]]
	YList=constructorObject[[6]]
	YVector=constructorObject[[7]]
	XList=constructorObject[[8]]
	XFrame=constructorObject[[9]]
	orderedTrt=constructorObject[[12]]
	
	
	optTrt_gem <- apply(cbind(1,XFrame),1,function(x){g(trainingObject[[3]][[2]],x)})
	PAB_gem <- sum(YVector[orderedTrt==optTrt_gem])/sum(orderedTrt==optTrt_gem)

	datt <- addingZ(dat, trainingObject[[1]])[[2]]
	
	es_gem <- abs(effectSize(datt$y,datt$trt,datt$Z))

	p <- ggplot(datt)+geom_point(aes(x=Z,y=y,color=factor(trt)))+
		 geom_smooth(aes(x=Z,y=y,color=factor(trt)),method="lm",fullrange=T)+
		 geom_vline(xintercept = trainingObject[[3]][[3]],colour="blue", linetype = "longdash")+
		 xlab("Combined Moderator")+ylab("")+ggtitle(paste("Method = ",trainingObject[[1]][[1]]))
	
	result <- list(
		"PAB_gem" = PAB_gem,
		"es_gem" = es_gem,
		"argumentedData" = datt,
		"combinedModeratorPlot" = p)
	return(result)
	
}

#' @rdname test
#' @export
test_outSample <- function(XFrame, y0, y1, trainingObject) #the cri is centered or not can be checked by the third element of criApObject
{

	
	optTrt_gem <- apply(cbind(1,XFrame),1,function(x){g(trainingObject[[3]][[2]],x)})
	PAB_gem <- sum(y0*(1-optTrt_gem)+y1*optTrt_gem)/length(optTrt_gem)

	result <- list(
		"PAB_gem" = PAB_gem)
	return(result)
	
}

