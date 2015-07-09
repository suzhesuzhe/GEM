
#' @import plyr
dataConstructor <- function(dat)
{
	colnames(dat)[1] <- "trt"    

	K <- length(unique(dat$trt))           #group categories
    p <- ncol(dat)-2            #p: the number of the predictors  
    N <- ddply(dat, .(trt), nrow)[,2]    #a vector recording the number of observation for each group

	centeredYList <- dlply(dat, .(trt), function(x)
    	{
    		apply(x[,2:ncol(dat)], 2, function(y){y-mean(y)})[,1]
    	}) 
	centeredYVector <- unlist(centeredYList) 
    
	uncenteredYList <- dlply(dat, .(trt),function(x){x[,2]})
	uncenteredYVector <- unlist(uncenteredYList)


    centeredXList <- dlply(dat, .(trt),function(x)
    	{
    		apply(x[,2:ncol(dat)], 2, function(y){y-mean(y)})[,2:(p+1)]
    	})  
	centeredXFrame <- ldply(centeredXList, function(x){x})[,-1]
    
	uncenteredXList <-dlply(dat, .(trt),function(x){x[,c(-1,-2)]})
    uncenteredXFrame <- ldply(uncenteredXList, function(x){x})[,-1]
	
    orderedTrt <- factor(ldply(centeredXList, function(x){x})[,1])
    
	result <- list(
		"K"=K,
		"p"=p,
		"N"=N,
		"centeredYList"=centeredYList,
		"centeredYVector"=centeredYVector,
		"uncenteredYList"=uncenteredYList,
		"uncenteredYVector"=uncenteredYVector,
		"centeredXList"=centeredXList,
		"centeredXFrame"= centeredXFrame,
		"uncenteredXList"=uncenteredXList,
		"uncenteredXFrame"=uncenteredXFrame,
		"orderedTrt"=orderedTrt)
	return(result)
}#here the ouput trt is a factor

