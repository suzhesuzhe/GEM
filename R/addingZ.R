
addingZ <- function(dat, alphaGemObject)
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
	
	datt <- cbind("trt"=factor(orderedTrt), "y"=YVector, XFrame, "Z"=as.matrix(XFrame) %*% alphaGemObject[[2]])
	
	result <- list(
		"alpha"=alphaGemObject[[2]],
		"argumentedData"=datt)
	return(result)
}