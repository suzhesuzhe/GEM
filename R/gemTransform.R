
gemTransform <- function(addingZObject)
{
	
    dat <- addingZObject[[2]]   #This selected uncentered designMatrix determine that this eta criteria is uncentered
    dat[,1] <- factor(dat[,1])
    mod <- lm(y~trt+Z+trt:Z,data=dat)
    #mod <- glm(y~trt+Z+trt:Z,data=dat,family = binomial)


	gamma_gem <- c(mod[[1]][2],mod[[1]][4])
    eta_gem <- c(gamma_gem[1],gamma_gem[2]*addingZObject[[1]])
	signZ <- Sign(mod[[1]][4])
    eta_gem <- eta_gem/sqrt(sum(eta_gem^2))
	intersectGem <- -mod[[1]][2]/mod[[1]][4]    ###warning, this constraint is indeed for the uncentered Z
	
	
    result <- list(
    	"gamma_gem" = gamma_gem,
    	"eta_gem" = eta_gem,
    	"intersectGem" = intersectGem)
	return(result)
}
