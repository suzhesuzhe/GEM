

alphaGem <- function(dat, type)
{
	constructorObject <- dataConstructor(dat)
	K=constructorObject[[1]]
	p=constructorObject[[2]]
	N=constructorObject[[3]]
	centeredYList=constructorObject[[4]]
	centeredYVector=constructorObject[[5]]
	centeredXList=constructorObject[[8]]
	centeredXFrame=constructorObject[[9]]
	orderedTrt=constructorObject[[12]]
	
    co <- cov(centeredXFrame)
    e <- Eigen(co)
    sqrtco <- e$vectors%*%diag(sqrt(e$values))%*%t(e$vectors)
    pro <- N/sum(N)
    pro1<-(N-1)/sum(N)
    Beta <- vector("list",K)
    for(i in 1:K)
    {
        Beta[[i]] <- as.matrix(lm(centeredYList[[i]]~centeredXList[[i]])$coeff[-1])
    }
    
    Beta.bar<-Reduce('+',lapply(1:K, function(j){pro[j]*Beta[[j]]}))
    B<-Reduce('+',lapply(1:K, function(j)
    	{
    		pro[j] * (Beta[[j]]- Beta.bar) %*% t(Beta[[j]] - Beta.bar)
    	}))
    # B is the matrix for the numerator method   
    D<-Reduce('+',lapply(1:K, function(j){pro1[j] * (Beta[[j]]) %*% t(Beta[[j]])}))
    # D is the matrix for the denominator method   
    A<-Reduce('+',lapply(1:K, function(j)
    	{
    		Reduce('*',(pro1[j]/(N[j]-1))*t(centeredYList[[j]])%*%centeredYList[[j]],diag(p))
    	}))-sqrtco%*%D%*%sqrtco
    # A is the matrix in the denominator for the F method  

    # maximizing the interaction term (nu)
    if(type=="nu")
    {                
        astar <- Eigen(sqrtco%*%B%*%sqrtco)$vector[,1]
        alpha_nu <- solve(sqrtco)%*%astar
        alpha_nu_2 <- Sign(Re(alpha_nu[1,1]))*Re(alpha_nu)
        alpha<-alpha_nu_2
    }
    
    # minimizing the error sum of squares (de)
    if(type=="de")
    {        
        astar <- Eigen(sqrtco%*%D%*%sqrtco)$vector[,1]
        alpha_de <- solve(sqrtco) %*% astar
        alpha_de_2 <- Sign(Re(alpha_de[1,1]))*Re(alpha_de)
        alpha<-alpha_de_2
    }
    
    # maximizing the F ratio (F)    
    if(type=="F")
    {
        astar<- Eigen(solve(A)%*%sqrtco%*%B%*%sqrtco)$vector[,1]
        alpha_F <- solve(sqrtco)%*%astar
        c <- sqrt(as.numeric(t(alpha_F)%*%co%*%alpha_F))
        alpha_F <- Re(alpha_F)/c
        
        alpha_F_2 <- Sign(Re(alpha_F[1,1])) * Re(alpha_F)
        alpha<-alpha_F_2
    }
    
	Z <- as.matrix(centeredXFrame) %*% alpha
	
	fit.full <- lm(centeredYVector~orderedTrt+Z+orderedTrt:Z)
	fit.reduced <- lm(centeredYVector~orderedTrt+Z)
	p_value <- anova(fit.full,fit.reduced)[2,6]
	
	if(K!=2) 
	{
		effect_size <- NA
	}else{
			effect_size <- abs(effectSize(centeredYVector,orderedTrt,Z))
	}
	results <- list("type"= type,      #1 The method for the gem implementation
                    "alpha"=alpha,
	  				"p_value" = p_value,
	  				"effect.size"=effect_size)  
	return(results)	
	
}
