#dhueduh
true_alpha <- function(co, #the p by p covariance matrix 
                   bet, #a list specifying the beta coefficients for each treatment group
                   pro, #pi k n_k/N#pi k star (n_k-1)/N
                   pro1=NULL, #pi k star (n_k-1)/N
                   r=NULL,
                   type# the % explained variance by the predictors (r^2), which is needed for \sigma_y
)
{
    p <- ncol(co)
    e <- Eigen(co)
    K <- length(bet)
    sqrtco <- e$vectors%*%diag(sqrt(e$values))%*%t(e$vectors)
    
    Beta_bar <-  Reduce('+',lapply(1:K, function(j){pro[j]*bet[[j]]}))
    B <- sqrtco %*% Reduce('+',lapply(1:K, function(j){pro[j]*(bet[[j]]- Beta_bar) %*% t(bet[[j]] - Beta_bar)})) %*% sqrtco
    
    if (type=="nu")
    {
        astar <- Eigen(B)
        alpha <- solve(sqrtco)%*%astar$vector[,1]
        alpha_2 <- Sign(Re(alpha)[1,1])*Re(alpha)
    }
    if (type=="de")
    {
        
        mm <- vector("list",K)
        for(i in 1:K) mm[[i]] <- pro1[i] * bet[[i]] %*% t(bet[[i]])
        D <- Reduce('+',lapply(1:K, function(j){mm[[j]]}))
        
        astar <- Eigen(sqrtco%*%D%*%sqrtco)$vector[,1]
        alpha <- solve(sqrtco) %*% astar
        alpha_2 <- Sign(Re(alpha[1,1]))*Re(alpha)
    }
    if (type=="F")
    {
        sigmay <- sapply(1:K, function(j) {(t(bet[[j]]) %*% co %*% bet[[j]])/r})
        
        A <- Reduce('+',lapply(1:K, function(j){Reduce('*',pro1[j]*sigmay[j],diag(p)) - pro1[j] * sqrtco%*%bet[[j]]%*%t(bet[[j]])%*%sqrtco}))
        
        astar<-Eigen(solve(A)%*%B)$vectors[,1]
        alpha <- solve(sqrtco)%*%astar
        # re-scale alhpa3 so that it satisfies the constraint alpha'*psix*alpha=1
        c <- sqrt(as.numeric(t(alpha)%*%co%*%alpha))
        #    alpha <- Re(alpha)/c
        alpha_2 <- Sign(Re(alpha[1,1])) * Re(alpha)
    }
        
    if (p==2) 
    {
        if(alpha_2[2,1]>=0) angle <- atan(alpha_2[2,1]/alpha_2[1,1])
        if(alpha_2[2,1]<0) angle <- atan(alpha_2[2,1]/alpha_2[1,1])+pi
    }
    if(p!=2) angle <- NA
    
    result <- list("true_alpha" = alpha_2,
                   "true_angle" = angle)
    return(result)
}


