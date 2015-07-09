
angle <- function (beta1, beta2)
{
	theta <- acos(sum(beta1*beta2)/sqrt(sum(beta1^2)*sum(beta2^2)))
	return(theta)
}

metric <- function(x, y)
{
    distance <- sqrt(sum((x-y)^2))
    return(distance) 
}

Sign <- function(x)
{
    if(x>=0) sign <- 1
    if(x < 0) sign <- -1
    return(sign)
}

Eigen <- function(dat)
{
    m1 <- eigen(dat)$values
    m2 <- eigen(dat)$vectors
    m2 <- apply(m2, 2, function(x){x * Sign(Re(x[1]))})
    results <- list("values"=m1,
                    "vectors"=m2)
    return(results)
}

g <- function (eta, x)
{
    if(t(eta) %*% x >=0) assignment <- 1
    else assignment <- 0
    return(assignment)
}
