#' Simulating data set under Gem scenario
#' @name gemData
#' @import MASS
#' 
#' @description These are two functions for simulating data set under Gem scenario.
#'
#' @param es A scalar indicating the effect size for each treatment group
#' @param r A scalar indicating the proportion of explained variance \eqn{R^2} for each treatment group
#' @param p A scalar indicating the numbder of predictors 
#' @param n A scalar indicating the number of observation for each treatment group
#' @param co A matrix of covariance structure for the predictors, if there is p predictors, 
#' this should be of dimension p by p, symmetric and positive definite.
#' @param beta1 A vector of length p recording coefficients of the predictors for the first treatment groups respectively
#' @param inter A vector of length 2 recording the intercept term for those two treatment groups respectively
#' 
#' @details These are two functions for simulating data set under Gem scenario, for both 
#' in-sample and out-sample types respectively. The data is simulated as a linear combination of 
#' certain predictors, which are multivariate Guassian distributed.  #The coefficients for the first 
#' treatment group is taken as an argument for the function while the coefficient of the second 
#' treamtment group is derived by the first one times a multiplier, which is calculated by \eqn{R^2} and effect size in the function. 
#' Since they are designed for simulating data to test with calculting gem, 
#' those data generators take the simplest setting for two treatment groups with same number 
#' of observations and same proportion of explained variance \eqn{R^2}. 
#' 
#' The \code{gemData_inSample} will generate only one reponse for each observation under certain treatment 
#' group, while \code{gemData_outSample} will generate two responses under both of the treatments. In general 
#' \code{gemData_inSample} is used to simulate a sample data to train the gem algorithm and \code{gemData_outSample}
#' is used for simulating the population.
#' 
#' As the gem algorithm is based on the linear model coefficients and the number of subjects
#' in each treatment group, see\url{http://en.wikipedia.org/wiki/Integer_overflow}, 
#' the function \code{gemData_inSample} calculates the theoretical value of the alpha under three gem methods. 
#
#' @return Output of these two functions for the data generator under Gem case is not the same,
#' and you can find explanation for each of them in the following:
#' 
#' For the function \code{nongemData_inSample}
#' \enumerate{
#' 		\item \code{dat} A dat set with the first column as treatment group index, the second 
#' 		column as response and the remaining columns as the predictors design matrix
#' 		\item \code{nu_true_alpha} Theoretical alpha value under numerator method
#' 		\item \code{de_true_alpha} Theoretical alpha value under denominator method
#' 		\item \code{F_true_alpha} Theoretical alpha value under F statistics method
#' 		\item \code{preset_co} The covariance matrix that used for simulating the data, mainly for reference
#' 		\item \code{beta} The beta coeffcients that used for simulating the data, mainly for reference
#' }
#' 
#' For the function \code{nongemData_outSample}
#' \enumerate{
#' 		\item \code{y0} Response for each observation under the first treatment assignment
#' 		\item \code{y1} Response for each observation under the second treatment assignment
#' 		\item \code{X} Design matrix for the predictors 
#' 		\item \code{optTrt} Optimal treatment assignment, for example, if for certain observation \code{y0>y1},
#' 		then the first treatment assignment is optimal
#' 		\item \code{oracle} Average of the response if every observation take the optimal treatment assignment
#' 		\item \code{oracle2} Average of the response if every observation not to take the optimal treatment assignment
#' }
#' 
#' @examples
#' bet <- c(1,1,1,1,1)
#' co <- matrix(0.2,5,5)
#' diag(co) <- 1
#' popu3 <- gemData_inSample(0.3,0.1,5,100,co,bet,c(0,1))
#' popu4 <- gemData_outSample(0.3,0.5,5,50000,co,bet,c(0,1))
#' @export
gemData_inSample <- function(es,  
                     r,
                     p,
                     n,
                     co, # here co is a matrix rather than a list, which assume that all treatment group's design matrix have the same covariance
                     beta1, #the coefficient of beta1
                     inter) #the intercept for each treatment group

    
{
    K <- ncol(co)
    pro <- rep(1/K,K)
    pro1 <- rep((n-1)/(n*K),K)
    bet <- vector("list",2)
    bet[[1]] <- bet1 <- beta1
	
    temp <- 2 * es^2 -r
    a <- (-r + sqrt(r^2-temp^2))/temp
    bet[[2]] <- bet2 <- bet1 * a
    
    X <- vector("list",2)#Here we create this list for the mod3 below
    
    X[[1]] <- x1 <- mvrnorm(n,rep(0,ncol(co)),co)
	X[[2]] <- x2 <- mvrnorm(n,rep(0,ncol(co)),co)
    
    SSR1 <- (t(bet1) %*% co %*% bet1) * (n-1) 
    SSE1 <- SSR1/r-SSR1
    s1 <- sqrt(SSE1/(n-2))
    
    SSR2 <- (t(bet2) %*% co %*% bet2) * (n-1)
    SSE2 <- SSR2/r-SSR2
    s2 <- sqrt(SSE2/(n-2))
    
    Y <- vector("list",2)#Here we create this list for the mod3 below
    Y[[1]] <- y1 <- x1 %*% bet1 + matrix(rnorm(n)*s1,n,1) + inter[1]
    Y[[2]] <- y2 <- x2 %*% bet2 + matrix(rnorm(n)*s2,n,1) + inter[2]
    
    
    trt1<- matrix(0,n,1)
    trt2<- matrix(1,n,1)
    
    dat1 <- as.data.frame(cbind(trt1,y1,x1))
    dat2 <- as.data.frame(cbind(trt2,y2,x2))
    dat <- rbind(dat1, dat2)
    colnames(dat)[1:2] <- c("trt","y")
    colnames(dat1)[1:2] <- c("trt","y")
    colnames(dat2)[1:2] <- c("trt","y")

    
    MM1 <- true_alpha(co, bet, pro1=NULL, pro, r=NULL, type="nu")
    MM2 <- true_alpha(co, bet, pro1, pro, r=NULL, type="de")
    MM3 <- true_alpha(co, bet, pro1, pro, r, type="F")
    
    results <- list("dat" = dat,                    #3
                    "nu_true_alpha" = MM1[[1]],     #4
                    "de_true_alpha" = MM2[[1]],     #6
                    "F_true_alpha" = MM3[[1]],      #8
                    "preset_co" = co,               #10
                    "bet" = bet)                    #11
    return(results)
}

#' @name gemData
#' @export
gemData_outSample <- function(es,  
                       		  r,
                              p,
                        	  n,
                              co, 
                              beta1, 
                              inter) 
    
{
    K <- ncol(co)
    pro <- rep(1/K,K)
    pro1 <- rep((n-1)/(n*K),K)
    bet <- vector("list",2)
    bet[[1]] <- bet1 <- beta1

    
    temp <- 2 * es^2 -r
    a <- (-r + sqrt(r^2-temp^2))/temp
    bet[[2]] <- bet2 <- bet1 * a
    
    X <- mvrnorm(n,rep(0,ncol(co)),co)
    
    SSR1 <- (t(bet1) %*% co %*% bet1) * (n-1) 
    SSE1 <- SSR1/r-SSR1
    s1 <- sqrt(SSE1/(n-2))
    
    SSR2 <- (t(bet2) %*% co %*% bet2) * (n-1)
    SSE2 <- SSR2/r-SSR2
    s2 <- sqrt(SSE2/(n-2))
    
    Y <- vector("list",2)#Here we create this list for the mod3 below
    Y[[1]] <- y1 <- X %*% bet1 + matrix(rnorm(n)*s1,n,1) + inter[1]
    Y[[2]] <- y2 <- X %*% bet2 + matrix(rnorm(n)*s2,n,1) + inter[2]
    optTrt <- as.numeric(Y[[1]]<=Y[[2]])
	oracle <- sum(y1*(1-optTrt)+y2*optTrt)/length(optTrt)
	oracle2 <- sum(y1*optTrt+y2*(1-optTrt))/length(optTrt)

	optTrt <- y1<y2
	pro0less1 <- sum(optTrt)/n


     results <- list("y0" = Y[[1]],                  #1
                     "y1" = Y[[2]],                  #2
                     "X" = X,
     				 "optTrt" = optTrt,
   					 "oracle"=oracle,
     				 "oracle2"=oracle2,
   				     "optTrt"=optTrt,
   					 "pro0less1"=pro0less1,
     				 "bet" = bet)                    #3
    return(results)
}


