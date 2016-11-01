#' Functions for Simulating Data 
#' @import MASS
#' 
#' @description When investigating the properties of GEM, the following three data generators are used in various simulations. 
#' They are designed to construct three specific types of data sets in the case of two treatment groups. See more detail in 
#' \cite{E Petkova, T Tarpey, Z Su, and RT Ogden. Generated effect modifiers (GEMs) in randomized clinical trials. Biostatistics, (First published online: July 27, 2016). doi: 10.1093/biostatistics/kxw035.}
#'
#' 
#' @param d A scalar indicating the effect size of the GEM when the data is generated under a GEM model 
#' @param R2 A scalar indicating the proportion of explained variance \eqn{R^2} for the entire data set
#' @param v2 A scalar indicating the proportion of explained variance \eqn{R^2} for the first treatment group
#' @param n A scalar indicating the number of observation in each treatment group, assumed to be the same.
#' @param co A \emph{p} by \emph{p} positive semidefinite matrix indicating the covariance matrix of the covariates
#' @param beta1 A vector of length \emph{p} giving the regression coefficients for the first treatment group
#' @param bet A list with two elements, each a vector of length \emph{p}, giving the regression coefficients for the two treatment groups respectively
#' @param inter A vector of length 2 recording the intercepts \eqn{\beta_{10},\beta_{20}} for the two treatment groups respectively
#' 
#' @details \code{data_generator1} is used to create data where the outcome is a linear function of the covariates \deqn{y_j = \beta_{j0} + X\beta_j + \epsilon, j = 1, 2, }  
#' and the coffcicients of covariates \eqn{\beta} are proportional between two treatment groups: \eqn{\beta_2 = b * \beta_1}. 
#' This type of data set matches perfectly with the motivation of GEM algorithm. \eqn{\beta_1} is set as an argument of the function while \eqn{\beta_2 = b * \beta_1}
#' is derived by controling \eqn{R^2} of the whole data and the effect size. See more detail in \cite{Kraemer, H. C. (2013). Discovering, comparing, and combining moderators of treatment 
#' on outcome after randomized clinical trials: a parametric approach. Statistics in medicine, 32(11), 1964-1973.}
#'
#' \code{data_generator2} is similar to the first one except that the coefficients of the covariates are not necessarily proportional. Hence two \eqn{\bold{\beta}}'s 
#' should be specified as arguments of the function.
#' 
#' \code{data_generator3} constructs a data set where the outcome under each treatment condition is given for all subjects. In addition, no error is added to the mean outcome.
#' This generator is useful for obtaining the "true" value of a treatment decision. This data generator is similar to data generator2 \deqn{y_j = \beta_{j0} + X\beta_j, j = 1,2.}
#'
#' @return The output from these functions are different:
#' 
#' For the function \code{data_generator1}
#' \enumerate{
#' 		\item \code{dat} A data frame with first and second column as treatment group index and outcome respectively,
#'   and each of the remaining columns as a covariate.
#' 		\item \code{bet} A list with two elements, each a vector of length \eqn{p}, giving the regression coefficients for the two treatment groups respectively
#' 		\item \code{error_12} A vector of length three represeting the standard deviation of \eqn{\epsilon}, the explained variance by the linear part for the first 
#'   and second treatment group respectively. 
#' }
#' 
#' For the function \code{data_generator2}
#' \enumerate{
#'   \item \code{dat} A data frame with first and second column as treatment group index and outcome respectively,
#'   and each of the remaining columns as a covariate.
#' 		\item \code{bet}  list with two elements, each a vector of length \eqn{p}, giving the regression coefficients for the two treatment groups respectively
#' 		\item \code{error} A scalar represeting the standard deviation of \eqn{\epsilon} 
#' }
#'
#' For the function \code{data_generator3}
#' \enumerate{
#' 		\item \code{y0} Outcome vector under the first treatment assignment
#' 		\item \code{y1} Outcome vector under the second treatment assignment
#' 		\item \code{X} Design matrix for the covariates 
#' 		\item \code{oracle} Average of the outcome if each subject takes the optimal treatment assignment
#' 		\item \code{invOracle} Average of the outcome if each subject does not take the optimal treatment assignment
#' } 
#' @examples
#' #constructing the covariance matrix
#' co <- matrix(0.2, 30, 30)
#' diag(co) <- 1
#' dataEx <- data_generator1(d = 0.3, R2 = 0.5, v2 = 1, n = 3000, 
#'                            co = co, beta1 = rep(1,30),inter = c(0,0))
#' #check the R squared of the simluated data set
#' dat <- dataEx[[1]]
#' summary(lm(V2~factor(trt)*(V3+V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14+V15+V16+
#' V17+V18+V19+V20+V21+V22+V23+V24+V25+V26+V27+V28+V29+V30+V31+V32),data=dat))
#' 
#' bigData <- data_generator3(n = 10000,co = co,bet =dataEx[[2]], inter = c(0,0))
#' @name data_generator
#' @export

data_generator1 <- function(d,    #effect size for the data set
                            R2,          # r square for the data set
					                       v2,          #the SSR for treatment group 1
                            n,           # number of observation for each treatment group
                            co,          # here co is a matrix rather than a list, which assume that all treatment group's design matrix have the same covariance
                            beta1,       #the coefficient of beta1
                            inter)       #the intercept for each treatment group

    
{
	    a2 = a0 = 2*d^2+(2-2*R2)*v2*d^2/R2 -1
	    a1 = 2

    	b  = (-a1 + sqrt(a1^2 - 4*a0*a2))/2/a2
    	s = 0.5*(1+b^2)*v2*(1-R2)/R2
    
    	temp <- v2/(t(beta1) %*% co %*% beta1)
	
    	bet <- vector("list",2)#Here we create this list for the mod3 below

    	beta1 <- beta1 * sqrt(temp)
	    beta2 <- beta1 * b
	    bet[[1]] <- beta1
    	bet[[2]] <- beta2
	
     X <- vector("list",2)#Here we create this list for the mod3 below
    
     X[[1]] <- x1 <- mvrnorm(n,rep(0,ncol(co)),co)
	    X[[2]] <- x2 <- mvrnorm(n,rep(0,ncol(co)),co)

     Y <- vector("list",2)#Here we create this list for the mod3 below
     Y[[1]] <- y1 <- x1 %*% beta1 + matrix(rnorm(n)*sqrt(s),n,1) + inter[1]
     Y[[2]] <- y2 <- x2 %*% beta2 + matrix(rnorm(n)*sqrt(s),n,1) + inter[2]
    
    
     trt1<- matrix(0,n,1)
     trt2<- matrix(1,n,1)
    
     dat1 <- as.matrix(cbind(y1,x1))
     dat2 <- as.matrix(cbind(y2,x2))
     dat <- as.data.frame(rbind(cbind(trt1,dat1), cbind(trt2,dat2)))
	    colnames(dat)[1] <- "trt"
	
     results <- list("dat" = dat,                    
                     "bet" = bet,
                     "error_12" = c(sqrt(s),t(beta1) %*% co %*% beta1,t(beta2) %*% co %*% beta2))   
     return(results)
}

#' @name data_generator
#' @export

data_generator2 <- function(n,    #A number specifying the number of observation for each group
                            co,   #co is the covariate covariance that is identical between treatment groups
                            R2,    #A scalar specifying the proportion of variance explained (R^2) 
                                  #by the predictors, same in all groups
                            bet,  #A two element list of legnth p vectors recording covariate coefficients for two treatment groups respectively
                            inter) #Vector of treatment groups' intercepts
{    
	bet1 <- as.matrix(bet[[1]])
	bet2 <- as.matrix(bet[[2]])
	
    X <- vector("list",2)#Here we create this list for the mod3 below
    
    X[[1]] <- x1 <- mvrnorm(n,rep(0,ncol(co)),co)
   	X[[2]] <- x2 <- mvrnorm(n,rep(0,ncol(co)),co)
    
   	s2 <- 0.5 * (1/R2-1) * (t(bet1) %*% co %*% bet1 + t(bet2) %*% co %*% bet2)
	
    y1 <- x1 %*% bet1 + matrix(rnorm(n)*sqrt(s2),n,1) + inter[1]
    y2 <- x2 %*% bet2 + matrix(rnorm(n)*sqrt(s2),n,1) + inter[2]
    
    
    trt1<- matrix(0,n,1)
    trt2<- matrix(1,n,1)
    
    dat1 <- as.matrix(cbind(y1,x1))
    dat2 <- as.matrix(cbind(y2,x2))
    dat <- as.data.frame(rbind(cbind(trt1,dat1), cbind(trt2,dat2)))
    colnames(dat)[1] <- "trt"
    results <- list("dat" = dat,                  
                    "bet" = bet,
             				   "error" = c(sqrt(s2))) 
   return(results)
}




#' @name data_generator
#' @export
data_generator3 <- function(n,     #the number of oberservations
                            co,    #co is the covariate covariance that is identical between treatment groups
                            bet,   #A two element list of legnth p vectors recording covariate coefficients for two treatment groups respectively
                            inter) #Vector of treatment groups' intercepts
{
    p <- ncol(co) #The number of predictors
    
	   bet0 <- bet[[1]]
   	bet1 <- bet[[2]]
	
   	X <- as.matrix(mvrnorm(n,rep(0,ncol(co)),co))

    y0 <- as.matrix(X %*% bet0 + inter[1])
    y1 <- as.matrix(X %*% bet1 + inter[2])
   	optTrt <- y0<=y1
   	oracle <- sum(y0*(1-optTrt)+y1*optTrt)/length(optTrt)
   	inv_oracle <- sum(y0*optTrt+y1*(1-optTrt))/length(optTrt)
    
    results <- list("y0" = y0,                  #1
                    "y1" = y1,                  #2
                    "X" = X,                    #3
               					"oracle"=oracle,            #4
   				            	"invOracle" = inv_oracle)   #5
    return(results)
}

 
