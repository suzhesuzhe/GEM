rm(list=ls())
require(devtools)

devtools::install_github(repo = "suzhesuzhe/GEM",force =T)

require(GEM)

packageDescription("GEM")

#covariance matrix, 10 predictors
co <- matrix(0.2, 10, 10)
diag(co) <- 1

#simulate gem type data
dataEx <- data_generator1(d = 0.3, R2 = 0.5, v2 = 1, n = 300, co = co, beta1 = rep(1,10),inter = c(0,0))

#extract the dataframe
dat <- dataEx[[1]]


#simluate another data of gigantic type
bigData <- data_generator3(n = 10000,co = co, bet =dataEx[[2]], inter = c(0,0))


#fit gem model with numerator method
model <- gem_fit(dat = dat, method = "nu")

#extract the augmented data with combined Z
aug <- model[[4]]

#calculated effect size, this is the same with model[[5]]
effectSize(response = aug[,2],treatment = aug[,1],moderator = aug$Z)                                                                                                                                    


#use the gem fit object to test itself
gem_test_sample(dat,model[[2]])

# use the gem fit object to test a dataset with observations under both treatment group ( gigantic data)
gem_test_simsample(bigData[[1]],bigData[[2]],bigData[[3]],model[[2]])


#calculate the permuted p value
permute_pvalue(dat = dat,method = "nu")


#to see all the functions in GEM package
ls(envir=as.environment("package:GEM"))


