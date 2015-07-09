
effectSize <- function(response, treatment, moderator)
{
    mod_scale <- scale(moderator, center=TRUE, scale=TRUE)
    treat_recode <- as.numeric(as.character(factor(treatment, labels=c(-0.5, 0.5))))
    mm <- lm(response~treat_recode + mod_scale + treat_recode * mod_scale)
    
    eff_size <- as.numeric(mm[[1]][4]/2/sqrt(summary(mm)[[6]]^2+mm[[1]][3]^2+mm[[1]][4]^2/4))
    
    return(eff_size)
    
}