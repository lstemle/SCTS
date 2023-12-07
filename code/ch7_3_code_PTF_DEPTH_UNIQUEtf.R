#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kery & J. Andy Royle
#
# Chapter 7. Modeling abundance using multinomial N-mixture models
# =========================================================================

library(AHMbook)
library(unmarked)

# 7.5 Example 1: Bird point counts based on removal sampling
# ==========================================================

#make data frame

PTF_data_R$unique <- as.factor(PTF_data_R$unique)#first make factor
PTFdata.list <- list(data = PTF_data_R[5:10], covariates = PTF_data_R[1:4])

#check it came out in 2 seperate ones with 408 X6 and 408 x4
View(PTFdata.list)



# 7.5.1 Setting up the data for analysis
# ------------------------------------------------------------------------
library(unmarked)

#this removes unique and keeps only double variables
#CTSFrame <- unmarkedFrameMPois(y = CTSdata.list$data,
#                               siteCovs = as.data.frame(scale(CTSdata.list$covariates[,-1])),
#                               type = "removal")

#trying to get the factor to be kept
PTFFrame <- unmarkedFrameMPois(y = PTFdata.list$data,
                               siteCovs = as.data.frame((PTFdata.list$covariates)),
                               type = "removal")
#Error in colMeans(x, na.rm = TRUE) : 'x' must be numeric, UNLESS REMOVE SCALE


# Fit models: multinomPois order of formulas: detection, abundance 
#16 total

f0tf <- multinomPois(~ 1 ~ 1, PTFFrame)
f1tf <- multinomPois(~ 1 ~ depth, PTFFrame)
f2tf <- multinomPois(~ 1 ~ unique, PTFFrame)
f3tf <- multinomPois(~ 1 ~ depth + unique, PTFFrame)
f4tf <- multinomPois(~unique ~ depth, PTFFrame)
f5tf <- multinomPois(~ depth ~ depth, PTFFrame)
f6tf <- multinomPois(~ depth ~ depth + unique, PTFFrame)
f7tf <- multinomPois(~ depth + unique ~ 1, PTFFrame)
f8tf <- multinomPois(~ unique ~1, PTFFrame)
f9tf <- multinomPois(~ depth ~ 1, PTFFrame) 
f10tf <- multinomPois(~ unique ~ unique, PTFFrame)
f11tf <- multinomPois(~ depth + unique ~ depth + unique, PTFFrame)
f12tf <- multinomPois(~ depth + unique ~ unique, PTFFrame)#gets a warning
f13tf <- multinomPois(~ depth + unique ~depth, PTFFrame)
f14tf <- multinomPois(~ unique ~ unique + depth, PTFFrame)
f15tf <- multinomPois(~ depth ~ unique, PTFFrame)



# Rank models by AIC
msPTF <- fitList(
  "lam(.)p(.)"                                = f0tf,
  "lam(depth)p(.)"                              = f1tf,
  "lam(unique)p(.)"                             = f2tf,
  "lam(depth+unique)p(.)"                         = f3tf,
  "lam(depth)p(unique)"                = f4tf,
  "lam(depth)p(depth)"                       = f5tf,
  "lam(depth+unique)p(depth)"              = f6tf,
  "lam(.)p(depth+unique)"              = f7tf,
  "lam(.)p(unique)"              = f8tf,
  "lam(.)p(depth)"                       = f9tf,
  "lam(unique)p(unique)"              = f10tf,
  "lam(depth+unique)p(depth+unique)"                       = f11tf,
  "lam(unique)p(depth+unique)"                       = f12tf,
  "lam(depth)p(depth+unique)"              = f13tf,
  "lam(depth+unique)p(unique)"                       = f14tf,
  "lam(unique)p(depth)"                       = f15tf)

(ms1PTF <- modSel(msPTF))
#f11 is the best bc lowest aic
f11tf

# Table with everything you could possibly need
coef(msPTF) # Only first 4 columns shown




output <- as(ms1CTS4, "data.frame")


msPTF.predicted.det <- predict(f11, type = "det") 
msPTF.predicted.det
warnings()

msPTF.predicted_unique_depth <- predict(f11tf, type = "state") 
msPTF.predicted_unique_depth


write.csv(msPTF.predicted.det, "lambda_dep_site_phi_depth_site_multiPos_PTF_det.csv")



# 7.5.3 Fitting models using function gmultmix
# ------------------------------------------------------------------------
ovenFrame <- unmarkedFrameGMM(ovendata.list$data,
                              siteCovs=as.data.frame(scale(ovendata.list$covariates[,-1])),
                              numPrimary=1,type = "removal")

fm0 <- gmultmix(lambdaformula = ~1, phiformula = ~1, pformula = ~1,
                data=ovenFrame)

# Fit Poisson models
fm1 <- gmultmix(~ ufc, ~ 1, ~  1, data = ovenFrame)
fm2 <- gmultmix(~ trba, ~ 1, ~ 1, data = ovenFrame)
fm3 <- gmultmix(~ ufc + trba, ~ 1, ~ 1, data = ovenFrame)
fm4 <- gmultmix(~ ufc + trba + ufc:trba, ~ 1, ~ 1, data = ovenFrame)
# Maybe p also depends on understory foliage?
fm5 <- gmultmix(~ ufc + trba, ~ 1, ~ ufc, data = ovenFrame)
fm6 <- gmultmix(~ ufc + trba + ufc:trba, ~ 1, ~ ufc, data = ovenFrame)

# Fit analogous NegBin models
fm0nb <- gmultmix(~ 1, ~ 1, ~ 1, mixture = "NB", data = ovenFrame)
fm1nb <- gmultmix(~ ufc, ~ 1, ~ 1, mixture = "NB", data = ovenFrame)
fm2nb <- gmultmix(~ trba, ~ 1, ~ 1, mixture = "NB", data = ovenFrame)
fm3nb <- gmultmix(~ ufc + trba , ~ 1, ~ 1, mixture = "NB", data = ovenFrame)
fm4nb <- gmultmix(~ ufc + trba + ufc:trba, ~ 1, ~ 1, mixture = "NB",
                  data = ovenFrame)
# maybe p also depends on understory foliage?
fm5nb <- gmultmix(~ ufc + trba, ~ 1, ~ ufc, mixture = "NB",
                  data = ovenFrame)
fm6nb <- gmultmix(~ ufc + trba + ufc:trba, ~ 1, ~ ufc, mixture = "NB",
                  data = ovenFrame)

# Rank models by AIC
gms <- fitList(
  "lam(.)p(.)"                                = fm0,
  "lam(ufc)p(.)"                              = fm1,
  "lam(trba)p(.)"                             = fm2,
  "lam(ufc+trba)p(.)"                         = fm3,
  "lam(ufc+trba+ufc:trba)p(.)"                = fm4,
  "lam(ufc+trba)p(ufc)"                       = fm5,
  "lam(ufc+trba+ufc:trba)p(ufc)"              = fm6,
  "NB,lam(.)p(.)"                             = fm0nb,
  "NB,lam(ufc)p(.)"                           = fm1nb,
  "NB,lam(trba)p(.)"                          = fm2nb,
  "NB,lam(ufc+trba)p(.)"                      = fm3nb,
  "NB,lam(ufc+trba+ufc:trba)p(.)"             = fm4nb,
  "NB,lam(ufc+trba)p(ufc)"                    = fm5nb,
  "NB,lam(ufc+trba+ufc:trba)p(ufc)"           = fm6nb)

(gms1 <- modSel(gms))

# Table with everything you could possibly need
output <- as(gms1, "data.frame")

# Summary results
gms1

fm2nb


print(coef(gms1), digits = 2)


# 7.5.3 Fitting models using function gmultmix
# ------------------------------------------------------------------------
CTSFrameG <- unmarkedFrameGMM(CTSdata.list$data,
                              siteCovs=as.data.frame(scale(CTSdata.list$covariates[,-1])),
                              numPrimary=1,type = "removal")

fmG0 <- gmultmix(lambdaformula = ~1, phiformula = ~1, pformula = ~1,
                data=CTSFrameG)

# Fit Poisson models
fmG1 <- gmultmix(~ depth, ~ 1, ~  1, data = CTSFrameG)
fmG2 <- gmultmix( ~ 1, ~ depth, ~ 1, data = CTSFrameG)
fmG3 <- gmultmix(~ depth, ~ depth, ~ 1, data = CTSFrameG)

# Rank models by AIC
gms <- fitList(
  "lam(.)p(.)"                                = fmG0,
  "lam(depth)p(.)"                            = fmG1,
  "lam(.)p(depth)"                             = fmG2,
  "lam(depth)p(depth)"                         = fmG3
  )

(gms1 <- modSel(gms))

# Table with everything you could possibly need
output <- as(gms1, "data.frame")

# Summary results
gms1

#nPars     AIC  delta    AICwt cumltvWt
#lam(depth)p(.)         3 1523.53   0.00  5.0e-01     0.50
#lam(depth)p(depth)     3 1523.53   0.00  5.0e-01     1.00
#lam(.)p(.)             2 1992.73 469.20 6.5e-103     1.00
#lam(.)p(depth)         2 1992.73 469.20 6.5e-103     1.00


print(coef(gms1), digits = 2)

#                   lambda(depth) lambda(Int) p(Int)
#lam(depth)p(.)              0.68       -0.39  -0.45
#lam(depth)p(depth)          0.68       -0.39  -0.45
#lam(.)p(.)                    NA        0.01  -0.45
#lam(.)p(depth)                NA        0.01  -0.45


##
GCTS.predicted <- predict(fmG0, type = "lambda") 
GCTS.predicted

write.csv(GCTS.predicted, "lambda_1_phi_1.csv")

GCTS.predicted2 <- predict(fmG1, type = "lambda") 
GCTS.predicted2
write.csv(GCTS.predicted2, "lambda_depth_phi_1.csv")

msCTS.predicted3 <- predict(fmG2, type = "lambda") 
msCTS.predicted3
write.csv(msCTS.predicted3, "lambda_1_phi_depth.csv")

msCTS.predicted4<- predict(fmG3, type = "lambda") 
msCTS.predicted4
write.csv(msCTS.predicted4, "lambda_depth_phi_depth2.csv")

msCTS.predicted5<- predict(fmG3, type = "standard") 
msCTS.predicted5
write.csv(msCTS.predicted4, "lambda_depth_phi_depth2.csv")

write.csv(msCTS.predicted4, "lambda_depth_phi_depth.csv")

write.csv(CTSdata, "CTSdata_R.csv")

library(tidyverse)
library(dplyr)


abund<- PTF_data_match %>% 
  group_by(.dots = c("unique")) %>% 
  summarize(sum_predicted_nmix_PTF = sum(Predicted))
write.csv(abund, "PTF_estimates.csv")


lm1<-lm(abund$sum_predicted_n~abund$sum_nl_estimate)
lm1
summary(lm1)

#pull into summary datasheet
write.csv(Pond_level_comparision_summary, "Pond_level_comparision_summ_final_transformed.csv")

shapiro.test(Pond_level_comparision_summary$sum_predicted_nmix_mp_final)#not normal
Pond_level_comparision_summary$log_nmix <- log(Pond_level_comparision_summary$sum_predicted_nmix_mp_final+0.0383)
shapiro.test(Pond_level_comparision_summary$log_nmix) #closer to normal now
hist(Pond_level_comparision_summary$log_nmix)
shapiro.test(Pond_level_comparision_summary$dendip)#not normal


Pond_level_comparision_summary$log_dendip <- log((Pond_level_comparision_summary$dendip+0.0625))
shapiro.test(Pond_level_comparision_summary$log_dendip)#normal now

lm2<-lm(Pond_level_comparision_summary$log_nmix~Pond_level_comparision_summary$log_dendip)
summary(lm2)
lm2
#Call:lm(formula = Pond_level_comparision_summary$log_nmix ~ Pond_level_comparision_summary$log_dendip)

#Residuals:
#  Min      1Q  Median      3Q     Max 
#-2.8428 -1.4572  0.2932  1.2988  2.1539 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                                 3.7421     0.6304   5.936 9.79e-05 ***
#  Pond_level_comparision_summary$log_dendip   2.0002     0.3866   5.174 0.000306 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 1.708 on 11 degrees of freedom
#Multiple R-squared:  0.7088,	Adjusted R-squared:  0.6823 
#F-statistic: 26.77 on 1 and 11 DF,  p-value: 0.0003064

plot(Pond_level_comparision_summary$log_nmix~Pond_level_comparision_summary$log_dendip)


cor.test(Pond_level_comparision_summary$log_nmix,Pond_level_comparision_summary$log_dendip,
         method = "spearman")

#S = 21.529, p-value = 1.672e-06 alternative hypothesis: true rho is not equal to 0
#  rho 0.9408537 

#redo
CTSFrameGS <- unmarkedFrameGMM(CTSdata.list$data,
                              siteCovs=as.data.frame((CTSdata.list$covariates)),
                              numPrimary=1,type = "removal")

fmG0 <- gmultmix(lambdaformula = ~1, phiformula = ~1, pformula = ~1,
                 data=CTSFrameGS)

# Fit Poisson models
fmG1 <- gmultmix(~ depth, ~ 1, ~  1, data = CTSFrameGS)
fmG2 <- gmultmix( ~ 1, ~ 1, ~ depth, data = CTSFrameGS)
fmG3 <- gmultmix(~ depth, ~ 1, ~ depth, data = CTSFrameGS)
fmG4 <- gmultmix(~ unique, ~ 1, ~ 1, data = CTSFrameGS)
fmG5 <- gmultmix(~ 1, ~ 1, ~unique, data = CTSFrameGS)
fmG6 <- gmultmix(~ unique + depth, ~ 1, ~ 1, data = CTSFrameGS)
#fmG7 <- gmultmix(~ unique, ~ 1, ~ unique, data = CTSFrameGS)
fmG8 <- gmultmix(~ 1, ~ 1, ~unique + depth, data = CTSFrameGS)
fmG9 <- gmultmix(~ unique + depth, ~ 1, ~ depth, data = CTSFrameGS)
#fmG10 <- gmultmix(~ unique, ~ 1, ~ unique +depth, data = CTSFrameGS)
fmG11 <- gmultmix(~ depth, ~ 1, ~unique + depth, data = CTSFrameGS)
fmG12 <- gmultmix(~ unique + depth, ~ 1, ~ depth +unique, data = CTSFrameGS)
fmG13 <- gmultmix(~ depth, ~ 1, ~ unique, data = CTSFrameGS)
fmG14 <- gmultmix(~ unique, ~ 1, ~ depth, data = CTSFrameGS)
fmG15 <- gmultmix(~ unique + depth, ~ 1, ~ unique, data = CTSFrameGS)


# Rank models by AIC
gmsG <- fitList(
  "lam(.)p(.)"                                = fmG0,
  "lam(depth)p(.)"                              = fmG1,
  "lam(.)p(depth)"                             = fmG2,
  "lam(depth)p(depth)"                         = fmG3,
  "lam(unique)p(.)"                = fmG4,
  "lam(.)p(unique)"                       = fmG5,
  "lam(unique+depth)p(.)"              = fmG6,
  "lam(unique)p(depth)"                             = fmG14,
  "lam(.)p(unique+depth)"                           = fmG8,
  "lam(unique+depth)p(depth)"                          = fmG9,
  "lam(unique+depth)p(unique)"                      = fmG15,
  "lam(depth)p(unique+depth)"             = fmG11,
  "lam(unique+depth)p(unique+depth)"                    = fmG12,
  "lam(depth)p(unique)"           = fmG13)

(gms1G <- modSel(gmsG))

CTSFrameG <- unmarkedFrameGMM(CTSdata.list$data,
                              siteCovs=as.data.frame(scale(CTSdata.list$covariates[,-1])),
                              numPrimary=1,type = "removal")

fmG0 <- gmultmix(lambdaformula = ~1, phiformula = ~1, pformula = ~1,
                 data=CTSFrameG)

# Fit Poisson models
fmG1 <- gmultmix(~ depth, ~ 1, ~  1, data = CTSFrameG)
fmG2 <- gmultmix( ~ 1, ~ depth, ~ 1, data = CTSFrameG)
fmG3 <- gmultmix(~ depth, ~ depth, ~ 1, data = CTSFrameG)

# Rank models by AIC
gms <- fitList(
  "lam(.)p(.)"                                = fmG0,
  "lam(depth)p(.)"                            = fmG1,
  "lam(.)p(depth)"                             = fmG2,
  "lam(depth)p(depth)"                         = fmG3
)
# Table with everything you could possibly need
output <- as(gms1, "data.frame")

# Summary results
gms1

fm2nb


print(coef(gms1), digits = 2)

# 7.5.4 Assessing model fit in unmarked
# ------------------------------------------------------------------------
set.seed(2015)
(gof <- parboot(f11, fitstats, nsim = 1000, report = 1)) #hmm doesn't look too good here

#t0 = 231312.7 17353.72 4672.322 
#Running in parallel on 3 cores. Bootstrapped statistics not reported.

#Call: parboot(object = f11, statistic = fitstats, nsim = 1000, report = 1)

#Parametric Bootstrap Statistics:
#             t0 mean(t0 - t_B) StdDev(t0 - t_B) Pr(t_B > t0)
#SSE          231313         214386            757.0            0
#Chisq         17354          14935             69.7            0
#freemanTukey   4672           3864             24.4            0

#t_B quantiles:
#  0%  2.5%   25%   50%   75% 97.5%  100%
#SSE          14720 15420 16385 16911 17487 18440 19649
#Chisq         2202  2290  2370  2417  2469  2555  2644
#reemanTukey   739   761   792   807   824   855   887

#t0 = Original statistic computed from data
#t_B = Vector of bootstrap samples

#the book says for a pvalues of .5, .1 and .5
#These results indicate that the best model appears to fit the data reasonably well with the bootstrap
#p-value not being extreme (not so close to 0 or 1) for any of the three fit statistics. 


set.seed(2015)#now9 is the best (same parameters) #doesn't look good either 
(gof <- parboot(fmG9, fitstats, nsim = 1000, report = 1))

#Call: parboot(object = fmG9, statistic = fitstats, nsim = 1000, report = 1)

#Parametric Bootstrap Statistics:
#           t0 mean(t0 - t_B) StdDev(t0 - t_B) Pr(t_B > t0)
#SSE          1167          786.4             39.4        0.000
#Chisq        2243          212.2            263.3        0.181
#freemanTukey  321           70.7             12.9        0.000

#t_B quantiles:
#  0% 2.5%  25%  50%  75% 97.5% 100%
#SSE           274  314  351  377  407   469  504
#Chisq        1419 1620 1854 1995 2178  2669 3256
#freemanTukey  209  226  241  250  259   275  292
#t0 = Original statistic computed from data
#t_B = Vector of bootstrap samples
