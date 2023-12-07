## Code for Cook et al.  by Arianne Messerman and Leyna Stemle
## Sonoma CTS larval density as predicted by pond attributes

setwd("C:/Users/leyna/Downloads")
#load("C:/Users/messe/Downloads/Sonoma_CTS_Stemle-Oct2020-Enviro-v4.RData")
dips<-read.csv("~/Sonoma_CTS_Larval_Dipnet_Data-v7-combined.csv")

#set up variables as correct type
str(dips)
dips$unique.id<-as.factor(dips$unique.id)
dips$surv.type<-as.factor(dips$surv.type)
dips$site<-as.factor(dips$site)
dips$date<-as.Date(dips$date, "%m/%d/%Y")
dips$pool.num<-as.factor(dips$pool.num)
dips$pool.type<-as.factor(dips$pool.type)
dips$depth.cm<-as.numeric(dips$depth.cm)
dips$surv.min<-as.integer(dips$surv.min)
dips$larv<-as.integer(dips$larv)
dips$pred<-as.factor(dips$pred)
dips$dry.date<-as.Date(dips$dry.date, "%m/%d/%Y")

#Examine response variable
hist(dips$larv)
hist(dips$dens)
hist(log10(dips$dens+1))
plot(density(na.omit(dips$dens)))
hist(sqrt(dips$dens+1))
hist(asin(dips$dens))
hist(asin(sqrt(dips$dens/100)))
#No typical transformations effective for normalizing larval density data

########################################################
#Create a randomization test using lmer() and a type III Anova with the car package
#Randomize density data 1000 times, run 1001 lmer with randomized data and true data
#Extract coefficients to compare groups, and X^2 test statistics to compare fixed effects
#Compare proportion randomized trials with more extreme coefficients or test statistics than true data as p-value
########################################################

library(lme4)
library(car)
library(lmerTest)

#Remove rows with NA for larval density
dips1<-dips[!is.na(dips$dens),]
sum(is.na(dips$dens))
str(dips)
str(dips1)
dips1$log.depth<-log(dips1$depth.cm+1)
dips1$log.max.depth<-log(dips1$max.depth)
dips1$log.area<-log(dips1$max.area+1)

#Scale and center covariates
for (i in 1:length(dips1$dens)) {
  dips1$std.yr[i] <- (dips1$yr[i]-mean(dips1$yr[], na.rm=TRUE))/sd(dips1$yr[], na.rm=TRUE)
  dips1$std.area[i] <- (dips1$log.area[i]-mean(dips1$log.area[], na.rm=TRUE))/sd(dips1$log.area[], na.rm=TRUE)
  dips1$std.depth[i] <- (dips1$log.depth[i]-mean(dips1$log.depth[], na.rm=TRUE))/sd(dips1$log.depth[], na.rm=TRUE)
  dips1$std.dry[i] <- (dips1$dry.date[i]-mean(dips1$dry.date[], na.rm=TRUE))/sd(dips1$dry.date[], na.rm=TRUE)
  dips1$std.max.depth[i] <- (dips1$log.max.depth[i]-mean(dips1$log.max.depth[], na.rm=TRUE))/sd(dips1$log.max.depth[], na.rm=TRUE)
}

#Examine each continuous predictor
hist(dips1$std.yr)#sampling year
hist(dips1$std.area)#maximum pond area
hist(dips1$std.depth)#pond depth at sampling period
hist(dips1$std.max.depth)#maximum pond depth overall
hist(dips1$std.dry)#date pond dried in 2007

library(corrplot)
M<-cor(dips1[,29:33], use="complete.obs")
corrplot(M, method="number")


modnatural <- lmer(dens ~ std.yr + dist_c_cbp + std.yr*dist_c_cbp +  (1|site) + (1|site:pool.num), data=dips1_naturalonly, na.action = na.exclude)
summary(modnatural)


#make main model looking at fixed and random effects on larval density
mod1<-lmer(dens ~ pool.type*std.yr + std.area + std.max.depth + 
        (1|site) + (1|site:pool.num), data=dips1, na.action = na.exclude) #std.dry removed due to correlation with std.depth
summary(mod1)
mod.aov<-Anova(mod1, type = "III")
summary(mod.aov)
mod.aov#Pre-randomization: pool.type, yr, depth, interaction effect significant (max area and dry date are not)
mod1

#Check pre-randomization pairwise comparisons to establish your expectations
library(emmeans)
emm<-emmeans(mod1, ~pool.type*std.yr + std.area + std.max.depth)
pairs(emm, simple="pool.type") #Compare pool types: all significant but constructed v. natural
emm1<-as.data.frame(emmeans(mod1, ~pool.type*std.yr + std.area + std.max.depth))
emm.true<-emm1$emmean

emm2<-emtrends(mod1, pairwise ~ pool.type, var = "std.yr")#Comparison of slopes of larval density by year given pool type
emm.contrast<-as.data.frame(emm2$contrasts)
emm.trends<-as.data.frame(emm2$emtrends)
emm.trends.true<-emm.trends$std.yr.trend

##Start randomization test
coeff.true<-fixef(mod1)#extract coefficients to compare groups
chisqr.true<-mod.aov$Chisq #extract X^2 values for fixed effects

#Build empty matrices for storage
rand.dens<-matrix(,nrow=length(dips1$dens), ncol=1000)
coeffs<-matrix(,nrow=1000, ncol=8)
chisqr<-matrix(,nrow=1000, ncol=6)
lsm<-matrix(,nrow=1000, ncol=3)
emm.tr<-matrix(,nrow=1000, ncol=3)

#Extract all data into vectors for ease of identifying model in for loop below
site1<-as.factor(c(dips1$site))
yr1<-c(dips1$std.yr)
pool.type1<-as.factor(c(dips1$pool.type))
pool.num1<-as.factor(c(dips1$pool.num))
area1<-c(dips1$std.area)
depth1<-c(dips1$std.max.depth)
dry1<-c(dips1$std.dry)

set.seed(8375)#So you can re-run the randomization later to get the same results

for (i in 1:1000){
  rand.dens[,i]<-sample(dips1$dens, length(dips1$dens), replace = FALSE)
  mod2<-lmer(rand.dens[,i] ~ pool.type1*yr1 + area1 + depth1 + (1|site1) + (1|site1:pool.num1))
  coeffs[i,]<-fixef(mod2)
  mod.aov2<-Anova(mod2, type = "III")
  chisqr[i,]<-mod.aov2$Chisq
  emm3<-as.data.frame(emmeans(mod2, ~pool.type1*yr1 + area1 + depth1))
  lsm[i,]<-emm3$emmean #Site EM means
  emm4<-emtrends(mod2, pairwise ~ pool.type1, var = "yr1")
  emm.trends2<-as.data.frame(emm4$emtrends)
  emm.tr[i,]<-emm.trends2$yr1.trend #Pool*yr slope EM means
}

View(coeffs)
View(chisqr)
View(lsm)
View(emm.tr)

coeff.prop<-c(1:length(coeff.true))
coeff.pval<-c(1:length(coeff.true))

for (i in 1:length(coeff.true)){
  coeff.prop[i]<-sum(ifelse(abs(coeffs[,i])>abs(coeff.true[i]), 1, 0)) #Number of 1000 coefficients greater than value of observed
  coeff.pval[i]<-2*(coeff.prop[i]/1000) #Calculate 2-sided p-value for each coefficient
}

chi.prop<-c(1:length(chisqr.true))
chi.pval<-c(1:length(chisqr.true))

for (i in 1:length(chisqr.true)){
  chi.prop[i]<-sum(ifelse(chisqr[,i]>chisqr.true[i], 1, 0)) #Number of 1000 Chi square values greater than observed
  chi.pval[i]<-2*(chi.prop[i]/1000) #Calculate 2-sided p-value for each Chi square value
}

#Investigate p-values
chi.pval 
coeff.pval 

#Compare differences of random group slopes for pool type fixed effect
alt.pairs<-matrix(,nrow=1000, ncol=3)
con.pairs<-matrix(,nrow=1000, ncol=3)#Matrix structure for 2 pairs with each site
nat.pairs<-matrix(,nrow=1000, ncol=3)

for (i in 1:3){
  alt.pairs[,i]<-lsm[,1]-lsm[,i]
  con.pairs[,i]<-lsm[,2]-lsm[,i]
  nat.pairs[,i]<-lsm[,3]-lsm[,i]
}

alt.true<- c(1:3)
con.true<- c(1:3)#Calculate differences in all slopes from observed data
nat.true<- c(1:3)

for (i in 1:3){
  alt.true[i]<-emm.true[1]-emm.true[i]
  con.true[i]<-emm.true[2]-emm.true[i]
  nat.true[i]<-emm.true[3]-emm.true[i]
}

alt.prop<-c(1:3)
alt.pval<-c(1:3)
con.prop<-c(1:3)
con.pval<-c(1:3) 
nat.prop<-c(1:3)
nat.pval<-c(1:3)

for (i in 1:3){
  alt.prop[i]<-sum(ifelse(abs(alt.pairs[,i])>abs(alt.true[i]), 1, 0))
  alt.pval[i]<-2*(alt.prop[i]/1000)
  con.prop[i]<-sum(ifelse(abs(con.pairs[,i])>abs(con.true[i]), 1, 0)) #Number of 1000 absolute slope value comparisons greater than absolute value of observed
  con.pval[i]<-2*(con.prop[i]/1000) #Calculate 2-sided p-value for each slope difference
  nat.prop[i]<-sum(ifelse(abs(nat.pairs[,i])>abs(nat.true[i]), 1, 0))
  nat.pval[i]<-2*(nat.prop[i]/1000) 
}

#mean values: natural>altered>constructed
#1=altered, 2=constructed, 3=natural
alt.pval 
con.pval 
nat.pval 
#Altered pools different from constructed and natural
#Constructed and natural did not differ from one another 
#altered=A, constructed=B, natural=B


#make plots to examine data
plot(as.factor(dips1$pool.type),dips1$dens, ylab="Larval density", xlab="Pool type")

library(dplyr)
str(dips1)
x.dens<-dips1 %>% 
  group_by(pool.type) %>% 
  summarize(means = mean(dens),
            sdd=sd(dens), 
            lower.ci = mean(dens) - qt(1 - (0.05 / 2), length(dens) - 1) * (sd(dens)/sqrt(length(dens))),
            upper.ci = mean(dens) + qt(1 - (0.05 / 2), length(dens) - 1) * (sd(dens)/sqrt(length(dens))))

x.dens<-x.dens[order(-x.dens$means),]#mean values: natural>altered>constructed

#Plot of observed mean with SD
#plot(x.dens$means, cex=3, pch=20, bg=1, ylab="Larval density", xlab="Poll type", ylim=c(-.5,1), xaxt="n", cex.axis=1.2, cex.lab=1.4)
#axis(1, 1:3, labels=x.dens$pool.type, cex.axis=1.2)
#mtext(side=3,line=-1.2, at = c(1:8), c("A", "AB", "B"), cex=1.1)
#segments(1:3, x.dens$means-x.dens$sdd, 1:3, x.dens$means+x.dens$sdd)

#Plot of observed mean with 95% CI
plot(x.dens$means, cex=3, pch=20, bg=1, ylab="Larval density", xlab="Pool type", xaxt="n", 
     ylim=c(0,0.3), cex.axis=1.2, cex.lab=1.4)
axis(1, 1:3, labels=x.dens$pool.type, cex.axis=1.2)
mtext(side=3,line=-1.2, at = c(1:3), c("A", "AB", "B"), cex=1.1)
segments(1:3, x.dens$lower.ci, 1:3, x.dens$upper.ci)

emm1.order<-emm1[order(-emm1$emmean),]
#Plot estimated marginal means with 95% CI
plot(emm1.order$emmean, cex=3, pch=20, bg=1, ylab="Larval density (larvae/min)", xlab="Pool type", xaxt="n", ylim=c(-0.5,.7), cex.axis=1.2, cex.lab=1.4)
axis(1, 1:3, labels=emm1.order$pool.type, cex.axis=1.2)
mtext(side=3,line=-1.2, at = c(1:8), c("  A", "A", "B"), cex=1.1)
segments(1:3, emm1.order$lower.CL, 1:3, emm1.order$upper.CL)


#Compare differences of randomized group slopes for pool.type*yr fixed interaction effects
#All pool.type*yr comparisons with altered*yr evident in coeff.pval[8:9]
#alt1.pairs<-coeffs[,8:9]#Randomized differences between each site and altered
alt1.pairs<-matrix(,nrow=1000, ncol=3)
con1.pairs<-matrix(,nrow=1000, ncol=3)#Matrix structure for 2 pairs with each pool.type*yr
nat1.pairs<-matrix(,nrow=1000, ncol=3)

for (i in 1:3){
  alt1.pairs[,i]<-emm.tr[,1]-emm.tr[,i]
  con1.pairs[,i]<-emm.tr[,2]-emm.tr[,i]
  nat1.pairs[,i]<-emm.tr[,3]-emm.tr[,i]
}

#alt1.true<- c(coeff.true[8:9])##All type comparisons with altered
alt1.true<- c(1:3)
con1.true<- c(1:3)#calculate differences in all slopes from observed data
nat1.true<- c(1:3)

for (i in 1:3){
  alt1.true[i]<-emm.trends.true[1]-emm.trends.true[i]
  con1.true[i]<-emm.trends.true[2]-emm.trends.true[i]
  nat1.true[i]<-emm.trends.true[3]-emm.trends.true[i]
}

alt1.prop<-c(1:3)
alt1.pval<-c(1:3) 
con1.prop<-c(1:3)
con1.pval<-c(1:3) 
nat1.prop<-c(1:3)
nat1.pval<-c(1:3)

for (i in 1:3){
  alt1.prop[i]<-sum(ifelse(abs(alt1.pairs[,i])>abs(alt1.true[i]), 1, 0))
  alt1.pval[i]<-2*(alt1.prop[i]/1000) 
  con1.prop[i]<-sum(ifelse(abs(con1.pairs[,i])>abs(con1.true[i]), 1, 0)) #Number of 1000 slope values greater than observed
  con1.pval[i]<-2*(con1.prop[i]/1000) #Calculate 2-sided p-value for each slope difference
  nat1.prop[i]<-sum(ifelse(abs(nat1.pairs[,i])>abs(nat1.true[i]), 1, 0))
  nat1.pval[i]<-2*(nat1.prop[i]/1000) 
}

#1=alt*yr, 2=con*yr, 3=nat*yr
alt1.pval
con1.pval 
nat1.pval 
#Altered slope different from constructed, not natural
#All other contrasts significant
#alt=A, con=B, nat=A


alt<-dips1[dips1$pool.type == 'altered',] 
con<-dips1[dips1$pool.type == 'constructed',] 
nat<-dips1[dips1$pool.type == 'natural',]

plot(dips1$yr,dips1$dens, col=as.factor(c(dips1$pool.type)), ylab="Larval density (larvae/min)", xlab="Year", pch=c(dips1$pool.type))
abline(summary(lm(dens~yr, data=alt)), col=1, lwd=3, lty=1)#slope=-0.020
abline(summary(lm(dens~yr, data=con)), col=2, lwd=3, lty=2)#slope=0.001
abline(summary(lm(dens~yr, data=nat)), col=3, lwd=3, lty=3)#slope=-0.020
legend(2016,7,legend=unique(dips1$pool.type), lwd=2, col=unique(dips1$pool.type), pch=unique(dips1$pool.type), lty=c(1:3))

library(ggplot2)
library(cowplot)

#Plot with 95% CI envelopes around trendlines
ggplot(dips1, aes(yr, dens, fill=factor(pool.type), color=factor(pool.type), linetype=factor(pool.type)))+
  geom_smooth(method = "lm", alpha=0.15)+
  theme_cowplot()+
  labs(y="Larval density (larvae/min)", x="Year")

##Identify pond depth and area that are predictive of CTS use
#Create a data sheet with 1 row per pond
library(dplyr)
str(dips1)
ponds<-dips1 %>% 
  group_by_(.dots=c("site","pool.num")) %>% 
  summarize(x.ln.area = mean(log(max.area)),
            x.ln.depth = mean(log.max.depth),
            max.depth = mean(max.depth),
            max.area = mean(max.area),
            yrs.occ = sum(cts.pres)
            )
ponds<-ponds[order(-ponds$yrs.occ),]
View(ponds)
write.csv(ponds,"ponds-table.csv")

m1<-lm(yrs.occ ~ x.ln.depth, data=ponds)
summary(m1)
anova(m1)

int<-exp(-18.55) #Backtransform from log
slope<-exp(6.38) #Backtransform from log

plot(ponds$max.depth, ponds$yrs.occ, xlab="Maximum pond depth (cm)", ylab="Number of years occupied")
plot(ponds$x.ln.depth, ponds$yrs.occ, xlab="Ln(maximum pond depth (cm))", ylab="Number of years occupied")

ggplot(ponds, aes(max.depth, yrs.occ))+
  geom_smooth(alpha=0.15)+
  geom_point()+
  theme_cowplot()+
  labs(y="Number of years occupied", x="Maximum pond depth (cm)")

ggplot(ponds, aes(x.ln.depth, yrs.occ))+
  geom_smooth(alpha=0.15)+
  geom_point()+
  theme_cowplot()+
  labs(y="Number of years occupied", x="Ln(maximum pond depth (cm))")

plot(ponds$x.ln.area, ponds$yrs.occ, xlab="Ln(maximum pond area (m^2))", ylab="Number of years occupied")
plot(ponds$max.area, ponds$yrs.occ, xlab="Maximum pond area (m^2)", ylab="Number of years occupied")

#SVL analyses and summary
library(dplyr)
str(CTS_Lengths_Sonoma)

#get mean length of SCTs by pond and year
length<-CTS_Lengths_Sonoma %>% 
  group_by(.dots=c("id", "yr")) %>% 
  summarize(x.tl= mean(tl)
  )

#get denisty in same format 
density_condensed<-dips1 %>% 
  group_by(.dots=c("unique.id", "yr")) %>% 
  summarize(dens= dens
  )


length$unique.id <- length$id
#match rows 


str(length)

#needs to be id and year combo
joined_df <- merge(length, dips1, by.x = "unique.id", 
                   by.y = "unique.id", all.x = FALSE, all.y = TRUE)

