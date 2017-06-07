##############################
#PREAMBLE TO SCRIPT FOR PAPER#
##############################

#Load necessary source code and libraries
#install.packages("nearfar") #package on CRAN; GPL-3 license
library(xtable)
library(nearfar)
library(AER)
library(MASS)
library(optmatch)

#For tables
source("https://raw.githubusercontent.com/joerigdon/Useful_Functions/master/Tables.R") #need rJava

#For figure
source("https://raw.githubusercontent.com/joerigdon/Useful_Functions/master/Figures.R")


##################
#SCRIPT FOR PAPER#
##################
set.seed(1212)

#2. An Illustrative Example
#Load illustrative data set (also attached with submission)
head(angrist)
df = angrist
df$IV = factor(df$IV) #make IV factor variable for table display

#Table 1: Variables in example summarized by quarter of birth
t1 = mktab(data=df, var.names=c("wage","educ","age","married","race"),ind.cat=c(0,0,0,1,1), group.name=c("qob"), cfn=describeMean, miss="always", pval=FALSE, tot="last", digit=1)
xtable(t1)

#Figure 2: Education and log weekly wage in example
pdf(width=8,height=5,"/Users/jrigdon/Box Sync/Rigdon/Sanjay/Paper1/Educ_wage.pdf")
boxplot(df$wage~df$educ, xlab="Years of education", ylab="Log weekly wage")
dev.off()

#Unadjusted OLS of log weekly wage on education to get first treatment effect estimate
unadj = lm(wage~educ,data=df)
summary(unadj)
confint(unadj)

#Check correlation of qob with residuals of unadj model above
cor.test(residuals(unadj),as.numeric(df$qob)) #0.004

#Adjusted OLS of log weekly wage on education to get second treatment effect estimate
adj = lm(wage~educ+age+married+race, data=df)
summary(adj)
confint(adj)


#Check correlation of qob with residuals of adj model above
cor.test(residuals(adj), as.numeric(df$qob)) #0.01

#OLS of education on IV to look at strength
l = lm(educ~factor(IV), data=df)
anova(l) #F of 1.97

#Same OLS with age, married, race added; look at partial F
m = lm(educ~age+married+race+factor(IV), data=df)
anova(m) #Partial F of 1.93

#2SLS regression to get third estimate of treatment effect in Table 3 (pre-near-far match)
m1i = ivreg(wage~age+married+race+educ|factor(IV)+age+married+race, data=df)
summary(m1i, diagnostics=TRUE)
confint(m1i)

#Near-far match
k = opt_nearfar(dta=df, trt="educ", covs=c("age", "married", "race"), iv="IV", trt.type="cont", imp.var=NA, tol.var=NA, adjust.IV=TRUE, max.time.seconds=300)
summary(k)

#Table 2: Variable summaries in example pre- and post-near-far match.  ASD = absolute standardized difference
#ASD pre-match
df$enc = ifelse(as.numeric(df$IV)<=2, 1, 0) #encouraged if qob is 3 or 4
match = matrix(NA, 515, 2)
match[,1] = which(df$enc==1)
match[1:485,2] = which(df$enc==0)
summ_matches(dta=df, covs=c("age", "married", "race"), iv="enc", match=match)

#ASD post-match
summary(k)

#Then filled in manually to latex file

#Table 3: Estimates (confidence intervals) for treatment effect of one-year increase in education on log weekly wages
#Apply same models for post-near-far match inference
df2 = df[as.numeric(k$match), ]

#Unadjusted OLS
unadj2 = lm(wage~educ, data=df2)
summary(unadj2)
confint(unadj2)

#Adjusted OLS
adj2 = lm(wage~educ+age+married+race, data=df2)
summary(adj2)
confint(adj2)

#2SLS
m2i = ivreg(wage~age+married+race+educ|IV+age+married+race, data=df2)
summary(m2i, diagnostics=TRUE)
confint(m2i)

#Compute effect ratio (not guaranteed to produce interval)
eff_ratio(dta=df, k$match, outc="wage", trt="educ", alpha=0.05) #0.29

#Figure 3: Kernel density plots of treatment effect estimates of education on earnings pre- and post-near-far match.  Red starts denot treatment effects summarized in Table 2

#Pre-match 2SLS estimates
get.coef = function(r1i) {
    s = summary(r1i)$coeff
    s[which(rownames(s)=="educ"),1]
}

all1 = c(
#Main effect of educ only
get.coef(ivreg(wage~educ|IV+age+married+race, data=df)),
#Educ plus one other
get.coef(ivreg(wage~age+educ|IV+age+married+race, data=df)),
get.coef(ivreg(wage~married+educ|IV+age+married+race, data=df)),
get.coef(ivreg(wage~race+educ|IV+age+married+race, data=df)),
#Educ plus two others
get.coef(ivreg(wage~age+married+educ|IV+age+married+race, data=df)),
get.coef(ivreg(wage~married+race+educ|IV+age+married+race, data=df)),
get.coef(ivreg(wage~race+age+educ|IV+age+married+race, data=df)),
#Educ plus all others
get.coef(ivreg(wage~age+married+race+educ|IV+age+married+race, data=df))
)

#Post-match 2SLS estimates
all2 = c(
#Main effect of educ only
get.coef(ivreg(wage~educ|IV+age+married+race, data=df2)),
#Educ plus one other
get.coef(ivreg(wage~age+educ|IV+age+married+race, data=df2)),
get.coef(ivreg(wage~married+educ|IV+age+married+race, data=df2)),
get.coef(ivreg(wage~race+educ|IV+age+married+race, data=df2)),
#Educ plus two others
get.coef(ivreg(wage~age+married+educ|IV+age+married+race, data=df2)),
get.coef(ivreg(wage~married+race+educ|IV+age+married+race, data=df2)),
get.coef(ivreg(wage~race+age+educ|IV+age+married+race, data=df2)),
#Educ plus all others
get.coef(ivreg(wage~age+married+race+educ|IV+age+married+race, data=df2))
)

#Get observed estimates in Table 3
df1 = approxfun(density(all1))
p1 = get.coef(ivreg(wage~age+married+race+educ|IV+age+married+race, data=df))

df11 = approxfun(density(all2))
p11 = get.coef(ivreg(wage~age+married+race+educ|IV+age+married+race, data=df2))

#Save figure
pdf("/Users/jrigdon/Box Sync/Rigdon/Sanjay/Paper1/Fig_2016-07-18.pdf",8,8)
densp(obj.list=list(all1,all2), y=7, xl=c(-0.3, 0.6), cols=c("black","black"), ltys=c(1,2), xtitle="Estimated 2SLS effect of education on earnings", mtitle="", leg=c("pre-match","post-match"), px=0.3, py=6)
points(p1, df1(p1),col=2,pch="*",cex=2)
points(p11, df11(p11),col=2,pch="*",cex=2)
dev.off()


#3. Simulated Example
#Simulated example
set.seed(172)
dta = mvrnorm(1000,c(10,10,10),matrix(c(1,-0.5,0.5,-0.5,1,0.5,0.5,0.5,1),3,3))
Zstar = dta[,1] #Part of Z that is correlated with unmeas conf
X.unmeas = dta[,2] #Unmeas conf
X.meas = dta[,3] #Meas conf
IV = rnorm(1000,10,1) #Instrumental variable
Z = 1+5*Zstar+3*X.meas+1*IV+rnorm(1000,0,10) #Observed treatment
Y = 1+1*Z+1*X.meas+5*X.unmeas+rnorm(1000,0,20) #Outcome
df.sim = data.frame(Y=Y,Z=Z,IV=IV,X=X.meas) #set up for near-far match

head(df.sim)

#How strong is IV
s = lm(Z~X+IV, data=df.sim)
anova(s)

#True model
a = lm(Y~Z+X.unmeas+X.meas)
summary(a)
confint(a)

#Model without unmeasured confounder
b = lm(Y~Z+X, data=df.sim)
summary(b)
confint(b)

#2SLS estimate using weak IV
c = ivreg(Y~X+Z|IV+X, data=df.sim)
summary(c, diagnostics=TRUE)
confint(c)

#Apply near-far matching
nf = opt_nearfar(dta=df.sim, trt="Z", covs="X", iv="IV", trt.type="cont", imp.var=NA, tol.var=NA, adjust.IV=TRUE, max.time.seconds=300)
summary(nf)

#Analytic cohort after near-far matching
df.sim2 = df.sim[as.numeric(nf$match), ]

#True model
a2 = lm(Y[as.numeric(nf$match)]~Z[as.numeric(nf$match)]+X.unmeas[as.numeric(nf$match)]+X.meas[as.numeric(nf$match)])
summary(a2)
confint(a2)

#Model without unmeasured confounder
b2 = lm(Y~Z+X, data=df.sim2)
summary(b2)
confint(b2)

#2SLS estimate using weak IV after near-far
d = ivreg(Y~X+Z|IV+X, data=df.sim2)
summary(d, diagnostics=TRUE)
confint(d)

eff_ratio(dta=df.sim, match=nf$match, outc="Y", trt="Z", alpha=0.05)

#Compare standard errors pre- and post-match
sqrt(vcov(c)[3,3]) #pre
sqrt(vcov(d)[3,3]) #post


#4. Key Analytical Choices
aa = angrist[, 5:7]
head(aa)
X2 = matrix(as.numeric(as.matrix(aa)),dim(aa)[1],dim(aa)[2])
jj = smahal(X=X2)
round(jj[1:5, 1:5], 2)

#Apply calipers
jj2 = calipers(distmat=jj, variable=aa$age, tolerance=0.2)
round(jj2[1:5, 1:5], 2)

#5. Special Scenarios
#5.1 Binary Treatment
df3 = data.frame(wage=df$wage, educ=ifelse(df$educ>12,1,0), IV=factor(df$IV), age=df$age, married=factor(df$married), race=factor(df$race))
head(df3)

#Pre-nearfar
m3 = glm(educ~age+married+race+IV, data=df3, family=binomial)
anova(m3)

#Propensity score match

#Pair matching within a propensity score caliper
ppty = glm(educ~age+married+race, family=binomial, data=df3)
#Compute MHD with caliper
mhd = match_on(educ~age+married+race, data=df3) + caliper(match_on(ppty), 2)
pm = pairmatch(mhd, data=df3)
#as.numeric(names(pm[!is.na(pm)])) #keep for analysis

#Near-far match
set.seed(44)
nf2 = opt_nearfar(dta=df3, trt="educ", covs=c("age", "married", "race"), iv="IV", trt.type="bin", imp.var=NA, tol.var=NA, adjust.IV=TRUE, max.time.seconds=300)
summary(nf2)

#Examine effect estimates from near-far and PSM
head(nf2$match)

df3a = df3[as.numeric(nf2$match), ]
m51a = lm(wage ~ educ + age + married + race, data=df3a)
summary(m51a)
confint(m51a)

df3b = df3[as.numeric(names(pm[!is.na(pm)])), ]
m51b = lm(wage ~ educ + age + married + race, data=df3b)
summary(m51b)
confint(m51b)

#Pre-match estimate
m51 = lm(wage ~ educ + age + married + race, data=df3)
summary(m51)
confint(m51)


#5.2 Prioritized Variables
set.seed(33)
nf3 = opt_nearfar(dta=df3, trt="educ", covs=c("age", "married", "race"), iv="IV", trt.type="bin", imp.var=c("age", "married"), tol.var=c(0.3,0.2), adjust.IV=TRUE, max.time.seconds=300)
summary(nf3)

#5.3 No Adjustment for Measured Confounders
set.seed(53)
nf4 = opt_nearfar(dta=df3, trt="educ", covs=c("age", "married", "race"), iv="IV", trt.type="bin", imp.var=NA, tol.var=NA, adjust.IV=FALSE, max.time.seconds=300)
summary(nf4)




