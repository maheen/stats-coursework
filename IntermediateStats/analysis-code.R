##Written by Julia Leonard & Maheen Shermohammed

require(pastecs)
require(nlme)
require(lattice)
require(betareg)
if (!require(scales)) {install.packages("scales"); require(scales)}
require(plotrix)
require(stats)
require(effects)
require(MASS)
if (!require(ggplot2)) install.packages("ggplot2"); require(ggplot2)
require(car)
require(coin)
require(pwr)

#import data
data <- read.csv("GATES_final.csv")
head(data)
attach(data)

freelunch_half <- factor(freelunch_half, labels = c("highSES", "lowSES"))
sex <- factor(sex, labels = c("female", "male"))

##### -- Does sex differ by SES group? - CHI SQUARED TEST #####

#make factors
freelunch_half <- factor(freelunch_half, labels = c("highSES", "lowSES"))
sex <- factor(sex, labels = c("female", "male"))

#look at table
chi_data<-table(sex,freelunch_half)
chi_data

#visualize data
op <- par(mfrow = c(1,2))
spineplot(chi_data, xlab = "sex", ylab = "lunch group", main = "sex by group")
barplot(t(chi_data), main = "sex by group", cex.names = 0.8, xlab = "Sex", 
        legend.text = colnames(chi_data), args.legend = list(x = "topleft"), beside = TRUE)
par(op)

#run chi squared test
resChi <- chisq.test(chi_data)
resChi  #sex is not significantly different by SES group

## permutation version: (in this case same as parametric chi^2, large frequencies)
chisq_test(as.table(chi_data), distribution = approximate(B=1000))  



##### -- Is working memory (WM) related to socioeconimic status (SES)? #####

###Let's start with descriptives
by(wm_span_05, freelunch_half, stat.desc, norm = TRUE) #WM not normally distributed for high SES

op <- par(mfrow = c(2,1))
hist(wm_span_05[freelunch_half == "highSES"], main = "Histogram High SES", xlab = "WM")
hist(wm_span_05[freelunch_half == "lowSES"], main = "Histogram Low SES", xlab = "WM")
par(op) #once again, WM not normally distributed for high SES

boxplot(wm_span_05 ~ freelunch_half, ylab = "WM", xlab = "SES")
points(tapply(wm_span_05, freelunch_half, mean, na.rm = TRUE), pch = 19, col = "red", type = "b", cex = 0.9) 

#Ok. Working memory is not normally distributed in the high SES group. For now let's pretend it is...

###Fit a model (which is just one-way anova / t-test for now)
fit_SES_WM1 <- lm(wm_span_05 ~ freelunch_half)
summary(fit_SES_WM1) 
#high SES has significantly greated WM than low SES (SES explains ~ 10% of the variance in WM)

#what if we account for gender?
#we will include this in the model as a random effect

###First visualize the individual regressions
xyplot(wm_span_05 ~ freelunch_half|sex, mainv = "Individual Regressions (Sex)", panel=function(x, y){
  panel.xyplot(x, y)
  panel.lmline(x, y, lty=2)
})

###Next we want to run and compare models
fitM1 <- gls(wm_span_05 ~ freelunch_half, method = "ML", na.action = na.exclude)
summary(fitM1)
#note RSE: 1.46383 
fitM2 <- lme(wm_span_05 ~ freelunch_half, random = ~1|sex, method = "ML", na.action = na.exclude)
summary(fitM2)
#RSE did not improve at all

#explicitly comparing models
anova(fitM1, fitM2) #there is no difference. Gender does not add to our model, we will drop it.

###What if we account for age?
hist(age)

fitM3 <- lme(wm_span_05 ~ freelunch_half, random = ~1|age, method = "ML", na.action = na.exclude)
summary(fitM3)
#RSE is a little better: 0.8490913

#explicitly compare models
anova(fitM1, fitM3) #difference is not significant

#I find this pretty surprising, as it suggests that there is not a significant relationship 
#between WM and age. I tested this explicitly, and found it to be true:
fit_age_WM <- lm(wm_span_05 ~ age)
summary(fit_age_WM)

#But, when we take another look at the distribution of the age variable, we find that the age
#range is very small, which puts the above results into better perspective. There might not
#be enough variation in the ages in this data set.
hist(age)


###WM is not normally distributed, though. In fact, it is bounded, so a beta regression is
###more appropriate.

#first we need to rescale WM to go from 0 to 1
wm_rescaled <- scales:::rescale(wm_span_05, to = c(0, 1))
hist(wm_rescaled)
boxplot(wm_rescaled ~ freelunch_half, ylab = "WM", xlab = "SES")

#can't have any NA's or exact 0's or 1's, so changing these values slightly
wm_rescaled <- replace(wm_rescaled, wm_rescaled==0, 0.01)
wm_rescaled <- replace(wm_rescaled, wm_rescaled==1, 0.99)
wm_rescaled_noNA <- wm_rescaled[!is.na(wm_rescaled)]

#look at beta distribution
bpar <- fitdistr(wm_rescaled_noNA, densfun = "beta", start = list(shape1 = 0.5, shape2 = 0.5))
xvals <- seq(0, 1, 0.01)
densvals <- dbeta(xvals, shape1 = bpar$estimate[1], shape2 = bpar$estimate[2])
hist(wm_rescaled, freq = FALSE)
lines(xvals, densvals, col = "red")

#fit the beta regression
fitBeta1 <- betareg(wm_rescaled ~ freelunch_half)
summary(fitBeta1)
#similar conclusion to previous analysis: SES is significantly related to WM

#data from our midterm suggests that the social variable opportunity_of_nurturance might
#also be related to WM
plot(wm_span_05 ~ opportunity_of_nurturance)

#one-way anova with covariate as response is not significant, which is good.
summary(aov(opportunity_of_nurturance ~ freelunch_half)) 

#we will include opportunity_of_nurturance in our model
fitBeta2 <- betareg(wm_rescaled ~ freelunch_half*opportunity_of_nurturance)
summary(fitBeta2)

#the effect got much smaller when we include this social variable, which is not significantly
#related to WM on its own:
fitBeta3 <- betareg(wm_rescaled ~ opportunity_of_nurturance)
summary(fitBeta3)

#could it be that opportunity_of_nurturance functions as a mediator? I don't know what the specific
#test for this is, but we can look at the graphs for each SES group to get a better sense of 
#what is going on.
plot(jitter(opportunity_of_nurturance[freelunch_half == "highSES"],1), jitter(wm_span_05[freelunch_half == "highSES"],1), xlab = "opportunity_of_nurturance", ylab="WM Span", main="GroupXOpportunity Interaction on WM", xlim = c(7,18), pch=19, col= adjustcolor("blue",alpha=0.5))
points(jitter(opportunity_of_nurturance[freelunch_half == "lowSES"],1), jitter(wm_span_05[freelunch_half == "lowSES"],1), xlab = "opportunity_of_nurturance", ylab="WM Span", pch=19, col=adjustcolor("red",alpha=0.5))
ablineclip(lm(wm_span_05[freelunch_half == "highSES"] ~ opportunity_of_nurturance[freelunch_half == "highSES"]), x1=8,x2=16, col = "blue", lwd = 2)
ablineclip(lm(wm_span_05[freelunch_half == "lowSES"] ~ opportunity_of_nurturance[freelunch_half == "lowSES"]), x1=8,x2=16, col = "red", lwd = 2)
par('usr') #figure out plot dimensions
legend(x=14.2,y=5.4,legend=c("lowSES","highSES"),fill=c("red","blue"),cex=.8,text.font=2)



##### -- Can we predict socioeconimic status (SES) with the social variables we collected? #####

#visualize social variables with histograms
social_variables <- c("high_contact_roles","number_people_social_network","embedded_network","attachment","social_integration","opportunity_of_nurturance","reassurance_of_worth","reliable_alliance","guidance")
op <- par(mfrow = c(3,3))
for (i in 1:length(social_variables)){
  hist(data[,social_variables[i]],xlab=social_variables[i],main=social_variables[i])
  n<-shapiro.test(data[,social_variables[i]])
  mtext(bquote(italic(p) == .(format(n$p.value, digits = 3))))
}
par(op)

##I need to fit a logistic regression

#have to remove all NA values (because I will eventually use 'step')
GATES_data = data.frame(freelunch_half,high_contact_roles,number_people_social_network,embedded_network,attachment,social_integration,opportunity_of_nurturance,reassurance_of_worth,reliable_alliance,guidance)
GATES_data<-GATES_data[complete.cases(GATES_data),]

#I can try to include all of the parameters, but there's too many and nothing comes out siginificant
fitSES <- glm(freelunch_half ~ ., family = binomial, data=GATES_data)
summary(fitSES)

plot(allEffects(fitSES))

#since I don't know what parameters to use, maybe a more principled way is to do a stepwise regression based on AIC
fitStep <- stats:::step(fitSES, direction = "backward")

#I find that the best model includes just high_contact_roles,embedded_network & guidance
#since removing any of these doesn't make the model better
fitSES_opt <- glm(freelunch_half ~ high_contact_roles + embedded_network + guidance, family = binomial, data=GATES_data)
summary(fitSES_opt)
#now embedded_network & guidance are significant

plot(allEffects(fitSES_opt))


##### -- What is the relationship with brain measure (amygdala volume)? #####

### ANCOVA###
# For the midterm I found an interaction with WM and SES on the R amyg and R basal amyg
# but I actually ran the wrong model. Now I know I should run an ANCOVA to do this
# becuase I have one categorical (SES) and one metric (WM) predictor.
#Also last time I didn't control for ICV which needs to be accounted for. 

#from the midterm:
anova<-aov(wm_span_05~freelunch_half*right_amyg_volume+estimatedtotalintracranialvol+sex)
summary(anova)
anova<-aov(wm_span_05~freelunch_half*R_basal_nucleus+estimatedtotalintracranialvol)
summary(anova)

data$lunch[data$freelunch_half=="0"] <- "high_SES"
data$lunch[data$freelunch_half=="1"] <- "low_SES"

#graph the interactions
ggplot(data, aes(x=data$R_basal_nucleus, y=data$wm_span_05, color=data$lunch)) + geom_point(shape=1) +
  scale_colour_hue(l=50) +# Use a slightly darker palette than normal
  ggtitle ("WM and R basal nucleus interaction by SES")+
  geom_smooth(method=lm,   # Add linear regression lines
              se=FALSE,    # Don't add shaded confidence region
              fullrange=T) # Extend regression lines

ggplot(data, aes(x=data$right_amyg_volume, y=data$wm_span_05, color=data$lunch)) + geom_point(shape=1) +
  scale_colour_hue(l=50) + # Use a slightly darker palette than normal
  ggtitle ("WM and R amygdala interaction by SES")+
  geom_smooth(method=lm,   # Add linear regression lines
              se=FALSE,    # Don't add shaded confidence region
              fullrange=T) # Extend regression lines

#hist of variables
par(mfrow=c(1,3))
hist(wm_span_05)
hist(right_amyg_volume)
hist(R_basal_nucleus)
par(op)

## boxplots with WM and lunch status as predictors
par(mfrow=c(1,2))
boxplot(right_amyg_volume~wm_span_05, ylab = "r amyg")
points(meantab[1,], col = "red", type = "b", pch = 19)
boxplot(R_basal_nucleus~wm_span_05 , ylab = "r basal amyg")
points(meantab[1,], col = "red", type = "b", pch = 19)
par(op)

par(mfrow=c(1,2))
boxplot(right_amyg_volume~data$lunch, ylab = "r amyg")
points(meantab[1,], col = "red", type = "b", pch = 19)
boxplot(R_basal_nucleus~data$lunch, ylab = "r basal amyg")
points(meantab[1,], col = "red", type = "b", pch = 19)
par(op)

## one-way ANOVA with covariate as response
summary(aov(wm_span_05~data$lunch))   ## significant. Not good, predictors are not independent...
#will go on...

#add predictors in model
summary(aov(right_amyg_volume~wm_span_05+data$lunch))
summary(aov(right_amyg_volume~data$lunch+wm_span_05))
summary(aov(R_basal_nucleus~wm_span_05+data$lunch))
summary(aov(R_basal_nucleus~data$lunch+wm_span_05)) #different order leads to different results - cant use type I

#But I'm most interested in the interaction between WM, SES, and the amygdala (controlling for ICV)

## check variance homogeneity assumptions, including for the interactions
leveneTest(wm_span_05~freelunch_half) 
leveneTest(right_amyg_volume, interaction(wm_span_05,data$lunch))   
leveneTest(R_basal_nucleus, interaction(wm_span_05,data$lunch))
#looks good

model1 <- aov(R_basal_nucleus ~ wm_span_05*freelunch_half)
summary(model1)
model2 <- aov(right_amyg_volume ~ wm_span_05*freelunch_half)
summary(model2)

## normality, patterns residuals
hist(residuals(model1))
qqnorm(residuals(model1))
qqline(residuals(model1))
shapiro.test(residuals(model1)) #technically, this violates assumptions

hist(residuals(model2))
qqnorm(residuals(model2))
qqline(residuals(model2))
shapiro.test(residuals(model2))

#want to add ICV to the model to control for brain size
model1<- aov(right_amyg_volume ~ wm_span_05*data$lunch +estimatedtotalintracranialvol)  #multilevel models
summary(model1)  #interaction still significant
model2<- aov(R_basal_nucleus ~ wm_span_05*data$lunch +estimatedtotalintracranialvol)  #multilevel models
summary(model2) #interaction still significant

#but dont have a balanced design so do ANOVA type II
## Type II ANOVA:
fitaov2 <- Anova(model1, type = "II")       ## Anova() from the car package
fitaov2  #interaction still sig
fitaov2 <- Anova(model2, type = "II")     
fitaov2  #interaction still sig

#significant so main effect is not interpretable, only interaction. Let's look at the plots again
ggplot(data, aes(x=data$wm_span_05, y=data$R_basal_nucleus, color=data$lunch)) + geom_point(shape=1) +
  scale_colour_hue(l=50) +# Use a slightly darker palette than normal
  ggtitle ("WM and R basal nucleus interaction by SES")+
  geom_smooth(method=lm,   # Add linear regression lines
              se=FALSE,    # Don't add shaded confidence region
              fullrange=T) # Extend regression lines

ggplot(data, aes(x=data$wm_span_05, y=data$right_amyg_volume, color=data$lunch)) + geom_point(shape=1) +
  scale_colour_hue(l=50) + # Use a slightly darker palette than normal
  ggtitle ("WM and R amygdala interaction by SES")+
  geom_smooth(method=lm,   # Add linear regression lines
              se=FALSE,    # Don't add shaded confidence region
              fullrange=T) # Extend regression lines

#I'm interested in explicitly testing these correlations (WM with amygdala volume) in each of
#the two groups. However, we've only done post-hoc t-tests. How would we do a post-hoc correlation
#test? Would I even need to correct for this? Just did a normal correlation for now:
free_lunch_data <- data[data$freelunch_half==1,]
paid_lunch_data <- data[data$freelunch_half==0,]

cor.test(free_lunch_data$wm_span_05,free_lunch_data$R_accessory_basal) #not sig
cor.test(free_lunch_data$wm_span_05,free_lunch_data$right_amyg_volume) #not sig
cor.test(paid_lunch_data$wm_span_05,paid_lunch_data$R_accessory_basal) #sig
cor.test(paid_lunch_data$wm_span_05,paid_lunch_data$right_amyg_volume) #sig
#the correlations b/w WM and amygdala volume are only significant in the high SES group

detach(data)

######### ----- Different Set of Questions Looking at Effort ----- #####

############Effort study part 1 ##################### 
data<-read.csv("effort_1_final.csv",header = TRUE, sep = ",")

head(data)
dim(data)
attach(data)

#look at data
table(Age, gender)  
table(Answer1)# Number of right answers on first trial = 13
table(answer2) # Number of right answers on second trial = 13
table(average) #10 kids got 50% correct, 8 got 100% correct

#one tailed and two tailed binomial tests
binom <-binom.test(sum(data$Answer1+data$answer2, na.rm=TRUE),36,alternative=c("greater")) # One-tailed test
binom 
binom.test(sum(data$Answer1+data$answer2,na.rm=TRUE),36) #Two-tailed test - kids scored above chance
binom 

#one and two tailed is above chance. This analysis looks at % of children's answers above chance, but does not look at 
# if kids got this 100% correct above chance. Both of these questions were asking about the same
# construct so it makes sense to see if they got them 100% above chance

#now lets look at kids who got both questions correct and see if they are above chance

binom.test(8,18,p = .25) #no, number of kids at ceiling is not above chance

#do this for q 1 and 2 separately to p =0.5
binom.test(sum(data$Answer1,na.rm=TRUE),18,.5)
binom.test(sum(data$answer2,na.rm=TRUE),18,.5)

#so number of kids who got both questions correct is not above chance. Number of kids who got question1 
# and question 2, tested individually, is also not above chance

detach(data)

############Effort study part 2 ##################### 
data<-read.csv("final_other_effort.csv",header = TRUE, sep = ",")

head(data)
attach(data)
summary(data$age)

#groups
effort <- data[data$condtion_num==1,]
no_effort <- data[data$condtion_num==2,]

#age and sex
table(effort$age, effort$gender) 
table(no_effort$age, no_effort$gender) 

#seconds_until_help - look at data histogram
x11()
op <- par(mfrow=c(2,1)) 
hist(effort$seconds_until_help, freq=TRUE, breaks = 10)
hist(no_effort$seconds_until_help,freq=TRUE, breaks = 10)
par(op)

ggplot(data, aes(x=seconds_until_help, fill=data$condtion)) + geom_density(alpha=.3)
ggplot(data, aes(x=seconds_until_help, fill=data$condtion)) + geom_histogram(alpha = 0.5, position = 'identity')

#normality
shapiro.test(effort$seconds_until_help)
shapiro.test(no_effort$seconds_until_help)
#not normal!!! so use medians!

#boxplot
x11()
op <- par(mfrow=c(1,2)) 
boxplot(effort$seconds_until_help, main ="effort")
boxplot(no_effort$seconds_until_help, main ="no effort")
par(op)

ggplot(na.omit(data[,c("seconds_until_help","condtion")]),aes(x=condtion, y=seconds_until_help, fill=condtion)) + geom_boxplot()

#ttests 
#seconds_until_help
t.test(data$seconds_until_help~data$condtion_num) #cant use this approach becuase data is not normally distributed

#ttest permutation test
oneway_test(seconds_until_help~condtion)

#not really normally distributed so utest to compare medians
wilcox.test(data$seconds_until_help ~ data$condtion_num)

#permutation test?
wilcox_test(seconds_until_help~condtion)

#I'm close to seeing a group difference, but the results arn't significant. Perhaps I need to 
#collect more data? Do a power analysis!

#power analysis
pwr.t.test(d = 0.8, power = 0.8) #So I want 26 kids/ group, so I need to collect more...

