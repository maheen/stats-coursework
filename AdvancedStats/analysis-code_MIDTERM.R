#Written by Maheen Shermohammed

require(car)
require(MASS)
require(effects)
require(nnet)
require(lme4)
require(lmerTest)

####################### SES-Schools Data Analysis #########################
#import data
data <- read.csv("GATES_final.csv")
head(data)
attach(data)

freelunch_half <- factor(freelunch_half, labels = c("highSES", "lowSES"))

### --- Is working memory (WM) related to socioeconimic status (SES)?
### WM is ordinal, so will have to use an ordinal regression

#visualize data
boxplot(wm_span_05 ~ freelunch_half,main="WM span by SES")

#prepare data for ordinal regression
wm_span_05_ord <- as.ordered(wm_span_05)

#fit a proportional odds model:
polrFit <- polr(wm_span_05_ord ~ freelunch_half)
summary(polrFit)
# The odds of having a higher working memory decrease by
exp(-1.26)
## for low socioeconomic status vs. high
#####is there a way I can test for the significance of this effect???

#we can look at the effect plots to visualize this effect
plot(allEffects(polrFit))
plot(allEffects(polrFit), style = "stacked") 
plot(allEffects(polrFit, latent = TRUE))     
#the large error bars for lowSES on this last graph are interesting. This might
#suggest a larger variability in this effect for the lowSES group, whereas the highSES group
#is stable. It's possible that this could be a meaningful psychological difference (e.g. low
#SES participants have a varying range of protective factors, whereas highSES participants peak
#out.) Alternatively, this could reflect a sampling bias- it's possible that the high SES subjects
#came from a less variable subject pool

#from work on this data last semester, we have reason to believe opportunity of nurturance
#might play a role in this effect. Let's include it in the model.
polrFit2 <- polr(wm_span_05_ord ~ freelunch_half*opportunity_of_nurturance)
summary(polrFit2)
#with opp of nurturance included, the effect of SES has completely changed. This model also 
#has a lower AIC. Could opportunity of nurturance function as a mediator?
plot(allEffects(polrFit2))
plot(allEffects(polrFit2), style = "stacked") 
plot(allEffects(polrFit2, latent = TRUE))     
#it looks like in the low SES group, subjects with lower WM spands tend to have a higher 
#opp. of nurt. (which means more people rely on you for their well being - this is a bad thing)
#how could I check for the siginificance of this interaction?


### --- Is there a relationship between school type and the non-cognitive trait GRIT
##school type is nominal, so will have to use multinomial regression
school_group_new <- factor(school_group_new)
levels(school_group_new)<-c("low_perf","charter","high_perf")
school_group_new <- relevel(school_group_new, ref = "charter")

hist(Grit_Average)

#fit multinomial regression
fit.multinom <- multinom(school_group_new ~ Grit_Average)
fit.multinom.sum <- summary(fit.multinom, Wald = TRUE)
fit.multinom.sum

#estimate p value using the normal approximation of the Wald test
p <- (1 - pnorm(abs(fit.multinom.sum$Wald.ratios), 0, 1)) * 2 
p
## the results are not at all significant. If they were, though, this is how we would interpret:
exp(coef(fit.multinom))
## The odds ratio for a one-unit increase in GRIT is 1.022 for being in charter vs. low performing public school
## The odds ratio for a one-unit increase in GRIT is 0.621 for being in charter vs. high performing public school

plot(allEffects(fit.multinom))

detach(data)

####################### GSR Data Analysis #########################

## --- Data Preparation

#import data
gsr_data<-read.csv("~/Documents/Stats3450/midterm/GSR_data.csv")

#only include data that has at least one of all 6 conditions
is_complete<-NULL
for(i in gsr_data$sub){
  bool<-length(gsr_data$conds_fin[gsr_data$sub==i])
  is_complete<-c(is_complete,bool)
}
gsr_data<-subset(gsr_data,is_complete==6)

#have no a priori hypotheses about GSR for positive images. I will not be looking at 
#those blocks going forward.
gsr_data <- subset(gsr_data, conds_fin!='PP' & conds_fin!='UP')
gsr_data$conds_fin<-levels(droplevels(gsr_data$conds_fin))
head(gsr_data)

attach(gsr_data)

#split conditions variable into appropriate factors reflecting certainty and affect
cond_factor <- as.factor(conds_fin)
pred_factor <- cond_factor
levels(pred_factor)[match(c("PN", "PX"),levels(pred_factor))] <- "predictable"
levels(pred_factor)[match(c("UN", "UX"),levels(pred_factor))] <- "unpredictable"

valence_factor <- cond_factor
levels(valence_factor)[match(c("PN", "UN"),levels(valence_factor))] <- "neutral"
levels(valence_factor)[match(c("PX", "UX"),levels(valence_factor))] <- "negative"

subs_factor<-as.factor(subs_fin)

## --- Descriptives

####visualizing
hist(peaks_fin)
qqP <- qqnorm(peaks_fin)
abline(lm(qqP$y ~ qqP$x))

#collapsing across subject
subjMeans <- tapply(peaks_fin, subs_fin, mean)
hist(subjMeans)
qqP <- qqnorm(subjMeans)
abline(lm(qqP$y ~ qqP$x))
boxplot(subjMeans, id.n = 2)

####look for influential points / outliers

#scatterplots, label most extreme values
scatterplot(peaks_fin ~ pred_factor, id.n = 3) 
scatterplot(peaks_fin ~ valence_factor, id.n = 3)
scatterplot(peaks_fin ~ age, id.n = 3)

#fit a linear regression and look at studentized residuals
fitlin <- lm(peaks_fin ~ valence_factor + pred_factor)
summary(fitlin)
qqPlot(fitlin, id.n = 3)
residualPlots(fitlin, type = "rstudent")

#run explicit outlier test
outlierTest(fitlin)

#a nice summary graph
influenceIndexPlot(fitlin, id.n = 3)
#the results of the above analysis suggest that we might benefit from statistical methods 
#that are robust to outliers. This is a rather large data set, so I'm not too worried about
#1 or 2 outliers skewing the distribution, but we will demonstrate these methods.

#also, I know that this should be a mixed effects model with subjects as a random effect.
#I go on to use that later, but since we're just demonstrating the robust regression, I pretend
#like it isn't for now.

## --- Robust Regression
fitrob <- rlm(peaks_fin ~ valence_factor + pred_factor)
summary(fitrob)

## compare with lm fit
summary(fitlin)

#can I compare directly like this?
anova(fitrob,fitlin)


## --- Polynomial Regression
#What I'm really interested in is how the effects of valence and predictability on GSR 
#differ across age. I study adolescence, so in addition to looking at how age increases
#linearly, I have reason to believe that something special might happen in adolescence (which 
#would look quadratic function) or might emerge at adolescence and then maintain into 
#adulthood (we will call this "adolescent-emergent")

#let's start by just looking at the effects of valence & predictability on GSR response
#using a two-way within subjects anova
fit_2way_within=aov(peaks_fin~(pred_factor*valence_factor)+
                      Error(subs_factor/(pred_factor*valence_factor)))
summary(fit_2way_within)

#make linear and quadratic age regressors to add to the model
lin <- scale(age, scale = FALSE)
quad <- (lin^2)

#linear regression
fitlin=aov(peaks_fin~(pred_factor*valence_factor*lin)+
             Error(subs_factor/(pred_factor*valence_factor)))
summary(fitlin)

#quadratic regression
fitquad <- aov(peaks_fin ~ (pred_factor*valence_factor*lin+pred_factor*valence_factor*quad)+
                 Error(subs_factor/(pred_factor*valence_factor))) 
summary(fitquad)
#suggests a significant interaction between predictability and quadratic age:
U_P<-peaks_fin[pred_factor=='unpredictable']-peaks_fin[pred_factor=='predictable'] #diff between GSRs for predictable and unpredictable
plot(U_P ~ age[pred_factor=='predictable'], xlab='age', main='Diff b/w Unpredictable & Predictable Blocks')
##how do I plot the appropriate fitted curve here?


#adolescent-emergent regression
#this one seems more complicated. We want a function that curves into adolescence, and
#then plateaus. Perhaps one way to make this is to take a cubic function and then artificially
#flatten out the second inflection? This is what my advisor does - is there a better way?
adol_emer <- (lin)^3
plot(lin,adol_emer)
adol_emer[lin>=0] <- median(adol_emer)
plot(lin,adol_emer)

fit_adol_emer=aov(peaks_fin~(pred_factor*valence_factor*adol_emer)+Error(subs_factor/(pred_factor*valence_factor)))
summary(fit_adol_emer)

#which of the 3 models is best? let's re-fit these using lme4 so we can compare
fitlin_lm <- lmer(peaks_fin ~ (pred_factor*valence_factor*lin)+
                    (1|subs_factor) + (1|valence_factor:subs_factor) + (1|valence_factor:subs_factor))
fitquad_lm <- lmer(peaks_fin ~ (pred_factor*valence_factor*lin+pred_factor*valence_factor*quad)+
                     (1|subs_factor) + (1|valence_factor:subs_factor) + (1|valence_factor:subs_factor))
fitadolEmer_lm <- lmer(peaks_fin ~ (pred_factor*valence_factor*adol_emer)+
                         (1|subs_factor) + (1|valence_factor:subs_factor) + (1|valence_factor:subs_factor))

anova(fitlin_lm,fitquad_lm,fitadolEmer_lm)
#looks like the model with the adolescent-emergent regressor is the best?

#it would be interesting to do a segmented regression with the spline at age ~13, 
#but I can't figure out how to do it with this complicated mixed effects model. 
#The code below does not work:
spline_age <- 13-mean(age) #13 mean centered
require(segmented)
fitseg <- segmented(fitlin_lm, seg.Z = ~ lin, psi = spline_age)

#But I can try the brute force solution, which would be to run 2 single linear regressions. 
fitlin_piece1 <- lmer(peaks_fin[age<13] ~ (pred_factor[age<13]*valence_factor[age<13]*lin[age<13])+
                    (1|subs_factor[age<13]) + (1|valence_factor[age<13]:subs_factor[age<13]) + (1|valence_factor[age<13]:subs_factor[age<13]))
int1 <- 0.53170 
beta1 <- -0.38811
fitlin_piece2 <- lmer(peaks_fin[age>=13] ~ (pred_factor[age>=13]*valence_factor[age>=13]*lin[age>=13])+
                        (1|subs_factor[age>=13]) + (1|valence_factor[age>=13]:subs_factor[age>=13]) + (1|valence_factor[age>=13]:subs_factor[age>=13]))
int2 <- 0.670687
beta2 <- -0.028695 

plot(lin, peaks_fin, xlab = "Mean-Centered Age", ylab="peak",main = "GSR Response", cex = 0.5)   ## scatterplot
lines(c(spline_age,-6), c(int1, int1 + spline_age*beta1), col='red', lwd=2)                ## first regression line
abline(v = spline_age, lty = 2)
lines(c(spline_age, 9), c(int2 + spline_age*beta2, int2 + 9*beta2), col='red', lwd=2)    ## second regression line
#Did I do that right???

detach(gsr_data)
