#Written by Maheen Shermohammed

if (!require(mclust)) {install.packages("mclust"); require(mclust)}
if (!require(cluster)) {install.packages("cluster"); require(cluster)}
if (!require(plot3d)) {install.packages("plot3d"); require(plot3d)}
require(rgl)
if (!require(flexmix)) {install.packages("flexmix"); require(flexmix)}
require(MASS)

#### --- Data Preparation -----

#import data
gsr_data<-read.csv("~/Documents/Stats3450/midterm/GSR_data.csv")
#gsr_data<-read.csv("~/Documents/Clocks_GSR/aggreg_event_data.csv")

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

#Make appropriate factors
gsr_data$conds_fin <- as.factor(gsr_data$conds_fin)
gsr_data$pred_factor <- gsr_data$conds_fin
levels(gsr_data$pred_factor)[match(c("PN", "PX"),levels(gsr_data$pred_factor))] <- "predictable"
levels(gsr_data$pred_factor)[match(c("UN", "UX"),levels(gsr_data$pred_factor))] <- "unpredictable"

gsr_data$valence_factor <- gsr_data$conds_fin
levels(gsr_data$valence_factor)[match(c("PN", "UN"),levels(gsr_data$valence_factor))] <- "neutral"
levels(gsr_data$valence_factor)[match(c("PX", "UX"),levels(gsr_data$valence_factor))] <- "negative"

gsr_data$subs_fin<-as.factor(gsr_data$subs_fin)

#Collapse across valence. Won't be looking at that today.
aggreg_gsr <- aggregate(cbind(gsr_data$peaks_fin,gsr_data$counts_fin), list(Subject = gsr_data$subs_fin, Pred = gsr_data$pred_factor), mean)
aggreg_gsr <- aggreg_gsr[with(aggreg_gsr, order(Subject)), ]
names(aggreg_gsr)[3:4] <- c("Peaks","Counts")
aggreg_gsr$Age <- gsr_data$age[gsr_data$valence_factor=='neutral']
head(aggreg_gsr)

#### --- Data Preparation Finished -----

#In the midterm, we found an effect of predictability. We will start the analysis with 
#the difference between unpredictable and predictable GSR response (U_P) as our dependent variable
aggreg_pred <- subset(aggreg_gsr,Pred=='predictable')
aggreg_unpred <- subset(aggreg_gsr,Pred=='unpredictable')
U_P<-aggreg_unpred$Peaks-aggreg_pred$Peaks #diff between GSRs for predictable and unpredictable
plot(U_P ~ aggreg_pred$Age)

## let's use the E-M method for 1-10 clusters by age and U_P
UP_bandclust <- Mclust(cbind(aggreg_pred$Age,U_P), G = 1:10)
UP_bandclust               ## it selects 2 clusters based on BIC 
UP_bandclust$BIC
plot(UP_bandclust, "BIC")   

#a 2 cluster solution...
UP_mclust <- Mclust(cbind(aggreg_pred$Age,U_P), G = 2)
summary(UP_mclust, parameters = TRUE)  ## the parameters from the mixture model
plot(UP_mclust, "classification", xlab = "Age")    ## classification
#it's pretty straightforward to see what's going on here. One of the clusters is primarily
#adults, the other is primarily kids

#a 3 cluster solution, just to see what it would look like..
UP_mclust_3 <- Mclust(cbind(aggreg_pred$Age,U_P), G = 3)
summary(UP_mclust_3, parameters = TRUE)
plot(UP_mclust_3, "classification")


#Now, what if we included scores on a scale of Intolerance to Uncertainty (IUS)?
IUS_data <- read.csv("~/Documents/Clocks_GSR/IUS.csv")
IUS_data <- subset(IUS_data,select = c(Subject_ID,Factor1,Factor2,TotalScore,Average))
IUS_data <- subset(IUS_data,Subject_ID %in% unique(aggreg_gsr$Subject))
hist(IUS_data$Average)
UP_IUS_data <- as.data.frame(na.omit(cbind(aggreg_pred$Age,U_P,IUS_data$Average)))
names(UP_IUS_data) <- c("Age","Unp_Pre","IUS_Avg")

UP_bandclust_IUS <- Mclust(UP_IUS_data, G = 1:10)
UP_bandclust_IUS               ## it selects 2 clusters based on BIC 
UP_bandclust_IUS$BIC
plot(UP_bandclust_IUS, "BIC")  

## visualize the clusters
UP_mclust_IUS <- Mclust(UP_IUS_data, G = 2)
table(UP_mclust_IUS$classification)
#...by first doing a PCA
clusplot(UP_IUS_data, UP_mclust_IUS$classification, main = "UP-IUS K-Means Solution", lines = 0, color = TRUE, labels = 5, cex = 0.8, cex.txt = 0.2)
#...by plotting them in 3d
plot3d(UP_IUS_data, col = UP_mclust_IUS$classification, pch = 19)
#...and in 2d, for comparison
plot(U_P~aggreg_pred$Age,col = UP_mclust_IUS$classification)

#Let's repeat the analysis with flexmix so we can use that going forward...
set.seed(8657)
fmstep1 <- stepFlexmix(~1, data = UP_IUS_data, k = 1:4, nrep = 10,
                       model=list(FLXMRglm(Unp_Pre ~ ., family = "gaussian"),
                                  FLXMRglm(Age ~ ., family = "gaussian"),
                                  FLXMRglm(IUS_Avg ~ ., family = "gaussian")))
plot(1:4, BIC(fmstep1), type = "b", pch = 19, xaxt = "n", xlab = "Number of Clusters", ylab = "BIC", main = "BIC Cluster Evaluation")
fm1 <- getModel(fmstep1, "BIC")               ## pull out the best model
fm1
#hmm, that didn't give me a 2 cluster solution. Am I doing something wrong?

#I can just move on with a 2 cluster solution going forward, for demonstration purposes..
set.seed(8697)
fmstep2 <- stepFlexmix(~1, data = UP_IUS_data, k = 2, nrep = 10,
                       model=list(FLXMRglm(Unp_Pre ~ ., family = "gaussian"),
                                  FLXMRglm(Age ~ ., family = "gaussian"),
                                  FLXMRglm(IUS_Avg ~ ., family = "gaussian")))
summary(fmstep2)
table(fmstep2@cluster)
plot3d(na.omit(UP_IUS_data), col = fmstep2@cluster, pch = 19)

#What would happen if I used IUS as a concomitant variable instead?
set.seed(8697)
fmstep3 <- stepFlexmix(~1, data = UP_IUS_data, k = 2, nrep = 10,concomitant=FLXMRglm(~IUS_Avg),
                       model=list(FLXMRglm(Unp_Pre ~ ., family = "gaussian"),
                                  FLXMRglm(Age ~ ., family = "gaussian")))
table(fmstep3@cluster)
plot(Unp_Pre ~ Age, data = UP_IUS_data, col = fmstep3@cluster, main = "IUS as concomitant") #this looks like a bad cluster solution, basically just showing the more extreme data points
plot(Unp_Pre ~ Age, data = UP_IUS_data, col = fmstep2@cluster, main = "IUS as part of clustering") #for comparison with the previous model


### Ok, U_P is a weird dependent variable. Ideally, we would like to have the actual GSR
### response as the DV, and have predictability as a repeated measure. We will build up to
### this model slowly.

#Let's start by building a mixture model with the GSR response collapsing across predictable
#and upredictable. 
allGSR <- aggregate(cbind(gsr_data$peaks_fin,gsr_data$counts_fin,gsr_data$age), list(Subject = gsr_data$subs_fin), mean)
names(allGSR)[-1] <- c("Peaks","Counts","Age")
head(allGSR)

## Mclust on Peaks and Age for 1-10 clusters
all_bandclust <- Mclust(allGSR[,c(4,2)], G = 1:10)
all_bandclust               ## it selects 2 clusters based on BIC 
all_bandclust$BIC
plot(all_bandclust, "BIC")   

#it selected a 2 cluster solution...
all_mclust <- Mclust(allGSR[,c(4,2)], G = 2)
summary(all_mclust, parameters = TRUE)  ## the parameters
plot(all_mclust, "classification") 
## once again, one of these clusters is primarily adults, and the other is primarily minors

#repeat with flexmix..
set.seed(45345)
all_fm1 <- stepFlexmix(~1, data = allGSR, k = 1:4, nrep = 10,
                       model=list(FLXMRglm(Peaks ~ ., family = "gaussian"),
                                  FLXMRglm(Age ~ ., family = "gaussian")))
plot(1:4, BIC(all_fm1), type = "b", pch = 19, xaxt = "n", xlab = "Number of Clusters", ylab = "BIC", main = "BIC Cluster Evaluation")
axis(1, at = 1:4, labels = 1:4)
fm_all <- getModel(all_fm1, "BIC")  ##best model is still 2 clusters!
fm_all
plot(Peaks ~ Age, data = allGSR, col = fm_all@cluster) #same as in mclust

#### Now we can also add in the *number* of responses a person had, which will require
#### a poisson distribution
allGSR$roundedCounts <- round(allGSR$Counts) #these are average counts, so I have to round
hist(allGSR$roundedCounts)

#First just checking to see how good of an approximation the poisson distribution is
mean(allGSR$roundedCounts)
var(allGSR$roundedCounts) #technically, mean should equal variance
lambda.c <- fitdistr(allGSR$roundedCounts, densfun = "Poisson")
xvals <- 0:max(allGSR$roundedCounts)
densvals <- dpois(xvals, lambda = lambda.c$estimate)   ## get the density values
hist(allGSR$roundedCounts, freq = FALSE, main = "GSR Counts Histogram", xlab = "Counts")
lines(xvals, densvals, col = "red") #not perfect, but it could be worse?

#fit the model
set.seed(78)
all_fm2 <- stepFlexmix(~1, data = allGSR, k = 1:4, nrep = 10,
                       model=list(FLXMRglm(Peaks ~ ., family = "gaussian"),
                                  FLXMRglm(Age ~ ., family = "gaussian"),
                                  FLXMRglm(roundedCounts ~ ., family = "poisson")))
plot(1:4, BIC(all_fm2), type = "b", pch = 19, xaxt = "n", xlab = "Number of Clusters", ylab = "BIC", main = "BIC Cluster Evaluation")
axis(1, at = 1:4, labels = 1:4) ##best model is 3 clusters
fm_all2 <- getModel(all_fm2, "BIC")  
fm_all2
plot(Peaks ~ Age, data = allGSR, col = fm_all2@cluster)
plot3d(allGSR[c(-1,-3)], col = fm_all2@cluster, pch = 19)

### Let's introduce the repeated measure, first without the counts data
head(aggreg_gsr)
set.seed(45346)
all_fm3 <- stepFlexmix(~Pred|Subject,data = aggreg_gsr, k = 1:5, nrep = 5,
                       model=list(FLXMRglm(Peaks ~ 1, family = "gaussian"),
                                  FLXMRglm(Age ~ 1, family = "gaussian")))
plot(1:5, BIC(all_fm3), type = "b", pch = 19, xaxt = "n", xlab = "Number of Clusters", ylab = "BIC", main = "BIC Cluster Evaluation")
axis(1, at = 1:5, labels = 1:5)
fm_all3 <- getModel(all_fm3, "BIC")  
fm_all3 #5 cluster solution
plot(Peaks ~ Age, data = aggreg_gsr, col = fm_all3@cluster)

#Building a full mixed effects mixture model, first without incorporating the influence of Age
aggreg_gsr$roundCounts <- round(aggreg_gsr$Counts)
Model_pks <- FLXMRglm(Peaks ~ ., family = "gaussian")
Model_cnt <- FLXMRglm(roundCounts ~ ., family = "poisson")
set.seed(45346)
mixreg1_step <- stepFlexmix(~ Pred|Subject, data = aggreg_gsr, k = 1:5, nrep=5, 
                            model = list(Model_pks,Model_cnt))
plot(1:5, BIC(mixreg1_step), type = "b", pch = 19, xaxt = "n", xlab = "Number of Clusters", ylab = "BIC", main = "BIC Cluster Evaluation")
getModel(mixreg1_step, "BIC")  #once again, 5 cluster solution
mixreg1 <- flexmix(~ Pred|Subject, data = aggreg_gsr, k = 5, model = list(Model_pks,Model_cnt))
mixreg1
summary(mixreg1)
parameters(mixreg1)  
mixreg2 <- modeltools::refit(mixreg1)           
summary(mixreg2, model = 1)
summary(mixreg2, model = 2)
plot(mixreg2, model = 1, main = "Gaussian Regression")
plot(mixreg2, model = 2, main = "Poisson Regression")

#Finally, we add in the influence of age to build the final model
mixreg3 <- flexmix(~ Age+Pred|Subject, data = aggreg_gsr, k = 5, model = list(Model_pks,Model_cnt))
mixreg3
summary(mixreg3)
parameters(mixreg3)
mixreg4 <- modeltools::refit(mixreg3)           
summary(mixreg4, model = 1)
summary(mixreg4, model = 2)
plot(mixreg4, model = 1, main = "Gaussian Regression")
plot(mixreg4, model = 2, main = "Poisson Regression")

op <- par(mfrow = c(1,2))
with(aggreg_gsr, plot(Age, Peaks, xlab = "Age", ylab = "Peaks", main = "Scatterplot Age and Peak (Gaussian)", col = mixreg3@cluster))
with(aggreg_gsr, plot(Age, roundCounts, xlab = "Age", ylab = "Counts", main = "Scatterplot Age and Counts (Poisson)", col = mixreg3@cluster))
par(op)
#these clusters don't look particularly meaningful to me

