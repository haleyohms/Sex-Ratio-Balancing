############################
# Estimating the proportion of migrants and residents in partially migratory species using sex ratio balancing
# H.A. Ohms, C.E. Jordan, A. Gitelman, D.A. Lytle
############################

############################
# The basic equation
# (equation 2 in the paper)
############################

# Assign migrant sex ratio
pi_g <- 0.8

# Assign resident sex ratio
pi_r <- 0.3

# Calculate proportion of migrants
rho_g <- ((0.5-pi_r)/(pi_g-pi_r))
rho_g


############################
# Computing the Bayesian HPD Interval
# (equation 6 in the paper)
############################

# Using the bat example:
# Migrants: 64:44  female:male 
# Residents: 1:126  female:male

# Posterior dist'n for migrants
likxprior_MigBats = function(p) dbeta(p,64,44)*2  #(p, females, males)
nc_MigBats = integrate(likxprior_MigBats,0.5,1)$value
pfun_MigBats = function(p) likxprior_MigBats(p)/nc_MigBats


# Posterior dist'n for residents 
likxprior_ResBats = function(p) dbeta(p,1,126)*2
nc_ResBats = integrate(likxprior_ResBats,0,0.5)$value
pfun_ResBats = function(p) likxprior_ResBats(p)/nc_ResBats


## Draw posterior samples

ss_MigBats <- seq(0.5,1,length=10000000)  # break the space into a million bits for numeric integration
pis_MigBats <- pfun_MigBats(ss_MigBats)

post_MigBats <- sample(ss_MigBats,10000,replace=T,prob=pis_MigBats)
hist(post_MigBats,prob=T,col="orange")


ss_ResBats <- seq(0,0.5,length=10000000)
pis_ResBats <- pfun_ResBats(ss_ResBats)

post_ResBats <- sample(ss_ResBats,10000,replace=T,prob=pis_ResBats)
hist(post_ResBats,prob=T,col="orange")


## Posterior proportion of bat migrants

out_Bats <- (0.5 - post_ResBats)/(post_MigBats-post_ResBats)
hist(out_Bats,prob=T,col="orange",nclass=25)

install.packages("TeachingDemos")
library(TeachingDemos)

emp.hpd(out_Bats)



###########
# SHEARWATERS
###########
# Migrants: 77:80  female:male 
# Residents: 1:14  female:male

## Migrants
likxprior_MigShear = function(p) dbeta(p,77,80)*2  #(p, females, males)
nc_MigShear = integrate(likxprior_MigShear,0.5,1)$value
pfun_MigShear = function(p) likxprior_MigShear(p)/nc_MigShear


## Residents 
likxprior_ResShear = function(p) dbeta(p,1,14)*2
nc_ResShear = integrate(likxprior_ResShear,0,0.5)$value
pfun_ResShear = function(p) likxprior_ResShear(p)/nc_ResShear


## Draw posterior samples

ss_MigShear <- seq(0.5,1,length=10000000)  # break the space into a million bits for numeric integration
pis_MigShear <- pfun_MigShear(ss_MigShear)

post_MigShear <- sample(ss_MigShear,10000,replace=T,prob=pis_MigShear)
hist(post_MigShear,prob=T,col="orange")


ss_ResShear <- seq(0,0.5,length=10000000)
pis_ResShear <- pfun_ResShear(ss_ResShear)

post_ResShear <- sample(ss_ResShear,10000,replace=T,prob=pis_ResShear)
hist(post_ResShear,prob=T,col="orange")


## Posterior proportion of migrant shearwaters

out_Shear <- (0.5 - post_ResShear)/(post_MigShear-post_ResShear)
hist(out_Shear,prob=T,col="orange",nclass=25)

meanShear <- mean(out_Shear)
medianShear <- median(out_Shear)

install.packages("TeachingDemos")
library(TeachingDemos)

emp.hpd(out_Shear)



########################
## Plots
########################


par(mfrow=c(2,3))
hist(post_MigShear,prob=T, xlim=c(0,1), col="yellow")
hist(post_ResShear,prob=T,xlim=c(0,1), col="blue")
hist(out_Shear,prob=T,xlim=c(0,1), col="green")

hist(post_MigBats,prob=T, xlim=c(0,1), col="lightblue")
hist(post_ResBats,prob=T,xlim=c(0,1), col="red")
hist(out_Bats,prob=T,xlim=c(0,1), col="purple")

# With different x-axis

par(mfrow=c(2,3), mgp=c(0,1.8,0))
hist(post_MigBats,prob=T, col="grey10", xlim=c(0.5, 0.8), cex.axis=1.5, yaxt='n', cex.main=1.5, main="Bat Migrant Sex Ratio", ylab="", xlab="")
hist(post_ResBats,prob=T, col="grey30", xlim=c(0, 0.12), yaxt='n', cex.axis=1.5, cex.main=1.5, main="Bat Resident Sex Ratio", ylab="", xlab="")
hist(out_Bats,prob=T, col="grey60", xlim=c(0.6, 1.0), yaxt='n', cex.axis=1.5, cex.main=1.5, main="Bat Proportion Migrant", ylab="", xlab="")
abline(h=0.8, col="black", lwd=3)

hist(post_MigShear,prob=T, col="grey10", xlim=c(0.5, 0.7), cex.axis=1.5, cex.main=1.5, yaxt='n', main="Shearwater Migrant Sex Ratio", ylab="", xlab="")
hist(post_ResShear,prob=T, col="grey30", xlim=c(0, 0.6), yaxt='n', cex.axis=1.5, cex.main=1.5, main="Shearwater Resident Sex Ratio", ylab="", xlab="")
hist(out_Shear,prob=T, col="grey60", breaks=20, xlim=c(0.5, 1), yaxt='n', cex.axis=1.5, cex.main=1.5, main="Shearwater Proportion Migrant", ylab="", xlab="")
abline(h=1.5, col="black", lwd=3)














# To make y-axis sum to 1
#h_MigShear <- hist(post_MigShear, breaks = 100, plot=FALSE)
#h_MigShear$counts=h_MigShear$counts/sum(h_MigShear$counts)
#plot(h_MigShear)


?hist()

hist(post_fs, col="orange", xlim=c(0,1), main="Posterior Pr(F|G)", xlab="")
hist(post1,prob=T,col="red", xlim=c(0,1), main="Bounded Priors: Post Pr(F|G)")

hist(post_fr, col="lightgreen", xlim=c(0,1), main="Posterior Pr(F|R)", xlab="")
hist(post2,prob=T,col="purple", xlim=c(0,1), main="Bounded Priors: Post Pr(F|R)" )

hist(ps, col="lightblue", xlim=c(-1, 1), breaks=1000, main="Posterior Pr(G)", xlab="")
hist(out,prob=T,col="yellow",nclass=25, xlim=c(0,1), main="Posterior Pr(G)")


h_MigShear <- hist(post_MigShear, breaks = 100, plot=FALSE)
h_MigShear$counts=h_MigShear$counts/sum(h_MigShear$counts)
plot(h_MigShear)







