## number of subjects ;
ns = 1000
## number of observations;
nis = 50

## population values
mu    = 3		#gamma00
beta1 = .3		#gamma10 or beta1
beta2 = .5		#gamma01
vare0 = -.5	#tau0
v1    = -.5		#alpha0
#v2    =  1.2		#?
#c12   =  -.35	#?
al   = 1				#alpha1		
ta1  = 1			#tau1
ta2  = .5			#tau2
set.seed  = 974657747

y <- array(NA, dim = c(ns*nis))			#empty y vector of length 1000*50
X1 <- array(NA, dim = c(ns*nis))		#xij
X2 <- array(NA, dim = c(ns*nis))		#wi
subj <- array(NA, dim = c(ns*nis))	#id
y.ni<- array(NA, dim = c(nis))			#empty ? vector of length 1000
x1.ni<- array(NA, dim = c(nis))			#?
sub.ni<- array(NA, dim = c(nis))		#?

for(n in 1:ns){										#for n in 1:1000 subjects
    sub1 = rnorm(1)								#sub1 variable on the zetaoi
    #sub2 = rnorm(1)							#?
    x2 = rnorm(1)    								#wi variable on gamma01 and alpha1
    s1  = sqrt(exp(v1 + x2*al))				#sigma_zeta0i
    #s12 = c12 / s1								#?
    #s2  = sqrt(v2 - (c12*c12 / s1*s1))	#?
    for(ni in 1:nis){								#for ni in 1:50 observations
        err =  rnorm(1)							#err variable on eij
        x1 = rnorm(1) 								#xij variable on beta1 and tau1
        g   = sqrt(exp(vare0 + x1*ta1 + x2*ta2))								#sigma_eij
        y.ni[ni]  =  mu + x1*beta1 + x2*beta2 + s1*sub1 + g*err		#response
        x1.ni[ni] <-  x1
        sub.ni[ni] <- n
    }
    y[((n-1)*nis+1):(n*nis)] <- y.ni
    X1[((n-1)*nis+1):(n*nis)] <- x1.ni
    X2[((n-1)*nis+1):(n*nis)] <- rep(x2, nis)
    subj[((n-1)*nis+1):(n*nis)] <- sub.ni
}

time = rep(c(1:50),1000)
scenario <- data.frame(y, X1, X2, subj, time)
head(scenario)

write.table(x=scenario, file='/Users/madelinecraft/Desktop/LSMsim.csv', sep = ",", row.names = FALSE, col.names = TRUE)

library(lattice)
xyplot(y ~ X1, groups = subj, type = 'p')

xyplot(y ~ X1|subj, groups = subj, type = 'l')
