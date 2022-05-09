

# NO COLOR INTERACTION
# NO BIRD INTERACTION
#prior
                          
#===============================================
mod <- 
"
model{
for (k in 1:Ntile) {
for (j in 1:Nspp)   {
for (i in 1:Nrep)  {
#==================================================
# LEVEL 1 # DETECTIONS FOLLOW A BERNOULLI
y[i, j, k] ~ dbern(p[i, j, k] * z[j, k])
# CREATE PREDICTED Y
ypred[i, j, k] ~ dbern(p[i, j, k] * z[j, k])
# LOGLIK

logit(p[i, j, k]) <- a0[j, k] + 
#-------  VISIT  ATTRIBUTES
a1[j, k] * duration[i, k]  +
a2[j, k] * numevent[i, k] +
#----- BIRD
a3[j, k] * size[j, 1] +
a4[j, k] * m.red[j,1] +
a5[j,k] * v.blue[j, 1] 
} # end nReps
#=====================================================
# LEVEL 2 # PRESENCE FOLLOW A BERNOULLI
# LINK LEVEL 1 TO 2 THRU PRESENCE, z

z[j, k] ~ dbern(w[j, k])

# MODEL PRESENCE, z
# model varies intercepts and slopes across species and site
logit(w[j, k]) <- b0[j, k] + 
# SITE ATTRIBUTES
b1[j,k] * elev[k,1] +
b2[j,k] * grass[k,1] +
# BIRD
b3[j,k] * d.vend[j, 1] +
b4[j,k] * d.nect[j, 1]  +
b5[j,k] * d.nect[j, 1] * size[j, 1] +
b6[j,k] * d.vend[j, 1] * d.inv[j,1]
 
             
# LEVEL 3: COEFFICIENTS FOLLOW A GAUSSIAN
# LINK LEVEL 3 TO 2 THRU INTERCEPTS AND COEFFICIENTS OF MODEL OF z         
a0[j,k] ~ dnorm(mu, tau); a1[j,k] ~ dnorm(mu, tau); a2[j,k] ~ dnorm(mu, tau)
a3[j,k] ~ dnorm(mu, tau)
a4[j,k] ~ dnorm(mu, tau)
a5[j,k] ~ dnorm(mu, tau)
a6[j,k] ~ dnorm(mu, tau)
#-----------------
b0[j,k] ~ dnorm(mu, tau); b1[j,k] ~ dnorm(mu, tau); b2[j,k] ~ dnorm(mu, tau);
b3[j,k] ~ dnorm(mu, tau); b4[j,k] ~ dnorm(mu, tau)
b5[j,k] ~ dnorm(mu, tau); b6[j,k] ~ dnorm(mu, tau); b7[j,k] ~ dnorm(mu, tau); 
} # end nSpp
} # end nSites

#===========================================================
# HYPERPRIORS
mu ~ dnorm(0, 10)
sd ~ dunif(0, 2)
tau <- 1/(sd * sd) #pow(sd, -2) 
}
"
