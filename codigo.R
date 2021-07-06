
pacman::p_load(ggplot2,CCA,dplyr)
select <- dplyr::select

p <- ggplot(data = d0)
p1 <- p + geom_histogram(aes(CHLA))

d1 <- d0 %>% mutate(logCHLA = log(CHLA))

p <- ggplot(data = d1)
p2 <- p + geom_histogram(aes(logCHLA))
windows()
p2

##################

######
Y <- select(d1, DO, logCHLA)

X <- select(d1, PC, NO23F, NH4F, TN, PO4F, PP, TP, WTEMP, PH)

#######
cc1 <- cc(X, Y)
cc2 <- comput(X, Y, cc1)
################
## Dimension contributions to variability
cc1$cor

###########
ev <- (1 - cc1$cor^2)

n <- dim(Y)[1]
p <- length(Y)
q <- length(X)
k <- min(p, q)
m <- n - 3/2 - (p + q)/2

w <- rev(cumprod(rev(ev)))

# initialize
d1 <- d2 <- f <- vector("numeric", k)

for (i in 1:k) {
    s <- sqrt((p^2 * q^2 - 4)/(p^2 + q^2 - 5))
    si <- 1/s
    d1[i] <- p * q
    d2[i] <- m * s - p * q/2 + 1
    r <- (1 - w[i]^si)/w[i]^si
    f[i] <- r * d2[i]/d1[i]
    p <- p - 1
    q <- q - 1
}

pv <- pf(f, d1, d2, lower.tail = FALSE)
(dmat <- cbind(WilksL = w, F = f, df1 = d1, df2 = d2, p = pv))

######################
# STANDARDIZED CANONICAL COEFFICIENTS

s1 <- diag(sqrt(diag(cov(Y))))
ycc <- s1 %*% cc1$ycoef
yname <- data.frame(names(Y))
ycc <- cbind(yname, ycc)
names(ycc) <- c('var', 'c1', 'c2')#, 'c3')

##
s2 <- diag(sqrt(diag(cov(X))))
xcc <- s2 %*% cc1$xcoef
xname <- data.frame(names(X))
xcc <- cbind(xname, xcc)
names(xcc) <- c('var', 'c1', 'c2')#, 'c3')

########################################
## STRUCTURE COEFFICIENTS
ysc <- cc2$corr.Y.yscores
ysc <- cbind(yname, ysc)
colnames(ysc) <- c('var','s1', 's2')#, 's3')

xsc <- cc2$corr.X.xscores
xsc <- cbind(xname, xsc)
colnames(xsc) <- c('var','s1', 's2')#, 's3')

##############
X2 <- left_join(xcc, xsc, by = 'var') %>%
arrange(desc(abs(s1)))

Y2 <- left_join(ycc, ysc, by = 'var') %>%
arrange(desc(abs(s1)))

Y2
X2


