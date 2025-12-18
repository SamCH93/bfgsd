## ----"main-setup", include = FALSE--------------------------------------------
## knitr options
library(knitr)
opts_chunk$set(fig.height = 4,
               echo = FALSE,
               warning = FALSE,
               message = FALSE,
               cache = FALSE,
               eval = TRUE,
               fig.align = "center")

## printed digits
options(scipen = 100000)

## should sessionInfo be printed at the end?
Reproducibility <- TRUE

## packages
## remotes::install_github(repo = "SamCH93/bfpwr", subdir = "package", ref = "gsd")
library(bfpwr)
library(xtable)
library(mvtnorm)
library(ggplot2)
library(ggpubr)
library(rpact)
library(scales)
library(dplyr)
library(tidyr)
library(ggrain)

## colors
col1 <- "#D55E00"
col2 <- "#0072B2"
col3 <- "#009E73"


## ----"cartoon-seqBF", fig.height = 4.5----------------------------------------
## drawing a carton of a sequential BF trajectory
set.seed(14)
n <- seq(1, 65, 1)
usd <- sqrt(2)
pm <- 0
psd <- 1
dpm <- pm
dpsd <- psd
m <- rnorm(n = 1, mean = dpm, sd = dpsd)
y <- rnorm(n = n, mean = m, sd = usd)
z <- sapply(X = n, FUN = function(ni) {
    est <- mean(y[1:ni])
    se <- usd/sqrt(ni)
    est/se
})
bf <- sapply(X = n, FUN = function(ni) {
    est <- mean(y[1:ni])
    se <- usd/sqrt(ni) #sd(y[1:ni])/sqrt(ni)
    bfpwr::bf01(estimate = est, se = se, null = 0, pm = pm, psd = psd)
})
bks <- c(1/30, 1/10, 1/3, 1, 3, 10)
labs <- c("1/30", "1/10", "1/3", "1", "3", "10")
greycol <- adjustcolor(col = 1, alpha.f = 0.5)
arrowcol <- adjustcolor(col = 1, alpha.f = 0.5)

k0 <- 3
k1 <- 1/3

layout(rbind(1, 2), heights = c(3.8, 5))
par(mar = c(0.5, 4, 0.5, 0.5))
plot(n, bf, type = "n", log = "y", ylab = bquote("Bayes factor BF"["01"]),
     xlab = "",
     xaxt = "n",
     yaxt = "n", ylim = c(1/30, 10))
abline(h = c(k0, k1), lty = c(2, 3), lwd = c(1, 1.5), col = greycol)
lines(n, bf, lwd = 1.5, type = "l")
axis(side = 2, at = bks, labels = labs, las = 1)
text(x = 4, y = 5, labels = bquote(italic(H[0])), col = greycol)
text(x = 4, y = 1/5, labels = bquote(italic(H[1])), col = greycol)

par(mar = c(4, 4, 0.5, 0.5))
zcrit0 <- sapply(X = n, FUN = function(ni) {
    est <- mean(y[1:ni])
    se <- usd/sqrt(ni)
    bfpwr:::zcrit(k = k0, se = se, mu = pm, tau = psd, type = "normal")
})
zcrit1 <- sapply(X = n, FUN = function(ni) {
    est <- mean(y[1:ni])
    se <- usd/sqrt(ni)
    bfpwr:::zcrit(k = k1, se = se, mu = pm, tau = psd, type = "normal")
})
plot(n, z, type = "n", ylab = bquote(italic(z) * "-statistic"),
     xlab = bquote("Sample size" ~ italic(n)), las = 1, ylim = c(-3, 3))
matlines(n, y = t(zcrit0), lty = 2, col = greycol)
matlines(n, y = t(zcrit1), lty = 3, col = greycol, lwd = 1.5)
lines(n, z, lwd = 1.5, type = "l")
text(x = 5, y = 2.8, labels = bquote(italic(H[1])), col = greycol)
text(x = 5, y = -2.8, labels = bquote(italic(H[1])), col = greycol)
text(x = 20, y = 0, labels = bquote(italic(H[0])), col = greycol)


## ----"critical-values-example"------------------------------------------------
## point null vs. directional alternative BF
bf0p <- function(z, se, mu, tau) {
    taup <- 1/sqrt(1/se^2 + 1/tau^2)
    mup <- (z/se + mu/tau^2)*taup^2
    sqrt(1 + tau^2/se^2)*exp(-0.5*(z^2 - (z - mu/se)^2/(1 + tau^2/se^2)))*
        (pnorm(mu/tau)/(pnorm(mup/taup)))
}
## ## check that correct BF
## z <- 2
## se <- 0.1
## mu <- 0.1
## tau <- 0.5
## bf0p(z = z, se = se, mu = mu, tau = tau)
## dnorm(z*se, 0, se)/integrate(f = function(x) {
##         dnorm(z*se, x, se)*dnorm(x, mu, tau)/(1 - pnorm(0, mu, tau))
##     }, lower = 0, upper = Inf)$value


n <- seq(50, 250, 50)
## get critical z-values from group sequential designs
zP <- getDesignGroupSequential(kMax = length(n), informationRates = n/max(n),
                               alpha = 0.025, sided = 1,
                               typeOfDesign = "P")$criticalValues
zOF <- getDesignGroupSequential(kMax = length(n), informationRates = n/max(n),
                                alpha = 0.025, sided = 1,
                                typeOfDesign = "OF")$criticalValues

## get critical z-values from sequential BF designs
k1 <- 1/10
k0 <- 3
se <- sqrt(2/n)
pm0 <- 0
pm1 <- 0.1
psd <- 1
pm2 <- 0.5
resnormal <- pbf01seq(n = n, k1 = k1, k0 = k0, se = sqrt(2/n), pm = pm0,
                      psd = psd, type = "normal", strict = FALSE)
resdirectional <- pbf01seq(n = n, k1 = k1, k0 = k0, se = sqrt(2/n), pm = 0,
                           psd = psd, type = "directional", strict = FALSE)
respoint1 <- pbf01seq(n = n, k1 = k1, k0 = k0, se = sqrt(2/n), pm = pm1,
                      psd = 0, type = "normal", strict = FALSE)
respoint2 <- pbf01seq(n = n, k1 = k1, k0 = k0, se = sqrt(2/n), pm = pm2,
                      psd = 0, type = "normal", strict = FALSE)
resmoment <- pbf01seq(n = n, k1 = k1, k0 = k0, se = sqrt(2/n), psd = pm2/sqrt(2),
                      dpm = pm, type = "moment", strict = FALSE)
resJZS <- ptbf01seq(n = n, k1 = k1, k0 = k0, plocation = pm0, pscale = 1/sqrt(2),
                    pdf = 1, dpm = pm, type = "two.sample",
                    alternative = "two.sided", strict = FALSE)
res0p <- sapply(X = se, FUN = function(se) {
    uniroot(f = function(z) {
        bf0p(z = z, se = se, mu = pm0, tau = psd) - k1
    }, interval = c(0, 10), tol = .Machine$double.eps)$root
})
zDF <- rbind(
    ## data.frame(n = n, z = resnormal$zk1[2,],
    ##            type = paste0("'BF:' ~ italic(H[0]) * ':' ~ theta == 0 ~ 'vs.' ~ italic(H[1]) * ':' ~ theta == 0 ~ 'with' ~ italic(k[1]) == 1/", 1/k1)),
    ## data.frame(n = n, z = resmoment$zk1[2,], type = "'BF: Moment'"),
    data.frame(n = n, z = resdirectional$zk1,
               type = paste0("'BF:' ~ italic(H[0]) * ':' ~ theta <= 0 ~ 'vs.' ~ italic(H[1]) * ':' ~ theta > 0 ~ 'with' ~  theta %~%  'N(' *", pm0, "* ',' ~ ",
                              psd, "* ') and' ~ italic(k[1]) == 1/", 1/k1)),
    data.frame(n = n, z = respoint1$zk1,
               type = paste0("'BF:' ~ italic(H[0]) * ':' ~ theta == 0 ~ 'vs.' ~ italic(H[1]) * ':' ~ theta ==", pm1, "~ 'with' ~ italic(k[1]) == 1/", 1/k1)),
    data.frame(n = n, z = res0p,
               type = paste0("'BF:' ~ italic(H[0]) * ':' ~ theta == 0 ~ 'vs.' ~ italic(H[1]) * ':' ~ theta > 0 ~ 'with' ~  theta ~ '|' ~ italic(H[1]) %~%  'N(' *", pm0, "* ',' ~ ",
                              psd, "* ')'[(0 * ',' ~ +infinity)] ~ 'and' ~ italic(k[1]) == 1/", 1/k1)),
    ## data.frame(n = n, z = respoint2$zk1,
    ##            type = paste0("'BF:' ~ italic(H[0]) * ':' ~ theta == 0 ~ 'vs.' ~ italic(H[1]) * ':' ~ theta ==", pm2, "~ 'with' ~ italic(k[1]) == 1/", 1/k1)),
    ## data.frame(n = n, z = resJZS$zk1[2,], type = "'BF: JZS'"),
    data.frame(n = n, z = zP,
               type = "'GSD: Pocock with' ~ alpha == 0.025"),
    data.frame(n = n, z = zOF,
               type = '"GSD: O\'Brien-Fleming with" ~ alpha == 0.025')
)


## ----"critical-values-plot", fig.height = 3.2---------------------------------
ggplot(data = zDF, aes(x = n, y = z, color = type)) +
    geom_line(alpha = 0.4, linewidth = 1) +
    geom_point(size = 2, alpha = 0.9) +
    scale_color_manual(values = palette.colors(n = length(unique(zDF$type)) + 1)[-5],
                       labels = scales::label_parse()) +
    guides(color = guide_legend(nrow = 3)) +
    ## expand_limits(y = 0) +
    labs(x = "Sample size per group",
         y = bquote("Critical" ~ italic(z) * "-value"), color = "") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "top")


## ----"plot-critical-z-values", fig.height = 6, fig.width = 11-----------------
## function to compute critical z-values for Bayes factors
zcrit <- function(k, i, mu, tau, type = "point") {
    if (type == "point") {
        if (tau == 0) {
            zcrit <- (i*mu^2 - log(k^2))/(2*sqrt(i)*mu)
        } else {
            A <- (mu^2/tau^2 + log(1 + tau^2*i) - log(k^2))*(1 + 1/(tau^2*i))
            M <- -mu/(sqrt(i)*tau^2)
            zcrit <- M + c(-1, 1)*sqrt(A)
        }
    } else if (type == "directional") {
        priorodds <- 1/pnorm(mu/tau) - 1
        zcrit <- (qnorm(1/(k*priorodds + 1))*sqrt(i + 1/tau^2) -
                  mu/tau^2)/sqrt(i)
    } else if (type == "moment") {
        Y <- (2*lamW::lambertW0(x = exp(0.5)*(1 + tau^2*i)^(3/2)/(2*k)) - 1)*(1 + 1/(i*tau^2))
        zcrit <- c(-1, 1)*sqrt(Y)
    } else {
        stop("supplied type not implemented")
    }
    return(zcrit)
}

z1seq <- z2seq <- seq(-4, 4, length.out = 10)
mua <- 0
n <- c(50, 100)
inf <- n/2

## colors for analysis-wise success regions
alph <- 0.4
col1 <- adjustcolor(col = "#0072B2", alpha.f = alph)
col2 <- adjustcolor(col = "#009E73", alpha.f = alph)
col3 <- adjustcolor(col = "#D55E00", alpha.f = alph)
col4 <- adjustcolor(col = "#CC79A7", alpha.f = alph)

par(mfrow = c(1, 2))
## only one critical value: point null vs. point alternative / directional testing
k0 <- 10
k1 <- 1/10
zcrit01 <- zcrit(k = k0, i = inf[1], mu = mua, tau = 0.5, type = "directional")
zcrit02 <- zcrit(k = k0, i = inf[2], mu = mua, tau = 0.5, type = "directional")
zcrit11 <- zcrit(k = k1, i = inf[1], mu = mua, tau = 0.5, type = "directional")
zcrit12 <- zcrit(k = k1, i = inf[2], mu = mua, tau = 0.5, type = "directional")
plot(z1seq, z2seq, type = "n", las = 1, pty = "s",
     xlab = bquote(italic(z)[1] ~ "(" * italic(z) * "-statistic at analysis 1)"),
     ylab = bquote(italic(z)[2] ~ "(" * italic(z) * "-statistic at analysis 2)"),
     xaxt = "n", yaxt = "n",
     main = "One critical value")
     ## main = "One critical value:\n point/directional null vs. point/directional alternative")
## Pr(BF01^1 < k1)
xpH11 <- c(zcrit11, 100, 100, zcrit11)
ypH11 <- c(-100, -100, 100, 100)
text(x = zcrit11 + 0.5, y = 0,
     labels = bquote("Pr(BF"["01"]^"1" < italic(k)[1] * ")"), adj = 0)
## Pr(k0 > BF01^1 > k1, BF01^2 < k1)
xpH12 <- c(zcrit01, zcrit11, zcrit11, zcrit01)
ypH12 <- c(zcrit12, zcrit12, 100, 100)
text(x = zcrit01 + 0.1, y = zcrit12 + 1.5,
     labels = bquote(atop("Pr(" * italic(k)[1] ~ "< BF"["01"]^"1" ~ "<" ~ italic(k)[0] * ",",
                          "BF"["01"]^"2" < italic(k)[1] * ")")),
     adj = 0)
## Pr(BF01^1 > k0)
xpH01 <- c(-100, zcrit01, zcrit01, -100)
ypH01 <- c(-100, -100, 100, 100)
text(x = zcrit01 - 0.5, y = 0,
     labels = bquote("Pr(BF"["01"]^"1" > italic(k)[0] * ")"), adj = 1)
## Pr(k0 > BF01^1 > k1, BF01^2 < k1)
xpH02 <- c(zcrit01, zcrit11, zcrit11, zcrit01)
ypH02 <- c(zcrit02, zcrit02, -100, -100)
text(x = zcrit01 + 0.1, y = zcrit02 - 1.5,
     labels = bquote(atop("Pr(" * italic(k)[1] ~ "< BF"["01"]^"1" ~ "<" ~ italic(k)[0] * ",",
                          "BF"["01"]^"2" > italic(k)[0] * ")")),
     adj = 0)
polygon(xpH11, ypH11, col = col1, border = FALSE)
polygon(xpH12, ypH12, col = col2, border = FALSE)
polygon(xpH01, ypH01, col = col3, border = FALSE)
polygon(xpH02, ypH02, col = col4, border = FALSE)
abline(v = c(zcrit01, zcrit11), lty = c(3, 2), lwd = 1.5)
abline(h = c(zcrit02, zcrit12), lty = c(3, 2), lwd = 1.5)
axis(side = 1, at = c(zcrit01, zcrit11),
     labels = c(expression(italic(z)["1,crit"](italic(k)[0])),
                expression(italic(z)["1,crit"](italic(k)[1]))))
axis(side = 2, at = c(zcrit02, zcrit12),
     labels = c(expression(italic(z)["2,crit"](italic(k)[0])),
                expression(italic(z)["2,crit"](italic(k)[1]))))

## two critical values: point null vs. normal alternative
k0 <- 2
zcrit01 <- zcrit(k = k0, i = inf[1], mu = mua, tau = 0.5, type = "point")
zcrit02 <- zcrit(k = k0, i = inf[2], mu = mua, tau = 0.5, type = "point")
zcrit11 <- zcrit(k = k1, i = inf[1], mu = mua, tau = 0.5, type = "point")
zcrit12 <- zcrit(k = k1, i = inf[2], mu = mua, tau = 0.5, type = "point")
plot(z1seq, z2seq, type = "n", las = 1, pty = "s",
     xlab = bquote(italic(z)[1] ~ "(" * italic(z) * "-statistic at analysis 1)"),
     ylab = bquote(italic(z)[2] ~ "(" * italic(z) * "-statistic at analysis 2)"),
     xaxt = "n", yaxt = "n",
     main = "Two critical values")
     ## main = "Two critical values:\n point null vs. normal alternative")
abline(v = c(zcrit01, zcrit11), lty = c(3, 3, 2, 2), lwd = 1.5)
abline(h = c(zcrit02, zcrit12), lty = c(3, 3, 2, 2), lwd = 1.5)
## Pr(BF01^1 < k1)
xpH11a <- c(zcrit11[1], -100, -100, zcrit11[1])
ypH11a <- c(-100, -100, 100, 100)
xpH11b <- c(zcrit11[2], 100, 100, zcrit11[2])
ypH11b <- c(-100, -100, 100, 100)
polygon(xpH11a, ypH11, col = col1, border = FALSE)
polygon(xpH11b, ypH11, col = col1, border = FALSE)
## Pr(k0 > BF01^1 > k1, BF01^2 < k1)
xpH12a <- c(zcrit01[1], zcrit11[1], zcrit11[1], zcrit01[1])
ypH12a <- c(zcrit12[1], zcrit12[1], -100, -100)
xpH12b <- c(zcrit01[1], zcrit11[1], zcrit11[1], zcrit01[1])
ypH12b <- c(zcrit12[2], zcrit12[2], 100, 100)
xpH12c <- c(zcrit01[2], zcrit11[2], zcrit11[2], zcrit01[2])
ypH12c <- c(zcrit12[2], zcrit12[2], 100, 100)
xpH12d <- c(zcrit01[2], zcrit11[2], zcrit11[2], zcrit01[2])
ypH12d <- c(zcrit12[1], zcrit12[1], -100, -100)
polygon(xpH12a, ypH12a, col = col2, border = FALSE)
polygon(xpH12b, ypH12b, col = col2, border = FALSE)
polygon(xpH12c, ypH12c, col = col2, border = FALSE)
polygon(xpH12d, ypH12d, col = col2, border = FALSE)
## Pr(BF01^1 > k0)
xpH01 <- c(zcrit01[1], zcrit01[1], zcrit01[2], zcrit01[2])
ypH01 <- c(100, -100, -100, 100)
polygon(xpH01, ypH01, col = col3, border = FALSE)
## Pr(k0 > BF01^1 > k1, BF01^2 < k1)
xpH02a <- c(zcrit01[1], zcrit11[1], zcrit11[1], zcrit01[1])
ypH02a <- c(zcrit02[1], zcrit02[1], zcrit02[2], zcrit02[2])
xpH02b <- c(zcrit01[2], zcrit11[2], zcrit11[2], zcrit01[2])
ypH02b <- c(zcrit02[1], zcrit02[1], zcrit02[2], zcrit02[2])
polygon(xpH02a, ypH02a, col = col4, border = FALSE)
polygon(xpH02b, ypH02b, col = col4, border = FALSE)
axis(side = 1, at = c(zcrit01, zcrit11),
     labels = c(expression(italic(z)["1,crit-"](italic(k)[0])),
                expression(italic(z)["1,crit+"](italic(k)[0])),
                expression(italic(z)["1,crit-"](italic(k)[1])),
                expression(italic(z)["1,crit+"](italic(k)[1]))))
axis(side = 2, at = c(zcrit02, zcrit12),
     labels = c(expression(italic(z)["2,crit-"](italic(k)[0])),
                expression(italic(z)["2,crit+"](italic(k)[0])),
                expression(italic(z)["2,crit-"](italic(k)[1])),
                expression(italic(z)["2,crit+"](italic(k)[1]))))


## ----"low-PV-trial-data"------------------------------------------------------
## https://doi.org/10.1016/S2352-3026(20)30373-2
## https://ars.els-cdn.com/content/image/1-s2.0-S2352302620303732-mmc1.pdf

## H1: p0 = 0.5 and p1 = 0.75
p0H1 <- 0.5
p1H1 <- 0.75
ORH1 <- (p1H1/(1 - p1H1)) / (p0H1/(1 - p0H1))

a1 <- 21
b1 <- 24 - 21
c1 <- 15
d1 <- 26 - 15
logOR1 <- log(a1*d1/(b1*c1))
selogOR1 <- sqrt(1/a1 + 1/b1 + 1/c1 + 1/d1)
z1 <- logOR1/selogOR1
## ## reported OR estimate 5.1 corresponds to this
OR1 <- exp(logOR1)
OR1ci <- exp(logOR1 + c(-1, 1)*selogOR1*1.96) # doesn't correspond to 1.1 to 32.5
## ## CI corresponds to this
## fisher.test(matrix(c(a1, b1, c1, d1), byrow = TRUE, ncol = 2))
## Epi::twoby2(matrix(c(a1, b1, c1, d1), byrow = TRUE, ncol = 2))
pm <- log(ORH1)
psd <- 0
bf1 <- bfpwr::bf01(estimate = logOR1, se = selogOR1, null = 0, pm = pm,
                   psd = psd)
a2 <- 42
b2 <- 50 - 42
c2 <- 30
d2 <- 50 - 30
logOR2 <- log(a2*d2/(b2*c2))
selogOR2 <- sqrt(1/a2 + 1/b2 + 1/c2 + 1/d2)
z2 <- logOR2/selogOR2
## ## reported OR estimate 3.5 corresponds to this
OR2 <- exp(logOR2)
OR2ci <- exp(logOR2 + c(-1, 1)*selogOR2*1.96) # doesn't correspond to 1.3 to 10.4
## ## CI corresponds to this
## fisher.test(matrix(c(a2, b2, c2, d2), byrow = TRUE, ncol = 2))
## Epi::twoby2(matrix(c(a2, b2, c2, d2), byrow = TRUE, ncol = 2))
bf2 <- bfpwr::bf01(estimate = logOR2, se = selogOR2, null = 0,  pm = pm,
                   psd = psd)

## OF stopping boundaries
zcritFut <- c(-3.7307, -2.5262, -1.9917)
zcritEff <- c(3.7307, 2.5262, 1.9917)
## rpact::getDesignGroupSequential(informationRates = c(0.33, 0.66, 1),
##                                 alpha = 0.025, typeOfDesign = "asOF")


## ----"Low-PV-trial-table", results = "asis"-----------------------------------
## BF sequential design based on logOR
k1 <- 1/10
k0 <- 10
n <- c(25, 50, 75)
se <- sqrt(1/(p0H1*(1 - p0H1)*n) + 1/(p1H1*(1 - p1H1)*n))
bfdesignH1 <- bfpwr::pbf01seq(k1 = k1, k0 = k0, se = se, n = n, pm = pm,
                              psd = psd, dpm = pm, dpsd = 0, type = "normal")
bfdesignH0 <- bfpwr::pbf01seq(k1 = k1, k0 = k0, se = se, n = n, pm = pm,
                              psd = psd, dpm = 0, dpsd = 0, type = "normal")
## bfdesignH1$zk0
## bfdesignH1$zk1
## plot(bfdesignH1)
zk1 <- bfpwr:::zcrit(k = k1, se = c(selogOR1, selogOR2), mu = pm, tau = psd,
                     type = "normal")
zk0 <- bfpwr:::zcrit(k = k0, se = c(selogOR1, selogOR2), mu = pm, tau = psd,
                     type = "normal")

## BF sequential design based on asin-sqrt difference
k1 <- 1/10
k0 <- 10
n <- c(25, 50, 75)
se <- sqrt(1/(p0H1*(1 - p0H1)*n) + 1/(p1H1*(1 - p1H1)*n))
se0 <- sqrt(1/(p0H1*(1 - p0H1)*n) + 1/(p0H1*(1 - p0H1)*n))
bfdesignH1 <- bfpwr::pbf01seq(k1 = k1, k0 = k0, se = se, n = n, pm = pm,
                              psd = psd, dpm = pm, dpsd = 0, type = "normal")
bfdesignH0 <- bfpwr::pbf01seq(k1 = k1, k0 = k0, se = se0, n = n, pm = pm,
                              psd = psd, dpm = 0, dpsd = 0, type = "normal")

formatround <- function(x, digits) {
    format(round(x, digits = digits), nsmall = digits)
}
restable <- data.frame(
    analysis = c("1", "2"),
    n0 = as.character(c(c1 + d1, c2 + d2)),
    n1 = as.character(c(a1 + b1, a2 + b2)),
    p0 = paste0(formatround(100*c(c1/(c1 + d1), c2/(c2 + d2)), 1), "\\%"),
    p1 = paste0(formatround(100*c(a1/(a1 + b1), a2/(a2 + b2)), 1), "\\%"),
    OR = paste0(formatround(c(OR1, OR2), 1), " (",
                formatround(c(OR1ci[1], OR2ci[1]), 1),
                " to ",
                formatround(c(OR1ci[2], OR2ci[2]), 1),
                ")"),
    z = formatround(c(z1, z2), 2),
    Decision = c("Continue", "Stop (efficacy)"),
    BF = paste0("1/", round(1/c(bf1, bf2), 1))
)

restablex <- xtable::xtable(restable)
colnames(restablex) <- c("Analysis", "$n_0$", "$n_1$", "$\\hat{\\pi}_0$",
                         "$\\hat{\\pi}_1$",
                         "$\\widehat{\\mathrm{OR}}$ (95\\% CI)", "$z$",
                         "Decision",
                         "$\\mathrm{BF}_{01}$")
print(restablex, floating = FALSE, include.rownames = FALSE,
      booktabs = TRUE,
      sanitize.text.function = function(x){x})


## ----"low-PV-design-plot", fig.height = 3.5-----------------------------------
pH1lab <- paste0("'Stop for' ~ italic(H[1]) ~ '(BF'['01'] <= 1/", 1/k1, " * ')'")
pH0lab <- paste0("'Stop for' ~ italic(H[0]) ~ '(BF'['01'] >= ", k0, " * ')'")
pInclab <- "'Inconclusive'"
cols <- c(2, 1, 4)
names(cols) <- c(pH0lab, pInclab, pH1lab)
dfh1 <- data.frame(plot(bfdesignH1, nullplot = FALSE, plot = FALSE)$pDF1,
                   dp = "H1")
dfh0 <- data.frame(plot(bfdesignH0, nullplot = FALSE, plot = FALSE)$pDF1,
                   dp = "H0")

H0lab <- "'under' ~ italic(H[0]) * ':' ~ theta == 0"
H1lab <- paste0("'under' ~ italic(H[1]) * ':' ~ theta == ", round(pm, 2))
plotdf <- rbind(dfh1, dfh0) |>
    tidyr::pivot_longer(cols = c("pH0", "pH1", "pInc"),
                        names_to = "type",
                        values_to = "probability") |>
    dplyr::mutate(setting = factor(dp, levels = c("H1", "H0"),
                                   labels = c(H1lab, H0lab)),
                  type = dplyr::case_when(type == "pH1" ~ pH1lab,
                                          type == "pH0" ~ pH0lab,
                                          type == "pInc" ~ pInclab),
                  type = factor(type, levels = c(pH1lab, pInclab, pH0lab)))
ggplot(data = plotdf, aes(x = n, y = probability, color = type)) +
    facet_wrap(~ setting, labeller = label_parsed) +
    geom_line(alpha = 0.3) +
    geom_point() +
    scale_y_continuous(breaks = seq(0, 1, 0.2),
                       labels = scales::percent, limits = c(0, 1)) +
    scale_x_continuous(breaks = seq(25, 100, 25)) +
    scale_color_manual(labels = scales::label_parse(),
                       values = cols) +
    labs(x = "Sample size (per group)", y = "Probability (by-analysis)", color = "") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "top",
          strip.background = element_rect(fill = "white"))


## ----"low-PV-design-sample-size"----------------------------------------------
## power function to use for root-finding maximum sample size
powerfun. <- function(nmax, m, H = "H1") {
    n <- seq(nmax/m, nmax, nmax/m)
    if (H == "H1") {
        se <- sqrt(1/(p0H1*(1 - p0H1)*n) + 1/(p1H1*(1 - p1H1)*n))
        dpm <- pm
    } else {
        se <- sqrt(1/(p0H1*(1 - p0H1)*n) + 1/(p0H1*(1 - p0H1)*n))
        dpm <- 0
    }
    res <- pbf01seq(k1 = k1, k0 = k0, se = se, n = n, pm = pm, psd = psd,
                    dpm = dpm, dpsd = 0, type = "normal")
    if (H == "H1") {
        p <- res$cumpH1[m]
    } else {
        p <- res$cumpH0[m]
    }
    return(p)
}
powerfun <- Vectorize(FUN = powerfun.)

## fix number of analyes and find maximum sample size to achieve target power
m <- 3
pow <- 0.9
rootfun <- function(nmax, m, H) {
    powerfun(nmax = nmax, m = m, H = H) - pow
}
nmax1 <- uniroot(f = rootfun, interval = c(10, 1000), m = m,
                 H = "H1")$root
nmax0 <- uniroot(f = rootfun, interval = c(10, 1000), m = m,
                 H = "H0")$root


## ----"low-PV-design-plot2", fig.height = 7.5----------------------------------
lowPVdesign <- function(nmax = 75, m = 3) {
    n <- seq(nmax/m, nmax, length.out = m)
    se1 <- sqrt(1/(p0H1*(1 - p0H1)*n) + 1/(p1H1*(1 - p1H1)*n))
    se0 <- sqrt(1/(p0H1*(1 - p0H1)*n) + 1/(p0H1*(1 - p0H1)*n))
    bfdesignH1 <- bfpwr::pbf01seq(k1 = k1, k0 = k0, se = se1, n = n, pm = pm,
                                  psd = psd, dpm = pm, dpsd = 0,
                                  type = "normal")
    bfdesignH0 <- bfpwr::pbf01seq(k1 = k1, k0 = k0, se = se0, n = n, pm = pm,
                                  psd = psd, dpm = 0, dpsd = 0, type = "normal")
    rbind(
        data.frame(nmax = nmax, m = m,
                   dp = "'under' ~ italic(H[0]) * ':' ~ theta == 0",
                   EN = bfdesignH0$EN,
                   SDN = sqrt(bfdesignH0$VarN),
                   p = c(tail(bfdesignH0$cumpH0, n = 1),
                         tail(bfdesignH0$cumpH1, n = 1)),
                   type = c("'Pr(Correct evidence)'", "'Pr(Misleading evidence)'")),
        data.frame(nmax = nmax, m = m,
                   dp = paste0("'under' ~ italic(H[1]) * ':' ~ theta == ",
                               round(pm, 2)),
                   EN = bfdesignH1$EN,
                   SDN = sqrt(bfdesignH1$VarN),
                   p = c(tail(bfdesignH1$cumpH1, n = 1),
                         tail(bfdesignH1$cumpH0, n = 1)),
                   type = c("'Pr(Correct evidence)'", "'Pr(Misleading evidence)'"))
          )
}
grid <- expand.grid(nmax = seq(50, 150, 5),
                    m = seq(1, 10, 1))
plotDF <- do.call("rbind", lapply(X = seq(1, nrow(grid)), FUN = function(i) {
    lowPVdesign(nmax = grid$nmax[i], m = grid$m[i])
}))

probplot <- ggplot(plotDF,
                   aes(x = nmax, y = p,
                       color = factor(m, ordered = TRUE))) +
    facet_grid(type ~ dp, scales = "free_y",
               labeller = label_parsed, switch = "y") +
    geom_line(alpha = 0.9) +
    ## geom_point() +
    scale_y_continuous(labels = scales::percent) +
    ## expand_limits(y = 0) +
    guides(color = guide_legend(nrow = 1)) +
    labs(x = "Maximum sample size (per group)", y = "",
         color = "Number of analyses") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "top",
          axis.text.x = element_blank(),
          axis.title.x.bottom = element_blank(),
          axis.ticks.x.bottom = element_blank(),
          strip.background.y = element_blank(),
          strip.background.x = element_rect(fill = "white"),
          strip.placement = "outside",
          strip.text.y = element_text(size = 9))

Nplotdf1 <- subset(plotDF, type == "'Pr(Correct evidence)'") |>
    dplyr::mutate(yvar = EN, type = "'E(sample size)'")
Nplotdf2 <- subset(plotDF, type == "'Pr(Correct evidence)'") |>
    dplyr::mutate(yvar = SDN, type = "'SD(sample size)'")
    ## dplyr::mutate(yvar = SDN/EN, type = "'COV(sample size)'")
Nplotdf <- rbind(Nplotdf1, Nplotdf2)
Nplot <- ggplot(Nplotdf,
                aes(x = nmax, y = yvar, color = factor(m, ordered = TRUE))) +
    facet_grid(type ~ dp, scales = "free_y",
               labeller = label_parsed, switch = "y") +
    geom_line(alpha = 0.9, show.legend = FALSE) +
    labs(x = "Maximum sample size (per group)", y = "") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "top",
          strip.background = element_blank(),
          strip.text.x = element_text(colour = NA, margin = margin(t = 0, b = 0)),
          strip.placement = "outside",
          strip.text.y = element_text(size = 9))
ggpubr::ggarrange(plotlist = list(probplot, Nplot), ncol = 1, align = "v",
                  heights = c(1, 1))


## ----"preclinical-data", fig.height = 5---------------------------------------
## library(metaDigitise)
## ## click on "Control" as 1 and "Low" as 2
## kangdat <- metaDigitise(dir = "kang-plots/", summary = FALSE)
## kangdat <- res$`data-kangetal.png`
## kangdat$weight_loss <- kangdat$y
## kangdat$dose <- ifelse(kangdat$x < 1.5, "control",
##                 ifelse(kangdat$x < 2.5, "low",
##                 ifelse(kandat$x < 3.5, "medium", "high")))
## dat <- kangdat[,c("dose", "weight_loss")]
## write.csv(x = dat, file = "kangdat.csv", row.names = FALSE)

## plot data
kangdat <- read.csv("kangdat.csv")
doselvls <- c("control", "low", "medium", "high")
doselabs <- c("Control", "Low", "Medium", "High")
doselabs2 <- c(paste0("Control (n = ", sum(kangdat$dose == "control"), ")"),
               paste0("Low (n = ", sum(kangdat$dose == "low"), ")"),
               paste0("Medium (n = ", sum(kangdat$dose == "medium"), ")"),
               paste0("High (n = ", sum(kangdat$dose == "high"), ")"))
kangdat$dosefac <- factor(kangdat$dose, levels = doselvls, labels = doselabs)
kangdat$dosefac2 <- factor(kangdat$dose, levels = doselvls, labels = doselabs2)
kangplot1 <- ggplot(kangdat, aes(x = dosefac2, y = weight_loss, fill = dosefac2)) +
    geom_hline(yintercept = 0, lty = 2, alpha = 0.3) +
    ggrain::geom_rain(show.legend = FALSE, alpha = 0.5,
                      point.args = list(alpha = 1, pch = 1, size = 2,
                                        show.legend = FALSE)) +
    scale_fill_manual(values = palette.colors(n = length(doselvls))) +
    labs(x = "Dose level", y = "Weight loss (grams)") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())

analyzekangdat <- function(control, treatment, pm, psd) {
    nseq <- seq(2, pmax(length(control), length(treatment)), 1)
    do.call("rbind", lapply(X = nseq, FUN = function(i) {
        controli <- na.omit(control[1:i])
        treatmenti <- na.omit(treatment[1:i])
        esti <- mean(treatmenti) - mean(controli)
        sei <- sqrt(var(treatmenti)/length(treatmenti) +
                   var(controli)/length(controli))

        ## calculate z-test BF
        if (is.finite(sei)) {
            ## bf01 <- bfpwr::dirbf01(estimate = esti, se = sei, null = 0, pm = pm,
            ##                        psd = psd)
            bf01 <- bfpwr::bf01(estimate = esti, se = sei, null = 0, pm = pm,
                                psd = psd)
        } else {
            bf01 <- NaN
        }

        ## return results
        data.frame(ncontrol = length(controli), ntreatment = length(treatmenti),
                   MD = esti, MDse = sei, bf01 = bf01)
    }))
}

## mean and standard deviation of normal prior assigned to mean difference
pm <- 5
psd <- 0

## evidence thresholds
k1 <- 1/10
k0 <- 10

set.seed(11)
control <- sample(x = subset(kangdat, dose == "control")$weight_loss,
                  replace = FALSE)
treat <- list("Low" = sample(x = subset(kangdat, dose == "low")$weight_loss,
                             replace = FALSE),
              "Medium" = sample(x = subset(kangdat, dose == "medium")$weight_loss,
                                replace = FALSE),
              "High" = sample(x = subset(kangdat, dose == "high")$weight_loss,
                              replace = FALSE))
kanglist <- lapply(X = seq_along(treat), FUN = function(i) {
    res <- analyzekangdat(control = control, treatment = treat[[i]], pm = pm,
                          psd = psd)
    decision <- "indecisive"
    for (ni in seq(1, nrow(res))) {
        ntreat <- res$ntreatment[ni]
        ncontrol <- res$ncontrol[ni]
        if (res$bf01[ni] > k0) {
            decision <- "H0"
            break
        }
        if (res$bf01[ni] < k1) {
            decision <- "H1"
            break
        }
    }
    ## how many rats are saved?
    ntreatsaved <- res$ntreatment[nrow(res)] - ntreat
    ncontrolsaved <- res$ncontrol[nrow(res)] - ncontrol
    list(df = data.frame(group = names(treat)[i], res),
         decision = data.frame(group = names(treat)[i], decision = decision,
                               ntreat = ntreat, ncontrol = ncontrol,
                               ntreatsaved = ntreatsaved,
                               ncontrolsaved = ncontrolsaved))
})
kangres <- do.call("rbind", lapply(X = kanglist, FUN = function(x) x$df))
kangres$group <- factor(kangres$group, levels = c("Low", "Medium", "High"),
                        ordered = TRUE)
bfbks <- c(1/10^6, 1/10^5, 1/10^4, 1/10^3, 1/10^2, 1/10, 1, 10, 100)
bflabs <- c("1/10^6", "1/10^5", "1/10^4", "1/10^3", "1/10^2", "1/10", "1", "10",
            "10^2")
kangplot2 <- ggplot(data = kangres,
                    aes(x = ncontrol + ntreatment, y = bf01, color = group)) +
    geom_hline(yintercept = 1, alpha = 0.5, lty = 2) +
    geom_line() +
    geom_point() +
    scale_y_log10(breaks = bfbks, labels = parse(text = bflabs)) +
    expand_limits(y = 100) +
    scale_x_continuous(breaks = seq(4, max(kangres$ncontrol + kangres$ntreatment), 2)) +
    labs(x = "Total sample size (control and treatment)",
         y = bquote("BF"["01"])) +
    guides(color = "none") +
    scale_color_manual(values = palette.colors(n = length(doselvls))[-1]) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
ggpubr::ggarrange(plotlist = list(kangplot1, kangplot2), ncol = 1, align = "hv")

## how many rats are saved? (not double-counting control rats as they contribute to all comparisons)
ratssaved <- do.call("rbind", lapply(X = kanglist, FUN = function(x) x$decision)) |>
    summarise(treatsaved = sum(ntreatsaved),
              controlsaved = min(ncontrolsaved),
              totalsaved = sum(ntreatsaved) + min(ncontrolsaved))


## ----"kang-design", fig.height = 3.5------------------------------------------
k1 <- 1/10
k0 <- 10
designkangdat <- function(control, treatment, pm, psd, setting) {
    nseq <- seq(2, 30, 2)
    MD <- mean(treatment) - mean(control)
    MDse <- sqrt(var(treatment)/(length(treatment)/2) +
                 var(control)/(length(control)/2))
    seseq <- sqrt(var(treatment)/(nseq/2) + var(control)/(nseq/2))
    res <- pbf01seq(k1 = k1, k0 = k0, se = seseq, n = nseq, pm = pm, psd = psd,
                    dpm = MD, dpsd = MDse, type = "normal")
    plotdf <- plot(res, plot = FALSE)
    rbind(data.frame(plotdf$pDF1, H = "Effect from data", setting = setting,
                     MD = MD, MDse = MDse),
          data.frame(plotdf$pDF2, H = "No effect", setting = setting, MD = 0,
                     MDse = 0))
}
treat2 <- treat
kangdesignres <- do.call("rbind", lapply(X = seq_along(treat2), FUN = function(i) {
    designkangdat(control = control, treatment = treat2[[i]], pm = pm, psd = psd,
                  setting = names(treat2)[i])
}))

pH1lab <- paste0("'Stop for' ~ italic(H[1]) ~ '(BF'['01'] <= 1/", 1/k1, " * ')'")
pH0lab <- paste0("'Stop for' ~ italic(H[0]) ~ '(BF'['01'] >= ", k0, " * ')'")
pInclab <- "'Inconclusive'"
cols <- c(2, 1, 4)
names(cols) <- c(pH0lab, pInclab, pH1lab)
kangdesignreslong <- kangdesignres |>
    tidyr::pivot_longer(cols = c("pH0", "pH1", "pInc"),
                        names_to = "type",
                        values_to = "probability") |>
    dplyr::mutate(setting = factor(setting, levels = c("Low", "Medium", "High")),
                  type = dplyr::case_when(type == "pH1" ~ pH1lab,
                                          type == "pH0" ~ pH0lab,
                                          type == "pInc" ~ pInclab),
                  type = factor(type, levels = c(pH1lab, pInclab, pH0lab)))

ggplot(subset(kangdesignreslong, MD != 0),
       aes(x = n, y = probability, color = type)) +
    ## facet_grid(H ~ setting) +
    facet_wrap(~ setting + MD + MDse,
               labeller = label_bquote(.(as.character(setting)) ~
                                           "dose:" ~
                                           theta %~% "N(" * .(round(MD, 1)) *
                                            "," ~ .(round(MDse, 1))^2 *
                                                ")")) +
    geom_line(alpha = 0.3) +
    geom_point() +
    scale_y_continuous(breaks = seq(0, 1, 0.2),
                       labels = scales::percent, limits = c(0, 1)) +
    scale_color_manual(values = cols, labels = scales::label_parse()) +
    labs(x = "Total sample size (control and treatment)",
         y = "Probability (by-analysis)", color = "") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "top",
          strip.background = element_rect(fill = "white"))


## ----"schoenbrodt-wagenmakers-example"----------------------------------------
## sample size and number of looks
nmin <- 40
nmax <- 100
step <- 1
n <- seq(nmin, nmax, step)
type <- "two.sample"

## one-sided JZS analysis prior
plocation <- 0
pscale <- 1/sqrt(2)
pdf <- 1
alternative <- "greater"

## design prior
dpm <- 0.5
dpsd <- 0.1

## BF thresholds
k0 <- 6
k1 <- 1/30

startbfpwr <- Sys.time()
res1 <- ptbf01seq(k1 = k1, k0 = k0, n = n, plocation = plocation,
                  pscale = pscale, pdf = pdf, dpm = dpm, dpsd = dpsd,
                  type = type, alternative = alternative)
res0 <- ptbf01seq(k1 = k1, k0 = k0, n = n, plocation = plocation,
                  pscale = pscale, pdf = pdf, dpm = 0, dpsd = 0, type = type,
                  alternative = alternative)
endbfpwr <- Sys.time()
timebfpwr <- difftime(endbfpwr, startbfpwr, units = "secs")

## startBFDA <- Sys.time()
## ## remotes::install_github("nicebread/BFDA", subdir = "package")
## library(BFDA)
## nMC <- 1000
## BFDAsim <- BFDA.sim(expected.ES = rnorm(n = nMC, mean = dpm, sd = dpsd),
##                     type = "t.between",
##                     alternative = alternative,
##                     prior = list("Cauchy",
##                                  list(prior.location = plocation,
##                                       prior.scale = pscale)),
##                     design = "sequential", stepsize = step, n.min = nmin,
##                     n.max = nmax, B = nMC, seed = 42, cores = 10)
## plot(BFDAsim, boundary = c(1/k1, 1/k0), n.trajectories = 60)
## BFDA.analyze(BFDAsim, design = "sequential", n.min = nmin, n.max = nmax,
##              boundary = c(1/k1, 1/k0))
## BFDAsim0 <- BFDA.sim(expected.ES = 0, type = "t.between",
##                      alternative = alternative,
##                      prior = list("Cauchy",
##                                  list(prior.location = plocation,
##                                       prior.scale = pscale)),
##                      design = "sequential", stepsize = step, n.min = nmin,
##                      n.max = nmax, B = nMC, seed = 43, cores = 10)
## plot(BFDAsim0, boundary = c(1/k1, 1/k0), n.trajectories = 60)
## BFDA.analyze(BFDAsim0, design = "sequential", n.min = nmin, n.max = nmax,
##              boundary = c(1/k1, 1/k0))
## endBFDA <- Sys.time()
## timeBFDA <- difftime(endBFDA, startBFDA, units = "mins")


## ----"schoenbrodt-wagenmakers-example-figure"---------------------------------
plotres <- plot(res1, plot = FALSE)
H0lab <- "'under' ~ theta == 0"
H1lab <- paste0("'under' ~ theta %~% 'N(", round(dpm, 2), ", ", round(dpsd, 2), "'^2 * ')'")
pH1lab <- paste0("'Stop for' ~ italic(H[1]) ~ '(BF'['01'] <= 1/", 1/k1, " * ')'")
pH0lab <- paste0("'Stop for' ~ italic(H[0]) ~ '(BF'['01'] >= ", k0, " * ')'")
pInclab <- "'Inconclusive'"
cols <- c(2, 1, 4)
names(cols) <- c(pH0lab, pInclab, pH1lab)
plotdf <- rbind(data.frame(plotres$pDF1, dp = H1lab),
                data.frame(plotres$pDF2, dp = H0lab)) |>
    tidyr::pivot_longer(cols = c("pH0", "pH1", "pInc"),
                        names_to = "type",
                        values_to = "probability") |>
    dplyr::mutate(type = dplyr::case_when(type == "pH1" ~ pH1lab,
                                          type == "pH0" ~ pH0lab,
                                          type == "pInc" ~ pInclab),
                  type = factor(type, levels = c(pH1lab, pInclab, pH0lab)))

ggplot(data = plotdf, aes(x = n, y = probability, color = type)) +
    facet_wrap(~ dp, ncol = 2,
               labeller = label_parsed) +
    geom_line(alpha = 0.3) +
    geom_point(size = 0.5) +
    labs(x = "Sample size (per group)", y = "Probability (by-analysis)", color = "") +
    scale_y_continuous(breaks = seq(0, 1, 0.2),
                       labels = scales::percent, limits = c(0, 1)) +
    scale_color_manual(values = cols, labels = scales::label_parse()) +
    guides(color = guide_legend(override.aes = list(size = 1))) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "top",
          strip.background = element_rect(fill = "white"))


## ----"appendix-package", fig.height = 6, echo = TRUE, size = "small"----------
## group sequential design features not in CRAN version yet
remotes::install_github(repo = "SamCH93/bfpwr", subdir = "package", ref = "gsd")
library(bfpwr) # load package
## set up sequential t-test Bayes factor design
design <- ptbf01seq(
    k1 = 1/10, # Bayes factor threshold for H1
    k0 = 6, # Bayes factor threshold for H0
    type = "two.sample", # two-sample t-test
    n = seq(20, 100, 20), # per-group sample sizes at analyses
    ## specify one-sided Jeffreys-Zellner-Siow analysis prior
    plocation = 0, pscale = 1/sqrt(2), pdf = 1, alternative = "greater",
    ## specify normal design prior around SMD = 0.5 with small stand. deviation
    dpm = 0.5, dpsd = 0.05
)
design # print design summary
plot(design) # plot design under design prior (top) and under H0 (bottom)


## ----"kang-additional-plot", fig.height = 3-----------------------------------
set.seed(42)
nsim <- 1000
k1 <- 1/10
k0 <- 10
kangsimres <- do.call("rbind", lapply(X = seq(1, nsim), FUN = function(j) {
    ## permute data order
    control <- sample(x = subset(kangdat, dose == "control")$weight_loss,
                      replace = FALSE)
    treat <- list("Low" = sample(x = subset(kangdat, dose == "low")$weight_loss,
                                 replace = FALSE),
                  "Medium" = sample(x = subset(kangdat, dose == "medium")$weight_loss,
                                    replace = FALSE),
                  "High" = sample(x = subset(kangdat, dose == "high")$weight_loss,
                                  replace = FALSE))
    ## compute sequential BFs
    do.call("rbind", lapply(X = seq_along(treat), FUN = function(i) {
        res <- analyzekangdat(control = control, treatment = treat[[i]],
                              pm = pm, psd = psd)
        decision <- "indecisive"
        for (ni in seq(1, nrow(res))) {
            ntreat <- res$ntreatment[ni]
            ncontrol <- res$ncontrol[ni]
            if (res$bf01[ni] > k0) {
                decision <- "H0"
                break
            }
            if (res$bf01[ni] < k1) {
                decision <- "H1"
                break
            }
        }
        ntreatsaved <- res$ntreatment[nrow(res)] - ntreat
        ncontrolsaved <- res$ncontrol[nrow(res)] - ncontrol
        data.frame(simulation = j, group = names(treat)[i], decision = decision,
                   ntreat = ntreat, ncontrol = ncontrol,
                   ntreatsaved = ntreatsaved, ncontrolsaved = ncontrolsaved)
    }))
}))

kangsimsummaries <- kangsimres |>
    group_by(simulation) |>
    summarise(ntreatsaved = sum(ntreatsaved),
              ncontrolsaved = min(ncontrolsaved)) |>
    mutate(nsavedtotal = ntreatsaved + ncontrolsaved)

ggplot(data = kangsimsummaries, aes(x = nsavedtotal)) +
    geom_bar(color = 1, fill = "lightgrey", width = 0.5) +
    scale_x_continuous(breaks = seq(0, 50, 1)) +
    labs(x = "Number of rats saved", y = "Count") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position = "top",
          strip.background = element_rect(fill = "white"))


## ----"kang-additional-table", results = "asis"--------------------------------
## compute proportion of permutations that were stopped for H0/H1
kangsimtable <- kangsimres |>
    mutate(group = factor(group, levels = c("Low", "Medium", "High"))) |>
    group_by(group) |>
    summarise(pstopH0 = paste0(round(mean(decision == "H0")*100, 1), "\\%"),
              pstopH1 = paste0(round(mean(decision == "H1")*100, 1), "\\%"),
              pindecisive = paste0(round(mean(decision == "indecisive")*100, 1), "\\%"))
restablex <- xtable::xtable(kangsimtable)
colnames(restablex) <- c("Treatment group",
                         "Stop for $H_0$",
                         "Stop for $H_1$",
                         "Indecisive")
print(restablex, floating = FALSE, include.rownames = FALSE,
      booktabs = TRUE,
      sanitize.text.function = function(x){x})



## ----"sessionInfo2", echo = Reproducibility, results = Reproducibility--------
cat(paste(Sys.time(), Sys.timezone(), "\n"))
sessionInfo()

