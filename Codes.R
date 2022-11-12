# estimate zibb， p and nu
# then p*nu
library(readr)
library(ggplot2)
raw_data <- read_csv("raw_data.csv")
papers <- c(1,2,3,6,7, 10,13, 15:22)
length(papers) # 15
asym_data <- raw_data[papers, c(2,4,5,8,13,14)]

asym_data$rate <- asym_data$`No. of infections among asymptomatic contacts`/asym_data$`No. of contacts with asymptomatic infection`
mean(asym_data$rate==0)
mean(asym_data$rate[asym_data$rate>0])

# analysis
library(gamlss)
mod1<-gamlss(asym_data$rate~1,sigma.formula=~1, nu.formula=~1, family=BEZI)
fitted(mod1)[1]
summary(mod1)

# expected value of the response
res <- c(fitted(mod1,"nu")[1], fitted(mod1,"mu")[1], fitted(mod1,"sigma")[1], meanBEZI(mod1)[1])
names(res) <- c('nu', 'mu', 'phi', 'estimation')
res

# qBEZI(0.975, mu = res[2], sigma = res[3], nu = res[1], lower.tail = TRUE, log.p = FALSE)

# simulation to see the accuracy
{
  N <- c(5, 10, 15, 20, 30)
  Mu <- Nu <- Phi <- Est <- matrix(NA, 100, length(N))
  for(i in 1:length(N)) {
    for(j in 1:100){
      dat<-rBEZI(N[i], mu=res[2], sigma=res[3], nu=res[1])
      if(all(dat==0)) next
      mod1<-gamlss(dat~1,sigma.formula=~1, nu.formula=~1, family=BEZI)
      Mu[j,i] <- fitted(mod1,"mu")[1]
      Nu[j,i] <- fitted(mod1,"nu")[1]
      Phi[j,i] <- fitted(mod1,"sigma")[1]
      Est[j,i] <- meanBEZI(mod1)[1]
    }
  }
  Phi[Phi>1e+3] <- NA
  # evalute the estimation
  colnames(Mu) <- colnames(Nu) <- colnames(Phi) <- colnames(Est) <- as.character(N)
  par(mfrow=c(2,2))
  boxplot(Mu, main=expression(mu)); abline(a=res[2], b=0, col=2)
  boxplot(Nu, main='nu'); abline(a=res[1], b=0, col=2)
  boxplot(Phi, main='phi'); abline(a=res[3], b=0, col=2)
  boxplot(Est, main='E(y)'); abline(a=res[4], b=0, col=2)
  
  # https://www.benjaminackerman.com/post/2019-03-08-equation_labels/
  dat1 <- data.frame(value=c(c(Nu), c(Mu), c(Phi), c(Est)),
                     N =rep(N, each=100), type=rep(c('pi', 'mu', 'phi','E(y)'), each=100*length(N)))
  dat1$N <- factor(dat1$N, levels = N)
  dat1$type <- factor(dat1$type, levels=c('pi', 'mu', 'phi','E(y)'))
  dat1$value[dat1$value>500] <- NA
  est <- data.frame(type=c('pi', 'mu', 'phi','E(y)'), value=res)
  est$type <- factor(est$type, levels = c('pi', 'mu', 'phi','E(y)'))
  g1 <- ggplot(data=dat1, aes(x=N, y=value, color=N, fill=N)) + 
    geom_boxplot(outlier.alpha = 0.1) + theme_classic() +
    facet_wrap(.~type, nrow=1, scales = 'free_y', label = "label_parsed") + 
    geom_hline(data=est, aes(yintercept=value),linetype="solid",color="black",size=0.7) +
    labs(x='', y='asymptomatic') + theme(legend.position = "none",
                                         strip.text = element_text(size=12),
                                         axis.title.y = element_text(size=16))
  
  
  
  estmeansd <- data.frame(N=factor(N), value=colMeans(Est, na.rm = T), 
                          sd=apply(Est, 2, function(x) sd(x, na.rm = T)))
  estmeansd$CI_lower=estmeansd$value -1.96 * estmeansd$sd
  estmeansd$CI_upper=estmeansd$value +1.96 * estmeansd$sd
  estmeansd$type <- '95% CI of E(y)'
  
  
  g2 <- ggplot(data=estmeansd, aes(x=N, y=value, color=N)) +
    geom_hline(yintercept = res[4], size=0.7) + theme_classic() +
    geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), width=.5) +
    geom_line() + theme(legend.position = "none", strip.text = element_text(size=12))  +
    geom_point(size=2) + labs(y='', x='') + facet_grid(.~type) 
  
  library(ggpubr)
  g12 <- ggarrange(g1, g2, nrow=1, ncol=2, widths = c(3.5,1)) 
  
  
}


# replicate 1000 times to calculate the sd for the estimation
mod1<-gamlss(asym_data$rate~1,sigma.formula=~1, nu.formula=~1, family=BEZI)
res <- c(fitted(mod1,"nu")[1], fitted(mod1,"mu")[1], fitted(mod1,"sigma")[1], meanBEZI(mod1)[1])
names(res) <- c('nu', 'mu', 'phi', 'estimation')
res
sum(asym_data$rate==0)
mean(asym_data$rate[asym_data$rate>0])

Est <-  c()
while(length(Est)<1001){
  dat<-rBEZI(15, mu=res[2], sigma=res[3], nu=res[1])
  # if(sum(dat==0)>12) next
  mod1<-gamlss(dat~1,sigma.formula=~1, nu.formula=~1, family=BEZI)
  Est <- c(Est,meanBEZI(mod1)[1])
}
Est <- Est[!is.na(Est)]; length(Est)
x <- c(mean(Est), sd(Est), mean(Est)-1.96*sd(Est), mean(Est)+1.96*sd(Est))

round(x*100, 4)


#----sub-population-----
sub <- asym_data[asym_data$Population=='Community', ]
sub$rate
mean(sub$rate[sub$rate>0])
mod1<-gamlss(sub$rate~1,sigma.formula=~1, nu.formula=~1, family=BEZI)
Est <-  c()
while(length(Est)<1001){
  dat<-rBEZI(nrow(sub), mu=fitted(mod1,"mu")[1], sigma=fitted(mod1,"sigma")[1], nu=fitted(mod1,"nu")[1])
  if(sum(dat==0)!=2) next
  mod2<-gamlss(dat~1,sigma.formula=~1, nu.formula=~1, family=BEZI)
  Est <-c(Est,meanBEZI(mod2)[1])
}
Est <- Est[!is.na(Est)]; length(Est)
x <- c(mean(Est), sd(Est), mean(Est)-1.96*sd(Est), mean(Est)+1.96*sd(Est))
round(x*100, 2)


sub <- asym_data[asym_data$Population!='Community', ]
sum(sub$rate==0)
mean(sub$rate[sub$rate>0]) # 0.03430434
mod1<-gamlss(sub$rate~1,sigma.formula=~1, nu.formula=~1, family=BEZI)
Est <-  c()
while(length(Est)<1001){
  dat<-rBEZI(nrow(sub), mu=fitted(mod1,"mu")[1], sigma=fitted(mod1,"sigma")[1], nu=fitted(mod1,"nu")[1])
  if(sum(dat==0)>6) next
  mod2<-gamlss(dat~1,sigma.formula=~1, nu.formula=~1, family=BEZI)
  Est <-c(Est,meanBEZI(mod2)[1])
}
Est <- Est[!is.na(Est)]; length(Est)
x <- c(mean(Est), sd(Est), mean(Est)-1.96*sd(Est), mean(Est)+1.96*sd(Est))
round(x*100, 2)


# sub-household

house <- c(2/45, 0/196, 5/42, 0/14, 0/4)
nohouse <- c(1/61, 0/259, 2/80, 0/38)
mean(house[house>0])
mean(nohouse[nohouse>0])

mod1<-gamlss(house~1,sigma.formula=~1, nu.formula=~1, family=BEZI)
Est <-  c()
while(length(Est)<1001){
  dat<-rBEZI(length(house), mu=fitted(mod1,"mu")[1], sigma=fitted(mod1,"sigma")[1], nu=fitted(mod1,"nu")[1])
  if(sum(dat==0)!=2) next
  mod2<-gamlss(dat~1,sigma.formula=~1, nu.formula=~1, family=BEZI)
  Est <- c(Est, meanBEZI(mod2)[1])
}
Est <- Est[!is.na(Est)]; length(Est)
x <- c(mean(Est), sd(Est), mean(Est)-1.96*sd(Est), mean(Est)+1.96*sd(Est))
round(x*100, 2)


# country
china <- asym_data$rate[c(1,4,5,8:12,14,15)]
others <- asym_data$rate[-c(1,4,5,8:12,14,15)]
length(china)
sum(china==0)
mean(china[china>0])
mod1<-gamlss(china~1,sigma.formula=~1, nu.formula=~1, family=BEZI)
Est <-  c()
while(length(Est)<1001){
  dat<-rBEZI(length(china), mu=fitted(mod1,"mu")[1], sigma=fitted(mod1,"sigma")[1], nu=fitted(mod1,"nu")[1])
  # if(sum(dat==0)>3) next
  mod2<-gamlss(dat~1,sigma.formula=~1, nu.formula=~1, family=BEZI)
  Est <- c(Est, meanBEZI(mod2)[1])
}
Est <- Est[!is.na(Est)]; length(Est)
x <- c(mean(Est), sd(Est), mean(Est)-1.96*sd(Est), mean(Est)+1.96*sd(Est))
round(x*100, 2)


# Pre-symptomatic
asy <- asym_data$rate
pre <- c(12/585, 7/88, 23/236, 0/11, 47/922, 11/250)
mean(pre[pre>0])

mod1<-gamlss(pre~1, sigma.formula=~1, nu.formula=~1, family=BEZI)

Mu <- Nu <- Phi <- Est <-  c()
while(length(Est)<1001){
  dat<-rBEZI(length(pre), mu=fitted(mod1,"mu")[1], sigma=fitted(mod1,"sigma")[1], nu=fitted(mod1,"nu")[1])
  if(sum(dat==0)>2) next
  mod2<-gamlss(dat~1,sigma.formula=~1, nu.formula=~1, family=BEZI)
  Mu <- c(Mu, fitted(mod2,"mu")[1])
  Nu <- c(Nu, fitted(mod2,"nu")[1])
  Phi <- c(Phi, fitted(mod2,"sigma")[1])
  Est <- c(Est, meanBEZI(mod2)[1])
}
Est <- Est[!is.na(Est)]; length(Est)
x <- c(mean(Est), sd(Est), mean(Est)-1.96*sd(Est), mean(Est)+1.96*sd(Est))
round(x*100, 2) 

{
  # res <- c(fitted(mod1,"nu")[1], fitted(mod1,"mu")[1], fitted(mod1,"sigma")[1], meanBEZI(mod1)[1])
  # names(res) <- c('nu', 'mu', 'phi', 'estimation')
  # res
  res <- c(mean(Nu), mean(Mu), mean(Phi), mean(Est))
  
  N <- c(10, 15, 20, 30)
  Mu <- Nu <- Phi <- Est <- matrix(NA, 100, length(N))
  for(i in 1:length(N)) {
    for(j in 1:100){
      dat<-rBEZI(N[i], mu=res[2], sigma=res[3], nu=res[1])
      if(all(dat==0)) next
      mod1<-gamlss(dat~1,sigma.formula=~1, nu.formula=~1, family=BEZI)
      Mu[j,i] <- fitted(mod1,"mu")[1]
      Nu[j,i] <- fitted(mod1,"nu")[1]
      Phi[j,i] <- fitted(mod1,"sigma")[1]
      Est[j,i] <- meanBEZI(mod1)[1]
    }
  }
  dat1 <- data.frame(value=c(c(Nu), c(Mu), c(Phi), c(Est)),
                     N =rep(N, each=100), type=rep(c('pi', 'mu', 'phi','E(y)'), each=100*length(N)))
  dat1$N <- factor(dat1$N, levels = N)
  dat1$type <- factor(dat1$type, levels=c('pi', 'mu', 'phi','E(y)'))
  dat1$value[dat1$value>500] <- NA
  est <- data.frame(type=c('pi', 'mu', 'phi','E(y)'), value=res)
  est$type <- factor(est$type, levels = c('pi', 'mu', 'phi','E(y)'))
  g3 <- ggplot(data=dat1, aes(x=N, y=value, color=N, fill=N)) + 
    geom_boxplot(outlier.alpha = 0.1) + theme_classic() +
    facet_wrap(.~type, nrow=1, scales = 'free_y', label = "label_parsed") + 
    geom_hline(data=est, aes(yintercept=value),linetype="solid",color="black",size=0.7) +
    labs(x='', y='presymptomatic') + theme(legend.position = "none",
                                           strip.text = element_text(size=12),
                                           axis.title.y = element_text(size=16))
  
  estmeansd <- data.frame(N=factor(N), value=colMeans(Est, na.rm = T), 
                          sd=apply(Est, 2, function(x) sd(x, na.rm = T)))
  estmeansd$CI_lower=estmeansd$value -1.96 * estmeansd$sd
  estmeansd$CI_upper=estmeansd$value +1.96 * estmeansd$sd
  estmeansd$type <- '95% CI of E(y)'
  
  g4 <- ggplot(data=estmeansd, aes(x=N, y=value, color=N)) +
    geom_hline(yintercept = res[4], size=0.7) + theme_classic() +
    geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), width=.5) +
    geom_line() + theme(legend.position = "none", strip.text = element_text(size=12))  +
    geom_point(size=2) + labs(y='', x='') + facet_grid(.~type) 
  
  g34 <- ggarrange(g3, g4, nrow=1, ncol=2, widths = c(3.5,1)) 
  
}

# symptomatic
sym <- c(28/1010,126/2001, 0/52, 5/130, 22/2644,117/2305,34/210, 20/471)
mean(sym[sym>0])
mod1<-gamlss(sym~1,sigma.formula=~1, nu.formula=~1, family=BEZI)
Mu <- Nu <- Phi <- Est <-  c()
while(length(Est)<1001){
  dat<-rBEZI(length(sym), mu=fitted(mod1,"mu")[1], sigma=fitted(mod1,"sigma")[1], nu=fitted(mod1,"nu")[1])
  if(sum(dat==0)>1) next
  mod2<-gamlss(dat~1,sigma.formula=~1, nu.formula=~1, family=BEZI)
  Mu <- c(Mu, fitted(mod2,"mu")[1])
  Nu <- c(Nu, fitted(mod2,"nu")[1])
  Phi <- c(Phi, fitted(mod2,"sigma")[1])
  Est <- c(Est, meanBEZI(mod2)[1])
}
Est <- Est[!is.na(Est)]; length(Est)
x <- c(mean(Est), sd(Est), mean(Est)-1.96*sd(Est), mean(Est)+1.96*sd(Est))
round(x*100, 2)

{
  # res <- c(fitted(mod1,"nu")[1], fitted(mod1,"mu")[1], fitted(mod1,"sigma")[1], meanBEZI(mod1)[1])
  # names(res) <- c('nu', 'mu', 'phi', 'estimation')
  # res
  res <- c(mean(Nu), mean(Mu), mean(Phi), mean(Est))
  N <- c(10, 15, 20, 30)
  Mu <- Nu <- Phi <- Est <- matrix(NA, 100, length(N))
  for(i in 1:length(N)) {
    for(j in 1:100){
      dat<-rBEZI(N[i], mu=res[2], sigma=res[3], nu=res[1])
      if(all(dat==0)) next
      mod1<-gamlss(dat~1,sigma.formula=~1, nu.formula=~1, family=BEZI)
      Mu[j,i] <- fitted(mod1,"mu")[1]
      Nu[j,i] <- fitted(mod1,"nu")[1]
      Phi[j,i] <- fitted(mod1,"sigma")[1]
      Est[j,i] <- meanBEZI(mod1)[1]
    }
  }
  dat1 <- data.frame(value=c(c(Nu), c(Mu), c(Phi), c(Est)),
                     N =rep(N, each=100), type=rep(c('pi', 'mu', 'phi','E(y)'), each=100*length(N)))
  dat1$N <- factor(dat1$N, levels = N)
  dat1$type <- factor(dat1$type, levels=c('pi', 'mu', 'phi','E(y)'))
  dat1$value[dat1$value>500] <- NA
  est <- data.frame(type=c('pi', 'mu', 'phi','E(y)'), value=res)
  est$type <- factor(est$type, levels = c('pi', 'mu', 'phi','E(y)'))
  g5 <- ggplot(data=dat1, aes(x=N, y=value, color=N, fill=N)) + 
    geom_boxplot(outlier.alpha = 0.1) + theme_classic() +
    facet_wrap(.~type, nrow=1, scales = 'free_y', label = "label_parsed") + 
    geom_hline(data=est, aes(yintercept=value),linetype="solid",color="black",size=0.7) +
    labs(x='', y='symptomatic') + theme(legend.position = "none",
                                        strip.text = element_text(size=12),
                                        axis.title.y = element_text(size=16))
  
  estmeansd <- data.frame(N=factor(N), value=colMeans(Est, na.rm = T), 
                          sd=apply(Est, 2, function(x) sd(x, na.rm = T)))
  estmeansd$CI_lower=estmeansd$value -1.96 * estmeansd$sd
  estmeansd$CI_upper=estmeansd$value +1.96 * estmeansd$sd
  estmeansd$type <- '95% CI of E(y)'
  
  g6 <- ggplot(data=estmeansd, aes(x=N, y=value, color=N)) +
    geom_hline(yintercept = res[4], size=0.7) + theme_classic() +
    geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), width=.5) +
    geom_line() + theme(legend.position = "none", strip.text = element_text(size=12))  +
    geom_point(size=2) + labs(y='', x='') + facet_grid(.~type) 
  
  g56 <- ggarrange(g5, g6, nrow=1, ncol=2, widths = c(3.5,1)) 
  
}

{
  g16 <- ggarrange(g12, g34, g56, nrow=3, ncol=1, labels = c('(a)', '(b)','(c)'))
  g16 <- annotate_figure(g16, bottom = text_grob('sample size', 
                                                 size=16, color='black'))
  ggsave('eval1.pdf', g16, width = 12, height = 8)
}
# sensitivity analyses
# asymptomatic
res_asy <- matrix(NA, 15, 6)
for(i in 1:15){
  x <- asym_data$`No. of infections among asymptomatic contacts`[-i]
  y <- asym_data$`No. of contacts with asymptomatic infection`[-i]
  tmp1 <- x/y
  mod1<-gamlss(tmp1~1,sigma.formula=~1, nu.formula=~1, family=BEZI)
  Est <-  c()
  while(length(Est)<1001){
    dat<-rBEZI(length(tmp1), mu=fitted(mod1,"mu")[1], sigma=fitted(mod1,"sigma")[1], nu=fitted(mod1,"nu")[1])
    mod2<-gamlss(dat~1,sigma.formula=~1, nu.formula=~1, family=BEZI)
    Est <- c(Est, meanBEZI(mod2)[1])
  }
  Est <- Est[!is.na(Est)]; length(Est)
  z <- c(mean(Est), sd(Est), mean(Est)-1.96*sd(Est), mean(Est)+1.96*sd(Est))
  z <- round(z*100, 2)
  res_asy[i, ] <- c(as.character(asym_data[i, 1]), sum(x), sum(y), z[1], z[3], z[4])
}
View(res_asy)
write.csv(res_asy, 's_asy.csv')
c
# Pre-symptomatic
pre_a <- c('Chaw 2020', 'Han 2020', 'Liu 2020', 'Park 2020', 'Shi 2020', 'Zhang 2020' )
pre_x <- c(12, 7, 23, 0, 47, 11)
pre_y <- c(585, 88, 236, 11, 922,250)
res_pre <- matrix(NA, length(pre_x), 6)
for(i in 1:length(pre_y)){
  x <- pre_x[-i]
  y <- pre_y[-i]
  tmp1 <- x/y
  mod1<-gamlss(tmp1~1,sigma.formula=~1, nu.formula=~1, family=BEZI)
  Est <-  c()
  while(length(Est)<1001){
    dat<-rBEZI(length(tmp1), mu=fitted(mod1,"mu")[1], sigma=fitted(mod1,"sigma")[1], nu=fitted(mod1,"nu")[1])
    if(sum(dat==0)>3) next
    mod2<-gamlss(dat~1,sigma.formula=~1, nu.formula=~1, family=BEZI)
    Est <- c(Est, meanBEZI(mod2)[1])
  }
  Est <- Est[!is.na(Est)]; length(Est)
  z <- c(mean(Est), sd(Est), mean(Est)-1.96*sd(Est), mean(Est)+1.96*sd(Est))
  z <- round(z*100, 2)
  res_pre[i, ] <- c(pre_a[i], sum(x), sum(y), z[1], z[3], z[4])
}
View(res_pre)
write.csv(res_pre, 's_pre.csv')

# symptomatic
sym_a <- c('Chaw 2020', 'Yin 2020', 'Han 2020', 'Jiang X 2020', 'Cheng 2020', 'Luo 2020', 'Park 2020', 'Shi 2020')
sym_x <- c(8, 126, 0, 5, 22,117,34, 20)
sym_y <- c(1010,2001, 52, 130, 2644,2305,210, 471)
res_sym <- matrix(NA, length(sym_x), 6)
for(i in 1:length(sym_a)){
  x <- sym_x[-i]
  y <- sym_y[-i]
  tmp1 <- x/y
  mod1<-gamlss(tmp1~1,sigma.formula=~1, nu.formula=~1, family=BEZI)
  Est <-  c()
  while(length(Est)<1001){
    dat<-rBEZI(length(tmp1), mu=fitted(mod1,"mu")[1], sigma=fitted(mod1,"sigma")[1], nu=fitted(mod1,"nu")[1])
    if(sum(dat==0)>2) next
    mod2<-gamlss(dat~1,sigma.formula=~1, nu.formula=~1, family=BEZI)
    Est <- c(Est, meanBEZI(mod2)[1])
  }
  Est <- Est[!is.na(Est)]; length(Est)
  z <- c(mean(Est), sd(Est), mean(Est)-1.96*sd(Est), mean(Est)+1.96*sd(Est))
  z <- round(z*100, 2)
  res_sym[i, ] <- c(sym_a[i], sum(x), sum(y), z[1], z[3], z[4])
}
View(res_sym)
write.csv(res_sym, 's_sym.csv')


## significant test
{
  # asy vs pre
  xa <- asym_data$`No. of infections among asymptomatic contacts`
  ya <- asym_data$`No. of contacts with asymptomatic infection`
  
  xp <- c(12, 7, 23, 0, 47, 11)
  yp <- c(585, 88, 236, 11, 47, 922, 250)
  
  prop.test(x=c(sum(xa), sum(xp)), n=c(sum(ya), sum(yp))) # 3.158e-13
  prop.test(x=c(sum(xa), sum(xp)), n=c(sum(ya), sum(yp)), alternative = 'less', correct = T) # 1.579e-13
  
  # asy vs sym
  xs <- c(28, 126, 0, 5, 22, 117, 34, 20)
  ys <- c(1010, 2001, 52, 130, 2644, 2305, 210, 471)
  
  prop.test(x=c(sum(xa), sum(xs)), n=c(sum(ya), sum(ys))) # 3.709e-13
  
  # pre vs sym
  prop.test(x=c(sum(xp), sum(xs)), n=c(sum(yp), sum(ys))) #  0.1707
  
  
  # house vs nohouse
  xh <- c(2, 0, 5, 0, 0)
  yh <- c(45, 196, 42, 14, 4)
  
  xnh <- c(1, 0, 2, 0)
  ynh <- c(61, 259, 80, 38)
  
  prop.test(x=c(sum(xh), sum(xnh)), n=c(sum(yh), sum(ynh)), alternative = 'greater', correct = F)
  # 0.02893
  
  # china vs others
  xc <- xa[c(1,4,5,8:12,14,15)]
  yc <- ya[c(1,4,5,8:12,14,15)]
  
  xo <- xa[-c(1,4,5,8:12,14,15)]
  yo <- ya[-c(1,4,5,8:12,14,15)]
  
  prop.test(x=c(sum(xc), sum(xo)), n=c(sum(yc), sum(yo)), alternative = 'less') # 0.01427
}













