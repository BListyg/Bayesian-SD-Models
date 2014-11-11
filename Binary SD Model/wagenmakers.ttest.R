#generate binary data and analyze it
rm(list=ls())
require(rjags)
require(polspline)
require(ggplot2)
source('sim.data.markov.R')
source('multiplot.R')
source('HDIofMCMC.R')
#----simuluated data set parameters----
p.rep01 = .8
p.rep11 = .8
n.subj = 20
n.sims = 50 * n.subj

trans.mat = matrix(c(1-p.rep01,p.rep01,
                     1-p.rep11,p.rep11),2,2,byrow=T)
#----mcmc parameters----
parameters = c('delta')     
adapt.steps = 500               
burnin.steps = 3000
n.chains = 3 
num.saved.steps=50000
thin.steps = 1
check.conv=F
wd = getwd()
model.file.name = '/models/oneSample_t-test_binomial.txt'
model.string = paste(wd,model.file.name,'sep'='')

#----simulate data set---- 

mcmc = sim.data.markov(trans.mat,n.sims)
subj.num = rep(1:n.subj,each=(n.sims/n.subj))
binary.frame = data.frame(response=mcmc)

#----recode data----
#1=target 0=foil, 1=yes, 0=no
binary.frame = binary.frame - 1
attach(binary.frame)
binary.frame$n1r <- c(rep(0,1),head(response,-1))
attach(binary.frame)
binary.frame$YY[n1r==1 & response==1]=1
binary.frame$YY[n1r==1 & response==0]=0
binary.frame$NY[n1r==0 & response==1]=1
binary.frame$NY[n1r==0 & response==0]=0
binary.frame$YN = as.numeric(!binary.frame$YY)
binary.frame$NN = as.numeric(!binary.frame$NY)
binary.frame$Y_ = binary.frame$YY + binary.frame$YN
binary.frame$N_ = binary.frame$NY + binary.frame$NN

#----aggregate data----
attach(binary.frame)
agg.data = aggregate(binary.frame, by = list(subj.num), FUN = sum, na.rm=TRUE)
detach(binary.frame)
attach(agg.data)
##----run the model----

data.list = list(YY=YY,NY=NY,
                 Y_=Y_,N_=N_,
                 NSUBJ=n.subj)

n.iter = ceiling( ( num.saved.steps * thin.steps ) / n.chains )
jags.model = jags.model(model.string, data=data.list,  
                        n.chains=n.chains, n.adapt=adapt.steps, inits= )
cat( "Burning in the MCMC chain...\n" )
update(jags.model, n.iter=burnin.steps )
cat( "Sampling final MCMC chain...\n" )
coda.samples = coda.samples(jags.model, variable.names=parameters, 
                            n.iter=n.iter, thin=thin.steps )
#----check convergence----
if(check.conv){
  windows()
  gelman.plot(coda.samples)
  
  windows()
  xyplot(coda.samples)
  
  windows()
  densityplot(coda.samples)
}

#---examine the results----
mcmc.chain = as.matrix( coda.samples )

delta.alpha = mcmc.chain[,paste('delta')]
HDIofMCMC(delta.alpha)
#----BFs for delta.alpha----
fit.posterior.alpha = logspline(delta.alpha)

# 95% confidence interval:
x0.alpha=qlogspline(0.025,fit.posterior.alpha)
x1.alpha=qlogspline(0.975,fit.posterior.alpha)

posterior.alpha     = dlogspline(0, fit.posterior.alpha) # this gives the pdf at point delta = 0
prior.alpha         = dnorm(0)                     #height of the prior  
BF01.alpha          = posterior.alpha/prior.alpha

delta.alpha.density = density(delta.alpha) #density object
delta.alpha.density.y = delta.alpha.density$y 
delta.alpha.density.x = delta.alpha.density$x
delta.alpha.density.frame <- data.frame(densityType = rep(c('Posterior'),length(delta.alpha.density.x)),x=delta.alpha.density.x,y=delta.alpha.density.y)
x = seq(from = -4,to=4,length=10000)
y = dnorm(x)
dnorm.frame <- data.frame(densityType = rep('Prior',length(x)),x=x,y=y)
delta.alpha.frame <- rbind(delta.alpha.density.frame,dnorm.frame)
delta.alpha.points.frame <- data.frame(densityType = c("Posterior","Prior"),x=c(0,0),y=c(posterior.alpha,prior.alpha))


#----plot the results----
p1 <- ggplot(data=delta.alpha.frame,aes(x=x,y=y,group=densityType,color=densityType)) +
  geom_line() +
  geom_abline(intercept=0,slope=0) +
  scale_linetype_manual(values=c("solid", "longdash")) +
  geom_point(data=delta.alpha.points.frame,aes(x=x,y=y),show_guide=F,size=5) +
  xlab(expression(delta)) +
  ylab("Density") +  
  ggtitle("Bias, No Sequential Dependencies") +
  labs(linetype="") +
  scale_x_continuous(breaks = seq(from = -4,to=4,by=1),limits=c(-4,4)) +
  theme_bw() + 
  theme(axis.line = element_line(color = "black"),
        axis.title = element_text(size=20),
        axis.title.y = element_text(vjust=.25),
        axis.title.x = element_text(vjust=.25),
        axis.text = element_text(size=15),
        axis.ticks = element_line(color='black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.justification = c(1,1),
        legend.position = c(1,1),
        legend.text = element_text(size=15),
        legend.title = element_blank()
  )



png(file="sd.ttest.plot.png",width=800,height=800,res=160)
p1
dev.off()