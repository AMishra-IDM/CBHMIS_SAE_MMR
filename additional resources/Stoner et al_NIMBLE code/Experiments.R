#########################################
##     Correcting Under-Reporting      ##
##       Simulation Experiments        ##   
#########################################

# Part One: Observed values versus prior information.
sim_code_1=nimbleCode({
  for(i in 1:N){
    pi[i] <- obs[i]+(1-obs[i])*ilogit(b0+b1*w[i])
    lambda[i] <- exp(a0+a1*x[i]+phi[i])
    z[i] ~ dpois(lambda[i]*pi[i])
  }
  phi[1:N] ~ dcar_normal(adj=adj[1:l_adj],num=n_adj[1:N],tau=tau,zero_mean = 1)
  a0 ~ dnorm(0,sd=10)
  a1 ~ dnorm(0,sd=10)
  b0 ~ dnorm(0,sd=prior_sd)
  b1 ~ dnorm(0,sd=10)
  nu ~ T(dnorm(0,1),0,1)
  tau <- 1/nu^2
})

proportion_observed_1=sort(rep(c(0,0.05,0.1,0.15,0.2,0.25,0.3),6)) # Proportion of y completely observed.
prior_sd_1=rep(c(1.2,1,0.8,0.6,0.4,0.2),7) # Standard deviation of prior for beta_0.

# Setup NIMBLE.
sim_constants_1=list(N=sim$N,x=sim$x,w=sim$w,adj=sim$adj,l_adj=length(sim$adj),n_adj=sim$n_adj)
sim_data_1=list(z=sim$z,prior_sd=prior_sd_1[1],obs=rep(0,sim$N))
sim_inits_1=list(a0=0,a1=0,b0=0,b1=0,nu=1,phi=rep(0,sim$N))
sim_model_1 <- nimbleModel(sim_code_1, sim_constants_1,sim_data_1,sim_inits_1)
sim_compiled_model_1<-compileNimble(sim_model_1,resetFunctions = TRUE)

sim_mcmc_conf_1 <- configureMCMC(sim_model_1,monitors=c('a0','a1','b0','b1','tau','lambda','pi'),useConjugacy = FALSE)
sim_mcmc_conf_1$removeSamplers(c('a0','a1','b0','b1','nu'))
# Add automated-factor slice sampler.
sim_mcmc_conf_1$addSampler(target=c('a0','a1','b0','b1','nu'),type='AF_slice',control=list(adaptInterval=1000))

sim_mcmc_1<-buildMCMC(sim_mcmc_conf_1)
sim_compiled_mcmc_1<-compileNimble(sim_mcmc_1, project = sim_model_1,resetFunctions = TRUE)

sim_mcmc_list_1=list()

sim_observed_1=sample(1:sim$N,replace=FALSE,size=max(proportion_observed_1)*sim$N)

# Run the model varying the proportion of y observed 
# and the standard deviation of the prior for beta_0.
for(i in 1:42){
  sim_compiled_model_1$prior_sd=prior_sd_1[i]
  sim_compiled_model_1$obs=rep(0,sim$N)
  sim_compiled_model_1$z=sim$z
  sim_compiled_model_1$obs[sim_observed_1[0:(proportion_observed_1[i]*sim$N)]]=1
  sim_compiled_model_1$z[sim_observed_1[0:(proportion_observed_1[i]*sim$N)]]=
    sim$y[sim_observed_1[0:(proportion_observed_1[i]*sim$N)]]
  sim_samples<-runMCMC(sim_compiled_mcmc_1,inits=sim_inits_1,
                       nchains = 1, nburnin=100000,niter = 200000,samplesAsCodaMCMC = TRUE,
                       summary = FALSE, WAIC = FALSE,thin=thin_multiplier*10,setSeed=seed)
  sim_mcmc_list_1[[i]]=as.mcmc(as.mcmc.list(sim_samples))
}

# Produce y samples.
sim_y_1=list()
for(i in 1:42){
  lambda=sim_mcmc_list_1[[i]][,5:(4+sim$N)]
  pi=sim_mcmc_list_1[[i]][,(5+sim$N):(4+2*sim$N)]
  z=sim$z
  z[sim_observed_1[0:(proportion_observed_1[i]*sim$N)]]=
    sim$y[sim_observed_1[0:(proportion_observed_1[i]*sim$N)]]
  sim_y_1[[i]]=matrix(nrow=dim(lambda)[1],ncol=sim$N)
  for(j in 1:sim$N){
    sim_y_1[[i]][,j]=z[j]+rpois(dim(lambda)[1],lambda[,j]*(1-pi[,j]))
  }
}

# Predictive uncertainty.
error_1=unlist(lapply(sim_y_1,function(x){
  lmse=apply(x,1,function(x){
    log(mean((x-sim$y)^2))
  })
  return(mean(lmse))
}))


ggplot(data.frame(x=proportion_observed_1,y=prior_sd_1,e=error_1))+
  geom_tile(aes(x=x,y=y,fill=e),alpha=1)+
  labs(
    y=expression('S.D. of Prior for '*beta[0]),
    x=expression('Proportion of '*y[s]*' Observed'),
    title=expression(log('Mean Squared Error'))
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    plot.title = element_text(size= 14, hjust=0.01, color = "#4e4d47", margin = margin(b = 0.1, t = 0.4, l = 2, unit = "cm"))
  )+
  scale_fill_viridis(begin=0.3,end=0.7,name='')
ggsave('sim_lmse.pdf',device='pdf',width=4.5,height=3)

# Part Two: Sensitivity to inofmative prior and
# under-reporting covariate.

sim_code_2=nimbleCode({
  for(i in 1:N){
    lambda[i] <- exp(a0+a1*x[i]+phi[i]) 
    pi[i] <- ilogit(b0+b1*v[i]+gamma[i]) 
    z[i] ~ dpois(lambda[i]*pi[i])
    gamma[i] ~ dnorm(0,sd=epsilon)
  }
  phi[1:N] ~ dcar_normal(adj=adj[1:l_adj],num=n_adj[1:N],tau=tau,zero_mean = 1)
  a0 ~ dnorm(0,sd=10)
  a1 ~ dnorm(0,sd=10)
  b0 ~ dnorm(prior_mean,sd=prior_sd)
  b1 ~ dnorm(0,sd=10)  
  nu ~ T(dnorm(0,1),0,1)
  tau <- 1/nu^2
  epsilon ~ T(dnorm(0,sd=1),0,)
})

prior_mean_2=b0+rep(sort(rep(c(-1.8,-1.2,-0.6,0,0.6,1.2,1.8),6)),6) # Vary the mean of the prior for beta_0.
prior_sd_2=rep(rep(c(1.2,1,0.8,0.6,0.4,0.2),7),6) # Vary the standard deviation of the prior for beta_0.
v_index_2=sort(rep(1:6,42)) # Vary the strength of the under-reporting covariate.

# Setup NIMBLE.

sim_constants_2=list(N=sim$N,x=sim$x,adj=sim$adj,l_adj=length(sim$adj),n_adj=sim$n_adj)
sim_data_2=list(z=sim$z,prior_sd=prior_sd_2[1],prior_mean=prior_mean_2[1],v=sim$v[,v_index_2[1]])
sim_inits_2=list(a0=0,a1=0,b0=0,b1=0,nu=1,epsilon=1,phi=rep(0,sim$N)) 

sim_model_2 <- nimbleModel(sim_code_2, sim_constants_2,sim_data_2,sim_inits_2)
sim_compiled_model_2<-compileNimble(sim_model_2,resetFunctions = TRUE)

sim_mcmc_conf_2 <- configureMCMC(sim_model_2,monitors=c('a0','a1','b0','b1','tau','epsilon','lambda','pi'),useConjugacy = FALSE)
sim_mcmc_conf_2$removeSamplers(c('a0','a1','b0','b1','nu','epsilon')) 
sim_mcmc_conf_2$addSampler(target=c('a0','a1','b0','b1','nu','epsilon'),type='AF_slice',control=list(adaptInterval=1000))

sim_mcmc_2<-buildMCMC(sim_mcmc_conf_2)
sim_compiled_mcmc_2<-compileNimble(sim_mcmc_2, project = sim_model_2,resetFunctions = TRUE)

sim_mcmc_list_2=list()

# The following loop can take around 8 hours.
for(i in 1:252){
  sim_compiled_model_2$prior_mean=prior_mean_2[i]
  sim_compiled_model_2$prior_sd=prior_sd_2[i]
  sim_compiled_model_2$v=sim$v[,v_index_2[i]]
  sim_samples_2<-runMCMC(sim_compiled_mcmc_2,inits=sim_inits_2,
                         nchains = 1, nburnin=100000,niter = 200000,samplesAsCodaMCMC = TRUE,
                         summary = FALSE, WAIC = FALSE,thin=thin_multiplier*10,setSeed=seed)
  sim_mcmc_list_2[[i]]=as.mcmc(as.mcmc.list(sim_samples_2))
}

# Produce y samples.
sim_y_2=lapply(sim_mcmc_list_2,function(x){
  lambda=x[,6:(5+sim$N)]
  pi=x[,(6+sim$N):(5+2*sim$N)]
  y=matrix(nrow=dim(lambda)[1],ncol=sim$N)
  for(j in 1:sim$N){
    y[,j]=sim$z[j]+rpois(dim(lambda)[1],lambda[,j]*(1-pi[,j]))
  }
  return(y)
})

# Compute coverage of 95% prediction intervals.
coverage_2=unlist(lapply(sim_y_2,function(x){
  covered=numeric(sim$N)
  for(j in 1:sim$N){
    covered[j]=sim$y[j]>=quantile(x[,j],0.025)&sim$y[j]<=quantile(x[,j],0.975)
  }
  return(mean(covered))
}))

# Compute mean error of log(lambda).
lambda_bias_2=unlist(lapply(sim_mcmc_list_2,function(x){
  log_lambda_means=apply(log(x[,6:(5+sim$N)]),2,mean)
  return(mean(log_lambda_means-log(sim$lambda)))
}))

# Compute root mean squared error of log(lambda).
lambda_error_2=unlist(lapply(sim_mcmc_list_2,function(x){
  log_lambda_means=apply(log(x[,6:(5+sim$N)]),2,mean)
  return(sqrt(mean((log_lambda_means-log(sim$lambda))^2)))
}))

ggdata_2=data.frame(x=prior_mean_2,y=prior_sd_2,z=v_index_2,cov=cor(sim$v)[v_index_2],
                  c=coverage_2,l=lambda_error_2,b=lambda_bias_2)

# Coverage plot when correlation is 0.6.
ggplot(filter(ggdata_2,z==3))+
  geom_tile(aes(x=x,y=y,fill=c))+
  geom_label(aes(x=x,y=y,label=sprintf("%0.2f", round(c, digits = 2)),colour=c),label.r=unit(0.25, "lines"))+
  labs(
    y=expression('S.D. of Prior for '*beta[0]),
    x=expression('Mean of Prior for '*beta[0]),
    title='Correlation 0.6'
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    plot.title = element_text(size= 14, hjust=0.01, color = "#4e4d47", margin = margin(b = 0.1, t = 0.4, l = 2, unit = "cm"))
  )+guides(colour=FALSE,fill=FALSE)+
  scale_fill_viridis(begin=0.1,end=0.5)+scale_colour_viridis(begin=0.1,end=0.5)
ggsave('sim_coverage.pdf',device='pdf',width=4.5,height=3)

# Comparison of variation in prediction interval coverage over the different covariates.
c1=ggplot(filter(ggdata_2,z==1))+geom_label(aes(x=x,y=y,label=c,colour=c))+ggtitle('Correlation 1')+
  xlab(expression('Mean of Prior for '*beta[0]))+ylab(expression('S.D. of Prior for '*beta[0]))+
  guides(colour=FALSE)+ scale_color_viridis(begin=0.1,end=0.5)+
  scale_x_continuous(limits=c(-2,2))+scale_y_continuous(limit=c(0.15,1.25))
c2=ggplot(filter(ggdata_2,z==2))+geom_label(aes(x=x,y=y,label=c,colour=c))+ggtitle('Correlation 0.8')+
  xlab(expression('Mean of Prior for '*beta[0]))+ylab(expression('S.D. of Prior for '*beta[0]))+
  guides(colour=FALSE)+ scale_color_viridis(begin=0.1,end=0.5)+
  scale_x_continuous(limits=c(-2,2))+scale_y_continuous(limit=c(0.15,1.25))
c3=ggplot(filter(ggdata_2,z==3))+geom_label(aes(x=x,y=y,label=c,colour=c))+ggtitle('Correlation 0.6')+
  xlab(expression('Mean of Prior for '*beta[0]))+ylab(expression('S.D. of Prior for '*beta[0]))+
  guides(colour=FALSE)+ scale_color_viridis(begin=0.1,end=0.5)+
  scale_x_continuous(limits=c(-2,2))+scale_y_continuous(limit=c(0.15,1.25))
c4=ggplot(filter(ggdata_2,z==4))+geom_label(aes(x=x,y=y,label=c,colour=c))+ggtitle('Correlation 0.4')+
  xlab(expression('Mean of Prior for '*beta[0]))+ylab(expression('S.D. of Prior for '*beta[0]))+
  guides(colour=FALSE)+ scale_color_viridis(begin=0.1,end=0.5)+
  scale_x_continuous(limits=c(-2,2))+scale_y_continuous(limit=c(0.15,1.25))
c5=ggplot(filter(ggdata_2,z==5))+geom_label(aes(x=x,y=y,label=c,colour=c))+ggtitle('Correlation 0.2')+
  xlab(expression('Mean of Prior for '*beta[0]))+ylab(expression('S.D. of Prior for '*beta[0]))+
  guides(colour=FALSE)+scale_color_viridis(begin=0.1,end=0.5)+
  scale_x_continuous(limits=c(-2,2))+scale_y_continuous(limit=c(0.15,1.25))
c6=ggplot(filter(ggdata_2,z==6))+geom_label(aes(x=x,y=y,label=c,colour=c))+ggtitle('Correlation 0')+
  xlab(expression('Mean of Prior for '*beta[0]))+ylab(expression('S.D. of Prior for '*beta[0]))+
  guides(colour=FALSE)+ scale_color_viridis(begin=0.1,end=0.5)+
  scale_x_discrete(limits=c(-2,2))+scale_y_continuous(limit=c(0.15,1.25))

pdf(file='full_coverage.pdf',width=9,height=9)
multiplot(c1,c2,c3,c4,c5,c6,cols=2)
dev.off()

# Comparison of variation in lambda prediction error over the different covariates.
l1=ggplot(filter(ggdata_2,z==1))+geom_label(aes(x=x,y=y,label=round(l,2),colour=round(l,2)))+ggtitle('Correlation 1')+
  xlab(expression('Mean of Prior for '*beta[0]))+ylab(expression('S.D. of Prior for '*beta[0]))+
  guides(colour=FALSE)+scale_color_viridis(begin=0.1,end=0.5)+
  scale_x_continuous(limits=c(-2,2))+scale_y_continuous(limit=c(0.15,1.25))
l2=ggplot(filter(ggdata_2,z==2))+geom_label(aes(x=x,y=y,label=round(l,2),colour=round(l,2)))+ggtitle('Correlation 0.8')+
  xlab(expression('Mean of Prior for '*beta[0]))+ylab(expression('S.D. of Prior for '*beta[0]))+
  guides(colour=FALSE)+scale_color_viridis(begin=0.1,end=0.5)+
  scale_x_continuous(limits=c(-2,2))+scale_y_continuous(limit=c(0.15,1.25))
l3=ggplot(filter(ggdata_2,z==3))+geom_label(aes(x=x,y=y,label=round(l,2),colour=round(l,2)))+ggtitle('Correlation 0.6')+
  xlab(expression('Mean of Prior for '*beta[0]))+ylab(expression('S.D. of Prior for '*beta[0]))+
  guides(colour=FALSE)+scale_color_viridis(begin=0.1,end=0.5)+
  scale_x_continuous(limits=c(-2,2))+scale_y_continuous(limit=c(0.15,1.25))
l4=ggplot(filter(ggdata_2,z==4))+geom_label(aes(x=x,y=y,label=round(l,2),colour=round(l,2)))+ggtitle('Correlation 0.4')+
  xlab(expression('Mean of Prior for '*beta[0]))+ylab(expression('S.D. of Prior for '*beta[0]))+
  guides(colour=FALSE)+scale_color_viridis(begin=0.1,end=0.5)+
  scale_x_continuous(limits=c(-2,2))+scale_y_continuous(limit=c(0.15,1.25))
l5=ggplot(filter(ggdata_2,z==5))+geom_label(aes(x=x,y=y,label=round(l,2),colour=round(l,2)))+ggtitle('Correlation 0.2')+
  xlab(expression('Mean of Prior for '*beta[0]))+ylab(expression('S.D. of Prior for '*beta[0]))+
  guides(colour=FALSE)+scale_color_viridis(begin=0.1,end=0.5)+
  scale_x_continuous(limits=c(-2,2))+scale_y_continuous(limit=c(0.15,1.25))
l6=ggplot(filter(ggdata_2,z==6))+geom_label(aes(x=x,y=y,label=round(l,2),colour=round(l,2)))+ggtitle('Correlation 0')+
  xlab(expression('Mean of Prior for '*beta[0]))+ylab(expression('S.D. of Prior for '*beta[0]))+
  guides(colour=FALSE)+scale_color_viridis(begin=0.1,end=0.5)+
  scale_x_continuous(limits=c(-2,2))+scale_y_continuous(limit=c(0.15,1.25))

pdf(file='lambda_error.pdf',width=9,height=9)
multiplot(l1,l2,l3,l4,l5,l6,cols=2)
dev.off()

# Coverage, mean error of log(lambda) and R.M.S.E. of log(lambda).
cov1=ggplot(ggdata_2)+
  geom_point(aes(x=cov,y=c),col=vp[7])+
  geom_smooth(aes(x=cov,y=c),se=TRUE,colour=vp[7],fill=vp[7])+
  labs(
    y=expression('P.I. Coverage for '*y[s]),
    x=expression('Correlation with '*w[s])
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA)
  )+guides(colour=FALSE,fill=FALSE)
e1=ggplot(ggdata_2)+
  geom_point(aes(x=cov,y=b),col=vp[9])+
  geom_smooth(aes(x=cov,y=b),se=TRUE,colour=vp[9],fill=vp[9])+
  labs(
    y=expression('Mean Error of '*log(lambda[s])),
    x=expression('Correlation with '*w[s])
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA)
  )+guides(colour=FALSE,fill=FALSE)
e2=ggplot(ggdata_2)+
  geom_point(aes(x=cov,y=l),col=vp[11])+
  geom_smooth(aes(x=cov,y=l),se=TRUE,colour=vp[11],fill=vp[11])+
  labs(
    y=expression('R.M.S.E. of '*log(lambda[s])),
    x=expression('Correlation with '*w[s])
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA)
  )+guides(colour=FALSE,fill=FALSE)

pdf(file='sim_strength.pdf',width=9,height=2.5)
multiplot(cov1,e1,e2,cols=3)
dev.off()

# Examine results for covariate correlation 0.6, prior mean 0.6 
# and prior standard deviation 0.6.

index_2=112
lambda_2=sim_mcmc_list_2[[index_2]][,6:(5+sim$N)]
pi_2=sim_mcmc_list_2[[index_2]][,(6+sim$N):(5+2*sim$N)]
y_2=matrix(nrow=dim(lambda_2)[1],ncol=sim$N)
for(j in 1:sim$N){
  y_2[,j]=sim$z[j]+rpois(dim(lambda_2)[1],lambda_2[,j]*(1-pi_2[,j]))
}

pi_sorted_2=expit(sim_mcmc_list_2[[index_2]][,3]+sim_mcmc_list_2[[index_2]][,4]%*%t(sort(sim$v[,3])))
phi_2=log(lambda_2)-sim_mcmc_list_2[[index_2]][,1]-sim_mcmc_list_2[[index_2]][,2]%*%t(sim$x)

post_1=ggplot(data=data.frame(a0=as.numeric(sim_mcmc_list_2[[index_2]][,1]),
                             x=seq(3.2,4.8,length=1000),prior=dnorm(seq(3,5,length=1000),0,10)))+
  geom_line(aes(x=x,y=prior),colour='#22211d')+
  stat_density(mapping=aes(x=a0),adjust=3,fill=vp[6],alpha=0.5)+
  labs(
    y='Density',
    x=expression(alpha[0])
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
    )+geom_vline(xintercept = a0,colour='#22211d')
post_2=ggplot(data=data.frame(a1=as.numeric(sim_mcmc_list_2[[index_2]][,2]),
                              x=seq(0.6,1.6,length=1000),prior=dnorm(seq(0.25,1.75,length=1000),0,10)))+
  geom_line(aes(x=x,y=prior),colour='#22211d')+
  stat_density(mapping=aes(x=a1),adjust=3,fill=vp[7],alpha=0.5)+
  labs(
    y='Density',
    x=expression(alpha[1])
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )+geom_vline(xintercept = a1,colour='#22211d')
post_3=ggplot(data=data.frame(b0=as.numeric(sim_mcmc_list_2[[index_2]][,3]),
                              x=seq(-1.7,3,length=1000),
                  prior=dnorm(seq(-1.7,3,length=1000),prior_mean_2[index_2],prior_sd_2[index_2])))+
  geom_line(aes(x=x,y=prior),colour='#22211d')+
  stat_density(mapping=aes(x=b0),adjust=3,fill=vp[9],alpha=0.5)+
  labs(
    y='Density',
    x=expression(beta[0])
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )+geom_vline(xintercept = b0,colour='#22211d')
post_4=ggplot(data=data.frame(p=sim$pi[sort.int(sim$v[,3],index.return=TRUE)$ix],x=sort(sim$v[,3]),
                              l=apply(pi_sorted_2,2,quantile,probs=0.025),
                              m=apply(pi_sorted_2,2,quantile,probs=0.5),
                              u=apply(pi_sorted_2,2,quantile,probs=0.975)))+
  geom_point(mapping=aes(x=x,y=p),colour=vp[10])+
  geom_ribbon(aes(x=x,ymin=l,ymax=u),alpha=0.5,fill=vp[10])+
  geom_line(aes(x=x,y=m),colour=vp[10],size=1)+
  labs(
    y=expression('Reporting Probability ('*pi[s]*')'),
    x=expression('Under-Reporting Covariate ('*v['s,3']*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )
post_5=ggplot(data=data.frame(x=sim$phi,y=apply(phi_2,2,mean)))+
  geom_abline(slope=1,intercept=0,colour="#22211d")+
  geom_point(mapping=aes(x=x,y=y),colour=vp[8])+
  labs(
    y=expression('True Spatial Effect ('*phi[s]*')'),
    x=expression('Mean Predicted Spatial Effect ('*phi[s]*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )
pred_1=ggplot(data.frame(x=sim$y,l=apply(y_2,2,quantile,0.025),u=apply(y_2,2,quantile,0.975),
                         m=apply(y_2,2,mean)))+geom_abline(slope=1,intercept=0,colour="#22211d")+
  geom_point(mapping=aes(x=x,y=l),colour=vp[8])+
  geom_point(mapping=aes(x=x,y=u),colour=vp[10])+
  labs(
    y=expression('Predicted Count ('*y[s]*')'),
    x=expression('True Count ('*y[s]*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )


pdf(file='sim_post.pdf',width=9,height=4)
multiplot(post_1,post_3,post_2,post_4,post_5,pred_1,cols=3)
dev.off()

# Part Three: Classification of covariates into process covariate(s)
# and under-reporting covariate.

sim_code_3=nimbleCode({
  for(i in 1:N){
    pi[i] <- ilogit(b0+b1*w[i]+gamma[i])
    lambda[i] <- exp(a0+a1*x[i]+phi[i]+theta[i])
    z[i] ~ dpois(lambda[i]*pi[i])
    gamma[i] ~ dnorm(0,sd=epsilon)
    theta[i] ~ dnorm(0,sd=sigma)
  }
  phi[1:N]~dcar_normal(adj=adj[1:l_adj],num=n_adj[1:N],tau=tau,zero_mean = 1)
  a0 ~ dnorm(0,sd=10)
  a1 ~ dnorm(0,sd=10)
  b0 ~ dnorm(0,sd=0.6)
  b1 ~ dnorm(0,sd=10)
  epsilon ~ T(dnorm(0,1),0,)
  sigma ~ T(dnorm(0,1),0,)
  nu ~ T(dnorm(0,1),0,)
  tau <- 1/nu^2
})

# Setup NIMBLE.
sim_constants_3=list(N=sim$N,adj=sim$adj,l_adj=length(sim$adj),n_adj=sim$n_adj)
sim_data_3=list(z=sim$z,x=sim$x,w=sim$w)
sim_inits_3=list(a0=0,a1=0,b0=0,b1=0,nu=1,phi=rep(0,sim$N),epsilon=0.1,sigma=0.1)
sim_model_3 <- nimbleModel(sim_code_3, sim_constants_3,sim_data_3,sim_inits_3)
sim_compiled_model_3<-compileNimble(sim_model_3,resetFunctions = TRUE)

sim_mcmc_conf_3 <- configureMCMC(sim_model_3,monitors=c('a0','a1','b0','b1','epsilon','tau','lambda','pi','sigma'),useConjugacy = FALSE)
sim_mcmc_conf_3$removeSamplers(c('a0','a1','b0','b1','nu','epsilon','sigma'))
# Add automated-factor slice sampler.
sim_mcmc_conf_3$addSampler(target=c('a0','a1','b0','b1','nu','epsilon','sigma'),type='AF_slice',control=list(adaptInterval=1000))

sim_mcmc_3<-buildMCMC(sim_mcmc_conf_3)
sim_compiled_mcmc_3<-compileNimble(sim_mcmc_3, project = sim_model_3,resetFunctions = TRUE)

# First run the model with the covariates correctly classified.
sim_samples_3_1<-as.mcmc(as.mcmc.list(runMCMC(sim_compiled_mcmc_3,inits=sim_inits_3,
                     nchains = 1, nburnin=100000,niter = 200000,samplesAsCodaMCMC = TRUE,
                     summary = FALSE, WAIC = FALSE,thin=thin_multiplier*10,setSeed=seed)))

lambda_3_1=sim_samples_3_1[,6:(5+sim$N)]
pi_3_1=sim_samples_3_1[,(6+sim$N):(5+2*sim$N)]
y_3_1=matrix(nrow=dim(lambda_3_1)[1],ncol=sim$N)
for(j in 1:sim$N){
  y_3_1[,j]=sim$z[j]+rpois(dim(lambda_3_1)[1],lambda_3_1[,j]*(1-pi_3_1[,j]))
}

class_1=ggplot(data.frame(x=sim$y,m=apply(y_3_1,2,median)))+
  geom_abline(slope=1,intercept=0,colour="#22211d")+xlab(expression('True Count ('*y[s]*')'))+
  geom_point(mapping=aes(x=x,y=m),colour=vp[7])+
  labs(
    y=expression('Median Predicted Count'),
    x=expression('True Count ('*y[s]*')'),
    title='Correct Classification'
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.title = element_text(size= 14, hjust=0.01, color = "#4e4d47", margin = margin(b = 0.1, t = 0.4, l = 2, unit = "cm")),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA)
  )

# Now swap the covariates around.
sim_compiled_model_3$w=sim$x
sim_compiled_model_3$x=sim$w

sim_samples_3_2<-as.mcmc(as.mcmc.list(runMCMC(sim_compiled_mcmc_3,inits=sim_inits_3,
                                              nchains = 1, nburnin=100000,niter = 200000,samplesAsCodaMCMC = TRUE,
                                              summary = FALSE, WAIC = FALSE,thin=thin_multiplier*10,setSeed=seed)))

lambda_3_2=sim_samples_3_2[,6:(5+sim$N)]
pi_3_2=sim_samples_3_2[,(6+sim$N):(5+2*sim$N)]
y_3_2=matrix(nrow=dim(lambda_3_2)[1],ncol=sim$N)
for(j in 1:sim$N){
  y_3_2[,j]=sim$z[j]+rpois(dim(lambda_3_2)[1],lambda_3_2[,j]*(1-pi_3_2[,j]))
}

class_2=ggplot(data.frame(x=sim$y,m=apply(y_3_2,2,median)))+
  geom_abline(slope=1,intercept=0,colour="#22211d")+xlab(expression('True Count ('*y[s]*')'))+
  geom_point(mapping=aes(x=x,y=m),colour=vp[11])+
  labs(
    y=expression('Median Predicted Count'),
    x=expression('True Count ('*y[s]*')'),
    title='Incorrect Classification'
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.title = element_text(size= 14, hjust=0.01, color = "#4e4d47", margin = margin(b = 0.1, t = 0.4, l = 2, unit = "cm")),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA)
  )

# Now remove both covariates.
sim_compiled_model_3$w=rep(0,sim$N)
sim_compiled_model_3$x=rep(0,sim$N)

sim_samples_3_3<-as.mcmc(as.mcmc.list(runMCMC(sim_compiled_mcmc_3,inits=sim_inits_3,
                                              nchains = 1, nburnin=100000,niter = 200000,samplesAsCodaMCMC = TRUE,
                                              summary = FALSE, WAIC = FALSE,thin=thin_multiplier*10,setSeed=seed)))

lambda_3_3=sim_samples_3_3[,6:(5+sim$N)]
pi_3_3=sim_samples_3_3[,(6+sim$N):(5+2*sim$N)]
y_3_3=matrix(nrow=dim(lambda_3_3)[1],ncol=sim$N)
for(j in 1:sim$N){
  y_3_3[,j]=sim$z[j]+rpois(dim(lambda_3_3)[1],lambda_3_3[,j]*(1-pi_3_3[,j]))
}

class_3=ggplot(data.frame(x=sim$y,m=apply(y_3_3,2,median)))+
  geom_abline(slope=1,intercept=0,colour="#22211d")+xlab(expression('True Count ('*y[s]*')'))+
  geom_point(mapping=aes(x=x,y=m),colour=vp[9])+
  labs(
    y=expression('Median Predicted Count'),
    x=expression('True Count ('*y[s]*')'),
    title='No Covariates'
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.title = element_text(size= 14, hjust=0.01, color = "#4e4d47", margin = margin(b = 0.1, t = 0.4, l = 2, unit = "cm")),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA)
  )

pdf(file='sim_class.pdf',width=9,height=2.5)
multiplot(class_1,class_3,class_2,cols=3)
dev.off()

