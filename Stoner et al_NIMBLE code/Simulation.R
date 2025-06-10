#########################################
##   Correcting for Under-Reporting    ##
##            Simulation               ##
#########################################

set.seed(seed)

sim=list()
sim$N=100 # Number of observations.

# Set-up adjacency structure for ICAR model.
sim$A=adjacency.matrix(10)
sim$n_adj=rowSums(sim$A)
sim$D=diag(sim$n_adj)
sim$adj=as.carAdjacency(sim$A)$adj

tau=4 # Spatial effect precision parameter.

# Note: the function to simulate ICAR fields calls eigen
# which we found produced numerically different results
# on different machines, so to reproduce the plots in the
# paper we include our simulated phi.
sim$phi <- ricar_simple(tau*(sim$D-sim$A)) 
sim$phi <- simulated_phi # Overwrite.

# True covariates.
sim$x=runif(sim$N,-1,1)
sim$w=runif(sim$N,-1,1) # Uniform(-1,1) has variance 1/3

# Under-reporting covariate noise.
sim$gamma=rnorm(sim$N,0,1)
sim$gamma=sim$gamma-mean(sim$gamma)

# Imperfect under-reporting covariates.
sim$v=matrix(nrow=sim$N,ncol=6)
sim$v[,1]=sim$w
sim$v[,2]=0.8*sim$w + sqrt((1-0.8^2)/3)*sim$gamma 
sim$v[,3]=0.6*sim$w + sqrt((1-0.6^2)/3)*sim$gamma
sim$v[,4]=0.4*sim$w + sqrt((1-0.4^2)/3)*sim$gamma
sim$v[,5]=0.2*sim$w + sqrt((1-0.2^2)/3)*sim$gamma
sim$v[,6]=sqrt(1/3)*sim$gamma

cor(sim$v) # Check correlation is as desired.

a0=4
a1=1
b0=0
b1=2

sim$lambda=exp(a0+a1*sim$x+sim$phi) # Poisson means.

sim$y=rpois(sim$N,sim$lambda) # True counts.

sim$pi=expit(b0+b1*sim$v[,1]) # Reporting probabilities.

sim$z=rbinom(sim$N,sim$y,sim$pi) # Observed counts.

# Figure 1:
s1=ggplot(data=data.frame(x=sim$x,y=sim$y,z=sim$z,c=sim$v[,1],p=sim$pi),
          mapping=aes(x=x,y=y))+geom_point(colour=vp[7])+
  labs(
    y=expression('True Count ('*y[s]*')'),
    x=expression('Process Covariate ('*x[s]*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
    )
s2=ggplot(data=data.frame(x=sim$x,y=sim$y,z=sim$z,c=sim$v[,1],p=sim$pi),
          mapping=aes(x=c,y=y))+geom_point(colour=vp[9])+
  labs(
    y=expression('True Count ('*y[s]*')'),
    x=expression('Under-Reporting Covariate ('*w['s']*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )
s3=ggplot(data=data.frame(x=sim$x,y=sim$y,z=sim$z,c=sim$v[,1],p=sim$pi),
          mapping=aes(x=x,y=z))+geom_point(colour=vp[7])+
  labs(
    y=expression('Observed Count ('*z[s]*')'),
    x=expression('Process Covariate ('*x[s]*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )
s4=ggplot(data=data.frame(x=sim$x,y=sim$y,z=sim$z,c=sim$v[,1],p=sim$pi),
          mapping=aes(x=c,y=z))+geom_point(colour=vp[11])+
  labs(
    y=expression('Observed Count ('*z[s]*')'),
    x=expression('Under-Reporting Covariate ('*w['s']*')')
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA)
  )

pdf(file='sim_values.pdf',height=2.5,width=9)
multiplot(s1,s2,s4,cols=3)
dev.off()

# Spatial effect plot.

ggplot(data.frame(x=sort(rep(1:10,10)),y=rep(1:10,10),phi=sim$phi))+geom_tile(aes(x=x,y=y,fill=phi))+ggtitle('Spatial Effect')+
  guides(fill=guide_legend(title=expression(phi[s])))+scale_fill_viridis()
ggsave('sim_spatial.pdf',device='pdf',height=3.5,width=4.5)

# Under-reporting covariates plot.

s6=ggplot(data=data.frame(c=sim$v[,1],p=logit(sim$pi)),
          mapping=aes(x=c,y=p))+geom_point(colour=vp[11])+
  labs(
    y=expression('Reporting Probability ('*pi[s]*')'),
    x=expression('Under-Reporting Covariate ('*w['s']*')'),
    title='Correlation 1'
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.title = element_text(size= 14, hjust=0.01, color = "#4e4d47", margin = margin(b = 0.1, t = 0.4, l = 2, unit = "cm")),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA)
  )
s7=ggplot(data=data.frame(c=sim$v[,2],p=logit(sim$pi)),
          mapping=aes(x=c,y=p))+geom_point(colour=vp[10])+
  labs(
    y=expression('Reporting Probability ('*pi[s]*')'),
    x=expression('Under-Reporting Covariate ('*v['s,2']*')'),
    title='Correlation 0.8'
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.title = element_text(size= 14, hjust=0.01, color = "#4e4d47", margin = margin(b = 0.1, t = 0.4, l = 2, unit = "cm")),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA)
  )
s8=ggplot(data=data.frame(c=sim$v[,3],p=logit(sim$pi)),
          mapping=aes(x=c,y=p))+geom_point(colour=vp[9])+
  labs(
    y=expression('Reporting Probability ('*pi[s]*')'),
    x=expression('Under-Reporting Covariate ('*v['s,3']*')'),
    title='Correlation 0.6'
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.title = element_text(size= 14, hjust=0.01, color = "#4e4d47", margin = margin(b = 0.1, t = 0.4, l = 2, unit = "cm")),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA)
  )
s9=ggplot(data=data.frame(c=sim$v[,4],p=logit(sim$pi)),
          mapping=aes(x=c,y=p))+geom_point(colour=vp[8])+
  labs(
    y=expression('Reporting Probability ('*pi[s]*')'),
    x=expression('Under-Reporting Covariate ('*v['s,4']*')'),
    title='Correlation 0.4'
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.title = element_text(size= 14, hjust=0.01, color = "#4e4d47", margin = margin(b = 0.1, t = 0.4, l = 2, unit = "cm")),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA)
  )
s10=ggplot(data=data.frame(c=sim$v[,5],p=logit(sim$pi)),
          mapping=aes(x=c,y=p))+geom_point(colour=vp[7])+
  labs(
    y=expression('Reporting Probability ('*pi[s]*')'),
    x=expression('Under-Reporting Covariate ('*v['s,5']*')'),
    title='Correlation 0.2'
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.title = element_text(size= 14, hjust=0.01, color = "#4e4d47", margin = margin(b = 0.1, t = 0.4, l = 2, unit = "cm")),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA)
  )
s11=ggplot(data=data.frame(c=sim$v[,6],p=logit(sim$pi)),
          mapping=aes(x=c,y=p))+geom_point(colour=vp[6])+
  labs(
    y=expression('Reporting Probability ('*pi[s]*')'),
    x=expression('Under-Reporting Covariate ('*v['s,6']*')'),
    title='Correlation 0'
  )+
  theme(
    text = element_text(color = "#22211d"),
    plot.title = element_text(size= 14, hjust=0.01, color = "#4e4d47", margin = margin(b = 0.1, t = 0.4, l = 2, unit = "cm")),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA)
  )

pdf(file='sim_covariates.pdf',height=4.5,width=9)
multiplot(s6,s9,s7,s10,s8,s11,cols=3)
dev.off()


