library(MASS)
library(ggplot2)
mu <- c(1,1)    
Sigma <- matrix(c(.125, .01, .125, .125), 2)  # Covariance matrix
bivn <- mvrnorm(5000, mu = mu, Sigma = Sigma )
data<-data.frame(x=bivn[,1],y=bivn[,2],vowel='BACK_DIP',v='/o/')
mu <- c(1,0)
Sigma <- matrix(c(.125, .01, .125, .125), 2)  # Covariance matrix
bivn <- mvrnorm(5000, mu = mu, Sigma = Sigma )
data<-rbind(data,data.frame(x=bivn[,1],y=bivn[,2],vowel='BACK_MONO',v='/o/'))
mu <- c(0,0)
Sigma <- matrix(c(.125, .01, .125, .125), 2)  # Covariance matrix
bivn <- mvrnorm(5000, mu = mu, Sigma = Sigma )
data<-rbind(data,data.frame(x=bivn[,1],y=bivn[,2],vowel='FRONT_MONO',v='/o/'))
mu <- c(0,1)
Sigma <- matrix(c(.125, .01, .125, .125), 2)  # Covariance matrix
bivn <- mvrnorm(5000, mu = mu, Sigma = Sigma )
data<-rbind(data,data.frame(x=bivn[,1],y=bivn[,2],vowel='FRONT_DIP',v='/o/'))
mu <- c(mean(data$x),mean(data$y))
Sigma <- matrix(c(.5, .01, .5, .5), 2)  # Covariance matrix
bivn <- mvrnorm(5000, mu = mu, Sigma = Sigma )
data<-rbind(data,data.frame(x=bivn[,1],y=bivn[,2],vowel='u',v='/u/'))


mu <- c(mean(data$x),mean(data$y))
Sigma <- matrix(c(.5, .01, .5, .5), 2)  # Covariance matrix
bivn <- mvrnorm(5000, mu = mu, Sigma = Sigma )
categories<-rbind(data[data$vowel=='u',],data.frame(x=bivn[,1],y=bivn[,2],vowel='ALL',v='/o/'))


ggplot(data,aes(x=x,y=y,fill=vowel))+stat_ellipse(geom='polygon',alpha=0.35)+stat_ellipse(data=categories,aes(x=x,y=y),linetype='dotted')+theme_bw()+scale_fill_brewer(palette='Set1')+theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.position='none',panel.border=element_rect(size=1.5,color='black'),strip.background = element_blank(),strip.text.x = element_blank(),panel.grid.major = element_line(colour = "#808080"))+xlab('-F2')+ylab('Euclidean distance')+facet_wrap(~v)+ggsave('/Users/pplsuser/Desktop/PWL_submission/predictions.pdf',width=6,height=3,device=cairo_pdf)
