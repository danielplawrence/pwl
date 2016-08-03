Sys.setlocale("LC_ALL", "en_US.UTF-8") 
library(ggplot2)
library(plyr)
setwd('/Users/pplsuser/Desktop/Perception_data')
data<-read.csv('all_perception_data.csv')
names(data)<-c('N','Version','Date','Trial','Selected','Image_gender','Image_age','Image_class','Image_local','RT','TL','TR','ID','Trial_age','Trial_loc','Trial_class','Test','Sound','Speaker','Vowel','Variant','Image_selected','Image_age_rating','Image_edu_rating','Image_york_rating','Image_broad_rating','Image_posh_rating','Comments','Part_ID','Part_DOB','Part_gender','Part_edu','Parent_edu','Part_occ','Parent_occ','York_parents')
data$Part_DOB<-as.character(data$Part_DOB)
data$variable<-paste(data$Vowel,'_',data$Variant,sep="")
yob<-function(date){
	return(unlist(strsplit(date,'/'))[3])
}

data$YOB<-as.integer(unlist(lapply(data$Part_DOB,yob)))
get_age_2<-function(YOB){
return(ifelse(YOB<=1960,"O",ifelse(YOB<=1980,"M","Y")))
}
data$Age<-as.factor(sapply(data$YOB,'get_age_2'))
data$Edu_bin<-ifelse(data$Part_edu==1|data$Part_edu==2,'LW',ifelse(data$Part_edu==3|data$Part_edu==4,"W","M"))data<-data[!is.na(data$variable),]
bio<- read.csv('/Users/pplsuser/FAVE/FAVE-align-thesis/Biodata_summary_social.csv')
#############################################################
#We need to add MR1 MR2 and MR3 to the perception data
head(pdata)
#Try by unique combination of gender, yob, and educ
new_ID<-function(Part_ID){
	d<-data[data$Part_ID==Part_ID,]
	gen<-d$Part_gender[1]
	yob<-d$YOB[1]
	edu<-d$Part_edu[1]
	pedu<-d$Parent_edu[1]
	york<-d$York_parents[1]
return(subset(bio,Gender==gen&YOB==yob&Education==edu&Parents_Education==pedu&York==york)$New_ID[1])
}
levels(data$Part_ID)
new_ID(levels(data$Part_ID)[1])
data$New_ID<-sapply(data$Part_ID,'new_ID')
data[data$Part_ID=='Haze',]$New_ID<-'TO_M_1989_2'
#Great! -- now we pull out the social dimensions
get_dims<-function(Pseudo,dim){
	n<-length(bio[bio$New_ID==Pseudo,dim])
	if (n>0){
	return(as.numeric(bio[bio$New_ID==Pseudo,dim][1]))}
	else {return('NA')}
}
data$New_ID<-as.character(data$New_ID)
bio$New_ID<-as.character(bio$New_ID)
data$dim1<-sapply(data$New_ID,'get_dims','MR1')
data$dim2<-sapply(data$New_ID,'get_dims','MR2')
data$dim3<-sapply(data$New_ID,'get_dims','MR3')
data[data$YOB==1989,]$dim1<-rep(bio[bio$YOB==1989,]$MR1,384)
data[data$YOB==1989,]$dim2<-rep(bio[bio$YOB==1989,]$MR2,384)
data[data$YOB==1989,]$dim3<-rep(bio[bio$YOB==1989,]$MR3,384)
###############################################################
###############################################################
source(url('http://www.ling.upenn.edu/~joseff/scripts/recontrast.R'))
library(lme4)
library(lmerTest)
#################SOCIAL CLASS
#################OW
owsoc<-data[data$Vowel=='OW'&data$Test=='soc',]
owsoc$Image_class<-relevel(owsoc$Image_class,'M')
owsoc$Sound<-factor(owsoc$Sound)
owsoc$Variant<-factor(owsoc$Variant)
owsoc<-recontrast(owsoc)
owsoc$YOB<-scale(owsoc$YOB)
owsoc$dim1<-scale(owsoc$dim1)
owsoc$dim2<-scale(owsoc$dim2)
owsoc$dim3<-scale(owsoc$dim3)
owsoc$Variant<-mapply(get_ipa,vowel='GOAT',variant=owsoc$Variant)
null<-glmer(Image_class~1+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
set<-glmer(Image_class~1+Version+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
base<-glmer(Image_class~Variant+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
age<-glmer(Image_class~Variant*YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(null,base,age,dim3int,dim33way)
dim3int<-glmer(Image_class~Variant+dim3+(Variant:dim3)+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)




dim1<-glmer(Image_class~Variant+YOB+dim1+(Variant:YOB)+(Variant:dim1)+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(age,dim1)

dim2<-glmer(Image_class~Variant+YOB+dim2+(Variant:YOB)+(Variant:dim2)+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(age,dim2)

dim3int<-glmer(Image_class~Variant+YOB+dim3+(Variant:YOB)+(Variant:dim3)+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(age,dim3int)

dim33way<-glmer(Image_class~Variant*YOB*dim3+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(dim3int,dim33way)#NOPE

#Could be stable effects?
dim3int_steadyage<-glmer(Image_class~Variant+YOB+dim3+(Variant:dim3)+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(dim3int_steadyage,dim3int)
dim3age_steadyint<-glmer(Image_class~Variant+YOB+dim3+(Variant:YOB)+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(dim3age_steadyint,dim3int_steadyage)
#Stable effect of age
#I'll check gender again
dim3int_steadyage_steadygen<-glmer(Image_class~Variant+YOB+dim3+Part_gender+(Variant:dim3)+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(dim3int_steadyage,dim3int_steadyage_steadygen)
#Stable effect of gender,check for interaction
dim3int_genderint_steadyage<-glmer(Image_class~Variant+YOB+dim3+Part_gender+(Variant:Part_gender)+(Variant:dim3)+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(dim3int_steadyage,dim3int_genderint_steadyage)
#Marginal interaction with gender
#Test the three way
dim3int_gender3int_steadyage<-glmer(Image_class~Variant+YOB+Variant*dim3*Part_gender+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa',optCtrl=list(maxfun=1e5)),nAGQ=0)
anova(dim3int_gender3int_steadyage,dim3int_genderint_steadyage)
#Prefer the simple model -- interaction of variant and dim3 and gender with YOB as control
owmod<-dim3int_steadyage
#Let's plot the scores
owsoc$Variant<-as.factor(owsoc$Variant)
newdata<-expand.grid(YOB=min(owsoc$YOB):max(owsoc$YOB),Variant=levels(owsoc$Variant),dim3=-3:3)
newdata$Prob<-predict(owmod,newdata=newdata,re.form=NA,type='response')
ggplot(newdata,aes(x=dim3,y=Prob))+geom_smooth(aes(color=Variant),se=F)+theme_bw()+scale_color_brewer(palette='Set1')+facet_wrap(~Part_gender)
tempdat<-owmod@frame
tempdat$Prob<-fitted(owmod)
ggplot(tempdat[tempdat$Variant%in%levels(tempdat$Variant)[c(1,8,2,7)],],aes(x=dim3,y=Prob,color=factor(Variant,levels=levels(tempdat$Variant)[c(1,7,8,2,1)])))+geom_point(alpha=0.05)+geom_smooth(data=newdata[newdata$Variant%in%levels(newdata$Variant)[c(1,8,2,7)],],aes(y=Prob))+theme_bw()+scale_color_brewer(palette='Set1')+guides(color=guide_legend(title=NULL))+theme_bw()+theme(panel.grid.major = element_line(colour = "#808080"),legend.text = element_text(size=15),plot.margin = unit(c(0,0,0,0), "cm"),axis.title=element_text(size=12.5),axis.text.x=element_text(size=9))+scale_shape("Vowel")+xlab('Mobility index')+ylab("Probability of a 'WC' selection")+geom_hline(yintercept=0.5,color='black',linetype='dashed')
+ggsave('/Users/pplsuser/Desktop/Second_Year_Report_copy/Presentation/ow_sens_soc_LEL.pdf',height=2.5,width=3,device=cairo_pdf)
+facet_wrap(~Part_gender)
ggplot(newdata,aes(x=YOB,color=Variant))+geom_point(alpha=0.05)+geom_smooth(data=newdata,aes(y=Prob))+theme_bw()+scale_color_brewer(palette='Set1')+facet_wrap(~Part_gender)
ggplot(owmod@frame,aes(x=Variant,y=fitted(owmod),linetype=Part_gender,color=Variant))+geom_jitter(alpha=0.5)+geom_boxplot(aes(x=Variant,y=predict(owmod,type='response')))+theme_bw()+scale_color_brewer(palette='Set1')


#Let's plot the scores
write.csv(newdata,'predictions_ow_genderdim3.csv')
save(dim3gen,file='ow_genderdim3mod.R')
#Plots for Murcia
tempdat<-owmod@frame
tempdat$Prob<-fitted(owmod)
ggplot(tempdat[tempdat$Variant%in%levels(tempdat$Variant)[c(1,8,2,7)],],aes(x=dim3,y=Prob,color=factor(Variant,levels=levels(tempdat$Variant)[c(1,7,8,2,1)])))+geom_point(alpha=0.05)+geom_smooth(data=newdata[newdata$Variant%in%levels(newdata$Variant)[c(1,8,2,7)],],aes(y=Prob))+theme_bw()+scale_color_brewer(palette='Set1')+guides(color=guide_legend(title=NULL))+theme_bw()+theme(panel.grid.major = element_line(colour = "#808080"),legend.text = element_text(size=15),plot.margin = unit(c(0,0,0,0), "cm"),axis.title=element_text(size=12.5),axis.text.x=element_text(size=9))+scale_shape("Vowel")+xlab('Mobility index')+ylab("Probability of a 'WC' selection")+geom_hline(yintercept=0.5,color='black',linetype='dashed')
ggplot(tempdat,aes(x=factor(Variant,levels=levels(tempdat$Variant)[c(5,8,6,7,2,1,3,4)]),y=Prob))+geom_violin()+geom_jitter(aes(group=New_ID),alpha=0.005)+stat_summary(fun.y=mean,fun.ymin=min,fun.ymax=max, geom="pointrange", color="blue",size=0.25)+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5), axis.text.y = element_blank())+xlab('Variant')+ylab('')+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+ggtitle('/o/')
+ggsave('/Users/pplsuser/Desktop/Murcia/o_fronting_perception_class.pdf',width=2.25,height=2.5,device=cairo_pdf)
ggplot(tempdat,aes(x=factor(Variant,levels=levels(tempdat$Variant)[c(5,8,6,7,2,1,3,4)]),y=Prob,group=New_ID))+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('')+ylab("Probability of a WC selection")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+ggtitle('/o/')+stat_summary(fun.y=mean,geom='point',color='black')+stat_summary(fun.data=mean_se,geom='line',color='blue')+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('')+ylab("P (WC)")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+ggtitle('/o/')+facet_wrap(~New_ID,ncol=8)+ggsave('/Users/pplsuser/Desktop/Murcia/o_ind.pdf',device=cairo_pdf,width=11,height=7.5)
ggplot(tempdat,aes(x=factor(Variant,levels=levels(tempdat$Variant)[c(5,8,6,7,2,1,3,4)]),y=Prob,group=New_ID))+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('')+ylab("Probability of a WC selection")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+ggtitle('/o/')+stat_summary(fun.y=mean,geom='point',color='black')+stat_summary(fun.data=mean_se,geom='line',color='blue')+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('')+ylab("P (WC)")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+ggtitle('/o/')+geom_rect(data = subset(tempdat,New_ID%in%c('JA_M_1997_3','MH_M_1958_2','SP_M_1997_1','SC_M_1959_5','HK_M_1998_1')),aes(fill = 'pink'),xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf,alpha = 0.003)+theme(legend.position='NULL')+geom_rect(data = subset(tempdat,New_ID%in%c('KI_F_1997_3','MA_M_1954_2','SA_M_1987_6','LY_F_1955_7','SC_M_1948_6','YB_F_1994_5')),aes(fill = ''),xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf,alpha = 0.003)+facet_wrap(~New_ID,ncol=8)
+ggsave('/Users/pplsuser/Desktop/Murcia/o_ind_highlight.pdf',device=cairo_pdf,width=10,height=7.5)
+ggsave('/Users/pplsuser/Desktop/Murcia/o_ind.pdf',device=cairo_pdf,width=9,height=7.5)
ggplot(tempdat,aes(x=factor(Variant,levels=levels(tempdat$Variant)[c(5,8,6,7,2,1,3,4)]),y=Prob,group=New_ID))+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('')+ylab("Probability of a WC selection")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+ggtitle('/o/')+stat_summary(fun.y=mean,geom='point',color='black')+stat_summary(fun.data=mean_se,geom='line',color='blue')+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5,color='black'))+xlab('')+ylab("P (WC)")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+theme(panel.border=element_rect(size=1.5),strip.background = element_blank(),strip.text.x=element_text(color='black'))+xlab('')+ggtitle('/o/')+facet_wrap(~New_ID,ncol=4)+ggsave('/Users/pplsuser/Desktop/Labphon Posters/o_ind.pdf',device=cairo_pdf,width=6.75,height=12)
ggplot(newdata,aes(x=dim3,y=Prob))+geom_smooth(aes(color=factor(Variant,levels=levels(newdata$Variant)[c(5,8,6,7,1,2,4,3)])))+theme_bw()+scale_color_brewer(palette='Set1')+guides(color=guide_legend(title=NULL))+theme_bw()+theme(panel.grid.major = element_line(colour = "#808080"))+ylab('P (WC)')+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+geom_hline(yintercept=0.5,linetype='dotted')+xlab('Mobility index')+theme(panel.border=element_rect(size=1.5,color='black'))+ggtitle('/o/')+ggsave('/Users/pplsuser/Desktop/Murcia/o_perception_mobility.pdf',width=4,height=4,device=cairo_pdf)

newdat<-subset(tempdat,New_ID%in%c('MH_M_1958_2','MA_M_1954_2','AR_F_1988_3'))
newdat$New_ID<-factor(newdat$New_ID,levels=c('AR_F_1988_3','MH_M_1958_2','MA_M_1954_2'))
levels(newdat$New_ID)<-c('Pattern 1','Pattern 2','Pattern 3')
ggplot(newdat,aes(x=factor(Variant,levels=levels(newdat$Variant)[c(5,8,6,7,2,1,3,4)]),y=Prob,group=New_ID))+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('')+ylab("Probability of a WC selection")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+ggtitle('/o/')+stat_summary(fun.y=mean,geom='point',color='black')+stat_summary(fun.data=mean_se,geom='line',color='blue')+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('')+ylab("P (WC)")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+ggtitle('/o/')+theme(legend.position='NULL')+facet_grid(.~New_ID)+geom_rect(data = subset(newdat,New_ID=='Pattern 2'),aes(fill = 'purple'),xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf,alpha = 0.003)+theme(panel.border=element_rect(size=1.5))+xlab('')+ylab("P (WC)")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+ggtitle('/o/')+facet_grid(.~New_ID)+geom_rect(data = subset(newdat,New_ID=='Pattern 3'),aes(fill = 'pink'),xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf,alpha = 0.003)+geom_rect(data = subset(newdat,New_ID=='Pattern 3'),fill = 'pink',xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf,alpha = 0.003)+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('')+ylab("P (WC)")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+theme(legend.position='NULL')+ggtitle('/o/')+facet_grid(.~New_ID)+theme(legend.position='NULL')+ggsave('/Users/pplsuser/Desktop/Murcia/o_ind_examples.pdf',device=cairo_pdf,height=4,width=5)
	)
[c(5,8,6,7,2,1,3,4)]

uwsoc$Variant<-factor(uwsoc$Variant)
uwnewdata<-expand.grid(Image_class=levels(uwsoc$Image_class),Variant=levels(uwsoc$Variant),dim3=min(uwsoc$dim3):max(uwsoc$dim3),YOB=min(uwsoc$YOB):max(uwsoc$YOB))
uwnewdata$Prob<-predict(dim3,newdata=uwnewdata,re.form=NA,type='response')
uwnewdata$Net<-cut(uwnewdata$dim3,breaks=4)
uwnewdata$Net<-as.factor(uwnewdata$Net)
uwmod<-dim3threeway
uwtempdat<-uwmod@frame
uwtempdat$Prob<-fitted(uwmod)
ggplot(uwtempdat,aes(x=factor(Variant,levels=levels(uwtempdat$Variant)[c(6,5,4,3,2,1)]),y=Prob))+geom_violin()+geom_jitter(aes(group=New_ID),alpha=0.005)+stat_summary(fun.y=mean,fun.ymax=max,fun.ymin=min, geom="pointrange", color="blue",size=0.25)+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('Variant')+ylab("Probability of a WC selection")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))
+ggtitle('/u/')+ggsave('/Users/pplsuser/Desktop/Murcia/u_fronting_perception_class.pdf',width=2.5,height=2.5,device=cairo_pdf)
ggplot(uwtempdat,aes(x=factor(Variant,levels=levels(uwtempdat$Variant)[c(6,5,4,3,2,1)]),y=Prob,group=New_ID))+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('Variant')+ylab("Probability of a WC selection")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+ggtitle('/u/')+stat_summary(fun.y=mean,geom='point',color='black')+stat_summary(fun.data=mean_se,geom='line',color='blue')+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('')+ylab("P (WC)")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+ggtitle('/u/')+geom_rect(data = subset(uwtempdat,New_ID%in%c('BE_F_1936_2','DO_F_1948_3','GR_M_1946_4','JE_F_1947_5','KE_M_1950_3','PS_F_1946_2','RP_F_1973_4','SL_F_1943_2','IH_F_1991_5')),aes(fill = 'skyblue'),xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf,alpha = 0.003)+geom_rect(data = subset(uwtempdat,New_ID%in%c('HE_F_1966_2','HK_M_1998_1','IL_M_1947_4','JA_M_1997_3','JC_M_1997_3','JC_M_1952_3','JO_M_1999_2','KB_F_1955_2','LI_F_1997_3','LY_F_1955_7','MC_F_1950_3','ME_M_1992_2','MH_M_1958_2','SC_M_1948_6','SC_M_1959_5','SP_M_1997_1','VL_F_1975_2')),aes(fill = 'purple'),xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf,alpha = 0.003)+facet_wrap(~New_ID,ncol=8)+theme(legend.position='NULL')+ggsave('/Users/pplsuser/Desktop/Murcia/u_ind_highlight3.pdf',device=cairo_pdf,width=9,height=7.5)
ggplot(uwtempdat,aes(x=factor(Variant,levels=levels(uwtempdat$Variant)[c(6,5,4,3,2,1)]),y=Prob,group=New_ID))+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('')+ylab("Probability of a WC selection")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+ggtitle('/o/')+stat_summary(fun.y=mean,geom='point',color='black')+stat_summary(fun.data=mean_se,geom='line',color='blue')+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5,color='black'))+xlab('')+ylab("P (WC)")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+theme(panel.border=element_rect(size=1.5),strip.background = element_blank(),strip.text.x=element_text(color='black'))+xlab('')+ggtitle('/u/')+facet_wrap(~New_ID,ncol=4)+ggsave('/Users/pplsuser/Desktop/Labphon Posters/u_ind.pdf',device=cairo_pdf,width=6.5,height=12)

+ggsave('/Users/pplsuser/Desktop/Murcia/u_fronting_perception_class.pdf',width=2.5,height=2.5,device=cairo_pdf)
uwnewdat<-subset(uwtempdat,New_ID%in%c('RK_F_1993_3','HE_F_1966_2','BE_F_1936_2'))
uwnewdat$New_ID<-factor(uwnewdat$New_ID,levels=c('RK_F_1993_3','HE_F_1966_2','BE_F_1936_2'))
levels(uwnewdat$New_ID)<-c('Pattern 1','Pattern 2','Pattern 3')

ggplot(uwnewdat,aes(x=factor(Variant,levels=levels(uwtempdat$Variant)[c(6,5,4,3,2,1)]),y=Prob,group=New_ID))+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('Variant')+ylab("Probability of a WC selection")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+ggtitle('/u/')+stat_summary(fun.y=mean,geom='point',color='black')+stat_summary(fun.data=mean_se,geom='line',color='blue')+geom_hline(yintercept=0.5,linetype='dotted')+geom_rect(data = subset(uwnewdat,New_ID=='Pattern 2'),aes(fill = 'purple'),xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf,alpha = 0.003)+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('')+ylab("P (WC)")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+ggtitle('/u/')+facet_grid(.~New_ID)+geom_rect(data = subset(uwnewdat,New_ID=='Pattern 3'),aes(fill = 'pink'),xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf,alpha = 0.003)+geom_rect(data = subset(uwnewdat,New_ID=='Pattern 3'),fill = 'pink',xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf,alpha = 0.003)+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('')+ylab("P (WC)")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+theme(legend.position='NULL')+ggtitle('/u/')+facet_grid(.~New_ID)+theme(legend.position='NULL')+ggsave('/Users/pplsuser/Desktop/Murcia/u_ind_examples.pdf',device=cairo_pdf,width=5,height=4)

#An experiment with clustering
library(Rmisc)
library(reshape2)
d<-summarySE(data=uwtempdat,measurevar='Prob',groupvars=c('New_ID','Variant'))
mat<-reshape(timevar='Variant',idvar='New_ID',data=d[,c(1,2,4)],direction='wide')
row.names(mat)<-mat$New_ID
dist <- dist(as.matrix(mat)) 
clust<-hclust(dist)
clusters<-cutree(clust,k=4)
clusters<-melt(clusters)
clusters$New_ID<-row.names(clusters)

fit <- cascadeKM(mat[,2:7], 1, 10, iter = 1000)
calinski.best <- as.numeric(which.max(fit$results[2,]))
calinski.best
clusters<-unlist(fit$partition[,calinski.best])
clusters<-melt(clusters)
clusters$New_ID<-row.names(clusters)


get_clust<-function(ID){
	return(clusters[clusters$New_ID==ID,]$value)
}
uwtempdat$Clust<-sapply(uwtempdat$New_ID,'get_clust')
uwtempdat$Clust<-as.factor(uwtempdat$Clust)
ggplot(uwtempdat,aes(x=factor(Variant,levels=levels(uwtempdat$Variant)[c(6,5,4,3,2,1)]),y=Prob,group=Clust))+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('')+ylab("Probability of a WC selection")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+ggtitle('/o/')+stat_summary(fun.y=mean,geom='point',color='black')+stat_summary(fun.data=mean_se,geom='line',color='black')+stat_summary(fun.y=mean,geom='point',color='black',alpha=0.05,aes(group=New_ID))+stat_summary(fun.data=mean_se,geom='line',alpha=0.1,aes(color=Clust,group=New_ID))+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5,color='black'))+xlab('')+ylab("P (WC)")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+theme(panel.border=element_rect(size=1.5),strip.background = element_blank(),strip.text.x=element_blank())+xlab('')+theme(legend.position='NULL')+theme(plot.margin=unit(c(0,0,0,0),'cm'))+ggtitle('/u/')+facet_wrap(~Clust,ncol=4)+scale_color_brewer(palette='Set1')+ggsave('/Users/pplsuser/Desktop/Labphon Posters/u_clust.pdf',device=cairo_pdf,width=1.5,height=1.5)
+ggsave('/Users/pplsuser/Desktop/Labphon Posters/u_ind.pdf',device=cairo_pdf,width=6.5,height=12)


BIC<-mclustBIC(mat[,2:7])
fit<-Mclust(mat[,2:7],BIC=BIC)
clusters<-melt(fit$classification)
clusters$New_ID<-row.names(clusters)


d<-summarySE(data=tempdat,measurevar='Prob',groupvars=c('New_ID','Variant'))
mat<-reshape(timevar='Variant',idvar='New_ID',data=d[,c(1,2,4)],direction='wide')
row.names(mat)<-mat$New_ID
dist <- dist(as.matrix(mat)) 
clust<-hclust(dist)
clusters<-cutree(clust,k=6)
clusters<-melt(clusters)
clusters$New_ID<-row.names(clusters)
get_clust<-function(ID){
	return(clusters[clusters$New_ID==ID,]$value[1])
}
tempdat$Clust<-sapply(tempdat$New_ID,'get_clust')
tempdat$Clust<-as.factor(tempdat$Clust)
ggplot(tempdat,aes(x=factor(Variant,levels=levels(tempdat$Variant)[c(5,8,6,7,2,1,3,4)]),y=Prob,group=Clust))+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('')+ylab("Probability of a WC selection")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+ggtitle('/o/')+stat_summary(fun.y=mean,geom='point',color='black')+stat_summary(fun.data=mean_se,geom='line',color='black')+stat_summary(fun.y=mean,geom='point',color='black',alpha=0.05,aes(group=New_ID))+stat_summary(fun.data=mean_se,geom='line',alpha=0.1,aes(color=Clust,group=New_ID))+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5,color='black'))+xlab('')+ylab("P (WC)")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+theme(panel.border=element_rect(size=1.5),strip.background = element_blank(),strip.text.x=element_blank())+xlab('')+theme(legend.position='NULL')+ggtitle('/o/')+facet_wrap(~Clust,ncol=4)+scale_color_brewer(palette='Set1')+theme(plot.margin= unit(c(0,0,0,0), "cm"))+ggsave('/Users/pplsuser/Desktop/Labphon Posters/o_clust.pdf',device=cairo_pdf,width=4,height=2.5)

library(mclust)
library(pvclust)
mat<-na.omit(mat)
fit <- pvclust(t(mat[,2:9]),method.dist="euclidean")
clusters<-pvpick(fit,alpha=.99)$clusters
n.obs <- sapply(clusters, length)
seq.max <- seq_len(max(n.obs))
clu <- (sapply(clusters, "[", i = seq.max))
clusters<-melt(clu)
names(clusters)<-c('id','value','New_ID')
clusters$New_ID<-as.character(clusters$New_ID)
clusters<-na.omit(clusters)

BIC<-mclustBIC(mat[,2:9])
fit<-Mclust(mat[,2:9],BIC=BIC)
clusters<-melt(fit$classification)
clusters$New_ID<-row.names(clusters)


#Trying on model parameters
mat<-ranef(dim3int)$New_ID
BIC<-mclustBIC(mat[,1:8])
fit<-Mclust(mat[,1:8],BIC=BIC)
clusters<-melt(fit$classification)
clusters$New_ID<-row.names(clusters)
fit <- pvclust(t(mat[,1:8]),method.dist="euclidean")

fit <- cascadeKM(mat[,2:9], 1, 10, iter = 1000)
calinski.best <- as.numeric(which.max(fit$results[2,]))
calinski.best
clusters<-unlist(fit$partition[,calinski.best])
clusters<-melt(clusters)
clusters$New_ID<-row.names(clusters)

wss <- (nrow(mat[,1:8])-1)*sum(apply(mat[,1:8,2,var))
  for (i in 2:15) wss[i] <- sum(kmeans(mydata,
                                       centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

clusters<-kmeans(mat[,2:9],centers=2)[1]$cluster
clusters<-melt(clusters)
clusters$New_ID<-row.names(clusters)


#Demonstration of possibilities for variability in indexical perception
owdemo<-owmod@frame
owdemo$Prob<-fitted(owmod)
owdemo<-owdemo[owdemo$New_ID=='MI_M_1961_7',]
owdemo2<-owdemo[owdemo$New_ID=='MI_M_1961_7',]
owdemo3<-owdemo[owdemo$New_ID=='MI_M_1961_7',]
owdemo4<-owdemo[owdemo$New_ID=='MI_M_1961_7',]
owdemo$New_ID<-'A'
owdemo2$New_ID<-'B'
owdemo3$New_ID<-'C'
owdemo4$New_ID<-'D'
owdemo<-rbind(owdemo,owdemo2,owdemo3,owdemo3,owdemo4)
ggplot(owdemo,aes(x=factor(Variant,levels=levels(owdemo$Variant)[c(5,8,6,7,2,1,3,4)]),y=Prob,group=New_ID))+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('Variant')+ylab("Probability of a WC selection")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+stat_summary(fun.y=mean,geom='line',color='blue')+stat_summary(fun.y=mean,geom='point',color='black')+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('')+ylab("P (WC)")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+facet_wrap(~New_ID,ncol=4)+ggsave('/Users/pplsuser/Desktop/Murcia/uniform.pdf',width=7,height=2,device=cairo_pdf)

owdemo<-owmod@frame
owdemo$Prob<-fitted(owmod)
owdemo<-tempdat[tempdat$New_ID=='MI_M_1961_7',]
owdemo2<-tempdat[tempdat$New_ID=='SC_M_1959_5',]
owdemo3<-tempdat[tempdat$New_ID=='PA_M_2000_1',]
owdemo4<-tempdat[tempdat$New_ID=='TO_M_1996_3',]
owdemo$New_ID<-'A'
owdemo2$New_ID<-'B'
owdemo3$New_ID<-'C'
owdemo4$New_ID<-'D'
owdemo<-rbind(owdemo,owdemo2,owdemo3,owdemo3,owdemo4)
ggplot(owdemo,aes(x=factor(Variant,levels=levels(owdemo$Variant)[c(5,8,6,7,2,1,3,4)]),y=Prob,group=New_ID))+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('Variant')+ylab("Probability of a WC selection")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+stat_summary(fun.data=mean_se,geom='line',color='blue')+stat_summary(fun.y=mean,geom='point')+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+ylab("P (WC)")+xlab('')+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+facet_wrap(~New_ID,ncol=4)+ggsave('/Users/pplsuser/Desktop/Murcia/quant_diffs.pdf',width=7,height=2,device=cairo_pdf)


owdemo<-owmod@frame
owdemo$Prob<-fitted(owmod)
owdemo<-tempdat[tempdat$New_ID=='MI_M_1961_7',]
owdemo2<-tempdat[tempdat$New_ID=='SC_M_1959_5',]
owdemo3<-tempdat[tempdat$New_ID=='PA_M_2000_1',]
owdemo4<-tempdat[tempdat$New_ID=='TO_M_1996_3',]
owdemo$New_ID<-'A'
owdemo2$New_ID<-'B'
owdemo3$New_ID<-'C'
owdemo4$New_ID<-'D'
owdemo<-rbind(owdemo,owdemo2,owdemo3,owdemo3,owdemo4)
ggplot(owdemo,aes(x=factor(Variant,levels=levels(owdemo$Variant)[c(5,8,6,7,2,1,3,4)]),y=Prob,group=New_ID))+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('')+ylab("P (WC)")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+stat_summary(fun.data=mean_se,geom='line',color='blue')+stat_summary(fun.y=mean,geom='point')+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('')+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+facet_wrap(~New_ID,ncol=4)+ggsave('/Users/pplsuser/Desktop/Murcia/qual_diffs.pdf',width=7,height=2,device=cairo_pdf)


owdemo<-owmod@frame
owdemo$Prob<-fitted(owmod)
owdemo<-tempdat[tempdat$New_ID=='MI_M_1961_7',]
owdemo2<-tempdat[tempdat$New_ID=='SP_M_1997_1',]
owdemo3<-tempdat[tempdat$New_ID=='MA_M_1954_2',]
owdemo4<-tempdat[tempdat$New_ID=='MI_F_1974_7',]
owdemo$New_ID<-'A'
owdemo2$New_ID<-'B'
owdemo3$New_ID<-'C'
owdemo4$New_ID<-'D'
owdemo4$Variant<-factor(owdemo4$Variant,levels=levels(owdemo4$Variant)[c(5,8,6,7,2,1,3,4)])
owdemo4$Variant<-as.character(owdemo4$Variant)
owdemo4$Variant[owdemo4$Variant=="o:"]<-'MIDMONO'
owdemo4$Variant[owdemo4$Variant=="ɘy"]<-"o:"
owdemo4$Variant[owdemo4$Variant=='MIDMONO']<-"ɘy"
owdemo4$Variant[owdemo4$Variant=="ɵ:"]<-'MIDMONO'
owdemo4$Variant[owdemo4$Variant=="ɘʊ"]<-"ɵ:"
owdemo4$Variant[owdemo4$Variant=='MIDMONO']<-"ɘʊ"
owdemo4$Variant[owdemo4$Variant=="ø:"]<-'MIDMONO'
owdemo4$Variant[owdemo4$Variant=="əʉ"]<-"ø:"
owdemo4$Variant[owdemo4$Variant=='MIDMONO']<-"əʉ"
owdemo4$Variant[owdemo4$Variant=="oʊ"]<-'MIDMONO'
owdemo4$Variant[owdemo4$Variant=="əʊ"]<-"oʊ"
owdemo4$Variant[owdemo4$Variant=='MIDMONO']<-"əʊ"
owdemo4$New_ID<-'D'
owdemo<-rbind(owdemo,owdemo2,owdemo3,owdemo3,owdemo4)
ggplot(owdemo,aes(x=factor(Variant,levels=levels(owdemo$Variant)[c(5,8,6,7,2,1,3,4)]),y=Prob,group=New_ID))+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('Variant')+ylab("Probability of a WC selection")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+stat_summary(fun.data=mean_se,geom='line',color='blue')+stat_summary(fun.y=mean,geom='point')+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('Variant')+ylab("Probability of a WC selection")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+facet_wrap(~New_ID,ncol=4)




#lmerTest gives me significant main effects for all 
#significant interactions for back dip, front dip, mid, mono, marginal on front mono and front onset
#################UW
uwsoc<-data[data$Vowel=='UW'&data$Test=='soc',]
uwsoc$Image_class<-relevel(uwsoc$Image_class,'M')
uwsoc$Sound<-factor(uwsoc$Sound)
uwsoc$Variant<-factor(uwsoc$Variant)
uwsoc<-recontrast(uwsoc)
uwsoc$YOB<-scale(uwsoc$YOB)
uwsoc$dim1<-scale(uwsoc$dim1)
uwsoc$dim2<-scale(uwsoc$dim2)
uwsoc$dim3<-scale(uwsoc$dim3)
uwsoc$Variant<-mapply(get_ipa,vowel='GOOSE',variant=uwsoc$Variant)
uwsoc$Variant<-factor(uwsoc$Variant)
base<-glmer(Image_class~Variant+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
age<-glmer(Image_class~Variant*YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(base,age)
#age is an improvement over baseline. now I'll add dim3 as it was sig last time
dim1<-glmer(Image_class~Variant*YOB+dim1+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
dim1int<-glmer(Image_class~Variant*YOB*dim1+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(base,age)
#just age
dim2<-glmer(Image_class~Variant*YOB*dim2+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
dim3<-glmer(Image_class~Variant*YOB*dim3+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
dim3noint<-glmer(Image_class~Variant*YOB+dim3+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
dim3agecontrol<-glmer(Image_class~Variant*dim3+YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
dim3twowayint<-glmer(Image_class~Variant+YOB+(Variant:YOB)+(Variant:dim3)+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
dim3threeway<-glmer(Image_class~Variant*YOB*dim3+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
dim3<-dim3threeway
dim3<-dim3agecontrol
#best model is with three-way interaction of dim3 and YOB
#plots¡
#Awkward to visualise so split into quartiles

uwnewdata<-expand.grid(Image_class=levels(uwsoc$Image_class),Variant=levels(uwsoc$Variant),dim3=-3:3,YOB=min(uwsoc$YOB):max(uwsoc$YOB))
uwnewdata$Prob<-predict(dim3,newdata=uwnewdata,re.form=NA,type='response')
uwnewdata$Net<-cut(uwnewdata$dim3,breaks=3)
uwnewdata$Net<-as.factor(uwnewdata$Net,levels=c('L'))
ggplot(uwnewdata,aes(x=YOB,y=Prob))+geom_smooth(aes(color=Variant))+theme_bw()+scale_color_brewer(palette='Set1')+facet_wrap(~Net)

ggplot(uwnewdata,aes(x=YOB,y=Prob))+geom_smooth(aes(color=factor(Variant,levels=levels(uwsoc$Variant)[c(3,6,2,1,5,4)]),se=F))+theme_bw()+scale_color_brewer(palette='Set1')+guides(color=guide_legend(title=NULL))+theme_bw()+theme(panel.grid.major = element_line(colour = "#808080"))+ylab('P (WC)')+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+geom_hline(yintercept=0.5,linetype='dotted')+xlab('Year of birth')+ggtitle('/u/')+theme(panel.border=element_rect(size=1.5,color='black'))+ggsave('/Users/pplsuser/Desktop/Murcia/u_perception_age.pdf',width=4,height=4,device=cairo_pdf)
ggplot(uwnewdata,aes(x=dim3,y=Prob))+geom_smooth(aes(color=factor(Variant,levels=levels(uwsoc$Variant)[c(3,6,2,1,5,4)]),se=F))+theme_bw()+scale_color_brewer(palette='Set1')+guides(color=guide_legend(title=NULL))+theme_bw()+theme(panel.grid.major = element_line(colour = "#808080"))+ylab('P (WC)')+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+geom_hline(yintercept=0.5,linetype='dotted')+xlab('Mobility index')+theme(legend.position='NULL')+theme(panel.border=element_rect(size=1.5,color='black'))+ggtitle('/u/')+ggsave('/Users/pplsuser/Desktop/Murcia/u_perception_mobility.pdf',width=3.25,height=4,device=cairo_pdf)




mod<-dim3@frame
mod$Prob<-fitted(dim3)
mod$YOB<-mod$YOB * attr(mod$YOB, 'scaled:scale') + attr(mod$YOB, 'scaled:center')
newdata$YOB<-newdata$YOB * attr(newdata$YOB, 'scaled:scale') + attr(newdata$YOB, 'scaled:center')
ggplot(mod,aes(x=YOB,y=Prob,color=Variant))+geom_point(alpha=0.05)+geom_smooth(data=uwnewdata,aes(y=Prob),se=F)+theme_bw()+scale_color_brewer(palette='Set1')+facet_wrap(~Net)

#################EY
eysoc<-data[data$Vowel=='EY'&data$Test=='soc',]
eysoc$Image_class<-relevel(eysoc$Image_class,'M')
eysoc$Sound<-factor(eysoc$Sound)
eysoc$Variant<-factor(eysoc$Variant)
eysoc<-recontrast(eysoc)
eysoc$YOB<-scale(eysoc$YOB)
eysoc$dim1<-scale(eysoc$dim1)
eysoc$dim2<-scale(eysoc$dim2)
eysoc$dim3<-scale(eysoc$dim3)
base<-glmer(Image_class~Variant+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=eysoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
age<-glmer(Image_class~Variant*YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=eysoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
stableage<-glmer(Image_class~Variant+YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=eysoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
stableagestabledim<-glmer(Image_class~Variant+YOB+dim2+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=eysoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
dim2int<-glmer(Image_class~Variant+YOB+(Variant:dim2)+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=eysoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
stableagedim2int<-glmer(Image_class~Variant*dim2+YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=eysoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
jsutdim2<-glmer(Image_class~Variant*dim2+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=eysoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(dim2int,dim2intgen)
dim2int<-stableagedim2int
#Best model is with stable YOB and interaction of variant and dim2
#It's weird that 2 comes out here
dim3int<-glmer(Image_class~Variant+YOB+(Variant:dim3)+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=eysoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
#Basically buy it...try gender
dim2intgen<-glmer(Image_class~Variant+YOB*Part_gender+(Variant:dim2)+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=eysoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
#Basically buy it...try gender
eynewdata<-expand.grid(Image_class=levels(eysoc$Image_class),Variant=levels(eysoc$Variant),dim2=-3:3,YOB=1935:2000)
eynewdata$Prob<-predict(dim2int,newdata=eynewdata,re.form=NA,type='response')
ggplot(eynewdata,aes(x=dim2,y=Prob,color=Variant))+geom_point()
ggplot(dim2int@frame,aes(x=dim2,y=fitted(dim2int)))+geom_point()
geom_smooth(aes(color=Variant),se=F)+theme_bw()+scale_color_brewer(palette='Set1')
+facet_wrap(~Net)
ind<-ranef(dim2int)$New_ID

get_year<-function(Part_ID){
return(strsplit(Part_ID,split='_')[[1]][3])}
ind$Part_ID<-as.character(row.names(ind))
ind<-melt(ind)
ind<-cbind(ind,ldply(lapply(ind$Part_ID,get_year),'data.frame'))
ind$data.frame=as.numeric(as.character(ind$data.frame))
ggplot(ind,aes(x=data.frame,y=value,color=variable))+geom_point(alpha=0.5)+geom_line(aes(y=eynewdata$Prob))
+geom_line(y=predict(dim2int))

ggplot(dim2int@frame,aes(x=dim2,y=fitted(dim2int),color=Variant))+geom_point()+geom_smooth(data=eynewdata,aes(y=Prob))+theme_bw()+scale_color_brewer(palette='Set1')
ggplot(dim2int@frame,aes(x=YOB,y=fitted(dim2int),color=Variant))+geom_point()+geom_smooth(data=eynewdata,aes(y=Prob))+theme_bw()+scale_color_brewer(palette='Set1')
#So it looks like the effect of social class on ey perception primarily targets the diphthong
ggplot(dim2int@frame,aes(x=dim2,y=fitted(dim2int),color=Variant))+geom_jitter()+geom_line(color='black',aes(group=New_ID),alpha=0.25)+geom_smooth(data=eynewdata,aes(y=Prob))+theme_bw()+scale_color_brewer(palette='Set1')
+geom_smooth(data=eynewdata,aes(y=Prob))+theme_bw()+scale_color_brewer(palette='Set1')+geom_line(color='black',linetype='dotted')+geom_smooth()
###########Urban/rural
#################OW
owloc<-data[data$Vowel=='OW'&data$Test=='loc',]
owloc$Image_local<-relevel(owloc$Image_local,'N')
owloc$Sound<-factor(owloc$Sound)
owloc$Variant<-factor(owloc$Variant)
owloc<-recontrast(owloc)
base<-glmer(Image_local~Variant+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
age<-glmer(Image_local~Variant*YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(base,age)
#no improvement, try dim1
dim2<-glmer(Image_local~Variant*dim1+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(base,dim2)
dim2<-glmer(Image_local~Variant*dim2+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(base,dim2)
dim3<-glmer(Image_local~Variant*dim3+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(base,dim3)
#no effect of anything!
#Check the variant!
null<-glmer(Image_local~1+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
base<-glmer(Image_local~Variant+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(null,base)
#Prefer base
yob<-glmer(Image_local~Variant*YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(yob,base)
#Okay!
ownewdata<-expand.grid(Image_class=levels(owloc$Image_class),Variant=levels(owloc$Variant),dim2=-3:3,YOB=1935:2000)
ownewdata$Prob<-predict(base,newdata=ownewdata,re.form=NA,type='response')
ggplot(yob@frame,aes(x=YOB,y=fitted(base),color=Variant))+geom_point()+geom_smooth(data=ownewdata,aes(y=Prob))+theme_bw()+scale_color_brewer(palette='Set1')
#################UW
uwloc<-data[data$Vowel=='UW'&data$Test=='loc',]
uwloc$Image_local<-relevel(uwloc$Image_local,'N')
uwloc$Sound<-factor(uwloc$Sound)
uwloc$Variant<-factor(uwloc$Variant)
uwloc<-recontrast(uwloc)
null<-glmer(Image_local~1+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
#Null model is base for uw on the urban rural dimension
#################EY
eyloc<-data[data$Vowel=='EY'&data$Test=='loc',]
eyloc$Image_class<-relevel(eyloc$Image_class,'M')
eyloc$Sound<-factor(eyloc$Sound)
eyloc$Variant<-factor(eyloc$Variant)
eyloc<-recontrast(eysoc)
base<-glmer(Image_local~Variant+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=eyloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
age<-glmer(Image_local~Variant*YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=eyloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(base,age)
#no improvement, try dim1
dim1<-glmer(Image_local~Variant*dim1+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=eyloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(base,dim1)
dim2<-glmer(Image_local~Variant*dim2+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=eyloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(base,dim2)
dim3<-glmer(Image_local~Variant*dim3+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=eyloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(base,dim3)
###########Age
#################OW
owage<-data[data$Vowel=='OW'&data$Test=='age',]
owage$Image_age<-relevel(owage$Image_age,'Y')
owage$Sound<-factor(owage$Sound)
owage$Variant<-factor(owage$Variant)
owage<-recontrast(owage)
null<-glmer(Image_age~1+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owage,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
base<-glmer(Image_age~Variant+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owage,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
age<-glmer(Image_age~Variant*YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owage,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(null,base)
anova(age,base)
#no improvement, try dim1
dim1<-glmer(Image_age~Variant*dim1+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owage,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(base,dim1)
dim2<-glmer(Image_age~Variant*dim2+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owage,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(base,dim2)
dim3<-glmer(Image_age~Variant*dim3+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owage,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(base,dim3)
#Okay!
ownewdata<-expand.grid(Variant=levels(owage$Variant))
ownewdata$Prob<-predict(base,newdata=ownewdata,re.form=NA,type='response')
ggplot(base@frame,aes(x=Variant,y=fitted(base),color=Variant))+geom_boxplot()
+geom_smooth(data=ownewdata,aes(y=Prob))+theme_bw()+scale_color_brewer(palette='Set1')
#################UW
uwage<-data[data$Vowel=='UW'&data$Test=='age',]
uwage$Image_age<-relevel(uwage$Image_age,'Y')
uwage$Sound<-factor(uwage$Sound)
uwage$Variant<-factor(uwage$Variant)
uwage<-recontrast(uwage)
null<-glmer(Image_age~1+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwage,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
base<-glmer(Image_age~Variant+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwage,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
age<-glmer(Image_age~Variant*YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwage,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(null,base)
#################EY
eyage<-data[data$Vowel=='EY'&data$Test=='age',]
eyage$Image_age<-relevel(eyage$Image_age,'Y')
eyage$Sound<-factor(eyage$Sound)
eyage$Variant<-factor(eyage$Variant)
eyage<-recontrast(eyage)
null<-glmer(Image_age~1+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=eyage,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
base<-glmer(Image_age~Variant+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=eyage,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
age<-glmer(Image_age~Variant*YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=eyage,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(null,base)
#no improvement, try dim1
dim1<-glmer(Image_age~Variant*dim1+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=eyage,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(base,dim1)
dim2<-glmer(Image_age~Variant*dim2+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=eyage,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(base,dim2)
dim3<-glmer(Image_age~Variant*dim3+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=eyage,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(base,dim3)
dim2mod<-dim2
#dim2 -- could it be a stable effect?
dim2stable<-glmer(Image_age~Variant+dim2+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=eyage,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
anova(dim2,dim2stable)
#it's a stable effect -- more Yorkshire you are, the more biased toward an 'O' selection you are
eynewdata<-expand.grid(Image_class=levels(eysoc$Image_class),Variant=levels(eysoc$Variant),dim2=-3:3,YOB=1935:2000)
eynewdata$Prob<-predict(dim2,newdata=eynewdata,re.form=NA,type='response')
ggplot(eynewdata,aes(x=dim2,y=Prob,color=Variant))+geom_point()
ggplot(dim2mod@frame,aes(x=dim2,y=fitted(dim2mod),color=Variant))+geom_point()+geom_smooth()



get_ipa<-function(vowel,variant){
  target<-paste(vowel,variant,sep="_")
  GOAT_BACK_DIP='o\u028A'
  GOAT_MID_DIP='\u0259\u0289'
  GOAT_FRONT_DIP='\u0258y'
  GOAT_MID_ONSET='\u0259\u028A'
  GOAT_FRONT_ONSET='\u0258\u028A'
  GOAT_BACK_MONO='o:'
  GOAT_MID_MONO='\u0275:'
  GOAT_FRONT_MONO='\u00F8:'
  GOOSE_BACK='\u028Au'
  GOOSE_BACK_LOW='\u0264u'
  GOOSE_MID='\u0289u'
  GOOSE_MID_LOW='\u0258u'
  GOOSE_FRONT='\u026Au'
  GOOSE_FRONT_LOW='eu'
  FACE_MONO='\u025B'
  FACE_DIP='\u025B\u026A'
  ipa=get(target)
  return(ipa)
}

#21st June. I'm checking the models because I noticed that the 
#random intercepts were set to 0 -- this could be changing how my clustering comes out
library(Rmisc)
library(vegan)
library(mclust)
library(reshape2)

dim3int_steadyage<-glmer(Image_class~Variant+YOB+dim3+(Variant:dim3)+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
owmod<-dim3int_steadyage
ranef(owmod)
#This one was fine
tempdat<-owmod@frame
tempdat$Prob<-fitted(owmod)

#Trying on fitted values first but will probably want to use coeffs and exclude intercepts
d<-summarySE(data=tempdat,measurevar='Prob',groupvars=c('New_ID','Variant'))
mat<-reshape(timevar='Variant',idvar='New_ID',data=d[,c(1,2,4)],direction='wide')
row.names(mat)<-mat$New_ID

fit <- cascadeKM(mat[,2:9], 1, 10, iter = 1000)
calinski.best <- as.numeric(which.max(fit$results[2,]))
calinski.best
clusters<-unlist(fit$partition[,calinski.best])
clusters<-melt(clusters)
clusters$New_ID<-row.names(clusters)
get_clust<-function(ID){
	return(clusters[clusters$New_ID==ID,]$value[1])
}
tempdat$Clust<-sapply(tempdat$New_ID,'get_clust')
tempdat$Clust<-as.factor(tempdat$Clust)
ggplot(tempdat,aes(x=factor(Variant,levels=levels(tempdat$Variant)[c(5,8,6,7,2,1,3,4)]),y=Prob,group=Clust))+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('')+ylab("Probability of a WC selection")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+ggtitle('/o/')+stat_summary(fun.y=mean,geom='point',color='black')+stat_summary(fun.data=mean_se,geom='line',color='black')+stat_summary(fun.y=mean,geom='point',color='black',alpha=0.05,aes(group=New_ID))+stat_summary(fun.data=mean_se,geom='line',alpha=0.1,aes(color=Clust,group=New_ID))+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5,color='black'))+xlab('')+ylab("P (WC)")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+theme(panel.border=element_rect(size=1.5),strip.background = element_blank(),strip.text.x=element_text(color='black'))+xlab('')+theme(legend.position='NULL')+ggtitle('/o/')+facet_wrap(~Clust,ncol=4)+scale_color_brewer(palette='Set1')
+ggsave('/Users/pplsuser/Desktop/Labphon Posters/o_clust.pdf',device=cairo_pdf,width=6.85,height=5)

#Trying on random slopes

mat<-ranef(owmod)$New_ID
fit <- cascadeKM(mat[,1:8], 1, 10, iter = 1000)
calinski.best <- as.numeric(which.max(fit$results[2,]))
calinski.best
clusters<-unlist(fit$partition[,calinski.best])
clusters<-melt(clusters)
clusters$New_ID<-row.names(clusters)
tempdat$Clust<-sapply(tempdat$New_ID,'get_clust')
tempdat$Clust<-as.factor(tempdat$Clust)
ggplot(tempdat,aes(x=factor(Variant,levels=levels(tempdat$Variant)[c(5,8,6,7,2,1,3,4)]),y=Prob,group=Clust))+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('')+ylab("Probability of a WC selection")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+ggtitle('/o/')+stat_summary(fun.y=mean,geom='point',color='black')+stat_summary(fun.data=mean_se,geom='line',color='black')+stat_summary(fun.y=mean,geom='point',color='black',alpha=0.05,aes(group=New_ID))+stat_summary(fun.data=mean_se,geom='line',alpha=0.1,aes(color=Clust,group=New_ID))+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5,color='black'))+xlab('')+ylab("P (WC)")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+theme(panel.border=element_rect(size=1.5),strip.background = element_blank(),strip.text.x=element_text(color='black'))+xlab('')+theme(legend.position='NULL')+ggtitle('/o/')+facet_wrap(~Clust,ncol=4)+scale_color_brewer(palette='Set1')
#add main effects
owsoc$Variant<-factor(owsoc$Variant)
ndat<-expand.grid(YOB=mean(owsoc$YOB),dim3=mean(owsoc$dim3),Variant=levels(owsoc$Variant))
main<-predict(owmod,newdata=expand.grid(YOB=mean(owsoc$YOB),dim3=mean(owsoc$dim3),Variant=levels(owsoc$Variant)),re.form=NA,type='response')
ndat$Prob<-main
ndat2<-ndat 
ndat$Clust<-'1'
ndat2$Clust<-'2'
ndat<-rbind(ndat,ndat2)
ggplot(tempdat,aes(x=factor(Variant,levels=levels(tempdat$Variant)[c(5,8,6,7,2,1,3,4)]),y=Prob,group=Clust))+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('')+ylab("Probability of a WC selection")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+ggtitle('/o/')+stat_summary(fun.y=mean,geom='point',color='black')+stat_summary(fun.data=mean_se,geom='line',color='black')+stat_summary(fun.y=mean,geom='point',color='black',alpha=0.05,aes(group=New_ID))+stat_summary(fun.data=mean_se,geom='line',alpha=0.1,aes(color=Clust,group=New_ID))+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5,color='black'))+xlab('')+ylab("P (WC)")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+theme(panel.border=element_rect(size=1.5),strip.background = element_blank(),strip.text.x=element_text(color='black'))+xlab('')+theme(legend.position='NULL')+ggtitle('/o/')+facet_wrap(~Clust,ncol=4)+scale_color_brewer(palette='Set1')
#add main effects
#Plot coeffs maybe it will show the patterns better
mat$New_ID<-row.names(clusters)
mat$Clust<-sapply(mat$New_ID,'get_clust')
mat<-melt(mat,id.vars=c('New_ID','Clust'))
names(mat)<-c("New_ID" ,  "Clust"  ,  "Variant", "Prob")
mat$Clust<-factor(mat$Clust)
ggplot(mat,aes(x=factor(Variant,levels=levels(mat$Variant)[c(5,8,6,7,2,1,3,4,9)]),y=Prob,group=Clust))+geom_hline(yintercept=0,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+stat_summary(fun.y=mean,geom='point',color='black')+stat_summary(fun.data=mean_se,geom='line',color='black')+stat_summary(fun.y=mean,geom='point',color='black',alpha=0.05,aes(group=New_ID))+stat_summary(fun.data=mean_se,geom='line',alpha=0.1,aes(color=Clust,group=New_ID))+facet_wrap(~Clust)
xlab('')+ylab("Probability of a WC selection")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+ggtitle('/o/')+stat_summary(fun.y=mean,geom='point',color='black')+stat_summary(fun.data=meanse,geom='line',color='black')+stat_summary(fun.y=mean,geom='point',color='black',alpha=0.05,aes(group=New_ID))+stat_summary(fun.data=mean_se,geom='line',alpha=0.1,aes(color=Clust,group=New_ID))+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5,color='black'))+xlab('')+ylab("P (WC)")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+theme(panel.border=element_rect(size=1.5),strip.background = element_blank(),strip.text.x=element_text(color='black'))+xlab('')+theme(legend.position='NULL')+ggtitle('/o/')+facet_wrap(~Clust,ncol=4)+scale_color_brewer(palette='Set1')
#Interesting -- I wonder if I could add the main effect
main<-melt(fixef(owmod))
main$Variant<-row.names(main)
main<-main[1:8,]
main$New_ID<-'Main_effect'
main$Clust<-'1'
main2<-main
main2$Clust<-'2'
main<-rbind(main,main2)
names(main)<-c('Prob','Variant','New_ID','Clust')
main<-main[,c('New_ID','Clust','Variant','Prob')]
ggplot(mat,aes(x=factor(Variant,levels=levels(mat$Variant)[c(5,8,6,7,2,1,3,4,9)]),y=Prob,group=Clust))+geom_hline(yintercept=0,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+stat_summary(fun.y=mean,geom='point',color='black')+stat_summary(fun.data=mean_se,geom='line',color='black',linetype='dotted')+stat_summary(fun.y=mean,geom='point',color='black',alpha=0.05,aes(group=New_ID))+stat_summary(fun.data=mean_se,geom='line',alpha=0.1,aes(color=Clust,group=New_ID))+stat_summary(fun.data=mean_se,geom='line',color='black',data=main)+facet_wrap(~Clust)
#Pretty cool
#Now /u/ 

uwsoc<-data[data$Vowel=='UW'&data$Test=='soc',]
uwsoc$Image_class<-relevel(uwsoc$Image_class,'M')
uwsoc$Sound<-factor(uwsoc$Sound)
uwsoc$Variant<-factor(uwsoc$Variant)

uwsoc$Variant<-mapply(get_ipa,vowel='GOOSE',variant=uwsoc$Variant)
uwsoc$Variant<-factor(uwsoc$Variant)
uwsoc$YOB<-scale(uwsoc$YOB,center=TRUE)
dim3threeway<-glmer(Image_class~Variant*YOB*dim3+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
uwmod<-dim3threeway

uwtempdat<-uwmod@frame
uwtempdat$Prob<-fitted(uwmod)

uwsoc$Variant<-factor(uwsoc$Variant)
uwnewdata<-expand.grid(Image_class=levels(uwsoc$Image_class),Variant=levels(uwsoc$Variant),dim3=min(uwsoc$dim3):max(uwsoc$dim3),YOB=min(uwsoc$YOB):max(uwsoc$YOB))
uwnewdata$Prob<-predict(dim3,newdata=uwnewdata,re.form=NA,type='response')
uwnewdata$Net<-cut(uwnewdata$dim3,breaks=4)
uwnewdata$Net<-as.factor(uwnewdata$Net)
uwmod<-dim3threeway
uwtempdat<-uwmod@frame
uwtempdat$Prob<-fitted(uwmod)


mat<-ranef(uwmod)$New_ID
fit <- cascadeKM(mat[,1:6], 1, 10, iter = 1000)
calinski.best <- as.numeric(which.max(fit$results[2,]))
calinski.best
clusters<-unlist(fit$partition[,calinski.best])
clusters<-melt(clusters)
clusters$New_ID<-row.names(clusters)
uwtempdat$Clust<-sapply(uwtempdat$New_ID,'get_clust')
uwtempdat$Clust<-as.factor(uwtempdat$Clust)
ggplot(uwtempdat,aes(x=factor(Variant,levels=levels(uwtempdat$Variant)[c(6,5,4,3,2,1)]),y=Prob,group=Clust))+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+xlab('')+ylab("Probability of a WC selection")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+ggtitle('/o/')+stat_summary(fun.y=mean,geom='point',color='black')+stat_summary(fun.data=mean_se,geom='line',color='black')+stat_summary(fun.y=mean,geom='point',color='black',alpha=0.05,aes(group=New_ID))+stat_summary(fun.data=mean_se,geom='line',alpha=0.1,aes(color=Clust,group=New_ID))+geom_hline(yintercept=0.5,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5,color='black'))+xlab('')+ylab("P (WC)")+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+theme(panel.border=element_rect(size=1.5),strip.background = element_blank(),strip.text.x=element_text(color='black'))+xlab('')+theme(legend.position='NULL')+ggtitle('/u/')+facet_wrap(~Clust,ncol=4)+scale_color_brewer(palette='Set1')+ggsave('/Users/pplsuser/Desktop/Labphon Posters/u_clust.pdf',device=cairo_pdf,width=10,height=5)

mat$New_ID<-row.names(clusters)
mat$Clust<-sapply(mat$New_ID,'get_clust')
mat<-melt(mat,id.vars=c('New_ID','Clust'))
names(mat)<-c("New_ID" ,  "Clust"  ,  "Variant", "Prob")
mat$Clust<-factor(mat$Clust)


#Interesting -- I wonder if I could add the main effect
main<-melt(fixef(uwmod))
main$Variant<-row.names(main)
main<-main[1:8,]
main$New_ID<-'Main_effect'
main$Clust<-'1'
main2<-main
main2$Clust<-'2'
main<-rbind(main,main2)
names(main)<-c('Prob','Variant','New_ID','Clust')
main<-main[,c('New_ID','Clust','Variant','Prob')]
ggplot(mat,aes(x=factor(Variant,levels=levels(mat$Variant)[c(6,5,4,3,2,1)]),y=Prob,group=Clust))+geom_hline(yintercept=0,linetype='dotted')+theme_bw()+theme(panel.border=element_rect(size=1.5))+stat_summary(fun.y=mean,geom='point',color='black')+stat_summary(fun.data=mean_se,geom='line',color='black',linetype='dotted')+stat_summary(fun.y=mean,geom='point',color='black',alpha=0.05,aes(group=New_ID))+stat_summary(fun.data=mean_se,geom='line',alpha=0.1,aes(color=Clust,group=New_ID))+stat_summary(fun.data=mean_se,geom='line',color='black',data=main)+facet_wrap(~Clust)
#Pretty cool


#Now I want to plot confidence intervals for the interactions
owsoc$Variant<-factor(owsoc$Variant)
owmod<-glmer(Image_class~Variant+YOB+dim3+(Variant:dim3)+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
newdata<-expand.grid(YOB=min(owsoc$YOB):max(owsoc$YOB),Variant=levels(owsoc$Variant),dim3=-3:3)
bootfit <- bootMer(owmod, FUN=function(x)predict(x, newdata, re.form=NA,type='response'),nsim=1000)
newdata$lci<-apply(bootfit$t, 2, quantile,0.025)
newdata$uci<-apply(bootfit$t, 2, quantile,0.975)
newdata$med<-apply(bootfit$t, 2, quantile,0.5)
newdata$se<-apply(bootfit$t, 2, sd)
save(newdata,file='owboot.R')
save(owmod,file='owmod.R')
save(bootfit,file='bootfit.R')
newdata$Prob<-predict(owmod,newdata=newdata,type='response',re.form=NA)
ggplot(newdata,aes(x=dim3,y=Prob))+geom_line(aes(group=Variant),size=1.5,linetype='dotdash')+geom_ribbon(aes(ymax=Prob+(se),ymin=Prob-(se),fill=factor(Variant,levels=levels(Variant)[c(5,8,6,7,2,1,3,4)])),alpha=0.25)+theme_bw()+scale_fill_brewer(palette='Set1')+guides(fill=guide_legend(title=NULL))+theme_bw()+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+theme(panel.grid.major = element_line(colour = "#808080"))+ylab('P (WC)')+geom_hline(yintercept=0.5,linetype='dotted')+xlab('Mobility index')+ggtitle('/o/')+theme(panel.border=element_rect(size=1.5,color='black'))+ggsave('/Users/pplsuser/Desktop/Labphon_Posters/o_perception_dim3_sd.pdf',width=3.3,height=3.3,device=cairo_pdf)
se<-function(x,n){
return(sd(x)/(sqrt(n)))
}

uwmod<-glmer(Image_class~Variant*YOB*dim3+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
uwnewdata<-expand.grid(YOB=min(uwsoc$YOB):max(uwsoc$YOB),Variant=levels(uwsoc$Variant),dim3=-3:3)
uwbootfit <- bootMer(uwmod, FUN=function(x)predict(x, uwnewdata, re.form=NA,type='response'),nsim=1000,verbose=T,parallel='multicore',ncpus=4)
uwnewdata$lci<-apply(uwbootfit$t, 2, quantile,0.025)
uwnewdata$uci<-apply(uwbootfit$t, 2, quantile,0.975)
uwnewdata$med<-apply(uwbootfit$t, 2, quantile,0.5)
uwnewdata$se<-apply(uwbootfit$t, 2, sd)
write.csv(uwnewdata,file='uwboot.csv')
save(uwmod,file='uwmod.R')
save(uwbootfit,file='uwbootfit.R')

uwnewdata$Prob<-predict(uwmod,newdata=uwnewdata,type='response',re.form=NA)
ggplot(uwnewdata[uwnewdata$dim3==mean(uwnewdata$dim3),],aes(x=YOB,y=Prob))+geom_line(aes(group=Variant),size=1.5,linetype='dotdash')+geom_ribbon(aes(ymax=Prob+(se),ymin=Prob-(se),fill=factor(Variant,levels=levels(Variant)[c(3,6,2,1,5,4)])),alpha=0.25)+theme_bw()+scale_fill_brewer(palette='Set1')+guides(fill=guide_legend(title=NULL))+theme_bw()+theme(panel.grid.major = element_line(colour = "#808080"),legend.position='NULL')+ylab('P (WC)')+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+geom_hline(yintercept=0.5,linetype='dotted')+xlab('Year of birth')+ggtitle('/u/')+theme(panel.border=element_rect(size=1.5,color='black'))+ggsave('/Users/pplsuser/Desktop/Labphon_Posters/u_perception_age_sd.pdf',width=2.6,height=3,device=cairo_pdf)
ggplot(uwnewdata[uwnewdata$YOB==mean(uwnewdata$YOB),],aes(x=dim3,y=Prob))+geom_line(aes(group=Variant),size=1.5,linetype='dotdash')+geom_ribbon(aes(ymax=Prob+(se),ymin=Prob-(se),fill=factor(Variant,levels=levels(Variant)[c(3,6,2,1,5,4)])),alpha=0.25)+theme_bw()+scale_fill_brewer(palette='Set1')+guides(fill=guide_legend(title=NULL))+theme_bw()+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+theme(panel.grid.major = element_line(colour = "#808080"),axis.text.y = element_blank(),axis.ticks.y = element_blank())+ylab('')+geom_hline(yintercept=0.5,linetype='dotted')+xlab('Mobility index')+ggtitle('/u/')+theme(panel.border=element_rect(size=1.5,color='black'))+ggsave('/Users/pplsuser/Desktop/Labphon_Posters/u_perception_dim3_sd.pdf',width=3,height=3,device=cairo_pdf)




####Model comparison plots

#################SOCIAL CLASS
#################OW
owsoc<-data[data$Vowel=='OW'&data$Test=='soc',]
owsoc$Image_class<-relevel(owsoc$Image_class,'M')
owsoc$Sound<-factor(owsoc$Sound)
owsoc$Variant<-factor(owsoc$Variant)
owsoc<-recontrast(owsoc)
owsoc$YOB<-scale(owsoc$YOB)
owsoc$dim1<-scale(owsoc$dim1)
owsoc$dim2<-scale(owsoc$dim2)
owsoc$dim3<-scale(owsoc$dim3)

owsoc$Variant<-mapply(get_ipa,vowel='GOAT',variant=owsoc$Variant)
ownull<-glmer(Image_class~1+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
owbase<-glmer(Image_class~Variant+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
owage<-glmer(Image_class~Variant+YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
owageint<-glmer(Image_class~Variant*YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
owdim3<-glmer(Image_class~Variant+dim3+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
owdim3int<-glmer(Image_class~Variant*dim3+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
owbothints<-glmer(Image_class~Variant+(YOB*Variant)+(dim3*Variant)+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
owdim33way<-glmer(Image_class~Variant*dim3*YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)

uwsoc<-data[data$Vowel=='UW'&data$Test=='soc',]
uwsoc$Image_class<-relevel(uwsoc$Image_class,'M')
uwsoc$Sound<-factor(uwsoc$Sound)
uwsoc$Variant<-factor(uwsoc$Variant)
uwsoc<-recontrast(uwsoc)
uwsoc$YOB<-scale(uwsoc$YOB)
uwsoc$dim1<-scale(uwsoc$dim1)
uwsoc$dim2<-scale(uwsoc$dim2)
uwsoc$dim3<-scale(uwsoc$dim3)
uwsoc$Variant<-mapply(get_ipa,vowel='GOOSE',variant=uwsoc$Variant)

uwnull<-glmer(Image_class~1+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
uwset<-glmer(Image_class~1+Version+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
uwbase<-glmer(Image_class~Variant+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
uwage<-glmer(Image_class~Variant+YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
uwageint<-glmer(Image_class~Variant*YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
uwdim3<-glmer(Image_class~Variant+dim3+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
uwdim3int<-glmer(Image_class~Variant*dim3+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
uwbothints<-glmer(Image_class~Variant+(Variant*dim3)+(Variant*YOB)+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
uwdim33way<-glmer(Image_class~Variant*dim3*YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)

#Chav
data$ChavTrial="No"
data$ChavTrial[grep('M_Y_WC_NL',data$TL)]="Yes"
data$ChavTrial[grep('M_Y_WC_NL',data$TR)]="Yes"
data$Chav="Nonchav"
data$Chav[grep('M_Y_WC_NL',data$Image_selected)]="Chav"
chav<-subset(data,ChavTrial=="Yes")

owchav<-chav[chav$Vowel=='OW',]
owchav$Sound<-factor(owchav$Sound)
owchav$Variant<-factor(owchav$Variant)
owchav<-recontrast(owchav)
owchav$YOB<-scale(owchav$YOB)
owchav$dim1<-scale(owchav$dim1)
owchav$dim2<-scale(owchav$dim2)
owchav$dim3<-scale(owchav$dim3)
owchav$Variant<-mapply(get_ipa,vowel='GOAT',variant=owchav$Variant)
owchav$Variant<-factor(owchav$Variant)
owchav<-recontrast(owchav)
owchav$Chav<-factor(owchav$Chav)
owchav$Chav<-relevel(owchav$Chav,'Nonchav')

owchavnull<-glmer(Chav~1+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owchav,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
owchavbase<-glmer(Chav~Variant+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owchav,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
owchavage<-glmer(Chav~Variant+YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owchav,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
owchavageint<-glmer(Chav~Variant*YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owchav,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
owchavdim3<-glmer(Chav~Variant+dim3+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owchav,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
owchavdim3int<-glmer(Chav~Variant*dim3+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owchav,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
owchavbothints<-glmer(Chav~Variant+(YOB*Variant)+(dim3*Variant)+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owchav,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
owchavdim33way<-glmer(Chav~Variant*dim3*YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owchav,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)

uwchav<-chav[chav$Vowel=='UW',]
uwchav$Image_class<-relevel(uwchav$Image_class,'M')
uwchav$Sound<-factor(uwchav$Sound)
uwchav$Variant<-factor(uwchav$Variant)
uwchav$YOB<-scale(uwchav$YOB)
uwchav$dim1<-scale(uwchav$dim1)
uwchav$dim2<-scale(uwchav$dim2)
uwchav$dim3<-scale(uwchav$dim3)
uwchav$Variant<-mapply(get_ipa,vowel='GOOSE',variant=uwchav$Variant)
uwchav$Variant<-factor(uwchav$Variant)
uwchav$Chav<-factor(uwchav$Chav)
uwchav$Chav<-relevel(uwchav$Chav,'Nonchav')


uwchavnull<-glmer(Chav~1+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwchav,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
uwchavbase<-glmer(Chav~Variant+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwchav,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
uwchavage<-glmer(Chav~Variant+YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwchav,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
uwchavageint<-glmer(Chav~Variant*YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwchav,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
uwchavdim3<-glmer(Chav~Variant+dim3+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwchav,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
uwchavdim3int<-glmer(Chav~Variant*dim3+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwchav,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
uwchavbothints<-glmer(Chav~Variant+(Variant*dim3)+(Variant*YOB)+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwchav,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
uwchavdim33way<-glmer(Chav~Variant*dim3*YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwchav,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)

#Local
uwloc<-data[data$Vowel=='UW'&data$Test=='loc',]
uwloc$Image_local<-relevel(uwloc$Image_local,'N')
uwloc$Sound<-factor(uwloc$Sound)
uwloc$Variant<-factor(uwloc$Variant)
uwloc<-recontrast(uwloc)
uwloc$YOB<-scale(uwloc$YOB)
uwloc$dim1<-scale(uwloc$dim1)
uwloc$dim2<-scale(uwloc$dim2)
uwloc$dim3<-scale(uwloc$dim3)
uwloc$Variant<-mapply(get_ipa,vowel='GOOSE',variant=uwloc$Variant)

uwlocnull<-glmer(Image_local~1+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
uwlocbase<-glmer(Image_local~Variant+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
uwlocage<-glmer(Image_local~Variant+YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
uwlocageint<-glmer(Image_local~Variant*YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
uwlocdim3<-glmer(Image_local~Variant+dim3+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
uwlocdim3int<-glmer(Image_local~Variant*dim3+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
uwlocbothints<-glmer(Image_local~Variant+(Variant*dim3)+(Variant*YOB)+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
uwlocdim33way<-glmer(Image_local~Variant*dim3*YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)

owloc<-data[data$Vowel=='OW'&data$Test=='loc',]
owloc$Image_local<-relevel(owloc$Image_local,'N')
owloc$Sound<-factor(owloc$Sound)
owloc$Variant<-factor(owloc$Variant)
owloc<-recontrast(owloc)
owloc$YOB<-scale(owloc$YOB)
owloc$dim1<-scale(owloc$dim1)
owloc$dim2<-scale(owloc$dim2)
owloc$dim3<-scale(owloc$dim3)
owloc$Variant<-mapply(get_ipa,vowel='GOAT',variant=owloc$Variant)

owlocnull<-glmer(Image_local~1+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
owlocbase<-glmer(Image_local~Variant+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
owlocage<-glmer(Image_local~Variant+YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
owlocageint<-glmer(Image_local~Variant*YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
owlocdim3<-glmer(Image_local~Variant+dim3+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
owlocdim3int<-glmer(Image_local~Variant*dim3+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
owlocbothints<-glmer(Image_local~Variant+(YOB*Variant)+(dim3*Variant)+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)
owlocdim33way<-glmer(Image_local~Variant*dim3*YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)

save.image('all_perception_models.R')
#Bootstraps?
library(boot)
bootfun<-function(fit){
return(AIC(logLik(fit)))
}

ownullresults <- bootMer(ownull,FUN=bootfun ,nsim=100,verbose=T)
owbaseresults <- bootMer(owbase,FUN=bootfun ,nsim=100,verbose=T)
owageresults <- bootMer(owage,FUN=bootfun ,nsim=100,verbose=T)
owageintresults <- bootMer(owageint,FUN=bootfun ,nsim=100,verbose=T)
owdim3results <- bootMer(owdim3,FUN=bootfun ,nsim=100,verbose=T)
owdim3intresults <- bootMer(owdim3int,FUN=bootfun ,nsim=100,verbose=T)
owbothintsresults <- bootMer(owbothints,FUN=bootfun ,nsim=100,verbose=T)
ow3wayresults <- bootMer(owdim33way,FUN=bootfun ,nsim=100,verbose=T)


uwnullresults <- bootMer(uwnull,FUN=bootfun ,nsim=100,verbose=T)
uwbaseresults <- bootMer(uwbase,FUN=bootfun ,nsim=100,verbose=T)
uwageresults <- bootMer(uwage,FUN=bootfun ,nsim=100,verbose=T)
uwageintresults <- bootMer(uwageint,FUN=bootfun ,nsim=100,verbose=T)
uwdim3results <- bootMer(uwdim3,FUN=bootfun ,nsim=100,verbose=T)
uwdim3intresults <- bootMer(uwdim3int,FUN=bootfun ,nsim=100,verbose=T)
uwbothintsresults <- bootMer(uwbothints,FUN=bootfun ,nsim=100,verbose=T)
uw3wayresults <- bootMer(uwdim33way,FUN=bootfun ,nsim=100,verbose=T)


owchavnullresults <- bootMer(owchavnull,FUN=bootfun ,nsim=100,verbose=T)
owchavbaseresults <- bootMer(owchavbase,FUN=bootfun ,nsim=100,verbose=T)
owchavageresults <- bootMer(owchavage,FUN=bootfun ,nsim=100,verbose=T)
owchavageintresults <- bootMer(owchavageint,FUN=bootfun ,nsim=100,verbose=T)
owchavdim3results <- bootMer(owchavdim3,FUN=bootfun ,nsim=100,verbose=T)
owchavdim3intresults <- bootMer(owchavdim3int,FUN=bootfun ,nsim=100,verbose=T)
owchavbothintsresults <- bootMer(owchavbothints,FUN=bootfun ,nsim=100,verbose=T)
owchav3wayresults <- bootMer(owchavdim33way,FUN=bootfun ,nsim=100,verbose=T)

uwchavnullresults <- bootMer(uwchavnull,FUN=bootfun ,nsim=100,verbose=T)
uwchavbaseresults <- bootMer(uwchavbase,FUN=bootfun ,nsim=100,verbose=T)
uwchavageresults <- bootMer(uwchavage,FUN=bootfun ,nsim=100,verbose=T)
uwchavageintresults <- bootMer(uwchavageint,FUN=bootfun ,nsim=100,verbose=T)
uwchavdim3results <- bootMer(uwchavdim3,FUN=bootfun ,nsim=100,verbose=T)
uwchavdim3intresults <- bootMer(uwchavdim3int,FUN=bootfun ,nsim=100,verbose=T)
uwchavbothintsresults <- bootMer(uwchavbothints,FUN=bootfun ,nsim=100,verbose=T)
uwchav3wayresults <- bootMer(uwchavdim33way,FUN=bootfun ,nsim=100,verbose=T)

owlocnullresults <- bootMer(owlocnull,FUN=bootfun ,nsim=100,verbose=T)
owlocbaseresults <- bootMer(owlocbase,FUN=bootfun ,nsim=100,verbose=T)
owlocageresults <- bootMer(owlocage,FUN=bootfun ,nsim=100,verbose=T)
owlocageintresults <- bootMer(owlocageint,FUN=bootfun ,nsim=100,verbose=T)
owlocdim3results <- bootMer(owlocdim3,FUN=bootfun ,nsim=100,verbose=T)
owlocdim3intresults <- bootMer(owlocdim3int,FUN=bootfun ,nsim=100,verbose=T)
owlocbothintsresults <- bootMer(owlocbothints,FUN=bootfun ,nsim=100,verbose=T)
owloc3wayresults <- bootMer(owlocdim33way,FUN=bootfun ,nsim=100,verbose=T)

uwlocnullresults <- bootMer(uwlocnull,FUN=bootfun ,nsim=100,verbose=T)
uwlocbaseresults <- bootMer(uwlocbase,FUN=bootfun ,nsim=100,verbose=T)
uwlocageresults <- bootMer(uwlocage,FUN=bootfun ,nsim=100,verbose=T)
uwlocageintresults <- bootMer(uwlocageint,FUN=bootfun ,nsim=100,verbose=T)
uwlocdim3results <- bootMer(uwlocdim3,FUN=bootfun ,nsim=100,verbose=T)
uwlocdim3intresults <- bootMer(uwlocdim3int,FUN=bootfun ,nsim=100,verbose=T)
uwlocbothintsresults <- bootMer(uwlocbothints,FUN=bootfun ,nsim=100,verbose=T)
uwloc3wayresults <- bootMer(uwlocdim33way,FUN=bootfun ,nsim=100,verbose=T)

owaic<-data.frame(null=ownullresults$t,base=owbaseresults$t,age=owageresults$t,ageint=owageintresults$t,dim3=owdim3results$t,dim3int=owdim3intresults$t,bothints=owbothintsresults$t,threeway=ow3wayresults$t)
confints<-data.frame(t(apply(t(owaic),FUN=quantile,probs=c(0.025,0.5,0.975),MARGIN=1)))
names(confints)<-c('lo','mid','hi')
confints$Model<-row.names(confints)

owaic<-melt(owaic)
ggplot(owaic)+geom_violin(aes(x=variable,y=value))+geom_errorbar(data=confints,aes(x=Model,ymax=hi,ymin=lo))

#3rd august parametric bootstrap
library(pbkrtest)
ownull<-glmer(Image_class~1+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0,REML=F)
owbase<-glmer(Image_class~Variant+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0,REML=F)
PBmodcomp(ownull,owbase,nsim=100)
#Ridiculously slow

anova(owbase,owdim3int)

pbLRT<-function(modnull,mod,n.sim){
lrt.obs <- anova(mod, modnull)$Chisq[2]
lrt.sim <- numeric(n.sim)
dattemp<-mod@frame
pb <- txtProgressBar(min = 0, max = n.sim, style = 3)
for(i in 1:n.sim) {
        Sys.sleep(0.1) # progress bar related stuff
        ysim       <- unlist(simulate(modnull)) # simulate new observations from the null-model  
        modnullsim <- refit(modnull,ysim)
        modaltsim  <- refit(mod,ysim)
        lrt.sim[i] <- anova(modnullsim, modaltsim)$Chisq[2] # save the likelihood ratio test statistic
    setTxtProgressBar(pb, i)
    }
 return(lrt.sim)}

#3rd august model viz
#Better do some models with gender
owgender<-glmer(Image_class~Variant+Part_gender+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0,REML=F)
owgenderint<-glmer(Image_class~Variant*Part_gender+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0,REML=F)
owgenderdim3<-glmer(Image_class~(Variant*Part_gender)+(Variant*dim3)+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0,REML=F)
owgenderage<-glmer(Image_class~(Variant*Part_gender)+(Variant*YOB)+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0,REML=F)

owlocgender<-glmer(Image_class~Variant+Part_gender+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0,REML=F)
owlocgenderint<-glmer(Image_class~Variant*Part_gender+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0,REML=F)

owchavgender<-glmer(Chav~Variant+Part_gender+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owchav,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0,REML=F)
owchavgenderint<-glmer(Chav~Variant*Part_gender+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owchav,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0,REML=F)

uwgender<-glmer(Image_class~Variant+Part_gender+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0,REML=F)
uwgenderint<-glmer(Image_class~Variant*Part_gender+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0,REML=F)

uwlocgender<-glmer(Image_class~Variant+Part_gender+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0,REML=F)
uwlocgenderint<-glmer(Image_class~Variant*Part_gender+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0,REML=F)


uwchavgender<-glmer(Chav~Variant+Part_gender+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwchav,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0,REML=F)
uwchavgenderint<-glmer(Chav~Variant*Part_gender+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwchav,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0,REML=F)


aic<-function(mod){return(c(as.character(formula(mod))[3],AIC(logLik(mod))))}
mods<-list(ownull,owbase,owage,owageint,owdim3,owdim3int,owbothints,owdim33way)
mat<-sapply(mods,aic)
mat<-as.data.frame(t(mat))
names(mat)<-c('Model','AIC')
mat$Model<-as.character(mat$Model)
mat$Model<-sapply(strsplit(mat$Model,split='\\(1'),"[[", 1)
mat$Model<-factor(mat$Model)
mat$Model<-factor(mat$Model,levels=mat$Model)
levels(mat$Model) <- gsub(" ", "\n",mat$Model) 
mat$AIC<-as.numeric(as.character(mat$AIC))
ggplot(mat,aes(x=Model,y=AIC,group=0))+geom_point()+geom_line(linetype='dotted')
owmatclass<-mat
owmatclass$Test<-"'Working class'"

mods<-list(uwnull,uwbase,uwage,uwageint,uwdim3,uwdim3int,uwbothints,uwdim33way)
mat<-sapply(mods,aic)
mat<-as.data.frame(t(mat))
names(mat)<-c('Model','AIC')
mat$Model<-as.character(mat$Model)
mat$Model<-sapply(strsplit(mat$Model,split='\\(1'),"[[", 1)
mat$Model<-factor(mat$Model)
mat$Model<-factor(mat$Model,levels=mat$Model)
levels(mat$Model) <- gsub(" ", "\n",mat$Model) 
mat$AIC<-as.numeric(as.character(mat$AIC))
ggplot(mat,aes(x=Model,y=AIC,group=0))+geom_point()+geom_line(linetype='dotted')
uwmatclass<-mat
uwmatclass$Test<-"'Working class'"

mods<-list(owlocnull,owlocbase,owlocage,owlocageint,owlocdim3,owlocdim3int,owlocbothints,owlocdim33way)
mat<-sapply(mods,aic)
mat<-as.data.frame(t(mat))
names(mat)<-c('Model','AIC')
mat$Model<-as.character(mat$Model)
mat$Model<-sapply(strsplit(mat$Model,split='\\(1'),"[[", 1)
mat$Model<-factor(mat$Model)
mat$Model<-factor(mat$Model,levels=mat$Model)
levels(mat$Model) <- gsub(" ", "\n",mat$Model) 
mat$AIC<-as.numeric(as.character(mat$AIC))
ggplot(mat,aes(x=Model,y=AIC,group=0))+geom_point()+geom_line(linetype='dotted')
owmatloc<-mat
owmatloc$Test<-"'Rural'"


mods<-list(owchavnull,owchavbase,owchavage,owchavageint,owchavdim3,owchavdim3int,owchavbothints,owchavdim33way)
mat<-sapply(mods,aic)
mat<-as.data.frame(t(mat))
names(mat)<-c('Model','AIC')
mat$Model<-as.character(mat$Model)
mat$Model<-sapply(strsplit(mat$Model,split='\\(1'),"[[", 1)
mat$Model<-factor(mat$Model)
mat$Model<-factor(mat$Model,levels=mat$Model)
levels(mat$Model) <- gsub(" ", "\n",mat$Model) 
mat$AIC<-as.numeric(as.character(mat$AIC))
ggplot(mat,aes(x=Model,y=AIC,group=0))+geom_point()+geom_line(linetype='dotted')
owmatchav<-mat
owmatchav$Test<-"'Chav'"

mods<-list(uwlocnull,uwlocbase,uwlocage,uwlocageint,uwlocdim3,uwlocdim3int,uwlocbothints,uwlocdim33way)
mat<-sapply(mods,aic)
mat<-as.data.frame(t(mat))
names(mat)<-c('Model','AIC')
mat$Model<-as.character(mat$Model)
mat$Model<-sapply(strsplit(mat$Model,split='\\(1'),"[[", 1)
mat$Model<-factor(mat$Model)
mat$Model<-factor(mat$Model,levels=mat$Model)
levels(mat$Model) <- gsub(" ", "\n",mat$Model) 
mat$AIC<-as.numeric(as.character(mat$AIC))
ggplot(mat,aes(x=Model,y=AIC,group=0))+geom_point()+geom_line(linetype='dotted')
uwmatloc<-mat
uwmatloc$Test<-"'Rural'"

mods<-list(uwchavnull,uwchavbase,uwchavage,uwchavageint,uwchavdim3,uwchavdim3int,uwchavbothints,uwchavdim33way)
mat<-sapply(mods,aic)
mat<-as.data.frame(t(mat))
names(mat)<-c('Model','AIC')
mat$Model<-as.character(mat$Model)
mat$Model<-sapply(strsplit(mat$Model,split='\\(1'),"[[", 1)
mat$Model<-factor(mat$Model)
mat$Model<-factor(mat$Model,levels=mat$Model)
levels(mat$Model) <- gsub(" ", "\n",mat$Model) 
mat$AIC<-as.numeric(as.character(mat$AIC))
ggplot(mat,aes(x=Model,y=AIC,group=0))+geom_point()+geom_line(linetype='dotted')
uwmatchav<-mat
uwmatchav$Test<-"'Chav'"

owmat<-rbind(owmatclass,owmatloc,owmatchav)
levels(owmat$Model)<-c('Baseline','Variant','Variant + YOB','Variant * YOB','Variant + Mobility index','Variant * Mobility index','Both ints','Three-way int')
levels(owmat$Model) <- gsub(" ", "\n",owmat$Model) 
owmat$Test<-factor(owmat$Test)
owmat$Test<-factor(owmat$Test,levels=c("'Working class'","'Rural'","'Chav'"))
ggplot(owmat,aes(x=Model,y=AIC,group=0))+geom_point(size=0.7)+geom_line(linetype='dotted',color='#386cb0')+facet_wrap(~Test,ncol=1,scales='free_y')+scale_y_continuous(breaks = myBreaks)+ggtitle('Model selection: /o/')+theme_bw()+theme(panel.border=element_rect(size=1.5,color='black'),strip.background=element_blank(),strip.text=element_blank(),plot.title=element_text(size=10),axis.text=element_text(size=6),axis.title.x=element_text(size=10),axis.title.y=element_text(size=10))+ggsave('/Users/pplsuser/Desktop/PWL Submission/owaic.pdf',width=3.1,height=4.5,device=cairo_pdf)



uwmat<-rbind(uwmatclass,uwmatloc,uwmatchav)
levels(uwmat$Model)<-c('Baseline','Variant','Variant + YOB','Variant * YOB','Variant + Mobility index','Variant * Mobility index','Both ints','Three-way int')
levels(uwmat$Model) <- gsub(" ", "\n",uwmat$Model) 
uwmat$Test<-factor(uwmat$Test)
uwmat$Test<-factor(uwmat$Test,levels=c("'Working class'","'Rural'","'Chav'"))
ggplot(uwmat,aes(x=Model,y=AIC,group=0))+geom_point(size=0.7)+geom_line(linetype='dotted',color='#386cb0')+facet_wrap(~Test,ncol=1,scales='free_y')+scale_y_continuous(breaks = myBreaks)+ggtitle('Model selection: /u/')+theme_bw()+theme(panel.border=element_rect(size=1.5,color='black'),strip.background=element_blank(),strip.text=element_blank(),plot.title=element_text(size=10),axis.text=element_text(size=6),axis.title.x=element_text(size=10),axis.title.y=element_blank())+ggsave('/Users/pplsuser/Desktop/PWL Submission/uwaic.pdf',width=2.9,height=4.5,device=cairo_pdf)

myBreaks <- function(x){
    breaks <- quantile(x,probs=(c(0,0.25,0.5,0.75,1)))
    names(breaks)<-round(breaks)
    breaks
}

#2nd august making plots for PWPL, effects plots w/sampling from posterior
library(dplyr)
library(tidyr)
library(magrittr)
library(itsadug)
library(MASS)
library(reshape2)
library(mgcv)
library(coda)
library(ggplot2)
library(RColorBrewer)
#################OW
#################SOCIAL CLASS
owsoc<-data[data$Vowel=='OW'&data$Test=='soc',]
owsoc$Image_class<-relevel(owsoc$Image_class,'M')
owsoc$Sound<-factor(owsoc$Sound)
owsoc$Variant<-factor(owsoc$Variant)
owsoc$YOB<-scale(owsoc$YOB)
owsoc$dim1<-scale(owsoc$dim1)
owsoc$dim2<-scale(owsoc$dim2)
owsoc$dim3<-scale(owsoc$dim3)
owsoc$Variant<-mapply(get_ipa,vowel='GOAT',variant=owsoc$Variant)
owsoc$Variant<-factor(owsoc$Variant)
owsoc<-recontrast(owsoc)

ow<-glmer(Image_class~Variant*dim3+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)

newdat<-expand.grid(Variant=levels(owsoc$Variant),dim3=seq(-2,2,by=0.1))
newdat$Variant<-factor(newdat$Variant)
newdat<-recontrast(newdat)

mm<-model.matrix(~Variant*dim3,newdat)
beta <- fixef(ow)
Sigma <- vcov(ow)
posterior <- mvrnorm(n = 10000, mu = beta, Sigma = Sigma)
mm%*%posterior[1,]
mult_it <- function(x, y){
  x %*% y
}

sims <- apply(posterior, MARGIN= 1, FUN = mult_it, x = mm)

get_hpd <- function(x){
  x %>% as.mcmc() %>% HPDinterval()
}

hpd_interval <- plogis(t(apply(sims,MARGIN=1,FUN = get_hpd)))
sims<-data.frame(sims)
hpd_interval <-data.frame(hpd_interval)
hpd_interval<-cbind(newdat,hpd_interval)
names(hpd_interval)<-c('Variant','dim3','lo','hi')
sims<-cbind(newdat,sims)
sims<-melt(sims,id.vars=names(newdat))
sims$Variant<-factor(sims$Variant,levels=levels(sims$Variant)[c(5,8,6,7,2,1,3,4)])
hpd_interval$Variant<-factor(hpd_interval$Variant,levels=levels(hpd_interval$Variant)[c(5,8,6,7,2,1,3,4)])
ggplot(sims[sims$dim3==mean(sims$dim3),])+geom_violin(aes(x=Variant,y=plogis(value)),color='#386cb0')+stat_summary(aes(x=Variant,y=plogis(value)),size=0.15)+geom_errorbar(data=hpd_interval[hpd_interval$dim3==mean(hpd_interval$dim3),],aes(x=Variant,ymin=lo,ymax=hi),width=0.35)+theme_bw()+geom_hline(yintercept=0.5,linetype='dotted')+theme(panel.border=element_rect(size=1.5,color='black'),plot.title=element_text(size=10),axis.title.x=element_text(size=10),axis.title.y=element_text(size=10))+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+ylab('P(WC)')+ggtitle("/o/:'Working-class'")+ggsave('/Users/pplsuser/Desktop/PWL submission/ow_class.pdf',device=cairo_pdf,height=2.5,width=2.5)
ggplot(sims)+geom_smooth(aes(x=dim3,y=plogis(value),group=Variant,color=Variant),linetype='dotdash')+geom_ribbon(data=hpd_interval,aes(x=dim3,ymin=lo,ymax=hi,fill=Variant),alpha=0.05)+theme_bw()+theme(panel.border=element_rect(size=1.5,color='black'))+xlab('Mobility index')+ylab('P(WC)')+scale_fill_brewer(palette='Set1')+scale_color_brewer(palette='Set1')+ggsave('/Users/pplsuser/Desktop/PWL Submission/ow_class_dim3.pdf',device=cairo_pdf,width=3.25,height=3)
	+geom_smooth(data=hpd_interval,aes(y=lo/1.96,x=dim3,color=Variant),linetype='dotted')
ggplot(sims)+geom_smooth(aes(x=dim3,y=plogis(value),group=Variant,color=Variant),linetype='dotdash')+geom_ribbon(data=hpd_interval,aes(x=dim3,ymin=lo,ymax=hi,fill=Variant),alpha=0.1)+theme_bw()+theme(panel.border=element_rect(size=1.5,color='black'),plot.title=element_text(size=10),axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),legend.title=element_blank())+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+geom_hline(yintercept=0.5,color='black',linetype='dotted')+xlab('Mobility index')+ylab('P(WC)')+ggtitle("/o/:'Working class'") +scale_fill_brewer(palette='Set1')+scale_color_brewer(palette='Set1')+ggsave('/Users/pplsuser/Desktop/PWL Submission/ow_class_dim3.pdf',device=cairo_pdf,width=3.25,height=3)


#################LOCAL
owloc<-data[data$Vowel=='OW'&data$Test=='loc',]
owloc$Image_local<-relevel(owloc$Image_local,'N')
owloc$Sound<-factor(owloc$Sound)
owloc$Variant<-factor(owloc$Variant)
owloc$YOB<-scale(owloc$YOB)
owloc$dim1<-scale(owloc$dim1)
owloc$dim2<-scale(owloc$dim2)
owloc$dim3<-scale(owloc$dim3)
owloc$Variant<-mapply(get_ipa,vowel='GOAT',variant=owloc$Variant)
owloc$Variant<-factor(owloc$Variant)
owloc<-recontrast(owloc)

owlocmod<-glmer(Image_local~Variant+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)

newdat<-expand.grid(Variant=levels(owloc$Variant))
newdat$Variant<-factor(newdat$Variant)
newdat<-recontrast(newdat)

mm<-model.matrix(~Variant,newdat)
beta <- fixef(owlocmod)
Sigma <- vcov(owlocmod)
posterior <- mvrnorm(n = 10000, mu = beta, Sigma = Sigma)
mm%*%posterior[1,]
mult_it <- function(x, y){
  x %*% y
}

sims <- apply(posterior, MARGIN= 1, FUN = mult_it, x = mm)

get_hpd <- function(x){
  x %>% as.mcmc() %>% HPDinterval()
}

hpd_interval <- plogis(t(apply(sims,MARGIN=1,FUN = get_hpd)))
sims<-data.frame(sims)
hpd_interval <-data.frame(hpd_interval)
hpd_interval<-cbind(newdat,hpd_interval)
names(hpd_interval)<-c('Variant','lo','hi')
sims<-cbind(newdat,sims)
sims<-melt(sims,id.vars=names(newdat))
sims$Variant<-factor(sims$Variant,levels=levels(sims$Variant)[c(5,8,6,7,2,1,3,4)])
hpd_interval$Variant<-factor(hpd_interval$Variant,levels=levels(hpd_interval$Variant)[c(5,8,6,7,2,1,3,4)])
ggplot(sims)+geom_violin(aes(x=Variant,y=plogis(value)),color='#386cb0')+stat_summary(aes(x=Variant,y=plogis(value)),size=0.15)+geom_errorbar(data=hpd_interval,aes(x=Variant,ymin=lo,ymax=hi),width=0.35)+theme_bw()+geom_hline(yintercept=0.5,linetype='dotted')+theme(panel.border=element_rect(size=1.5,color='black'),plot.title=element_text(size=10),axis.title.x=element_text(size=10),axis.title.y=element_text(size=10))+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+ylab('P(Rural)')+ggtitle("/o/:'Rural'")+ggsave('/Users/pplsuser/Desktop/PWL submission/ow_local.pdf',device=cairo_pdf,height=2.5,width=2.5)


#Chav
data$ChavTrial="No"
data$ChavTrial[grep('M_Y_WC_NL',data$TL)]="Yes"
data$ChavTrial[grep('M_Y_WC_NL',data$TR)]="Yes"
data$Chav="Nonchav"
data$Chav[grep('M_Y_WC_NL',data$Image_selected)]="Chav"
chav<-subset(data,ChavTrial=="Yes")

owchav<-chav[chav$Vowel=='OW',]
owchav$Sound<-factor(owchav$Sound)
owchav$Variant<-factor(owchav$Variant)
owchav<-recontrast(owchav)
owchav$YOB<-scale(owchav$YOB)
owchav$dim1<-scale(owchav$dim1)
owchav$dim2<-scale(owchav$dim2)
owchav$dim3<-scale(owchav$dim3)
owchav$Chav<-factor(owchav$Chav)
owchav$Chav<-relevel(owchav$Chav,'Nonchav')
owchav$Variant<-mapply(get_ipa,vowel='GOAT',variant=owchav$Variant)
owchav$Variant<-factor(owchav$Variant)
owchav<-recontrast(owchav)

owchavmod<-glmer(Chav~Variant+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owchav,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)

newdat<-expand.grid(Variant=levels(owchav$Variant))
newdat$Variant<-factor(newdat$Variant)
newdat<-recontrast(newdat)

mm<-model.matrix(~Variant,newdat)
beta <- fixef(owchavmod)
Sigma <- vcov(owchavmod)
posterior <- mvrnorm(n = 10000, mu = beta, Sigma = Sigma)
mm%*%posterior[1,]
mult_it <- function(x, y){
  x %*% y
}

sims <- apply(posterior, MARGIN= 1, FUN = mult_it, x = mm)

get_hpd <- function(x){
  x %>% as.mcmc() %>% HPDinterval()
}

hpd_interval <- plogis(t(apply(sims,MARGIN=1,FUN = get_hpd)))
sims<-data.frame(sims)
hpd_interval <-data.frame(hpd_interval)
hpd_interval<-cbind(newdat,hpd_interval)
names(hpd_interval)<-c('Variant','lo','hi')
sims<-cbind(newdat,sims)
sims<-melt(sims,id.vars=names(newdat))
sims$Variant<-factor(sims$Variant,levels=levels(sims$Variant)[c(5,8,6,7,2,1,3,4)])
hpd_interval$Variant<-factor(hpd_interval$Variant,levels=levels(hpd_interval$Variant)[c(5,8,6,7,2,1,3,4)])
ggplot(sims)+geom_violin(aes(x=Variant,y=plogis(value)),color='#386cb0')+stat_summary(aes(x=Variant,y=plogis(value)),size=0.15)+geom_errorbar(data=hpd_interval,aes(x=Variant,ymin=lo,ymax=hi),width=0.35)+theme_bw()+geom_hline(yintercept=0.5,linetype='dotted')+theme(panel.border=element_rect(size=1.5,color='black'),plot.title=element_text(size=10),axis.title.x=element_text(size=10),axis.title.y=element_text(size=10))+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+ylab('P(Chav)')+ggtitle("/o/:'Chav'")+ggsave('/Users/pplsuser/Desktop/PWL submission/ow_chav.pdf',device=cairo_pdf,height=2.5,width=2.5)


#################UW
#################SOCIAL CLASS
uwsoc<-data[data$Vowel=='UW'&data$Test=='soc',]
uwsoc$Image_class<-relevel(uwsoc$Image_class,'M')
uwsoc$Sound<-factor(uwsoc$Sound)
uwsoc$Variant<-factor(uwsoc$Variant)
uwsoc$YOB<-scale(uwsoc$YOB)
uwsoc$dim1<-scale(uwsoc$dim1)
uwsoc$dim2<-scale(uwsoc$dim2)
uwsoc$dim3<-scale(uwsoc$dim3)
uwsoc$Variant<-mapply(get_ipa,vowel='GOOSE',variant=uwsoc$Variant)
uwsoc$Variant<-factor(uwsoc$Variant)
uwsoc<-recontrast(uwsoc)

uw<-glmer(Image_class~Variant+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)

newdat<-expand.grid(Variant=levels(uwsoc$Variant))
newdat$Variant<-factor(newdat$Variant)
newdat<-recontrast(newdat)

mm<-model.matrix(~Variant,newdat)
beta <- fixef(uw)
Sigma <- vcov(uw)
posterior <- mvrnorm(n = 10000, mu = beta, Sigma = Sigma)
mm%*%posterior[1,]
mult_it <- function(x, y){
  x %*% y
}

sims <- apply(posterior, MARGIN= 1, FUN = mult_it, x = mm)

get_hpd <- function(x){
  x %>% as.mcmc() %>% HPDinterval()
}

hpd_interval <- plogis(t(apply(sims,MARGIN=1,FUN = get_hpd)))
sims<-data.frame(sims)
hpd_interval <-data.frame(hpd_interval)
hpd_interval<-cbind(newdat,hpd_interval)
names(hpd_interval)<-c('Variant','lo','hi')
sims<-cbind(newdat,sims)
sims<-melt(sims,id.vars=names(newdat))
sims$Variant<-factor(sims$Variant,levels=levels(sims$Variant)[c(6,5,4,3,2,1)])
hpd_interval$Variant<-factor(hpd_interval$Variant,levels=levels(hpd_interval$Variant)[c(6,5,4,3,2,1)])
ggplot(sims)+geom_violin(aes(x=Variant,y=plogis(value)),color='#386cb0')+stat_summary(aes(x=Variant,y=plogis(value)),size=0.15)+geom_errorbar(data=hpd_interval,aes(x=Variant,ymin=lo,ymax=hi),width=0.35)+theme_bw()+geom_hline(yintercept=0.5,linetype='dotted')+theme(panel.border=element_rect(size=1.5,color='black'),plot.title=element_text(size=10),axis.title.x=element_text(size=10),axis.title.y=element_text(size=10))+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+ylab('P(WC)')+ggtitle("/u/:'Working-class'")+ggsave('/Users/pplsuser/Desktop/PWL submission/uw_class.pdf',device=cairo_pdf,height=2.5,width=2.5)

#################LOCAL
uwloc<-data[data$Vowel=='UW'&data$Test=='loc',]
uwloc$Image_local<-relevel(uwloc$Image_local,'N')
uwloc$Sound<-factor(uwloc$Sound)
uwloc$Variant<-factor(uwloc$Variant)
uwloc$YOB<-scale(uwloc$YOB)
uwloc$dim1<-scale(uwloc$dim1)
uwloc$dim2<-scale(uwloc$dim2)
uwloc$dim3<-scale(uwloc$dim3)
uwloc$Variant<-mapply(get_ipa,vowel='GOOSE',variant=uwloc$Variant)
uwloc$Variant<-factor(uwloc$Variant)
uwloc<-recontrast(uwloc)

uwlocal<-glmer(Image_local~Variant+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwloc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)

newdat<-expand.grid(Variant=levels(uwloc$Variant))
newdat$Variant<-factor(newdat$Variant)
newdat<-recontrast(newdat)

mm<-model.matrix(~Variant,newdat)
beta <- fixef(uwlocal)
Sigma <- vcov(uwlocal)
posterior <- mvrnorm(n = 10000, mu = beta, Sigma = Sigma)
mm%*%posterior[1,]
mult_it <- function(x, y){
  x %*% y
}

sims <- apply(posterior, MARGIN= 1, FUN = mult_it, x = mm)

get_hpd <- function(x){
  x %>% as.mcmc() %>% HPDinterval()
}

hpd_interval <- plogis(t(apply(sims,MARGIN=1,FUN = get_hpd)))
sims<-data.frame(sims)
hpd_interval <-data.frame(hpd_interval)
hpd_interval<-cbind(newdat,hpd_interval)
names(hpd_interval)<-c('Variant','lo','hi')
sims<-cbind(newdat,sims)
sims<-melt(sims,id.vars=names(newdat))
sims$Variant<-factor(sims$Variant,levels=levels(sims$Variant)[c(6,5,4,3,2,1)])
hpd_interval$Variant<-factor(hpd_interval$Variant,levels=levels(hpd_interval$Variant)[c(6,5,4,3,2,1)])
ggplot(sims)+geom_violin(aes(x=Variant,y=plogis(value)),color='#386cb0')+stat_summary(aes(x=Variant,y=plogis(value)),size=0.15)+geom_errorbar(data=hpd_interval,aes(x=Variant,ymin=lo,ymax=hi),width=0.35)+theme_bw()+geom_hline(yintercept=0.5,linetype='dotted')+theme(panel.border=element_rect(size=1.5,color='black'),plot.title=element_text(size=10),axis.title.x=element_text(size=10),axis.title.y=element_text(size=10))+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+ggtitle("/u/:'Rural'")+ylab('P(Rural)')+ggsave('/Users/pplsuser/Desktop/PWL submission/uw_local.pdf',device=cairo_pdf,height=2.5,width=2.5)

######Chav
data$ChavTrial="No"
data$ChavTrial[grep('M_Y_WC_NL',data$TL)]="Yes"
data$ChavTrial[grep('M_Y_WC_NL',data$TR)]="Yes"
data$Chav="Nonchav"
data$Chav[grep('M_Y_WC_NL',data$Image_selected)]="Chav"
chav<-subset(data,ChavTrial=="Yes")

uwchav<-chav[chav$Vowel=='UW',]
uwchav$Sound<-factor(uwchav$Sound)
uwchav$Variant<-factor(uwchav$Variant)
uwchav<-recontrast(uwchav)
uwchav$YOB<-scale(uwchav$YOB)
uwchav$dim1<-scale(uwchav$dim1)
uwchav$dim2<-scale(uwchav$dim2)
uwchav$dim3<-scale(uwchav$dim3)
uwchav$Chav<-factor(uwchav$Chav)
uwchav$Chav<-relevel(uwchav$Chav,'Nonchav')
uwchav$Variant<-mapply(get_ipa,vowel='GOOSE',variant=uwchav$Variant)
uwchav$Variant<-factor(uwchav$Variant)
uwchav<-recontrast(uwchav)

uwchavmod<-glmer(Chav~Variant+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwchav,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)

newdat<-expand.grid(Variant=levels(uwchav$Variant))
newdat$Variant<-factor(newdat$Variant)
newdat<-recontrast(newdat)

mm<-model.matrix(~Variant,newdat)
beta <- fixef(uwchavmod)
Sigma <- vcov(uwchavmod)
posterior <- mvrnorm(n = 10000, mu = beta, Sigma = Sigma)
mm%*%posterior[1,]
mult_it <- function(x, y){
  x %*% y
}

sims <- apply(posterior, MARGIN= 1, FUN = mult_it, x = mm)

get_hpd <- function(x){
  x %>% as.mcmc() %>% HPDinterval()
}

hpd_interval <- plogis(t(apply(sims,MARGIN=1,FUN = get_hpd)))
sims<-data.frame(sims)
hpd_interval <-data.frame(hpd_interval)
hpd_interval<-cbind(newdat,hpd_interval)
names(hpd_interval)<-c('Variant','lo','hi')
sims<-cbind(newdat,sims)
sims<-melt(sims,id.vars=names(newdat))
sims$Variant<-factor(sims$Variant,levels=levels(sims$Variant)[c(6,5,4,3,2,1)])
hpd_interval$Variant<-factor(hpd_interval$Variant,levels=levels(hpd_interval$Variant)[c(6,5,4,3,2,1)])
ggplot(sims)+geom_violin(aes(x=Variant,y=plogis(value)),color='#386cb0')+stat_summary(aes(x=Variant,y=plogis(value)),size=0.15)+geom_errorbar(data=hpd_interval,aes(x=Variant,ymin=lo,ymax=hi),width=0.35)+theme_bw()+geom_hline(yintercept=0.5,linetype='dotted')+theme(panel.border=element_rect(size=1.5,color='black'),plot.title=element_text(size=10),axis.title.x=element_text(size=10),axis.title.y=element_text(size=10))+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+ylab('P(Chav)')+ggtitle("/u/:'Chav'")+ggsave('/Users/pplsuser/Desktop/PWL submission/uw_chav.pdf',device=cairo_pdf,height=2.5,width=2.5)


#################interactions
#################OW CHAV Age

owchav<-chav[chav$Vowel=='OW',]
owchav$Sound<-factor(owchav$Sound)
owchav$Variant<-factor(owchav$Variant)
owchav<-recontrast(owchav)
owchav$dim1<-scale(owchav$dim1)
owchav$dim2<-scale(owchav$dim2)
owchav$dim3<-scale(owchav$dim3)
owchav$Chav<-factor(owchav$Chav)
owchav$Chav<-relevel(owchav$Chav,'Nonchav')
owchav$Variant<-mapply(get_ipa,vowel='GOAT',variant=owchav$Variant)
owchav$Variant<-factor(owchav$Variant)
owchav<-recontrast(owchav)

owchavmod<-glmer(Chav~Variant*YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=owchav,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)



newdat<-expand.grid(Variant=levels(owchav$Variant),YOB=1935:2001)
newdat$Variant<-factor(newdat$Variant)
newdat<-recontrast(newdat)

mm<-model.matrix(~Variant*YOB,newdat)
beta <- fixef(owchavmod)
Sigma <- vcov(owchavmod)
posterior <- mvrnorm(n = 10000, mu = beta, Sigma = Sigma)
mm%*%posterior[1,]
mult_it <- function(x, y){
  x %*% y
}

sims <- apply(posterior, MARGIN= 1, FUN = mult_it, x = mm)

get_hpd <- function(x){
  x %>% as.mcmc() %>% HPDinterval()
}

hpd_interval <- plogis(t(apply(sims,MARGIN=1,FUN = get_hpd)))
sims<-data.frame(sims)
hpd_interval <-data.frame(hpd_interval)
hpd_interval<-cbind(newdat,hpd_interval)
names(hpd_interval)<-c('Variant','YOB','lo','hi')
sims<-cbind(newdat,sims)
sims<-melt(sims,id.vars=names(newdat))
sims$Variant<-factor(sims$Variant,levels=levels(sims$Variant)[c(5,8,6,7,2,1,3,4)])
hpd_interval$Variant<-factor(hpd_interval$Variant,levels=levels(hpd_interval$Variant)[c(5,8,6,7,2,1,3,4)])
ggplot(sims)
ggplot(sims)+geom_smooth(aes(x=YOB,y=plogis(value),group=Variant,color=Variant),linetype='dotdash')+geom_ribbon(data=hpd_interval,aes(x=YOB,ymin=lo,ymax=hi,fill=Variant),alpha=0.1)+theme_bw()+theme(panel.border=element_rect(size=1.5,color='black'),plot.title=element_text(size=10),axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),legend.title=element_blank())+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+geom_hline(yintercept=0.5,color='black',linetype='dotted')+xlab('Year of birth')+ylab('P(Chav)')+ggtitle("/o/:'Chav'") +scale_fill_brewer(palette='Set1')+scale_color_brewer(palette='Set1')+ggsave('/Users/pplsuser/Desktop/PWL Submission/ow_chav_age.pdf',device=cairo_pdf,width=3.25,height=3)


#################UW
#################SOCIAL CLASS AGE
uwsoc<-data[data$Vowel=='UW'&data$Test=='soc',]
uwsoc$Image_class<-relevel(uwsoc$Image_class,'M')
uwsoc$Sound<-factor(uwsoc$Sound)
uwsoc$Variant<-factor(uwsoc$Variant)
uwsoc$dim1<-scale(uwsoc$dim1)
uwsoc$dim2<-scale(uwsoc$dim2)
uwsoc$dim3<-scale(uwsoc$dim3)
uwsoc$Variant<-mapply(get_ipa,vowel='GOOSE',variant=uwsoc$Variant)
uwsoc$Variant<-factor(uwsoc$Variant)
uwsoc<-recontrast(uwsoc)

uw<-glmer(Image_class~Variant*YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwsoc,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)

newdat<-expand.grid(Variant=levels(uwsoc$Variant),YOB=1935:2001)
newdat$Variant<-factor(newdat$Variant)
newdat<-recontrast(newdat)

mm<-model.matrix(~Variant*YOB,newdat)
beta <- fixef(uw)
Sigma <- vcov(uw)
posterior <- mvrnorm(n = 10000, mu = beta, Sigma = Sigma)
mm%*%posterior[1,]
mult_it <- function(x, y){
  x %*% y
}

sims <- apply(posterior, MARGIN= 1, FUN = mult_it, x = mm)

get_hpd <- function(x){
  x %>% as.mcmc() %>% HPDinterval()
}

hpd_interval <- plogis(t(apply(sims,MARGIN=1,FUN = get_hpd)))
sims<-data.frame(sims)
hpd_interval <-data.frame(hpd_interval)
hpd_interval<-cbind(newdat,hpd_interval)
names(hpd_interval)<-c('Variant','YOB','lo','hi')
sims<-cbind(newdat,sims)
sims<-melt(sims,id.vars=names(newdat))
sims$Variant<-factor(sims$Variant,levels=levels(sims$Variant)[c(6,5,4,3,2,1)])
hpd_interval$Variant<-factor(hpd_interval$Variant,levels=levels(hpd_interval$Variant)[c(6,5,4,3,2,1)])
ggplot(sims)+geom_smooth(aes(x=YOB,y=plogis(value),group=Variant,color=Variant),linetype='dotdash')+geom_ribbon(data=hpd_interval,aes(x=YOB,ymin=lo,ymax=hi,fill=Variant),alpha=0.1)+theme_bw()+theme(panel.border=element_rect(size=1.5,color='black'),plot.title=element_text(size=10),axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),legend.title=element_blank())+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+geom_hline(yintercept=0.5,color='black',linetype='dotted')+xlab('Year of birth')+ylab('P(WC)')+ggtitle("/u/:'Working class'") +scale_fill_brewer(palette='Set1')+scale_color_brewer(palette='Set1')+ggsave('/Users/pplsuser/Desktop/PWL Submission/uw_class_age.pdf',device=cairo_pdf,width=3.25,height=3)



uwchav<-chav[chav$Vowel=='UW',]
uwchav$Sound<-factor(uwchav$Sound)
uwchav$Variant<-factor(uwchav$Variant)
uwchav<-recontrast(uwchav)
uwchav$dim1<-scale(uwchav$dim1)
uwchav$dim2<-scale(uwchav$dim2)
uwchav$dim3<-scale(uwchav$dim3)
uwchav$Chav<-factor(uwchav$Chav)
uwchav$Chav<-relevel(uwchav$Chav,'Nonchav')
uwchav$Variant<-mapply(get_ipa,vowel='GOOSE',variant=uwchav$Variant)
uwchav$Variant<-factor(uwchav$Variant)
uwchav<-recontrast(uwchav)

uwchavmod<-glmer(Chav~Variant*YOB+(1|Sound/Variant)+(1|New_ID)+(0+Variant|New_ID),data=uwchav,family='binomial',verbose=1,control=glmerControl(optimizer='bobyqa'),nAGQ=0)



newdat<-expand.grid(Variant=levels(uwchav$Variant),YOB=1935:2001)
newdat$Variant<-factor(newdat$Variant)
newdat<-recontrast(newdat)

mm<-model.matrix(~Variant*YOB,newdat)
beta <- fixef(uwchavmod)
Sigma <- vcov(uwchavmod)
posterior <- mvrnorm(n = 10000, mu = beta, Sigma = Sigma)
mm%*%posterior[1,]
mult_it <- function(x, y){
  x %*% y
}

sims <- apply(posterior, MARGIN= 1, FUN = mult_it, x = mm)

get_hpd <- function(x){
  x %>% as.mcmc() %>% HPDinterval()
}

hpd_interval <- plogis(t(apply(sims,MARGIN=1,FUN = get_hpd)))
sims<-data.frame(sims)
hpd_interval <-data.frame(hpd_interval)
hpd_interval<-cbind(newdat,hpd_interval)
names(hpd_interval)<-c('Variant','YOB','lo','hi')
sims<-cbind(newdat,sims)
sims<-melt(sims,id.vars=names(newdat))
sims$Variant<-factor(sims$Variant,levels=levels(sims$Variant)[c(6,5,4,3,2,1)])
hpd_interval$Variant<-factor(hpd_interval$Variant,levels=levels(hpd_interval$Variant)[c(6,5,4,3,2,1)])
ggplot(sims)
ggplot(sims)+geom_smooth(aes(x=YOB,y=plogis(value),group=Variant,color=Variant),linetype='dotdash')+geom_ribbon(data=hpd_interval,aes(x=YOB,ymin=lo,ymax=hi,fill=Variant),alpha=0.1)+theme_bw()+theme(panel.border=element_rect(size=1.5,color='black'),plot.title=element_text(size=10),axis.title.x=element_text(size=10),axis.title.y=element_text(size=10),legend.title=element_blank())+scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25))+geom_hline(yintercept=0.5,color='black',linetype='dotted')+xlab('Year of birth')+ylab('P(Chav)')+ggtitle("/u/:'Chav'") +scale_fill_brewer(palette='Set1')+scale_color_brewer(palette='Set1')+ggsave('/Users/pplsuser/Desktop/PWL Submission/uw_chav_age.pdf',device=cairo_pdf,width=3.25,height=3)
