#load relevant packages
library(ggplot2)
library(data.table)
library(gtools)
library(RColorBrewer)
library(scales)
library(viridis)
library(gridExtra)
library(fmsb)
library(data.table)

# Function that returns Root Mean Squared Error
rmse <- function(error)
{
  sqrt(mean(error^2))
}


#download the data
download.file("https://www.epa.gov/sites/production/files/2017-02/nla12_keyvariables_data.csv",
              "D:/LimnologicalDiscretization/nla12_keyvariables_data.csv")

#read in the data
NLA<-fread("D:/LimnologicalDiscretization/nla12_keyvariables_data.csv",stringsAsFactors = TRUE)
NLA$INDEX_SITE_DEPTH.jit<-jitter(NLA$INDEX_SITE_DEPTH,0.0001)
NLA$CHLX_RESULT.jit<-jitter(NLA$CHLX_RESULT,0.0001)+0.0001
NLA<-NLA[is.na(INDEX_SITE_DEPTH)==F&
           is.na(CHLX_RESULT)==F]

NLA_discmod<-as.data.table((matrix(NA, nrow = 1124, ncol = 9)))

for (j in c(1:8)) {
  i<-(2^j)-1
  data<-data.frame(pred=as.factor(quantcut(log(NLA$INDEX_SITE_DEPTH.jit),(i+1),labels=FALSE)),
                   pred.lm=log(NLA$INDEX_SITE_DEPTH.jit),
                   resp=as.numeric(log(NLA$CHLX_RESULT.jit)))
  data<-as.data.table(data)
  data<-na.omit(data)
  
  #full model
  fit<-lm(data=data,resp~pred)
  summaryfit<-summary(fit)
  
  NLA_discmod[,j]<-(fit$fitted.values)
  
}

#now for the lm
data<-data.frame(pred=as.factor(quantcut(log(NLA$INDEX_SITE_DEPTH.jit),(i+1),labels=FALSE)),
                 pred.lm=log(NLA$INDEX_SITE_DEPTH.jit),
                 resp=as.numeric(log(NLA$CHLX_RESULT.jit)))
data<-as.data.table(data)
data<-na.omit(data)

#full model
fit<-lm(data=data,resp~pred.lm)
#summaryfit<-summary(fit)
NLA_discmod[,9]<-(fit$fitted.values)
NLA_discmod$INDEX_SITE_DEPTH.jit<-NLA$INDEX_SITE_DEPTH.jit
names(NLA_discmod)[c(1:9)]<-c("2 Discrete Groups",
                              "4 Discrete Groups",
                              "8 Discrete Groups",
                              "16 Discrete Groups",
                              "32 Discrete Groups",
                              "64 Discrete Groups",
                              "128 Discrete Groups",
                              "256 Discrete Groups",
                              "Continuous Gradient")

NLA.m<-melt(NLA_discmod,id="INDEX_SITE_DEPTH.jit")
names(NLA.m)<-c("INDEX_SITE_DEPTH.jit","discretization","CHLX_RESULT.mod")
NLA.m$category<-"Discrete"
NLA.m[discretization=="Continuous Gradient"]$category<-"Continuous"

mods<-data.frame(RMSE=numeric(0),
                 Pred_RMSE=numeric(0),
                 discretization=numeric(0),
                 RMSE.lm=numeric(0),
                 Pred_RMSE.lm=numeric(0)
                 )
summary(mods)
#i<-4

for (i in c(1:562)) {
  print(i)
  data<-data.frame(pred=as.factor(quantcut(log(NLA$INDEX_SITE_DEPTH.jit),(i+1),labels=FALSE)),pred.lm=log(NLA$INDEX_SITE_DEPTH.jit),resp=as.numeric(log(NLA$CHLX_RESULT.jit)))
  data<-as.data.table(data)
  data<-na.omit(data)
  
  #full model
  fit<-lm(data=data,resp~pred)
  summaryfit<-summary(fit)
  
  fit.lm<-lm(data=data,resp~pred.lm)
  summaryfit.lm<-summary(fit.lm)
  
  #cross validation
  
  RMSEvec<-rep(NA,10000)
  
  RMSEvec.lm<-rep(NA,10000)
  
  for(j in 1:10000){
  data_train<-data[,.SD[sample(x=.N, size = round((1124/(length(unique(data$pred))))/2,0))],by = pred]
  data_test<-subset(data, !(resp %in% data_train$resp))
  data_train2<-data_test[,.SD[sample(x=.N, size = round((1124/(length(unique(data$pred))))/2,0),replace=T)],by = pred]
  if(nrow(data_train)<562){data_train<-rbind(data_train,data_train2[sample(.N,(1124/2)-nrow(data_train))])}
  data_test<-subset(data, !(resp %in% data_train$resp))
  rm(data_train2)
 
  fit_train<-lm(data=data_train,resp~pred)
  data_test$fitted<-predict.lm(fit_train,data_test)
  data_test$residuals<-data_test$fitted-data_test$resp
  fit_test<-lm(data=data_test,fitted~resp)
  summaryfit_test<-summary(fit_test)
  
  RMSEi<-rmse(data_test$residuals)
  RMSEvec[j]<-RMSEi
  
  fit.lm_train<-lm(data=data_train,resp~pred.lm)
  data_test$fitted.lm<-predict.lm(fit.lm_train,data_test)
  data_test$residuals.lm<-data_test$fitted.lm-data_test$resp
  fit.lm_test<-lm(data=data_test,fitted.lm~resp)
  summaryfit.lm_test<-summary(fit_test)
  
  RMSEi.lm<-rmse(data_test$residuals.lm)
  RMSEvec.lm[j]<-RMSEi.lm
  
  }
  
  mods[i,]$RMSE<-rmse(fit$residuals)
  mods[i,]$Pred_RMSE<-mean(RMSEvec)
  mods[i,]$discretization<-i
  mods[i,]$RMSE.lm<-rmse(fit.lm$residuals)
  mods[i,]$Pred_RMSE.lm<-mean(RMSEvec.lm)
}

fwrite(mods,"D:/LimnologicalDiscretization/DATA_chlxModels_v4.csv")
mods<-fread("D:/LimnologicalDiscretization/DATA_chlxModels_v4.csv",stringsAsFactors = T)

####FIG 1: 
a<-ggplot()+
  geom_point(data=NLA, aes(y=(CHLX_RESULT),x=(INDEX_SITE_DEPTH)),colour="lightgrey",shape=16)+
  geom_line(data=NLA.m[#discretization!="V7"&
    discretization!="256 Discrete Groups"&
      discretization!="Continuous Gradient"
    ],aes(y=(exp(CHLX_RESULT.mod)),
          x=(INDEX_SITE_DEPTH.jit),
          group=as.factor(CHLX_RESULT.mod),
          colour=discretization),
    size=2)+
  scale_color_manual(values=c(rgb(127,191,123,maxColorValue=225),
                              rgb(127,191,123,maxColorValue=245),
                              rgb(127,191,123,maxColorValue=285),
                              rgb(127,191,123,maxColorValue=365),
                              rgb(127,191,123,maxColorValue=525),
                              rgb(127,191,123,maxColorValue=845),
                              rgb(127,191,123,maxColorValue=1485)))+
  geom_smooth(data=NLA.m[#discretization!="V7"&
    discretization=="Continuous Gradient"
    ],aes(y=(exp(CHLX_RESULT.mod)),
          x=(INDEX_SITE_DEPTH.jit)),
    size=2,method="lm",colour="#af8dc3")+
  scale_y_continuous(trans="log10",labels = function(x) format(x, scientific = FALSE))+
  scale_x_continuous(trans="log10")+
  facet_wrap(~discretization,nrow=2)+
  labs(y="Chla mg/m^3",x="Lake Depth (m)",colour=NULL)+
  theme_bw()+
  theme(legend.position = "none")

b<-ggplot()+
  geom_segment(data=mods[ParNo==2|
                           ParNo==4|
                           ParNo==8|
                           ParNo==16|
                           ParNo==32|
                           ParNo==64|
                           ParNo==128],aes(x=ParNo,xend=ParNo,y=RMSE+0.01,yend=RMSE-0.01),colour="grey",size=1)+
  geom_point(data=mods,aes(y=RMSE,x=ParNo,colour=..x..))+
  geom_smooth(data=mods,aes(y=RMSE,x=ParNo,colour=..x..),span=.1,se=F
  )+
  scale_color_continuous(low=rgb(127,191,123,maxColorValue=225),
                         high=rgb(127,191,123,maxColorValue=1485))+
  geom_smooth(data=mods,aes(y=mean(RMSE.lm),x=ParNo),colour="#762a83",linetype="dashed"#,se=F
  )+
  scale_x_continuous(limits=c(0,128))+
  labs(y="RMSE (log(Chla mg/m^3))",x="Number of Discrete Groups",colour=NULL)+
  theme_bw()+
  annotate("text", label = "Continuous Gradient", x = 90, y = 1.392,colour="#762a83") +
  annotate("text", label = "Discrete Groups", x = 90, y = 1.3,colour=rgb(127,191,123,maxColorValue=365)) +
  theme(legend.position = "none")

c<-ggplot()+
  geom_segment(data=mods[ParNo==2|
                           ParNo==4|
                           ParNo==8|
                           ParNo==16|
                           ParNo==32|
                           ParNo==64|
                           ParNo==128],aes(x=ParNo,xend=ParNo,y=Pred_RMSE+0.01,yend=Pred_RMSE-0.01),
               colour="darkgrey",size=1)+
  geom_point(data=mods,aes(y=Pred_RMSE,x=ParNo,colour=..x..))+
  geom_smooth(data=mods,aes(y=Pred_RMSE,x=ParNo,colour=..x..),span=.1,se=F
              )+
  scale_color_continuous(low=rgb(127,191,123,maxColorValue=225),
                         high=rgb(127,191,123,maxColorValue=1485))+
  geom_smooth(data=mods,aes(y=mean(Pred_RMSE.lm),x=ParNo),colour="#762a83",linetype="dashed"
              )+
  scale_x_continuous(limits=c(0,128))+
  scale_y_continuous(limits=c(1.36,1.58))+
  labs(y="Cross Validation RMSE (log(Chla mg/m^3))",x="Number of Discrete Groups",colour=NULL)+
  theme_bw()+
  annotate("text", label = "Continuous Gradient", x = 90, y = 1.4,colour="#762a83") +
  annotate("text", label = "Discrete Groups", x = 90, y = 1.542,colour=rgb(127,191,123,maxColorValue=365)) +
  theme(legend.position = "none")

tiff("D:/LimnologicalDiscretization/FIG2.tiff", units="in", width=8, height=7, res=300)
grid.arrange(a,b,c,layout_matrix=rbind(c(1),
                                      c(2,3)))
dev.off()

####Radar Plot for conceptual figure####
set.seed(154)
data <- as.data.frame(matrix((sample(c(1:1000),6, replace=TRUE)/1000), ncol=3))
colnames(data) <- c("Objectivity","Predictivity","Stability")
rownames(data) <- c("    Continuous 
    Gradient",
"    Discrete 
    Gradient")

# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each variable to show on the plot!
data <- rbind(c(0,0,0), c(1,1,1) , data)

# plot with default options:
radarchart(data)

# Color vector
colors_border=c( rgb(175,141,195,255,maxColorValue=255), rgb(127,191,123,255,maxColorValue=255) )
colors_in=c( rgb(175,141,195,100,maxColorValue=255), rgb(127,191,123,100,maxColorValue=255))
# plot with default options:
radarchart( data  , axistype=1 , 
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", caxislabels=seq(0,1,0.25), cglwd=1.2,
            #custom labels
            vlcex=1.2 
)

# Add a legend
legend(x=0.7, y=1, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "black", cex=1.2, pt.cex=3)

tiff("D:/LimnologicalDiscretization/FIG1.tiff", units="in", width=8, height=7, res=300)
radarchart( data  , axistype=1 , 
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", caxislabels=seq(0,1,0.25), cglwd=1.2,
            #custom labels
            vlcex=1.2 
)

# Add a legend
legend(x=0.7, y=1, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "black", cex=1.2, pt.cex=3)
dev.off()