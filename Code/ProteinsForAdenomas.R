library(openxlsx)
library(plyr)
library(dplyr)
library(data.table)
library(randomForest)
library(mixtools)
library(mvtnorm)
library(adabag)
library(MASS)
library(ggVennDiagram)
library(glmnet)
library(pROC)
library(ggrepel)

dir="/Data"







LimitsProt=read.xlsx(file.path(dir,'ProtAdenom1.xlsx'),rows=seq(2,4))
ProtAdenomas1=read.xlsx(file.path(dir,'ProtAdenom1.xlsx'),startRow=8)
ProtAdenomas1=ProtAdenomas1[grep('CRC',ProtAdenomas1$X2),]
ExtractBVNames=function(s){ind=regexpr('CRC',s);indCRC=regexpr('CRC',s);return(substr(s,indCRC,indCRC+13))}
ProtAdenomas1$X2=sapply(ProtAdenomas1$X2,ExtractBVNames)
ProteinNames=colnames(LimitsProt)[colnames(LimitsProt) %in% colnames(ProtAdenomas1)]
LimitsProt=LimitsProt[,ProteinNames]
ProtAdenomas1=ProtAdenomas1[,c('X2',colnames(LimitsProt))]
for (c in 2:ncol(ProtAdenomas1))
{Lower=grep("<",ProtAdenomas1[,c]);
ProtAdenomas1[Lower,c]=LimitsProt[1,c-1];
Upper=grep(">",ProtAdenomas1[,c]);
ProtAdenomas1[Upper,c]=LimitsProt[2,c-1];
ProtAdenomas1[,c]=as.numeric(ProtAdenomas1[,c])
}
colnames(ProtAdenomas1)[1]='bv.name'




LimitsProt=read.xlsx(file.path(dir,'ProtAdenom2.xlsx'),rows=seq(2,4))
ProtAdenomas2=read.xlsx(file.path(dir,'ProtAdenom2.xlsx'),startRow=8)
ProtAdenomas2=ProtAdenomas2[grep('CRC',ProtAdenomas2$X2),]
ExtractBVNames=function(s){ind=regexpr('CRC',s);indCRC=regexpr('CRC',s);return(substr(s,indCRC,indCRC+13))}
ProtAdenomas2$X2=sapply(ProtAdenomas2$X2,ExtractBVNames)
ProteinNames=colnames(LimitsProt)[colnames(LimitsProt) %in% colnames(ProtAdenomas2)]
LimitsProt=LimitsProt[,ProteinNames]
ProtAdenomas2=ProtAdenomas2[,c('X2',colnames(LimitsProt))]
for (c in 2:ncol(ProtAdenomas2))
{Lower=grep("<",ProtAdenomas2[,c]);
ProtAdenomas2[Lower,c]=LimitsProt[1,c-1];
Upper=grep(">",ProtAdenomas2[,c]);
ProtAdenomas2[Upper,c]=LimitsProt[2,c-1];
ProtAdenomas2[,c]=as.numeric(ProtAdenomas2[,c])
}
colnames(ProtAdenomas2)[1]='bv.name'


ProtAdenomas=rbind(ProtAdenomas1,ProtAdenomas2)





for (indFile in 2:4)
  
  
{
  

  
  LimitsProt=read.xlsx(file.path(dir,paste0('ProtAdenom',2*indFile-1,'.xlsx')),rows=seq(2,4))
  ProtAdenomas1=read.xlsx(file.path(dir,paste0('ProtAdenom',2*indFile-1,'.xlsx')),startRow=8)
  ProtAdenomas1=ProtAdenomas1[grep('CRC',ProtAdenomas1$X2),]
  ExtractBVNames=function(s){ind=regexpr('CRC',s);indCRC=regexpr('CRC',s);return(substr(s,indCRC,indCRC+13))}
  ProtAdenomas1$X2=sapply(ProtAdenomas1$X2,ExtractBVNames)
  ProteinNames=colnames(LimitsProt)[colnames(LimitsProt) %in% colnames(ProtAdenomas1)]
  LimitsProt=LimitsProt[,ProteinNames]
  LimitsProt=as.matrix(LimitsProt)
  ProtAdenomas1=ProtAdenomas1[,c('X2',ProteinNames)]
  for (c in 2:ncol(ProtAdenomas1))
  {Lower=grep("<",ProtAdenomas1[,c]);
  ProtAdenomas1[Lower,c]=LimitsProt[1,c-1];
  Upper=grep(">",ProtAdenomas1[,c]);
  ProtAdenomas1[Upper,c]=LimitsProt[2,c-1];
  ProtAdenomas1[,c]=as.numeric(ProtAdenomas1[,c])
  }
  colnames(ProtAdenomas1)[1]='bv.name'
  
  
  
  LimitsProt=read.xlsx(file.path(dir,paste0('ProtAdenom',2*indFile,'.xlsx')),rows=seq(2,4))
  ProtAdenomas2=read.xlsx(file.path(dir,paste0('ProtAdenom',2*indFile,'.xlsx')),startRow=8)
  ProtAdenomas2=ProtAdenomas2[grep('CRC',ProtAdenomas2$X2),]
  ExtractBVNames=function(s){ind=regexpr('CRC',s);indCRC=regexpr('CRC',s);return(substr(s,indCRC,indCRC+13))}
  ProtAdenomas2$X2=sapply(ProtAdenomas2$X2,ExtractBVNames)
  ProteinNames=colnames(LimitsProt)[colnames(LimitsProt) %in% colnames(ProtAdenomas2)]
  LimitsProt=LimitsProt[,ProteinNames]
  LimitsProt=as.matrix(LimitsProt)
  ProtAdenomas2=ProtAdenomas2[,c('X2',ProteinNames)]
  for (c in 2:ncol(ProtAdenomas2))
  {Lower=grep("<",ProtAdenomas2[,c]);
  ProtAdenomas2[Lower,c]=LimitsProt[1,c-1];
  Upper=grep(">",ProtAdenomas2[,c]);
  ProtAdenomas2[Upper,c]=LimitsProt[2,c-1];
  ProtAdenomas2[,c]=as.numeric(ProtAdenomas2[,c])
  }
  colnames(ProtAdenomas2)[1]='bv.name'
  
  
  
ProtAdenomas=join(ProtAdenomas,rbind(ProtAdenomas1,ProtAdenomas2),by='bv.name')
  
  
} 
  

  
  
NormalizePVals=function(x){r=(length(x)-rank(x)+1)/length(x);return(r)}
logCentered=function(x){return(abs(log(x/(1-x))))}



ProtAdenomasNormalized=apply(ProtAdenomas[,2:ncol(ProtAdenomas)],2,NormalizePVals)
ScoresProtAdenomas=apply(-log(ProtAdenomasNormalized),1,sum)

ProtAdenomasScores=data.frame(bv.name=ProtAdenomas$bv.name,ScoreForProt=ScoresProtAdenomas)



  

Adenomas=read.csv(file.path(dir,'AdenomasSummarySignaL.csv'))


Adenomas=filter(Adenomas,!is.na(Patient.type))
U=union(grep('10397',Adenomas$cd.name),grep('10398',Adenomas$cd.name))

Adenomas$Patient.type=as.character(Adenomas$Patient.type)
Adenomas$Signal=Adenomas$SignaL.Score


ExtractBVNamesFusion=function(s){ind=regexpr('CRC',s);indCRC=regexpr('CRC',s);return(substr(s,indCRC,indCRC+6))}
ProtAdenomasScores$bv.name=sapply(ProtAdenomasScores$bv.name,ExtractBVNamesFusion)
ProtAdenomasScores=ddply(ProtAdenomasScores,.(bv.name),summarise,ScoreForProt=mean(ScoreForProt))
Adenomas$bv.name=sapply(Adenomas$bv.name,ExtractBVNamesFusion)


TotScores=join(ProtAdenomasScores,Adenomas,by='bv.name')
TotScores=filter(TotScores,!is.na(Patient.type))
TotScores=filter(TotScores,!is.na(Signal))
sort(TotScores$ScoreForProt[TotScores$Patient.type=='case'])
sort(TotScores$ScoreForProt[TotScores$Patient.type=='control'])


PerfProt=sum(TotScores$ScoreForProt[TotScores$Patient.type=='case']>quantile(TotScores$ScoreForProt[TotScores$Patient.type=='control'],0.87))/length(TotScores$ScoreForProt[TotScores$Patient.type=='case'])
PerfSignaL=sum(TotScores$Signal[TotScores$Patient.type=='case']>quantile(TotScores$Signal[TotScores$Patient.type=='control'],0.87))/length(TotScores$Signal[TotScores$Patient.type=='case'])
TotScores$PerfCombined=(-log(NormalizePVals(TotScores$Signal))-log(NormalizePVals(TotScores$ScoreForProt)))


PerfCombined=sum(TotScores$PerfCombined[TotScores$Patient.type=='case']>quantile(TotScores$PerfCombined[TotScores$Patient.type=='control'],0.87))/length(TotScores$PerfCombined[TotScores$Patient.type=='case'])

print(paste('performance for proteins is',PerfProt))
print(paste('performance for SignaL is',PerfSignaL))

print(paste('performance for combined is',PerfCombined))


CatchProt=TotScores$bv.name[(TotScores$ScoreForProt>quantile(TotScores$ScoreForProt[TotScores$Patient.type=='control'],0.87)) & (TotScores$Patient.type=='case')]
CatchSignal=TotScores$bv.name[(TotScores$Signal>quantile(TotScores$Signal[TotScores$Patient.type=='control'],0.87)) & (TotScores$Patient.type=='case')]
CatchCombined=TotScores$bv.name[(TotScores$PerfCombined>quantile(TotScores$PerfCombined[TotScores$Patient.type=='control'],0.87)) & (TotScores$Patient.type=='case')]





x=list(Proteins=as.character(CatchProt),SignaL=as.character(CatchSignal), Combined=as.character(CatchCombined))
ggVennDiagram(x)+ggtitle('Venn diagram for fraction of TP by each marker. Sensitivity is 42.5% at Spec 87.5%')





ProtAdenomas$Status=rep('nn',nrow(ProtAdenomas))
ProtAdenomas$Signal=rep(0,nrow(ProtAdenomas))
for (ll in 1:nrow(TotScores)){inds=grep(TotScores$bv.name[ll],ProtAdenomas$bv.name);print((inds));ProtAdenomas$Status[inds]=TotScores$Patient.type[ll];ProtAdenomas$Signal[inds]=TotScores$Signal[ll]}
ProtAdenomas=filter(ProtAdenomas,Status !='nn')



Proteins.CancerSeek=read.xlsx(file.path(dir,'CancerSeekProt.xlsx'),sheet=6,startRow =3)
Proteins.CancerSeek.Colorectum=filter(Proteins.CancerSeek,Tumor.type=='Colorectum')
Proteins.CancerSeek.Normal=filter(Proteins.CancerSeek,Tumor.type=='Normal')

PCA.Proteins=prcomp(ProtAdenomas[,2:18],scale=TRUE)

PCA.Proteins=data.frame(Status=ProtAdenomas$Status,PCA.Dim.1=PCA.Proteins$x[,1],PCA.Dim.2=PCA.Proteins$x[,11])

PCA.Proteins %>% ggplot(.)+geom_point(aes(x=PCA.Dim.1,y=PCA.Dim.2,col=Status))


colnames(ProtAdenomas)[9]='CA-125'
FindPattern=function(x){return(grep(x,colnames(Proteins.CancerSeek)))}
unlist(sapply(colnames(ProtAdenomas[2:18]),FindPattern))






Col.To.Select=unlist(sapply(colnames(ProtAdenomas[2:18]),FindPattern))
colnames(Proteins.CancerSeek.Colorectum)[Col.To.Select]=colnames(ProtAdenomas)[2:18]
colnames(Proteins.CancerSeek.Normal)[Col.To.Select]=colnames(ProtAdenomas)[2:18]
Proteins.CancerSeek.Colorectum=Proteins.CancerSeek.Colorectum[,c(c(1,2,3),Col.To.Select)]
Proteins.CancerSeek.Normal=Proteins.CancerSeek.Normal[,c(c(1,2,3),Col.To.Select)]





for (j in 4:20)
{
Proteins.CancerSeek.Normal[,j]=as.numeric(gsub("\\*", "", Proteins.CancerSeek.Normal[,j]))
Proteins.CancerSeek.Colorectum[,j]=as.numeric(gsub("\\*", "", Proteins.CancerSeek.Colorectum[,j]))
}







DNormalTrain=Proteins.CancerSeek.Normal[,4:20]
DCancerTrain=Proteins.CancerSeek.Colorectum[,4:20]
DNormalTest=filter(ProtAdenomas,Status=='control')[,2:18]
DCancerTest=filter(ProtAdenomas,Status=='case')[,2:18]





ValidatedNormals = rep(0, nrow(Proteins.CancerSeek.Normal))
ValidatedCancers = rep(0, nrow(Proteins.CancerSeek.Colorectum))


for (r in 1:10)
  
{
  
  set.seed(10*r)
  foldsNormal = createFolds(1:nrow(Proteins.CancerSeek.Normal), k = 10, returnTrain = TRUE)
  foldsCancer = createFolds(1:nrow(Proteins.CancerSeek.Colorectum), k = 10, returnTrain = TRUE)
  
  for (s in 1:10)
    
  {
    
    
    fit=cv.glmnet(as.matrix(rbind(DNormalTrain[foldsNormal[[s]], ], DCancerTrain[foldsCancer[[s]], ])),as.factor(c(rep(0,nrow(DNormalTrain[foldsNormal[[s]], ])),rep(1,nrow(DCancerTrain[foldsCancer[[s]], ])))),family='binomial',alpha=1/2,nfolds=10,lower.limit=rep(0,ncol(DCancerTrain)))
    
    ValidatedNormals[-foldsNormal[[s]]]=ValidatedNormals[-foldsNormal[[s]]]+predict(fit,s='lambda.min',newx=as.matrix(DNormalTrain[-foldsNormal[[s]], ]),type='response')
    ValidatedCancers[-foldsCancer[[s]]]=ValidatedCancers[-foldsCancer[[s]]]+predict(fit,s='lambda.min',newx=as.matrix(DCancerTrain[-foldsCancer[[s]], ]),type='response')
   
    
  }
  
  
}

ValidatedNormals=ValidatedNormals/10
ValidatedCancers=ValidatedCancers/10

Proteins.CancerSeek.Normal$Protein_17_Validated=ValidatedNormals
Proteins.CancerSeek.Colorectum$Protein_17_Validated=ValidatedCancers












fit=cv.glmnet(as.matrix(rbind(((DNormalTrain)),((DCancerTrain)))),as.factor(c(rep(0,nrow(DNormalTrain)),rep(1,nrow(DCancerTrain)))),family='binomial',alpha=1/2,nfolds=10,lower.limit=rep(0,ncol(DCancerTrain)))
predNormalTrain=predict(fit,s='lambda.min',newx=as.matrix(DNormalTrain),type='response')
predCancerTrain=predict(fit,s='lambda.min',newx=as.matrix(DCancerTrain),type='response')
predNormalTest=predict(fit,s='lambda.min',newx=as.matrix(DNormalTest),type='response')
predCancerTest=predict(fit,s='lambda.min',newx=as.matrix(DCancerTest),type='response')
ProtAdenomas$ProtScore=predict(fit,s='lambda.min',newx=as.matrix(ProtAdenomas[,2:18]),type='response')
Proteins.CancerSeek.Normal$Protein_17_Training=as.vector(predict(fit,s='lambda.min',newx=as.matrix(DNormalTrain),type='response'))
Proteins.CancerSeek.Colorectum$Protein_17_Training=as.vector(predict(fit,s='lambda.min',newx=as.matrix(DCancerTrain),type='response'))







predCancerTest=(filter(ProtAdenomas,Status=='case'))$Combined
predNormalTest=(filter(ProtAdenomas,Status=='control'))$Combined


IndSignal=which((filter(ProtAdenomas,Status=='case'))$Signal>quantile((filter(ProtAdenomas,Status=='control'))$Signal,0.95))
IndProt=which((filter(ProtAdenomas,Status=='case'))$ProtScore>quantile((filter(ProtAdenomas,Status=='control'))$ProtScore,0.95))
IndSignalAlone=which((filter(ProtAdenomas,Status=='case'))$Signal>quantile((filter(ProtAdenomas,Status=='control'))$Signal,0.87))
IndProtAlone=which((filter(ProtAdenomas,Status=='case'))$ProtScore>quantile((filter(ProtAdenomas,Status=='control'))$ProtScore,0.87))
print(length(union(IndSignal,IndProt)))

IndSignalC=which((filter(ProtAdenomas,Status=='control'))$Signal>quantile((filter(ProtAdenomas,Status=='control'))$Signal,0.95))
IndProtC=which((filter(ProtAdenomas,Status=='control'))$ProtScore>quantile((filter(ProtAdenomas,Status=='control'))$ProtScore,0.95))
print(length(union(IndSignalC,IndProtC)))





Specificities=seq(0,1,0.01)

ROC=data.frame(FalsePositiveRate=rep(0,length(Specificities)),Sensitivity=rep(0,length(Specificities)))


for (s in 1:length(Specificities))

{
  
spec=Specificities[s]

IndSignalC=which((filter(ProtAdenomas,Status=='control'))$Signal>quantile((filter(ProtAdenomas,Status=='control'))$Signal,spec))
IndProtC=which((filter(ProtAdenomas,Status=='control'))$ProtScore>quantile((filter(ProtAdenomas,Status=='control'))$ProtScore,spec))
ROC[s,1]=(length(union(IndSignalC,IndProtC)))/32

IndSignal=which((filter(ProtAdenomas,Status=='case'))$Signal>quantile((filter(ProtAdenomas,Status=='control'))$Signal,spec))
IndProt=which((filter(ProtAdenomas,Status=='case'))$ProtScore>quantile((filter(ProtAdenomas,Status=='control'))$ProtScore,spec))
ROC[s,2]=(length(union(IndSignal,IndProt)))/40


}









ISorted = sort(
  ProtAdenomas$ProtScore,
  decreasing = TRUE,
  index.return = TRUE
)$ix
ProtAdenomas =  ProtAdenomas[ISorted, ]
pl1 =  ProtAdenomas %>% ggplot(.) + geom_point(aes(x = 1:72, y =
                                                     ProtScore))
pl1 + aes(colour = Status) + xlab('Patients') + ylab('Proteins score') +
  theme(plot.title = element_text(hjust = 0.5))


CatchProt=ProtAdenomas$bv.name[(IndProt)]
CatchSignal=ProtAdenomas$bv.name[(IndSignal)]
CatchCombined=ProtAdenomas$bv.name[union(IndProt,IndSignal)]

x=list(Proteins=as.character(CatchProt),SignaL=as.character(CatchSignal))
ggVennDiagram(x)

SignaLInfo=read.csv('SignaLResults.csv')
pl2 =  SignaLInfo %>% ggplot(.) + geom_point(aes(x = 1:72, y =
                                                     SignaL.Score))
pl2 + aes(colour = Status) + xlab('Patients') + ylab('SignaL score') +
  theme(plot.title = element_text(hjust = 0.5))

ROC$labs=paste0(ROC$FalsePositiveRate,'/',ROC$Sensitivity)
ROC=ROC[seq(101,1,-1),]
ROC %>% ggplot(.)+geom_step(aes(x=FalsePositiveRate,y=Sensitivity),linewidth=1)+geom_line(aes(x=rep(0.125,101),y=seq(0,1-1/101,1/101)),colour='red',size=1)+labs(title = "                ROC curve for SignaL + proteins", x = "False positive rate", y = "Sensitivity",size=2)+scale_y_continuous(breaks=c(seq(0,1,0.1)),limits=c(0,1))+scale_x_continuous(breaks=c(seq(0,1,0.25),0.125),limits=c(0,1))

