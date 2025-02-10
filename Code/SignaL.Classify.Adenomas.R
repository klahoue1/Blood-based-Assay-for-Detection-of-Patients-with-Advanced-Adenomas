library(openxlsx)
library(plyr)
library(dplyr)
library(data.table)
library(MASS)
library(nleqslv)
library(optimization)
library(ggplot2)
library(e1071)
library(randomForest)
library(caret)



dir="/Data"


#Data on patients
FullMetaData = read.csv(file.path(dir,'FullMetaDataSep21.csv'))
Unblind=read.xlsx(file.path(dir,'Unblind_AA.xlsx'))
colnames(Unblind)[1]='cd.name'


# Extract positions of arms
ChromArmsUn = read.table(file.path(dir,'ChromArms.txt'))[1:22, ]
ChromArms = data.frame(V4 = rep(0, 22))
for (i in 1:22) {
  ChromArms$V4[i] = ChromArmsUn$V4[ChromArmsUn$V1 == paste0('chr', i)]
}
AmpSizes = read.csv(file.path(dir,paste0('FinalAmpsFilteredOct3', 1, '.csv')))
ExtractPosition = function(s) {
  r = as.numeric(strsplit(gsub('[^0-9]', " ", s), " ")[[1]][1])
  return(r)
}

PArms = which(sapply(AmpSizes$V2, ExtractPosition) < ChromArms$V4[1])
QArms = which(sapply(AmpSizes$V2, ExtractPosition) > ChromArms$V4[1])



#For chromosome 1 arms P and Q, identify Signatures associated with long amplicons(InterpretFactor=-1), short amplicons (InterpretFactor=1),
# and neutral InterpretFactor=0)


WP = readRDS(file.path(dir,paste0('NonNegativeMatrixP', 1)))
Norms = (apply(WP, 2, sum))
Normalize = matrix(rep(Norms, nrow(WP)), nrow = nrow(WP), byrow = TRUE)
WP = WP / Normalize
MeanLengths = AmpSizes$V3[PArms] %*% WP
#Interpret signature for P arm
tail=2/3
InterpretFactorPArm = as.numeric(MeanLengths > quantile(MeanLengths,tail
)) - as.numeric(MeanLengths < quantile(MeanLengths,1-tail))

DistributionP=cbind(data.frame(Size=AmpSizes$V3[PArms]),WP)
colnames(DistributionP)[2:ncol(DistributionP)]=paste0('Sig',1:(ncol(DistributionP)-1))
SumDistP=ddply(DistributionP,.(Size),numcolwise(mean))
CorsSigP=(as.numeric(cor(SumDistP)[1,2:ncol(SumDistP)]))



WQ = readRDS(file.path(dir,paste0('NonNegativeMatrixQ', 1)))
Norms = (apply(WQ, 2, sum))
Normalize = matrix(rep(Norms, nrow(WQ)), nrow = nrow(WQ), byrow = TRUE)
WQ = WQ / Normalize
MeanLengths = AmpSizes$V3[QArms] %*% WQ
#Interpret signature for Q arm
InterpretFactorQArm = as.numeric(MeanLengths > quantile(MeanLengths,tail
)) - as.numeric(MeanLengths < quantile(MeanLengths, 1-tail))


DistributionQ=cbind(data.frame(Size=AmpSizes$V3[QArms]),WQ)
colnames(DistributionQ)[2:ncol(DistributionQ)]=paste0('Sig',1:(ncol(DistributionQ)-1))
SumDistQ=ddply(DistributionQ,.(Size),numcolwise(mean))
CorsSigQ=(as.numeric(cor(SumDistQ)[1,2:ncol(SumDistQ)]))










dimPArm = 20
dimQArm = 20



NMFTrainNormalPArm = read.csv(file.path(dir,paste0('NMFScoresNormalsTrain', 1)))[, c(1, 2:(dimPArm +
                                                                                     1))]
NMFTrainNormalQArm = read.csv(file.path(dir,paste0('NMFScoresNormalsTrain', 1)))[, c(1, (dimPArm +
                                                                                   2):(dimPArm + dimQArm + 1))]

NMFTrainCancerPArm = read.csv(file.path(dir,paste0('NMFScoresCancersTrain', 1)))[, c(1, 2:(dimPArm +
                                                                                     1))]

NMFTrainCancerQArm = read.csv(file.path(dir,paste0('NMFScoresCancersTrain', 1)))[, c(1, (dimPArm +
                                                                                   2):(dimPArm + dimQArm + 1))]


NMFTestNormalPArm = read.csv(file.path(dir,paste0('NMFScoresNormalsTest', 1)))[, c(1, 2:(dimPArm +
                                                                                   1))]

NMFTestNormalQArm = read.csv(file.path(dir,paste0('NMFScoresNormalsTest', 1)))[, c(1, (dimPArm +
                                                                                 2):(dimPArm + dimQArm + 1))]

NMFTestCancerPArm = read.csv(file.path(dir,paste0('NMFScoresCancersTest', 1)))[, c(1, 2:(dimPArm +
                                                                                   1))]

NMFTestCancerQArm = read.csv(file.path(dir,paste0('NMFScoresCancersTest', 1)))[, c(1, (dimPArm +
                                                                                 2):(dimPArm + dimQArm + 1))]
#Arms not present
SeqQ = c(13, 14, 15, 21, 22)
SeqP = c(-1)


#Repeating same procedure for other chromosomes


for (nchr in 2:22)
  
  
{
  print(paste('This is chromosome', nchr))
  
  
  AmpSizes = read.csv(file.path(dir,paste0('FinalAmpsFilteredOct3', nchr, '.csv')))
  
  
  
  
  if (!(nchr %in% union(SeqP, SeqQ)))
  {
    PArms = which(sapply(AmpSizes$V2, ExtractPosition) < ChromArms$V4[nchr])
    
    QArms = which(sapply(AmpSizes$V2, ExtractPosition) > ChromArms$V4[nchr])
    
  }
  
  if (nchr %in% SeqP)
  {
    PArms = which(sapply(AmpSizes$V2, ExtractPosition) < ChromArms$V4[nchr])
    QArms = PArms
  }
  
  if (nchr %in% SeqQ)
  {
    PArms = which(sapply(AmpSizes$V2, ExtractPosition) > ChromArms$V4[nchr])
    QArms = PArms
  }
  
  
  WP = readRDS(file.path(dir,paste0('NonNegativeMatrixP', nchr)))
  Norms = (apply(WP, 2, sum))
  Normalize = matrix(rep(Norms, nrow(WP)), nrow = nrow(WP), byrow = TRUE)
  WP = WP / Normalize
  MeanLengths = AmpSizes$V3[PArms] %*% WP
  InterpretFactorPArm = c(
    InterpretFactorPArm,
    as.numeric(MeanLengths > quantile(MeanLengths, tail)) - as.numeric(MeanLengths <
                                                                         quantile(MeanLengths, 1-tail))
  )
  
  DistributionP=cbind(data.frame(Size=AmpSizes$V3[PArms]),WP)
  colnames(DistributionP)[2:ncol(DistributionP)]=paste0('Sig',1:(ncol(DistributionP)-1))
  SumDistP=ddply(DistributionP,.(Size),numcolwise(mean))
  CorsSigP=c(CorsSigP,(as.numeric(cor(SumDistP)[1,2:ncol(SumDistP)])))
  
  
  WQ = readRDS(file.path(dir,paste0('NonNegativeMatrixQ', nchr)))
  Norms = (apply(WQ, 2, sum))
  Normalize = matrix(rep(Norms, nrow(WQ)), nrow = nrow(WQ), byrow = TRUE)
  WQ = WQ / Normalize
  MeanLengths = AmpSizes$V3[QArms] %*% WQ
  InterpretFactorQArm = c(
    InterpretFactorQArm,
    as.numeric(MeanLengths > quantile(MeanLengths, tail)) - as.numeric(MeanLengths <
                                                                         quantile(MeanLengths, 1-tail))
  )
  
  DistributionQ=cbind(data.frame(Size=AmpSizes$V3[QArms]),WQ)
  colnames(DistributionQ)[2:ncol(DistributionQ)]=paste0('Sig',1:(ncol(DistributionQ)-1))
  SumDistQ=ddply(DistributionQ,.(Size),numcolwise(mean))
  CorsSigQ=c(CorsSigQ,(as.numeric(cor(SumDistQ)[1,2:ncol(SumDistQ)])))
  
  
  
  NMFTrainNormalPArmInter = read.csv(file.path(dir,paste0('NMFScoresNormalsTrain', nchr)))[, c(1, 2:(dimPArm +
                                                                                               1))]
  
  NMFTrainNormalQArmInter = read.csv(file.path(dir,paste0('NMFScoresNormalsTrain', nchr)))[, c(1, (dimPArm +
                                                                                             2):(dimPArm + dimQArm + 1))]
  
  NMFTrainCancerPArmInter = read.csv(file.path(dir,paste0('NMFScoresCancersTrain', nchr)))[, c(1, 2:(dimPArm +
                                                                                               1))]
  
  NMFTrainCancerQArmInter = read.csv(file.path(dir,paste0('NMFScoresCancersTrain', nchr)))[, c(1, (dimPArm +
                                                                                             2):(dimPArm + dimQArm + 1))]
  
  
  NMFTestNormalPArmInter = read.csv(file.path(dir,paste0('NMFScoresNormalsTest', nchr)))[, c(1, 2:(dimPArm +
                                                                                             1))]
  
  NMFTestNormalQArmInter = read.csv(file.path(dir,paste0('NMFScoresNormalsTest', nchr)))[, c(1, (dimPArm +
                                                                                           2):(dimPArm + dimQArm + 1))]
  
  NMFTestCancerPArmInter = read.csv(file.path(dir,paste0('NMFScoresCancersTest', nchr)))[, c(1, 2:(dimPArm +
                                                                                             1))]
  
  NMFTestCancerQArmInter = read.csv(file.path(dir,paste0('NMFScoresCancersTest', nchr)))[, c(1, (dimPArm +
                                                                                           2):(dimPArm + dimQArm + 1))]
  
  
  
  NMFTrainNormalPArm = join(NMFTrainNormalPArm, NMFTrainNormalPArmInter, by =
                              'Name0')
  NMFTrainNormalQArm = join(NMFTrainNormalQArm, NMFTrainNormalQArmInter, by =
                              'Name0')
  
  print(dim(NMFTrainCancerPArm))
  print(dim(NMFTrainCancerQArm))
  
  NMFTrainCancerPArm = join(NMFTrainCancerPArm, NMFTrainCancerPArmInter, by =
                              'Name1')
  NMFTrainCancerQArm = join(NMFTrainCancerQArm, NMFTrainCancerQArmInter, by =
                              'Name1')
  
  NMFTestNormalPArm = join(NMFTestNormalPArm, NMFTestNormalPArmInter, by =
                             'Name0')
  NMFTestNormalQArm = join(NMFTestNormalQArm, NMFTestNormalQArmInter, by =
                             'Name0')
  
  NMFTestCancerPArm = join(NMFTestCancerPArm, NMFTestCancerPArmInter, by =
                             'Name1')
  NMFTestCancerQArm = join(NMFTestCancerQArm, NMFTestCancerQArmInter, by =
                             'Name1')
  
  
  
  
  
  
}

colnames(NMFTrainNormalPArm)[2:ncol(NMFTrainNormalPArm)] = paste0('score', seq(1, ncol(NMFTrainNormalPArm) -
                                                                                 1))
colnames(NMFTrainNormalQArm)[2:ncol(NMFTrainNormalQArm)] = paste0('score', seq(1, ncol(NMFTrainNormalQArm) -
                                                                                 1))

colnames(NMFTrainCancerPArm)[2:ncol(NMFTrainCancerPArm)] = paste0('score', seq(1, ncol(NMFTrainCancerPArm) -
                                                                                 1))
colnames(NMFTrainCancerQArm)[2:ncol(NMFTrainCancerQArm)] = paste0('score', seq(1, ncol(NMFTrainCancerQArm) -
                                                                                 1))

colnames(NMFTestNormalPArm)[2:ncol(NMFTestNormalPArm)] = paste0('score', seq(1, ncol(NMFTestNormalPArm) -
                                                                               1))
colnames(NMFTestNormalQArm)[2:ncol(NMFTestNormalQArm)] = paste0('score', seq(1, ncol(NMFTestNormalQArm) -
                                                                               1))

colnames(NMFTestCancerPArm)[2:ncol(NMFTestCancerPArm)] = paste0('score', seq(1, ncol(NMFTestCancerPArm) -
                                                                               1))
colnames(NMFTestCancerQArm)[2:ncol(NMFTestCancerQArm)] = paste0('score', seq(1, ncol(NMFTestCancerQArm) -
                                                                               1))

NMFTrainNormal = join(NMFTrainNormalPArm, NMFTrainNormalQArm, by = 'Name0')
NMFTrainCancer = join(NMFTrainCancerPArm, NMFTrainCancerQArm, by = 'Name1')

NMFTestNormal = join(NMFTestNormalPArm, NMFTestNormalQArm, by = 'Name0')
NMFTestCancer = join(NMFTestCancerPArm, NMFTestCancerQArm, by = 'Name1')


DNormalTrain = NMFTrainNormal[, 2:ncol(NMFTrainNormal)]
DCancerTrain = NMFTrainCancer[, 2:ncol(NMFTrainCancer)]


DNormalTest = NMFTestNormal[, 2:ncol(NMFTestNormal)]
DCancerTest = NMFTestCancer[, 2:ncol(NMFTestCancer)]

InterpretFactor = c(InterpretFactorPArm, InterpretFactorQArm)
CorsSig=c(CorsSigP,CorsSigQ)

I1 = intersect(which(apply(DNormalTrain, 2, median) < apply(DCancerTrain, 2, median)), which(InterpretFactor ==
                                                                                               -1))
I2 = intersect(which(apply(DNormalTrain, 2, median) > apply(DCancerTrain, 2, median)), which(InterpretFactor ==
                                                                                               1))

#Keeping only the consistent signatures
I = union(union(I1, I2), which(InterpretFactor==0))
#I=which(abs(CorsSig)>quantile(abs(CorsSig),0.))


DNormalTrain = DNormalTrain[, I]
DCancerTrain = DCancerTrain[, I]


DNormalTest = DNormalTest[, I]
DCancerTest = DCancerTest[, I]





#Train SVM

S = 0



ValidatedNormals = rep(0, nrow(DNormalTrain))
ValidatedCancers = rep(0, nrow(DCancerTrain))


for (r in 1:10)
  
{
  
  set.seed(10*r)
  foldsNormal = createFolds(1:nrow(DNormalTrain), k = 10, returnTrain = TRUE)
  foldsCancer = createFolds(1:nrow(DCancerTrain), k = 10, returnTrain = TRUE)
  
  for (s in 1:10)
    
  {
    i1 = sample(1:nrow(DNormalTrain), nrow(DNormalTrain) / 2)
    
    
    
    
    SVM.Model = svm(as.matrix(rbind(DNormalTrain[foldsNormal[[s]], ], DCancerTrain[foldsCancer[[s]], ])),
                    as.factor(c(rep(
                      0, nrow(DNormalTrain[foldsNormal[[s]], ])
                    ), rep(
                      1, nrow(DCancerTrain[foldsCancer[[s]], ])
                    ))),
                    cost = 1,
                    probability = TRUE)
    
    
    predCancerTest = attr(predict(SVM.Model, as.matrix(DCancerTest), probability =
                                    TRUE),
                          'probabilities')[, 2]
    predNormalTest = attr(predict(SVM.Model, as.matrix(DNormalTest), probability =
                                    TRUE),
                          'probabilities')[, 2]
    ValidatedCancers[-foldsCancer[[s]]] =  ValidatedCancers[-foldsCancer[[s]]]+attr(predict(SVM.Model, as.matrix(DCancerTrain[-foldsCancer[[s]], ]), probability =
                                                                                              TRUE),
                                                                                    'probabilities')[, 2]
    ValidatedNormals[-foldsNormal[[s]]] =ValidatedNormals[-foldsNormal[[s]]]+attr(predict(SVM.Model, as.matrix(DNormalTrain[-foldsNormal[[s]], ]), probability =
                                                                                            TRUE),
                                                                                  'probabilities')[, 2]
  }
  
  
}

ValidatedNormals=ValidatedNormals/10
ValidatedCancers=ValidatedCancers/10
































#Train SVM
SVM.Model = svm(as.matrix(rbind(DNormalTrain, DCancerTrain)),
                as.factor(c(rep(
                  0, nrow(DNormalTrain)
                ), rep(
                  1, nrow(DCancerTrain)
                ))),
                cost = 1,
                probability = TRUE)


predCancerTest = attr(predict(SVM.Model, as.matrix(DCancerTest), probability =
                                TRUE),
                      'probabilities')[, 2]
predNormalTest = attr(predict(SVM.Model, as.matrix(DNormalTest), probability =
                                TRUE),
                      'probabilities')[, 2]
predCancerTrain = attr(predict(SVM.Model, as.matrix(DCancerTrain), probability =
                                 TRUE),
                       'probabilities')[, 2]
predNormalTrain = attr(predict(SVM.Model, as.matrix(DNormalTrain), probability =
                                 TRUE),
                       'probabilities')[, 2]


ExtractBVNames=function(s){ind=regexpr('CRC',s);indCRC=regexpr('CRC',s);return(substr(s,indCRC,indCRC+6))}

#Extract Adenomas results on desired patients
AdenomasTot = read.xlsx(file.path(dir,'AdenomasGrant.xlsx'), sheet = 1)
AdenomasTot = filter(AdenomasTot, !is.na(Patient.type))
U = union(grep('10397', AdenomasTot$cd.name),
          grep('10398', AdenomasTot$cd.name))
Adenomas = AdenomasTot[-U, ]
AdenomasFirst=AdenomasTot[U, ]

AdenomasFirst$bv.name=sapply(AdenomasFirst$bv.name,ExtractBVNames)


#Table to match IDs
AdenomasMatchIds = read.xlsx(file.path(dir,'BrendaInfos.xlsx')) %>% filter(.,Previous.Hopkins.ID !='not applicable')
AdenomasMatchIds$Previous.Hopkins.ID=sapply(AdenomasMatchIds$Previous.Hopkins.ID,ExtractBVNames)



#Generating data frame with cancer test scores
Cancers = union((filter(FullMetaData, cancer == TRUE))$douville_name,NMFTestCancer$Name1)
CancerTestSignaL = data.frame(cd.name = NMFTestCancer$Name1[NMFTestCancer$Name1 %in% Cancers], SignaL.Score =
                                predCancerTest[NMFTestCancer$Name1 %in% Cancers])



#Generating data frame with normal test scores
Normals = (filter(FullMetaData, cancer == FALSE))$douville_name
NormalTestSignaL = data.frame(cd.name = NMFTestNormal$Name0[NMFTestNormal$Name0 %in% Normals], SignaL.Score =
                                predNormalTest[NMFTestNormal$Name0 %in% Normals])


CancerSignaLExact=CancerTestSignaL[grep('E',CancerTestSignaL$cd.name),]
CancerSignaLExact=join(CancerSignaLExact,Unblind,by='cd.name')

q1=quantile(NormalTestSignaL$SignaL.Score,0.99)
qPrime=quantile(NormalTestSignaL$SignaL.Score,0.985)
q2=quantile(NormalTestSignaL$SignaL.Score,0.95)
q3=quantile(NormalTestSignaL$SignaL.Score,0.9)

CancerSignaLExact$Pass.99Spec=CancerSignaLExact$SignaL.Score>q1
CancerSignaLExact$Pass.985Spec=CancerSignaLExact$SignaL.Score>qPrime
CancerSignaLExact$Pass.95Spec=CancerSignaLExact$SignaL.Score>q2
CancerSignaLExact$Pass.9Spec=CancerSignaLExact$SignaL.Score>q3




#Generating data frame with Adenomas scores
AdenomasResultsOnBlocks = data.frame(cd.name = NMFTestCancer$Name1[NMFTestCancer$Name1 %in% Adenomas$cd.name],
                                     SignaL.Score = predCancerTest[NMFTestCancer$Name1 %in% Adenomas$cd.name]) %>% join(Adenomas, ., by =
                                                                                                                          'cd.name')


AdenomasResultsOnBlocksFirst = data.frame(cd.name = NMFTestCancer$Name1[NMFTestCancer$Name1 %in% AdenomasFirst$cd.name],
                                          SignaL.Score = predCancerTest[NMFTestCancer$Name1 %in% AdenomasFirst$cd.name]) %>% join(AdenomasFirst, ., by =
                                                                                                                                    'cd.name')

AdenomasResultsOnBlocks$Previous.SignaL.Score=-1
AdenomasResultsOnBlocksFirst$Schoen_num=NA
for (p in 1:nrow(AdenomasResultsOnBlocks))
{
  
  IndMatch=which(AdenomasResultsOnBlocks$`Schoen.#`[p]==AdenomasMatchIds$`Sample.ID.(May.2021)`)
  w=which(AdenomasMatchIds$Previous.Hopkins.ID[IndMatch]==AdenomasResultsOnBlocksFirst$bv.name)
  if (length(IndMatch)>0 & length(w)>0)
  {
    print(w)
    AdenomasResultsOnBlocksFirst$Schoen_num[w]=AdenomasMatchIds$`Sample.ID.(May.2021)`[IndMatch]
    AdenomasResultsOnBlocks$Previous.SignaL.Score[p]=AdenomasResultsOnBlocksFirst$SignaL.Score[w]
  }
}