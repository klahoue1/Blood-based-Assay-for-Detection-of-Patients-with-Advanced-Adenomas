
library(openxlsx)
library(plyr)
library(dplyr)
library(data.table)
library(MASS)
library(nleqslv)
library(optimization)
library(NMF)
library(nnls)
library(feather)
library(feather)



dir="/Data"




FullMetaData = read.csv(file.path(dir,'FullMetaTraining.csv'))
TrainingNames=FullMetaData$douville_name


#ArmsPosition
ChromArmsUn = read.table(file.path(dir,'ChromArms.txt'))[1:22, ]
ChromArms = data.frame(V4 = rep(0, 22))
for (i in 1:22) {
  ChromArms$V4[i] = ChromArmsUn$V4[ChromArmsUn$V1 == paste0('chr', i)]
}
#For extracting position of amplicon from string with particular format
ExtractPosition = function(s) {
  r = as.numeric(strsplit(gsub('[^0-9]', " ", s), " ")[[1]][1])
  return(r)
}



#Arms not present
SeqQ = c(13, 14, 15, 21, 22)
#Just for checking
SeqP = c(-1)








for (nchr in 1:22)
  
{
  #Data on the chromosome
  AmpPrim = readRDS(file.path(dir,(paste0('FinalSetsFullAmp', nchr, '.rds'))))
  
  #The Cancers
  Cancers = FullMetaData$douville_name[FullMetaData$cancer == TRUE]
  Cancers = intersect(intersect(Cancers[!is.na(Cancers)], TrainingNames),BlockProts)
  
  #The Normals
  Normals = FullMetaData$douville_name[FullMetaData$cancer == FALSE]
  Normals = intersect(Normals[!is.na(Normals)], TrainingNames)
  
  #Split into cancers and normals
  Amp1 = filter(AmpPrim, douville_name %in% Cancers)
  Amp2 = filter(AmpPrim, douville_name %in% Normals)
  
  #Filtered amplicons
  AmpSizes = read.csv(paste0('FinalAmpsFilteredOct3', nchr, '.csv'))
  ColSelected = paste0('amplicon_id_', AmpSizes$X)
  A = which(colnames(Amp1)[2:ncol(Amp1)] %in% ColSelected)
  
  Amp1 = Amp1[, c(1, 1 + A)]
  Amp2 = Amp2[, c(1, 1 + A)]
  
  
  
  
  
  
  
  
  
  
  Amp = rbind(Amp1[, 2:ncol(Amp1)], Amp2[, 2:ncol(Amp2)])
  NormalizeFactor = apply(Amp, 1, sum)
  Amp = cbind(c(rep(TRUE, nrow(Amp1)), rep(FALSE, nrow(Amp2))), Amp)
  colnames(Amp)[1] = 'cancer'
  
  
  
  S0 = data.frame(Name0 = Amp2$douville_name)
  S1 = data.frame(Name1 = Amp1$douville_name)
  
  
  #This is in case we want to partition the chromosome (Not used in this version jump to 140)
  Clust = function(n) {
    stride = floor(n / 1)
    L = list()
    i = 1
    count = 1
    while (i < n) {
      L[[count]] = seq(i, min(i + stride, n))
      count = count + 1
      i = min(i + stride + 1, n)
    }
    if (i == n) {
      L[[count - 1]] = c(L[[count - 1]], n)
    }
    return(L)
  }
  
  
  L = Clust(ncol(Amp) - 1)
  
  
  
  
  
  
  for (chri in 1:length(L))
    
  {
    #AmpTot = Amp[, c(1, 1 + unique(L[[chri]]))]
    
    
    
    
    #PArms and QArms amplicons
    PArms = which(sapply(AmpSizes$V2, ExtractPosition) < ChromArms$V4[nchr])
    QArms = which(sapply(AmpSizes$V2, ExtractPosition) > ChromArms$V4[nchr])
    
    DistLengthP=as.numeric(table(AmpSizes[PArms,]$V3)/sum(table(AmpSizes[PArms,]$V3)))
    DistLengthQ=as.numeric(table(AmpSizes[QArms,]$V3)/sum(table(AmpSizes[QArms,]$V3)))
    
    
    
    
    
    
    #Keep only amplicon counts
    Data = Amp[, grep('amplicon', colnames(Amp))] %>% as.matrix(.)
    
    
    #Keeping only patients with with non zero read proportion >90%. I need to re-write this part, unnecessarily confusing (was meant for debugging)!!!
    # Zero counts proportion
    SumZeros = function(x) {
      return(sum(x == 0))
    }
    P = apply(Data, 1, SumZeros) / ncol(Data)
    Amp = Amp[P < 0.9 , ]
    Data = Data[P < 0.9, ]
    Y = 2 * as.numeric(Amp$cancer) - 1
    
    
    
    
    BV_Name=c(as.character(Amp1$douville_name),as.character(Amp2$douville_name))[P<0.9]
    BV_Normals = BV_Name[Amp$cancer == FALSE]
    BV_Cancers = BV_Name[Amp$cancer == TRUE]
    
    
    
    
    
    
    #Normals
    Data0 = Data[Y == -1, ]
    #Cancers
    Data1 = Data[Y == 1, ]
    
    
    
    #Reweihgting factors
    ProbasLengthsP=((table(AmpSizes[PArms,]$V3)/sum(table(AmpSizes[PArms,]$V3))))
    ProbasLengthsQ=((table(AmpSizes[QArms,]$V3)/sum(table(AmpSizes[QArms,]$V3))))
    LengthsP=as.numeric(names(ProbasLengthsP))
    LengthsQ=as.numeric(names(ProbasLengthsQ))
    ProbasLengthsP=as.numeric(ProbasLengthsP)
    ProbasLengthsQ=as.numeric(ProbasLengthsQ)
    AmpSizesP=AmpSizes[PArms,]
    AmpSizesQ=AmpSizes[QArms,]
    ProbasP=rep(0,nrow(AmpSizesP))
    ProbasQ=rep(0,nrow(AmpSizesQ))
    
    for (lp in 1:length(LengthsP))
    {
    ProbasP[AmpSizesP$V3==LengthsP[lp]]=1/sqrt(ProbasLengthsP[lp])
    }
    
    for (lq in 1:length(LengthsQ))
    {
      ProbasQ[AmpSizesQ$V3==LengthsQ[lq]]=1/sqrt(ProbasLengthsQ[lq])
    }
    
    ProbasP0=matrix(rep(ProbasP,nrow(Data0)),ncol=length(ProbasP),byrow=TRUE)
    ProbasP1=matrix(rep(ProbasP,nrow(Data1)),ncol=length(ProbasP),byrow=TRUE)
    
    ProbasQ0=matrix(rep(ProbasQ,nrow(Data0)),ncol=length(ProbasQ),byrow=TRUE)
    ProbasQ1=matrix(rep(ProbasQ,nrow(Data1)),ncol=length(ProbasQ),byrow=TRUE)
    
    
    
    
    #Splitting into arms and normalizing factors
    Data0Q = Data0[, QArms]*ProbasQ0
    if (nchr %in% SeqQ){Data0P = Data0[, PArms]} else {Data0P = Data0[, PArms]*ProbasP0}
    NormalizeFactor0P = apply(Data0P, 1, sum)
    NormalizeFactor0Q = apply(Data0Q, 1, sum)
    
    
    Data1Q = Data1[, QArms]*ProbasQ1
    if (nchr %in% SeqQ){Data1P = Data1[, PArms]} else {Data1P = Data1[, PArms]*ProbasP1}
    NormalizeFactor1P = apply(Data1P, 1, sum)
    NormalizeFactor1Q = apply(Data1Q, 1, sum)
    
    
    print(paste0('This is chromosome', nchr))
    print(dim(Data0))
    
    dimNMF = 20
    
    
    #NMF on the normalized counts: if one of the arms is not present in the chromosome, we duplicate the present arm.
    if (nchr %in% SeqP)
      
    {
      MP = rbind(Data0P / NormalizeFactor0P, Data1P / NormalizeFactor1P)
      MQ = rbind(Data0P / NormalizeFactor0P, Data1P / NormalizeFactor1P)
      
      
      WP = (nmf(t(MP), dimNMF))@fit@W
      WQ = WP
      
      
      HP = function(b) {
        return((nnls(WP, b)$x))
      }
      HQ = function(b) {
        return((nnls(WQ, b)$x))
      }
      
      
      
    }
    
    
    if (nchr %in% SeqQ)
      
    {
      MP = rbind(Data0Q / NormalizeFactor0Q, Data1Q / NormalizeFactor1Q)
      MQ = rbind(Data0Q / NormalizeFactor0Q, Data1Q / NormalizeFactor1Q)
      
      
      WP = (nmf(t(MP), dimNMF))@fit@W
      WQ = WP
      
      
      HP = function(b) {
        return((nnls(WP, b)$x))
      }
      HQ = function(b) {
        return((nnls(WQ, b)$x))
      }
      
      
      
    }
    
    
    
    if (!(nchr %in% union(SeqP, SeqQ)))
      
      
    {
      MP = rbind(Data0P / NormalizeFactor0P, Data1P / NormalizeFactor1P)
      MQ = rbind(Data0Q / NormalizeFactor0Q, Data1Q / NormalizeFactor1Q)
      
      
      
      
      
  
      
      #Loading the signatures from the training phase
      WP = readRDS(file.path(dir,paste0('NonNegativeMatrixP', nchr)))
      WQ = readRDS(file.path(dir,paste0('NonNegativeMatrixQ', nchr)))
      
      #Uncomment if you want to regenerate signatures 
      # WP = (nmf(t(MP), dimNMF))@fit@W
      # WQ = (nmf(t(MQ), dimNMF))@fit@W
      
      
      HP = function(b) {
        return((nnls(WP, b)$x))
      }
      HQ = function(b) {
        return((nnls(WQ, b)$x))
      }
      
    }
    
    
    
    
    SpaceP = t(apply(MP, 1, HP))
    SpaceQ = t(apply(MQ, 1, HQ))
    
    
    
    #Extracting the factors from the NMF (features)
    DataReg = cbind(SpaceP, SpaceQ)
    DataReg0 = DataReg[1:nrow(Data0), ]
    DataReg1 = DataReg[(nrow(Data0) + 1):nrow(DataReg), ]
    
    Scores0 = data.frame(Name0 = BV_Normals, score = DataReg0)
    Scores1 = data.frame(Name1 = BV_Cancers, score = DataReg1)
    
    S0 = join(S0, Scores0, by = 'Name0')
    S1 = join(S1, Scores1, by = 'Name1')
    
    
    
    
    
  }
  
  
  
  
  
  #Saving the signatures of the NMF (needed for test)
  #saveRDS(WP, file = paste0('NonNegativeMatrixP', nchr))
  #saveRDS(WQ, file = paste0('NonNegativeMatrixQ', nchr))
  
  S0 = filter(S0, !is.na(apply(S0[, 2:ncol(S0)], 1, sum)))
  S1 = filter(S1, !is.na(apply(S1[, 2:ncol(S1)], 1, sum)))
  write.csv(S0,
            file.path(dir,paste0('NMFScoresNormalsTrain', nchr)),
            row.names = FALSE)
  write.csv(S1,
            file.path(dir,paste0('NMFScoresCancersTrain', nchr)),
            row.names = FALSE)
  
  
  
  
  
  
  
  
  
  
}



saveRDS(BV_Normals, file.path(dir,'bvNormalsTrain'))
saveRDS(BV_Cancers, file.path(dir,'bvCancersTrain'))
