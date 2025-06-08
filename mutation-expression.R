load('.../data/FIs.Rdata')

preprocessing=function(patMutMatrix,FIs1,patOutMatrix)
{
  Mutrow=intersect(colnames(patMutMatrix),row.names(FIs1))
  Outcol=intersect(colnames(patOutMatrix),colnames(FIs1))
  pats=intersect(row.names(patMutMatrix),row.names(patOutMatrix))
  patMutMatrix=patMutMatrix[pats,Mutrow]
  patOutMatrix=patOutMatrix[pats,Outcol]
  index=which(rowSums(patMutMatrix)==0)
  if(length(index)>0)
  {
    patMutMatrix=patMutMatrix[-index,]
    patOutMatrix=patOutMatrix[-index,]
  }
  pats_M=patMutMatrix[,intersect(colnames(patMutMatrix),colnames(patOutMatrix))]
  pats_O=patOutMatrix[,intersect(colnames(patMutMatrix),colnames(patOutMatrix))]
  mutandout=pats_O+pats_M
  fs=c()
  for(i in 1:nrow(mutandout))
  {
    outs=which(mutandout[i,]==2)
    if(length(outs)>0)
    {
      outs_name=names(outs)
      rr=cbind(row.names(mutandout)[i],outs_name)
      fs=rbind(fs,rr)
    }
  }
  patOuts=patOutMatrix
  fs=as.data.frame(fs)
  if(length(fs)>0)
  {for(i in 1:nrow(fs))
  {
    patOuts[fs[i,1],fs[i,2]]=2
  }
  }
  patOuts[which(patOuts==1)]=2 
  patOuts[which(patOuts==2)]=0 
  return(list(patOuts=patOuts,patMutMatrix=patMutMatrix))
}

diffusion=function(lambda,PPI,patOuts,patMutMatrix)
{
  PPI1=PPI[match(row.names(FIs1),row.names(PPI)),match(colnames(FIs1),colnames(PPI))]
  mm=rowSums(PPI1)
  oo=colSums(PPI1)
  degree=mm%*%t(oo)
  degree1=degree^lambda
  Wij=(1-FIs1)/(1+exp(-degree1))
  FIs1=Wij
  diag(FIs1)=0
  turs1=0.5*patMutMatrix+0.5*patOuts%*%t(FIs1)
  turs2=0.5*patOuts+0.5*turs1%*%FIs1
  
  
  return(list(turs1=turs1,turs2=turs2))
}




Mutrow=intersect(colnames(patMutMatrix),row.names(PPI))
Outcol=intersect(colnames(patOutMatrix),colnames(PPI))
FIs=PPI[Mutrow,Outcol]
if(length(which(rowSums(FIs)==0))>0)
{
  FIs=FIs[-which(rowSums(FIs)==0),]
}
if(length(which(colSums(FIs)==0))>0)
{
  FIs=FIs[,-which(colSums(FIs)==0)]
}
Mutrow=intersect(colnames(patMutMatrix),row.names(FIs))
Outcol=intersect(colnames(patOutMatrix),colnames(FIs))
patMutMatrix=patMutMatrix[,Mutrow]
patOutMatrix=patOutMatrix[,Outcol]
pathways=read.table('.../data/BLCA.csv',sep = ',',header = TRUE)
ll=list()
for(i in 1:nrow(pathways))
{
  ll[[i]]=unlist(strsplit(pathways[i,8],'/'))
}
names(ll)=pathways[,2]
N=nrow(pathways)
FIs1=FIs
for(i in 1:nrow(FIs1))
{
  muts=row.names(FIs1)[i]
  m_p=length(which(muts==unlist(ll)))
  mp_n=names(which(muts==unlist(ll)))
  o_index=which(FIs1[i,]==1)
  for(j in 1:length(o_index))
  {
    outs=names(o_index)[j]
    o_p=length(which(outs==unlist(ll)))
    op_n=names(which(outs==unlist(ll)))
    mo_p=length(intersect(mp_n,op_n))
    R_value=phyper(mo_p,m_p,N-m_p,o_p)
    FIs1[muts,outs]=R_value
  }
}
FIs1[which(FIs1==1)]=0
save(FIs1,file='.../data/BLCA_FIs_ME.rdata')