load('.../data/FIs.Rdata')

preprocessing1=function(patMutMatrix1,FIs1m,patMethMatrix)
{
  Mutrow1=intersect(colnames(patMutMatrix1),row.names(FIs1m))
  Methcol=intersect(colnames(patMethMatrix),colnames(FIs1m))
  pats1=intersect(row.names(patMutMatrix1),row.names(patMethMatrix))
  patMutMatrix1=patMutMatrix1[pats1,Mutrow1]
  patMethMatrix=patMethMatrix[pats1,Methcol]
  index1=which(rowSums(patMutMatrix1)==0)
  if(length(index1)>0)
  {
    patMutMatrix1=patMutMatrix1[-index1,]
    patMethMatrix=patMethMatrix[-index1,]
  }
  pats_M1=patMutMatrix1[,intersect(colnames(patMutMatrix1),colnames(patMethMatrix))]
  pats_O1=patMethMatrix[,intersect(colnames(patMutMatrix1),colnames(patMethMatrix))]
  mutandout1=pats_O1+pats_M1
  fs1=c()
  for(i in 1:nrow(mutandout1))
  {
    outs1=which(mutandout1[i,]==2)
    if(length(outs1)>0)
    {
      outs_name1=names(outs1)
      rr1=cbind(row.names(mutandout1)[i],outs_name1)
      fs1=rbind(fs1,rr1)
    }
  }
  patMeths=patMethMatrix
  fs1=as.data.frame(fs1)
  if(length(fs1)>0)
  {for(i in 1:nrow(fs1))
  {
    patMeths[fs1[i,1],fs1[i,2]]=2
  }
  }
  patMeths[which(patMeths==1)]=2 
  patMeths[which(patMeths==2)]=0 
  return(list(patMeths=patMeths,patMutMatrix1=patMutMatrix1))
}

diffusion1=function(lambda1,PPI,patMeths,patMutMatrix1)
{
  PPI11=PPI[match(row.names(FIs1m),row.names(PPI)),match(colnames(FIs1m),colnames(PPI))]
  mm1=rowSums(PPI11)
  oo1=colSums(PPI11)
  degree1=mm1%*%t(oo1)
  degree11=degree1^lambda1
  Wij1=(1-FIs1m)/(1+exp(-degree11))
  FIs1m=Wij1
  diag(FIs1m)=0
  turs11=0.5*patMutMatrix1+0.5*patMeths%*%t(FIs1m)
  turs22=0.5*patMeths+0.5*turs11%*%FIs1m
  return(list(turs11=turs11,turs22=turs22))
}



Mutrow1=intersect(colnames(patMutMatrix),row.names(PPI))
Methcol=intersect(colnames(patMethMatrix),colnames(PPI))
FIsm=PPI[Mutrow1,Methcol]
if(length(which(rowSums(FIsm)==0))>0)
{
  FIsm=FIsm[-which(rowSums(FIsm)==0),]
}
if(length(which(colSums(FIsm)==0))>0)
{
  FIsm=FIsm[,-which(colSums(FIsm)==0)]
}
Mutrow1=intersect(colnames(patMutMatrix),row.names(FIsm))
Methcol=intersect(colnames(patMethMatrix),colnames(FIsm))
patMutMatrix1=patMutMatrix[,Mutrow1]
patMethMatrix=patMethMatrix[,Methcol]
pathways1=read.table('.../data/BLCA.csv',sep = ',',header = TRUE)
llm=list()
for(i in 1:nrow(pathways1))
{
  llm[[i]]=unlist(strsplit(pathways1[i,8],'/'))
}
names(llm)=pathways1[,2]
Nm=nrow(pathways1)
FIs1m=FIsm
for(i in 1:nrow(FIs1m))
{
  muts1=row.names(FIs1m)[i]
  m_p1=length(which(muts1==unlist(llm)))
  mp_n1=names(which(muts1==unlist(llm)))
  o_index1=which(FIs1m[i,]==1)
  for(j in 1:length(o_index1))
  {
    outs1=names(o_index1)[j]
    o_p1=length(which(outs1==unlist(llm)))
    op_n1=names(which(outs1==unlist(llm)))
    mo_p1=length(intersect(mp_n1,op_n1))
    R_value1=phyper(mo_p1,m_p1,Nm-m_p1,o_p1)
    FIs1m[muts1,outs1]=R_value1
  }
}
FIs1m[which(FIs1m==1)]=0
save(FIs1m,file='.../data/BLCA_FIs_MM.rdata')