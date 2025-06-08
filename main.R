load(".../data/BLCA.RData")
source('.../mutation-expression.R')
source('.../mutation-methylation.R')
source('.../calculation similarity.R')

l1=preprocessing(patMutMatrix,FIs1,patOutMatrix)
l2=diffusion(lambda=0.3,PPI,l1$patOuts,l1$patMutMatrix)

l11=preprocessing1(patMutMatrix,FIs1m,patMethMatrix)
l21=diffusion1(lambda1=0.3,PPI,l11$patMeths,l11$patMutMatrix1)

gene_lists= calculation_similarity(l2$turs1, l2$turs2, l21$turs11, l21$turs22, l1$patMutMatrix, l11$patMutMatrix1,mu = 0.8)

