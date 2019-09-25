ssrpairs=function(data1,data2,gid1,gid2){
	timestart<-Sys.time()
	data1=as.matrix(data1)
	data2=as.matrix(data2)
	#binom pvalue
	binom.p=function(i){
		gid11=gid1[-c(1:i)]
		pair1=cbind(gid1[i],gid11)
		coms1=data1[match(pair1[,1],gid1),,drop=F]-data1[match(pair1[,2],gid1),,drop=F]
		freq11=rowSums(coms1>0)
		freq12=rowSums(coms1<0)
		pvalue11=1-pbinom(freq11-1,ncol(coms1),0.5)
		pvalue12=1-pbinom(freq12-1,ncol(coms1),0.5)
		coms2=data2[match(pair1[,1],gid1),,drop=F]-data2[match(pair1[,2],gid1),,drop=F]
		freq21=rowSums(coms2>0)
		freq22=rowSums(coms2<0)
		pvalue21=1-pbinom(freq21-1,ncol(coms2),0.5)
		pvalue22=1-pbinom(freq22-1,ncol(coms2),0.5)
		pvalues=list(p11=pvalue11,p12=pvalue12,p21=pvalue21,p22=pvalue22)
		return(pvalues)
	}
	#所有pvalue合集
	library(doParallel)
	cores=detectCores() #检查可用核数
	cl <- makeCluster(cores-1)
	registerDoParallel(cl)
	result1=foreach(i=1:(length(gid1)-1),.combine='c') %dopar% binom.p(i)#-------------------Note
	stopCluster(cl)
	#分别提取4列p值
	p11=unlist(result1[seq(1,length(result1),4)])
	p12=unlist(result1[seq(2,length(result1),4)])
	p21=unlist(result1[seq(3,length(result1),4)])
	p22=unlist(result1[seq(4,length(result1),4)])
	pairs=t(combn(gid1,2))#耗时很长 提前跑出来 加载RData文件时间缩短很多
	#load("F://KRAS//results//result3//pairs_tcga.RData")
	index1=which(p11>p12)
	p=p11
	p[index1]=p12[index1]
	fdr1=p.adjust(p,method="BH")
	fdr11=fdr1
	fdr11[index1]=1
	fdr12=fdr1
	fdr12[setdiff(1:length(fdr12),index1)]=1
	index2=which(p21>p22)
	p=p21
	p[index2]=p22[index2]
	fdr2=p.adjust(p,method="BH")
	fdr21=fdr2
	fdr21[index2]=1
	fdr22=fdr2
	fdr22[setdiff(1:length(fdr22),index2)]=1
	reversepair1=pairs[fdr11<0.05&fdr22<0.05,]
	reversepair2=pairs[fdr12<0.05&fdr21<0.05,c(2,1)]
	reverse=rbind(reversepair1,reversepair2)
	timeend<-Sys.time()
	runningtime<-timeend-timestart
	get_ssrpairs=list(reverse=reverse,runningtime=runningtime)
	return(get_ssrpairs)
}
#expression
load("data_used.RData")
result1=list()
tongji=c()
for(i in 1:length(geneexp)){
	cat(i,"\n")
	exp1=geneexp[[i]]
	gid1=as.numeric(rownames(exp1))
	exp1=exp1[match(intergene,gid1),]
	cli1=cliinfo[[i]]
	variable=colnames(cli1)
	mut=exp1[,match(cli1[which(cli1$kras=="M"),"sample"],colnames(exp1))]
	if("braf" %in% variable){
		wt=exp1[,match(cli1[which(cli1$kras=="W"&cli1$braf=="W"),"sample"],colnames(exp1))]
	} else{wt=exp1[,match(cli1[which(cli1$kras=="W"),"sample"],colnames(exp1))]}
	result2=ssrpairs(mut,wt,intergene,intergene)
	tongji[i]=nrow(result2$reverse)
	result1[[i]]=result2
}
names(result1)=names(geneexp)
save(result1,tongji,file="reverse005.RData")
