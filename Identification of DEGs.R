#--------identify DE genes for RNAseq or microarray data
########data_type: RNAseq_counts, RNAseq_fpkm, microarray
IDEGs=function(case.exp,control.exp,geneid,data_type){
	if(data_type=="RNAseq_counts"){
		library(edgeR)
		colnames(control.exp)=paste("control",1:ncol(control.exp),sep="")
		colnames(case.exp)=paste("case",1:ncol(case.exp),sep="")
		#---kras mutation vs. mutation-like
		exp1=cbind(control.exp,case.exp)
		rownames(exp1)=geneid
		group=factor(c(rep("control",ncol(control.exp)),rep("case",ncol(case.exp))))
		#创建edgeR数据格式
		data1=DGEList(counts=exp1,genes=geneid,group=group)
		#过滤
		index1=rowSums(cpm(data1)>1)>=(ncol(exp1)/2)
		data1=data1[index1,]
		#标准化，默认为TMN
		data1=calcNormFactors(data1)
		data1=estimateCommonDisp(data1)
		data1=estimateTagwiseDisp(data1)
		et1=exactTest(data1)
		fdr=p.adjust(et1$table[,3],method="BH")
		et1$table[,1]=-et1$table[,1]
		geneid=as.numeric(rownames(et1$table))
		DEG1=cbind(geneid,et1$table,fdr)
		DEGlist=list(DEGs=DEG1,data_type=data_type,algorithm="edgeR")
	}
	if(data_type=="RNAseq_fpkm"|data_type=="microarray"){
		library(limma)
		labe1=c(rep(0,ncol(control.exp)),rep(1,ncol(case.exp)))
		exp1=cbind(control.exp,case.exp)
		rownames(exp1)=geneid
		design1=model.matrix(~0+factor(labe1))
		colnames(design1)=c("control","case")##构建专用的标签
		fit1=lmFit(exp1,design1)##对每个基因线性拟合
		contrast.matrix=makeContrasts(case-control,levels=design1)
		fit11=contrasts.fit(fit1,contrast.matrix)#估计两组间的系数和标准差
		fit12=eBayes(fit11)
		DEG1=topTable(fit12,adjust="BH",num=nrow(exp1))
		DEG1=cbind(as.numeric(rownames(DEG1)),DEG1);colnames(DEG1)[c(1,5,6)]=c("geneid","PValue","fdr")
		DEGlist=list(DEGs=DEG1,data_type=data_type,algorithm="limma")
	}
	return(DEGlist)
}

