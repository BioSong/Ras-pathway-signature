####-------identification of mutation-like samples 
Imlike=function(cancer.exp,cancer.cli,geneid){
	load("F://kras//results//result5//mlike_signatures.RData")
	mlikes=wlikes=list()
	for(i in 1:4){
		spair1=spairs[[i]]
		spair2=cbind(geneid[match(spair1[,1],geneid)],geneid[match(spair1[,2],geneid)])
		spair2=na.omit(spair2);th_half=nrow(spair2)/2
		com1=cancer.exp[match(spair2[,1],geneid),]-cancer.exp[match(spair2[,2],geneid),]
		score=colSums(com1>0);
		mlike=cancer.cli[score>=th_half,];wlike=cancer.cli[score<th_half,]
		mlikes[[i]]=mlike;wlikes[[i]]=wlike
	}
	mc.cli=unique(rbind(mlikes[[1]],mlikes[[2]],mlikes[[3]],mlikes[[4]]));
	wc.cli=unique(rbind(wlikes[[1]],wlikes[[2]],wlikes[[3]],wlikes[[4]]))
	wc.cli=wc.cli[match(setdiff(wc.cli[,1],mc.cli[,1]),wc.cli[,1]),]
	mc.exp=cancer.exp[,match(mc.cli$response1,colnames(cancer.exp))]
	wc.exp=cancer.exp[,match(wc.cli$response1,colnames(cancer.exp))]
	Imlike_result=list(mlike.cli=mc.cli,wlike.cli=wc.cli,mlike.exp=mc.exp,wlike.exp=wc.exp)
	return(Imlike_result)
}

