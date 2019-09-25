# go redundancy
go_redundancy=function(reverse11,fd1,method){
	if(method==0){
		reverse12=reverse11
	}else if(method==1){
		#----去冗余 method1---去掉的多
		reverse11=reverse1[order(fd1,decreasing=T),]
		fd1=fd1[order(fd1,decreasing=T)]
		regenes=unique(c(reverse11[,1],reverse11[,2]))
		index1=index2=list()
		for(i in 1:length(regenes)){
			g1=regenes[i]
			index11=which(reverse11[,1]==g1|reverse11[,2]==g1)
			index1[[i]]=index11[1]
			index2[[i]]=index11[-1]
		}
		index1=unlist(index1)
		index2=unlist(index2)
		index=setdiff(index1,index2)
		reverse12=reverse11[index,]
	}else if(method==2){
		#----去冗余 method2---去掉的相对少
		reverse11=reverse1[order(fd1,decreasing=T),]
		fd1=fd1[order(fd1,decreasing=T)]
		reverse12=reverse11
		for(i in 1:nrow(reverse11)){
			for(j in 1:2){
				g1=reverse11[i,j]
				index1=which(reverse12[,1]==g1|reverse12[,2]==g1)
				if(length(index1)<2){
					next}
				reverse12=reverse12[-index1[-1],]
			}
		}
	}
	return(reverse12)
}