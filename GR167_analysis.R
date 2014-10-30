#GR167
#まずはレシオデータ
data <- read.table("GR167_mas5_del.txt",row.names=1,header=T,sep="\t")
nsamp <- ncol(data)
ngene <- nrow(data)

data_samp <- read.table("SCADS_ChDB_140725_del.txt",header=T,sep="\t",quote="\"")
Compound <- as.character(data_samp[,3])
ref_col <- grep("Control",Compound)
#grepは部分一致した要素番号を返す
nref <- length(ref_col)
output <- matrix(0,ngene,nsamp-nref)
col_name <- c()

for(i in 1:nref){
	
	s <- ref_col[i]+1
	
	if(i == nref){		
		e <- nsamp		
	}else{		
		e <- ref_col[i+1]-1
	}
				
	ref <- ref_col[i]
	
	for(j in s:e){

		rat <- log(pmax(data[,j],50)/pmax(data[,ref],50))
		col_name <- c(col_name,substr(colnames(data[j]), 1, nchar(colnames(data[j]))-21))
		output[,j-i] <- rat
	
	}
	
}
colnames(output) <- col_name

#グループが書いてあるサンプルにレシオデータをしぼる
group_data <- read.table("140918_GR167_group.txt",header=T,sep="\t",quote="\"")
ratio_group <- output[,as.character(group_data[,1])]
write.table(ratio_group,"140918_GR167_ratio_for_pc_log_50_group.txt",append=F,row.names=F,col.names=F,sep="\t",quote=F)

#主成分分析
data <- read.table("140918_GR167_ratio_for_pc_log_50_group.txt",sep="\t")
res <- prcomp(t(data))
#sname <- substr(as.character(group_data[,1]), 9, nchar(as.character(group_data[,1])))
#プロットするときにサンプル名でプロットしたい場合はありにする

group_type <- group_data[,2]
type <- c("A","B","C","D","E")
ntype <- length(type)
type_num <- as.list(rep(NA,ntype))

for(i in 1:ntype){

type_num[[i]] <- grep(as.character(type[i]),group_type)	
}


#出力
#---------------140919_記号でプロットする場合は以下を使用--------------------------------------------------#


pdf("140918_GR167_plot_group_1_2.pdf")
for(i in 1:3){

plot(res$x[type_num[[i]],1],res$x[type_num[[i]],2],pch=i,cex=1.0,col=i,ann=T,xlim=c(-100,50),ylim=c(-50,100),xlab="res$x[,1]",ylab="res$y[,2]")
par(new=T)
	
}
plot(res$x[type_num[[5]],1],res$x[type_num[[5]],2],pch=18,cex=1.5,col=5,ann=T,xlim=c(-100,50),ylim=c(-50,100),xlab="res$x[,1]",ylab="res$y[,2]")
par(new=T)
plot(res$x[type_num[[4]],1],res$x[type_num[[4]],2],pch=4,cex=1.0,col="darkorange",ann=T,xlim=c(-100,50),ylim=c(-50,100),xlab="res$x[,1]",ylab="res$y[,2]")
par(new=T)
dev.off()

pdf("140918_GR167_plot_group_1_3.pdf")
for(i in 1:3){

plot(res$x[type_num[[i]],1],res$x[type_num[[i]],3],pch=i,cex=1.0,col=i,ann=T,xlim=c(-100,50),ylim=c(-80,50),xlab="res$x[,1]",ylab="res$y[,3]")
par(new=T)
	
}
plot(res$x[type_num[[5]],1],res$x[type_num[[5]],3],pch=18,cex=1.5,col=5,ann=T,xlim=c(-100,50),ylim=c(-80,50),xlab="res$x[,1]",ylab="res$y[,3]")
par(new=T)
plot(res$x[type_num[[4]],1],res$x[type_num[[4]],3],pch=4,cex=1.0,col="darkorange",ann=T,xlim=c(-100,50),ylim=c(-80,50),xlab="res$x[,1]",ylab="res$y[,3]")
par(new=T)
dev.off()

pdf("140918_GR167_plot_group_2_3.pdf")
for(i in 1:3){

plot(res$x[type_num[[i]],2],res$x[type_num[[i]],3],pch=i,cex=1.0,col=i,ann=T,xlim=c(-50,100),ylim=c(-80,50),xlab="res$x[,2]",ylab="res$y[,3]")
par(new=T)
	
}
plot(res$x[type_num[[5]],1],res$x[type_num[[5]],3],pch=18,cex=1.5,col=5,ann=T,xlim=c(-50,100),ylim=c(-80,50),xlab="res$x[,2]",ylab="res$y[,3]")
par(new=T)
plot(res$x[type_num[[4]],1],res$x[type_num[[4]],3],pch=4,cex=1.0,col="darkorange",ann=T,xlim=c(-50,100),ylim=c(-80,50),xlab="res$x[,2]",ylab="res$y[,3]")
par(new=T)
dev.off()


#----------------------------------------------------------------------------------------------------#





#---------------140919_サンプルネームでプロットする場合は以下を使用--------------------------------------------------#

pdf("140918_GR167_plot_group_1_2.pdf")

plot(res$x[,1],res$x[,2],type="n")

col <- c(1:3,5)
for(i in col){
	
	text(res$x[type_num[[i]],1],res$x[type_num[[i]],2],sname[type_num[[i]]],cex=0.8,col=i)
}

text(res$x[type_num[[4]],1],res$x[type_num[[4]],2],sname[type_num[[4]]],cex=0.8,col="darkorange")

#col4が水色でみにくいためオレンジに変更
dev.off()

pdf("140918_GR167_plot_group_1_3.pdf")
plot(res$x[,1],res$x[,3],type="n")

col <- c(1:3,5)
for(i in col){
	
	text(res$x[type_num[[i]],1],res$x[type_num[[i]],3],sname[type_num[[i]]],cex=0.8,col=i)

}
text(res$x[type_num[[4]],1],res$x[type_num[[4]],3],sname[type_num[[4]]],cex=0.8,col="darkorange")
#col4が水色でみにくいためオレンジに変更
dev.off()

pdf("140918_GR167_plot_group_2_3.pdf")
plot(res$x[,2],res$x[,3],type="n")

col <- c(1:3,5,6)
for(i in col){
	
	text(res$x[type_num[[i]],2],res$x[type_num[[i]],3],sname[type_num[[i]]],cex=0.8,col=i)

}
text(res$x[type_num[[4]],2],res$x[type_num[[4]],3],sname[type_num[[4]]],cex=0.8,col="darkorange")
#col4が水色でみにくいためオレンジに変更
dev.off()

#-----------------------------------------------------------------------------------------------------#









#0922主成分分析_3Dプロット　出力部分もちょっと改善

＃ーーーーーーーーーーーーーーーーーーーーデータ読み込みーーーーーーーーーーーーーーーーーーーーーーーーー
group_data <- read.table("140918_GR167_group.txt",header=T,sep="\t",quote="\"")
data <- read.table("140918_GR167_ratio_for_pc_log_50_group.txt",sep="\t")
＃ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー

res <- prcomp(t(data))
group_type <- group_data[,2]

＃ーーーーーーーーーーーーーーーーーーーーグループの種類を入力ーーーーーーーーーーーーーーーーーーーーーーーーーー
type <- c("A","B","C","D","E")
＃ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー

ntype <- length(type)
type_num <- as.list(rep(NA,ntype))
#5成分それぞれがlistになっている。listのlistみたいな感じで各成分の要素数は未設定でOK

for(i in 1:ntype){

	type_num[[i]] <- grep(as.character(type[i]),group_type)
	#grepは要素番号の集まりを返すので、type_num[[i]][[j]]のような指定や、type_num[[i]]<-c(type_num[[i]],grep.....)の必要ない

}

#出力

#install.packages("scatterplot3d")
#library(scatterplot3d)

＃ーーーーーーーーーーーーーーー3Dプロットに用いる3軸の主成分番号ーーーーーーーーーーーーーーーーーーーーーーーーー
comb <- list(c(1,2,3),c(1,3,2),c(2,1,3))
nplot <- 3
＃ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー

ncomb <- length(comb)

for(i in 1:ncomb){
	
	pdf(paste("140922_GR167_plot_group_",comb[[i]][[1]],"_",comb[[i]][[2]],"_",comb[[i]][[3]],".pdf",sep=""))

	max <- as.list(rep(NA,nplot))
	min <- as.list(rep(NA,nplot))
	interval <- as.list(rep(NA,nplot))
	
	for(k in 1:nplot){
		
			max[[k]] <- max(res$x[,comb[[i]][[k]]])
			min[[k]] <- min(res$x[,comb[[i]][[k]]])	
			interval[[k]] <- c(ceiling(max[[k]]),floor(min[[k]]))
			#ceilingはx未満でない最小の整数をかえす（くりあげ)、floorはx以上でない最大の整数を返す（くりさげ)	
	}

	
#ーーーーーーーー色とマークの特別指定ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー

	col <- list(rep(NA,ntype))
	mark <- c()


	for(j in 1:ntype){
	
			if(j == 4){
			
				col[[j]] <- "darkorange"
				
			}else{
				
				col[[j]] <- as.numeric(j)
				
			}
	}
				
	for(j in 1:ntype){
	
			if(j == 5){
			
				mark <- c(mark,18)
								
			}else{
				
				mark <- c(mark,j)
								
			}
	}
		
#ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー
			
	
	for(j in 1:ntype){
			
		scatterplot3d(x=res$x[type_num[[j]],comb[[i]][[1]]],y=res$x[type_num[[j]],comb[[i]][[2]]],z=res$x[type_num[[j]],comb[[i]][[3]]],pch=mark[[j]],color=col[[j]],xlim=interval[[1]],ylim=interval[[2]],zlim=interval[[3]],xlab=paste("THE_",comb[[i]][[1]],"_PC",sep=""),ylab=paste("THE_",comb[[i]][[2]],"_PC",sep=""),zlab=paste("THE_",comb[[i]][[3]],"_PC",sep=""))
		par(new=T)
	
	}
	
	dev.off()
}


＃ーーーーーーーーーーーーーーーーーーーー補足事項ーーーーーーーーーーーーーーーーーーーーーーーーーー
#補足事項
#xlim=eval(parse(text=paste("interval",com[1],sep="")))をモデルプログラムでは用いていた。
#意味は文字列をRの命令文として実行
#つまり、intervalcom[1]と書いてもどこまでが文字列かわからないのでそのままでは実行できない
#そういうときに便利なのがeval(parse)
#サンプルプログラムは以下

a <- c(10,20,30,40,50)

for(i in 1:5){
	eval(parse(text=paste("testestes",a[i],"<-i",sep="")))

}

> testestes10
[1] 1
> testestes20
[1] 2
> testestes30
[1] 3

＃ーーーーーーーーーーーーーーーーーーーー補足事項ーーーーーーーーーーーーーーーーーーーーーーーーーー

