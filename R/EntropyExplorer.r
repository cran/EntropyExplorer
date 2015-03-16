#input is a matrix
#nrows == gene count
#ncols == sample count
#
#code sample to produce such a matrix and call the functions
#data file should contain both the header and rownames
#>chick=read.table(file="Chicken_Metabolite.csv",header=TRUE, sep=',',row.names=1) #sep='\t'
#>chickm <- as.matrix(chick)		#both rownames(chick) and rownames(chickm) work fine
# chickm<-as.matrix(read.table(file="Chicken_Metabolite.csv",header=TRUE, sep='\t',row.names=1))

EntropyExplorer <- function(expm1,expm2,dmetric, otype,ntop,nperm,shift=c(0,0))
{
	if(missing(expm1))
	{
		stop("expm1 is a required input.")
	}
	if(missing(expm2))
	{
		stop("expm2 is a required input.")
	}
	
	if(missing(dmetric))
	{
		stop("dmetric is a required input.")
	}
	if(missing(otype))
	{
		stop("otype is a required input.")
	}
	if(tolower(dmetric)!="de" & tolower(dmetric)!="dcv" & tolower(dmetric)!="dse")
	{
		stop("supported dmetric: de, dcv, dse")
	}
	if(tolower(otype)!="v" & tolower(otype)!="p" & tolower(otype)!="bv" & tolower(otype)!="bp")
	{
		stop("supported otype: v, p, bv, bp")
	}
	if((tolower(dmetric)=="de" |  tolower(dmetric)=="dcv") & !missing(nperm))
	{
		stop("nperm need not be specified for differential expression analysis and differential CV analysis.")
	}
	if(tolower(dmetric)=="dse" &  tolower(otype)=="v" & !missing(nperm))
	{
		stop("nperm need not be specified for differential nse analysis with otype of v.")
	}
	
	min1 <- min(expm1,na.rm=TRUE)
	min2 <- min(expm2,na.rm=TRUE)
	if(min1 <=0 | min2 <=0)
	{
		if(missing(shift))
		{
			if(min1 <=0)
			{
				cat(paste("The minimum value of matrix 1 is", as.character(min1),"\n"),fill=TRUE)
			}
			if(min2 <=0)
			{
				cat(paste("The minimum value of matrix 2 is", as.character(min2),"\n"),fill=TRUE)
			}
			stop("matrix 1 or 2 has non-positive values. EntropyExplorer expects all data values to be positive. You may rerun this function using the shift parameter to shift the data upward, or adjust the data yourself.\n")
		}
		else#if(missing(shift))
		{
			if(length(shift)!=2)
			{
				stop("The shift parameter should have length 2.")
			}
			if(min1<=0)
			{
				if(tolower(shift[[1]])=="auto")
				{
					expm1<-expm1+(abs(min1)+0.001) #use 0.001 as a default shift
				}
				else if(is.na(suppressWarnings(as.numeric(shift[[1]]))))
				{
					stop("Your shift 1 is not \"auto\" or a valid number.")
				}
				else
				{
					s1<-as.numeric(shift[[1]])
					if(s1<=abs(min1))
					{
						stop(paste("The minimum value of matrix 1 is", as.character(min1), ". Your shift 1 is not large enough.\n"))
					}
					else
					{
						expm1<-expm1+s1
					}
				}
			}#if(min1<=0)
			else
			{
				cat(paste("All values of matrix 1 are positive. Your shift 1 is ignored.\n"),fill=TRUE)
			}
			if(min2<=0)
			{
				if(tolower(shift[[2]])=="auto")
				{
					expm2<-expm2+(abs(min2)+0.001) #use 0.001 as a default shift
				}
				else if(is.na(suppressWarnings(as.numeric(shift[[2]]))))
				{
					stop("Your shift 2 is not \"auto\" or a valid number.")
				}
				else
				{
					s2<-as.numeric(shift[[2]])
					if(s2<=abs(min2))
					{
						stop(paste("The minimum value of matrix 2 is", as.character(min2), ". Your shift 2 is not large enough.\n"))
					}
					else
					{
						expm2<-expm2+s2
					}
				}
			}#if(min2<=0)
			else
			{
				cat(paste("All values of matrix 2 are positive. Your shift 2 is ignored.\n"),fill=TRUE)
			}
		}
	}#if(min1 <=0 | min2 <=0)
	else
	{
		if(!missing(shift))
		{
			cat(paste("All values of the two matrices are positive. Your shifts are ignored.\n"),fill=TRUE)
		}
	}
	
	if(tolower(dmetric)=="de")#differential expression
	{
		if(tolower(otype)=="v")#differential value
		{
			return (diffexp1(expm1,expm2, ntop))
		}
		else if(tolower(otype)=="p")#differential p-value
		{
			return (diffexp2(expm1,expm2, ntop))
		}
		else if(tolower(otype)=="bp")
		{
			return (diffexp(expm1,expm2, ntop,2))
		}
		else if(tolower(otype)=="bv")
		{
			return (diffexp(expm1,expm2, ntop,1))
		}
	}
	else if(tolower(dmetric)=="dcv")#differential CV
	{
		if(tolower(otype)=="v")#differential value
		{
			return (diffcv(expm1,expm2, ntop))
		}
		else if(tolower(otype)=="p")#differential p-value
		{
			return (diffcvpFK(expm1,expm2, ntop))
		}
		else if(tolower(otype)=="bp")
		{
			return (diffcvall(expm1,expm2, ntop,2))
		}
		else if(tolower(otype)=="bv")
		{
			return (diffcvall(expm1,expm2, ntop,1))
		}
	}
	else if(tolower(dmetric)=="dse")#differential entropy
	{
		if(tolower(otype)=="v")#differential value
		{
			return (diffnse(expm1,expm2, ntop))
		}
		else if(tolower(otype)=="p")#differential p-value
		{
			return (diffnsepperm(expm1,expm2, ntop,nperm))
		}
		else if(tolower(otype)=="bp")
		{
			return (diffseall(expm1,expm2, ntop,nperm,2))
		}
		else if(tolower(otype)=="bv")
		{
			return (diffseall(expm1,expm2, ntop,nperm,1))
		}
	}
}#EntropyExplorer

#expm1 and expm2 are of class matrix; diffexp1 is function name
#returns a one-column matrix (absolute differences between means) with row number equal to ntop
diffexp1 <- function(expm1,expm2, ntop)
{
	dimv1 <- dim(expm1)
	dimv2 <- dim(expm2)
	stopifnot(dimv1[1] == dimv2[1])
	if(missing(ntop))
	{
		ntop<-dimv1[1]
	}
	stopifnot(ntop>0)
	stopifnot(ntop<=dimv1[1])
	genecount <- dimv1[1]
	de1 <- matrix(0,nrow=genecount,ncol=1,dimnames=list(rownames(expm1),"differential expression"))
	nacount <- 0
	for (i in 1:genecount)
	{
		x1 <- expm1[i,][!is.na(expm1[i,])]
		if (length(x1)<4)
		{
			de1[i,1] <- -1
			nacount <- nacount+1
			next
		}
		x2 <- expm2[i,][!is.na(expm2[i,])]
		if (length(x2)<4)
		{
			de1[i,1] <- -1
			nacount <- nacount+1
			next
		}
		e1 <- mean(x1)
		e2 <- mean(x2)
		de1[i,1] <- abs(e1-e2)
	}
	if(nacount>0)
	{
		warning(nacount, " gene(s)/probe ID(s) were excluded because of too few expression values in case or control.")
	}
	de1 <- de1[ order(-de1[,1]), ,drop=FALSE]#decreasing order
	if(ntop>dimv1[1]-nacount)
	{
		ntop <- dimv1[1]-nacount
	}
	return (de1[c(0:ntop), ,drop=FALSE]);
}#diffexp1

#expm1 and expm2 are of class matrix; diffexp2 is function name
#returns a one column matrix of p-values (t-test) with row-number equal to ntop
diffexp2 <- function(expm1,expm2,ntop)
{
	dimv1 <- dim(expm1)
	dimv2 <- dim(expm2)
	stopifnot(dimv1[1] == dimv2[1])
	if(missing(ntop))
	{
		ntop<-dimv1[1]
	}
	stopifnot(ntop>0)
	stopifnot(ntop<=dimv1[1])
	genecount <- dimv1[1]
	de1 <- matrix(0,nrow=genecount,ncol=1,dimnames=list(rownames(expm1),"p-value"))
	nacount <- 0
	for (i in 1:genecount)
	{
		x1 <- expm1[i,][!is.na(expm1[i,])]
		if (length(x1)<4)
		{
			de1[i,1] <- 2
			nacount <- nacount+1
			next
		}
		x2 <- expm2[i,][!is.na(expm2[i,])]
		if (length(x2)<4)
		{
			de1[i,1] <- 2
			nacount <- nacount+1
			next
		}
		myt=t.test(x1,x2)
		de1[i,1] <- myt[['p.value']]
	}
	if(nacount>0)
	{
		warning(nacount, " gene(s)/probe ID(s) were excluded because of too few expression values in case or control.")
	}
	de1 <- de1[order(de1[,1]), ,drop=FALSE]#increasing order
	if(ntop>dimv1[1]-nacount)
	{
		ntop <- dimv1[1]-nacount
	}
	return (de1[c(0:ntop), ,drop=FALSE]);
}#diffexp2

diffexp <- function(expm1,expm2, ntop,sorder)
{
	dimv1 <- dim(expm1)
	dimv2 <- dim(expm2)
	stopifnot(dimv1[1] == dimv2[1])
	if(missing(ntop))
	{
		ntop<-dimv1[1]
	}
	stopifnot(ntop>0)
	stopifnot(ntop<=dimv1[1])
	genecount <- dimv1[1]
	detwo <- matrix(0,nrow=genecount,ncol=2,dimnames=list(rownames(expm1),c("differential expression","p-value")))#col1: value; col2: p-value
	nacount <- 0
	for (i in 1:genecount)
	{
		x1 <- expm1[i,][!is.na(expm1[i,])]
		if (length(x1)<4)
		{
			detwo[i,1] <- -1
			detwo[i,2] <- 2
			nacount <- nacount+1
			next
		}
		x2 <- expm2[i,][!is.na(expm2[i,])]
		if (length(x2)<4)
		{
			detwo[i,1] <- -1
			detwo[i,2] <- 2
			nacount <- nacount+1
			next
		}
		e1 <- mean(x1)
		e2 <- mean(x2)
		detwo[i,1] <- abs(e1-e2)
		myt=t.test(x1,x2)
		detwo[i,2] <- myt[['p.value']]
	}
	if(nacount>0)
	{
		warning(nacount, " gene(s)/probe ID(s) were excluded because of too few expression values in case or control.")
	}
	if(sorder==2)
	{
		detwo<-detwo[ order(detwo[,sorder]), , drop=FALSE];
	}
	else
	{
		detwo<-detwo[ order(-detwo[,sorder]), , drop=FALSE];
	}
	if(ntop>dimv1[1]-nacount)
	{
		ntop <- dimv1[1]-nacount
	}
	return (detwo[c(0:ntop), ,drop=FALSE]);
}#diffexp

#expm1 and expm2 are of class matrix; diffcv is function name
diffcv <- function(expm1,expm2,ntop)
{
	dimv1 <- dim(expm1)
	dimv2 <- dim(expm2)
	stopifnot(dimv1[1] == dimv2[1])
	if(missing(ntop))
	{
		ntop<-dimv1[1]
	}
	stopifnot(ntop>0)
	stopifnot(ntop<=dimv1[1])
	genecount <- dimv1[1]
	dcv <- matrix(0,nrow=genecount,ncol=1,dimnames=list(rownames(expm1),"differential CV"))
	nacount <- 0
	for (i in 1:genecount)
	{
		x1 <- expm1[i,][!is.na(expm1[i,])]
		if (length(x1)<4)
		{
			dcv[i,1] <- -1
			nacount <- nacount+1
			next
		}
		x2 <- expm2[i,][!is.na(expm2[i,])]
		if (length(x2)<4)
		{
			dcv[i,1] <- -1
			nacount <- nacount+1
			next
		}
		cv1 <- sd(x1)/mean(x1)
		cv2 <- sd(x2)/mean(x2)
		dcv[i,1] <- abs(cv1-cv2)
	}
	if(nacount>0)
	{
		warning(nacount, " gene(s)/probe ID(s) were excluded because of too few expression values in case or control.")
	}
	dcv <- dcv[order(-dcv[,1]), ,drop=FALSE]
	if(ntop>dimv1[1]-nacount)
	{
		ntop <- dimv1[1]-nacount
	}
	return (dcv[c(0:ntop), ,drop=FALSE]);
}#diffcv

#expm1 and expm2 are of class matrix; diffcvpFK is function name
#returns a one-column matrix of pvalues
#based on Fligner-Killeen (median) test
diffcvpFK <- function(expm1,expm2,ntop)
{
	dimv1 <- dim(expm1)
	dimv2 <- dim(expm2)
	stopifnot(dimv1[1] == dimv2[1])
	if(missing(ntop))
	{
		ntop<-dimv1[1]
	}
	stopifnot(ntop>0)
	stopifnot(ntop<=dimv1[1])
	genecount <- dimv1[1]
	pcv <- matrix(0,nrow=genecount,ncol=1,dimnames=list(rownames(expm1),"p-value"))
	nacount<-0
	for (i in 1:genecount)
	{
		x1 <- expm1[i,][!is.na(expm1[i,])]
		if (length(x1)<4)
		{
			pcv[i,1] <- 2
			nacount<-nacount+1
			next
		}
		x2 <- expm2[i,][!is.na(expm2[i,])]
		if (length(x2)<4)
		{
			pcv[i,1] <- 2
			nacount<-nacount+1
			next
		}
		result <- fligner.test(list(log(x1),log(x2)))
		pcv[i,1] <- result$p.value
	}
	if(nacount>0)
	{
		warning(nacount, " gene(s)/probe ID(s) were excluded because of too few expression values in case or control.")
	}
	pcv <- pcv[order(pcv[,1]), ,drop=FALSE]
	if(ntop>dimv1[1]-nacount)
	{
		ntop <- dimv1[1]-nacount
	}
	return (pcv[c(0:ntop), ,drop=FALSE]);
}#diffcvpFK

diffcvall <- function(expm1,expm2, ntop, sorder)
{
	dimv1 <- dim(expm1)
	dimv2 <- dim(expm2)
	stopifnot(dimv1[1] == dimv2[1])
	if(missing(ntop))
	{
		ntop<-dimv1[1]
	}
	stopifnot(ntop>0)
	stopifnot(ntop<=dimv1[1])
	genecount <- dimv1[1]
	dcvtwo <- matrix(0,nrow=genecount,ncol=2,dimnames=list(rownames(expm1),c("differential CV","p-value")))#col1: value; col2: p-value
	nacount <- 0
	for (i in 1:genecount)
	{
		x1 <- expm1[i,][!is.na(expm1[i,])]
		if (length(x1)<4)
		{
			dcvtwo[i,1] <- -1
			dcvtwo[i,2] <- 2
			nacount <- nacount+1
			next
		}
		x2 <- expm2[i,][!is.na(expm2[i,])]
		if (length(x2)<4)
		{
			dcvtwo[i,1] <- -1
			dcvtwo[i,2] <- 2
			nacount <- nacount+1
			next
		}
		cv1 <- sd(x1)/mean(x1)
		cv2 <- sd(x2)/mean(x2)
		dcvtwo[i,1] <- abs(cv1-cv2)
		result <- fligner.test(list(log(x1),log(x2)))
		dcvtwo[i,2] <- result$p.value
	}
	if(nacount>0)
	{
		warning(nacount, " gene(s)/probe ID(s) were excluded because of too few expression values in case or control.")
	}
	if(sorder==2)
	{
		dcvtwo<-dcvtwo[ order(dcvtwo[,sorder]), , drop=FALSE];#ascending order
	}
	else
	{
		dcvtwo<-dcvtwo[ order(-dcvtwo[,sorder]), , drop=FALSE];#decending order
	}
	if(ntop>dimv1[1]-nacount)
	{
		ntop <- dimv1[1]-nacount
	}
	return (dcvtwo[c(0:ntop), ,drop=FALSE]);
}#diffcvall

#expm1 and expm2 are of class matrix; diffnse is function name
#returns a one-column matrix with row number equal to ntop
diffnse <- function(expm1,expm2,ntop)
{
	dimv1 <- dim(expm1)
	dimv2 <- dim(expm2)
	stopifnot(dimv1[1] == dimv2[1])
	if(missing(ntop))
	{
		ntop<-dimv1[1]
	}
	stopifnot(ntop>0)
	stopifnot(ntop<=dimv1[1])
	genecount <- dimv1[1]
	dnse <- matrix(0,nrow=genecount,ncol=1,dimnames=list(rownames(expm1),"differential entropy"))
	nacount<-0
	for (i in 1:genecount)
	{
		x1 <- expm1[i,][!is.na(expm1[i,])]
		if (length(x1)<4)
		{
			dnse[i,1] <- -1
			nacount<-nacount+1
			next
		}
		x2 <- expm2[i,][!is.na(expm2[i,])]
		if (length(x2)<4)
		{
			dnse[i,1] <- -1
			nacount<-nacount+1
			next
		}
		s <- sum(x1)
		t <- 0
		for (j in 1:length(x1))
		{
			t <- t+x1[j]/s*log2(x1[j]/s)
		}
		nse1 <- -t/log2(length(x1))
		
		s <- sum(x2)
		t <- 0
		for (j in 1:length(x2))
		{
			t <- t+x2[j]/s*log2(x2[j]/s)
		}
		nse2 <- -t/log2(length(x2))
		dnse[i,1] <- abs(nse1-nse2)
	}
	if(nacount>0)
	{
		warning(nacount, " gene(s)/probe ID(s) were excluded because of too few expression values in case or control.")
	}
	dnse <- dnse[order(-dnse[,1]), ,drop=FALSE]
	if(ntop>dimv1[1]-nacount)
	{
		ntop <- dimv1[1]-nacount
	}
	return (dnse[c(0:ntop), ,drop=FALSE]);
}#diffnse

#expm1 and expm2 are of class matrix; diffnsepperm is function name
#returns a one-column matrix of pvalues with row number equal to ntop
#uses the permutation test method
diffnsepperm <- function(expm1,expm2,ntop,nperm)
{
	dimv1 <- dim(expm1)
	dimv2 <- dim(expm2)
	stopifnot(dimv1[1] == dimv2[1])
	if(missing(ntop))
	{
		ntop<-dimv1[1]
	}
	if(missing(nperm))
	{
		nperm<-1000
	}
	stopifnot(ntop>0)
	stopifnot(nperm>=100)
	stopifnot(ntop<=dimv1[1])
	genecount <- dimv1[1]
	dnse <- matrix(0,nrow=genecount,ncol=1) #to store signed differential value
	nacount<-0
	
	for (i in 1:genecount)
	{
		x1 <- expm1[i,][!is.na(expm1[i,])]
		if (length(x1)<4)
		{
			dnse[i,1] <- 0
			nacount<-nacount+1
			next
		}
		x2 <- expm2[i,][!is.na(expm2[i,])]
		if (length(x2)<4)
		{
			dnse[i,1] <- 0
			nacount<-nacount+1
			next
		}
		s <- sum(x1)
		t <- 0
		for (j in 1:length(x1))
		{
			t <- t+x1[j]/s*log2(x1[j]/s)
		}
		nse1 <- -t/log2(length(x1))
		
		s <- sum(x2)
		t <- 0
		for (j in 1:length(x2))
		{
			t <- t+x2[j]/s*log2(x2[j]/s)
		}
		nse2 <- -t/log2(length(x2))
		dnse[i,1] <- nse1-nse2
	}
	if(nacount>0)
	{
		warning(nacount, " gene(s)/probe ID(s) were excluded because of too few expression values in case or control.")
	}

	pe <- matrix(0,nrow=genecount,ncol=1,dimnames=list(rownames(expm1),"p-value")) #to store p-value
	for (i in 1:genecount)
	{
		x1 <- expm1[i,][!is.na(expm1[i,])]
		x1Len <- length(x1)
		if (x1Len<4)
		{
			pe[i,1] <- 2
			next
		}
		x2 <- expm2[i,][!is.na(expm2[i,])]
		x2Len <- length(x2)
		if (x2Len<4)
		{
			pe[i,1] <- 2
			next
		}
		if(dnse[i,1]==0)
		{
			pe[i,1] <- 1
			next
		}
		xc <- c(x1,x2)
		xcLen <- x1Len + x2Len
		ecount<-0
		#Ncount <- 1000 #temporarily set 1000
		Ncount <- nperm
		for (j in 1:Ncount)
		{
			permvec <- sample.int(xcLen)
			xctemp <- xc[permvec] #sample(1:xcLen) produces a permutation of 1,2,...,xcLen
			x1temp <- xctemp[1:x1Len]
			s <- sum(x1temp)
			t <- 0
			for (k in 1:x1Len)
			{
				t <- t+x1temp[k]/s*log2(x1temp[k]/s)
			}
			nse1 <- -t/log2(x1Len)
			x2temp <- xctemp[(1+x1Len):xcLen]
			s <- sum(x2temp)
			t <- 0
			for (k in 1:x2Len)
			{
				t <- t+x2temp[k]/s*log2(x2temp[k]/s)
			}
			nse2 <- -t/log2(x2Len)
			if(dnse[i,1]>0)
			{
				if(nse1-nse2 >= dnse[i,1])
				{
					ecount <- ecount +1
				}
			}
			else #dnse[i,1]<0
			{
				if(nse1-nse2 <= dnse[i,1])
				{
					ecount <- ecount +1
				}
			}
		} # for j
		pe[i,1] <- ecount/Ncount
	}#for i

	pe <- pe[order(pe[,1]), ,drop=FALSE]
	if(ntop>dimv1[1]-nacount)
	{
		ntop <- dimv1[1]-nacount
	}
	return (pe[c(0:ntop), ,drop=FALSE]);
}#diffnsepperm; expand.grid? ReferenceClasses?

diffseall <- function(expm1,expm2, ntop,nperm,sorder)
{
	dimv1 <- dim(expm1)
	dimv2 <- dim(expm2)
	stopifnot(dimv1[1] == dimv2[1])
	if(missing(ntop))
	{
		ntop<-dimv1[1]
	}
	if(missing(nperm))
	{
		nperm<-1000
	}
	stopifnot(ntop>0)
	stopifnot(nperm>=100)
	stopifnot(ntop<=dimv1[1])
	genecount <- dimv1[1]
	dsetwo <- matrix(0,nrow=genecount,ncol=2,dimnames=list(rownames(expm1),c("differential entropy","p-value")))#col1: value; col2: p-value
	dnse <- matrix(0,nrow=genecount,ncol=1) #to store signed differential value
	nacount <- 0
	for (i in 1:genecount)
	{
		x1 <- expm1[i,][!is.na(expm1[i,])]
		if (length(x1)<4)
		{
			dsetwo[i,1] <- -1
			dnse[i,1] <- 0
			nacount <- nacount+1
			next
		}
		x2 <- expm2[i,][!is.na(expm2[i,])]
		if (length(x2)<4)
		{
			dsetwo[i,1] <- -1
			dnse[i,1] <- 0
			nacount <- nacount+1
			next
		}
		s <- sum(x1)
		t <- 0
		for (j in 1:length(x1))
		{
			t <- t+x1[j]/s*log2(x1[j]/s)
		}
		nse1 <- -t/log2(length(x1))
		
		s <- sum(x2)
		t <- 0
		for (j in 1:length(x2))
		{
			t <- t+x2[j]/s*log2(x2[j]/s)
		}
		nse2 <- -t/log2(length(x2))
		dsetwo[i,1] <- abs(nse1-nse2)
		dnse[i,1] <- nse1-nse2
	}
	if(nacount>0)
	{
		warning(nacount, " gene(s)/probe ID(s) were excluded because of too few expression values in case or control.")
	}
	for (i in 1:genecount)
	{
		x1 <- expm1[i,][!is.na(expm1[i,])]
		x1Len <- length(x1)
		if (x1Len<4)
		{
			dsetwo[i,2] <- 2
			next
		}
		x2 <- expm2[i,][!is.na(expm2[i,])]
		x2Len <- length(x2)
		if (x2Len<4)
		{
			dsetwo[i,2] <- 2
			next
		}
		if(dnse[i,1]==0)
		{
			dsetwo[i,2] <- 1
			next
		}
		xc <- c(x1,x2)
		xcLen <- x1Len + x2Len
		ecount<-0
		#Ncount <- 1000 #temporarily set 1000
		Ncount <- nperm
		for (j in 1:Ncount)
		{
			permvec <- sample.int(xcLen)
			xctemp <- xc[permvec] #sample(1:xcLen) produces a permutation of 1,2,...,xcLen
			x1temp <- xctemp[1:x1Len]
			s <- sum(x1temp)
			t <- 0
			for (k in 1:x1Len)
			{
				t <- t+x1temp[k]/s*log2(x1temp[k]/s)
			}
			nse1 <- -t/log2(x1Len)
			x2temp <- xctemp[(1+x1Len):xcLen]
			s <- sum(x2temp)
			t <- 0
			for (k in 1:x2Len)
			{
				t <- t+x2temp[k]/s*log2(x2temp[k]/s)
			}
			nse2 <- -t/log2(x2Len)
			if(dnse[i,1]>0)
			{
				if(nse1-nse2 >= dnse[i,1])
				{
					ecount <- ecount +1
				}
			}
			else #dnse[i,1]<0
			{
				if(nse1-nse2 <= dnse[i,1])
				{
					ecount <- ecount +1
				}
			}
		} # for j
		dsetwo[i,2] <- ecount/Ncount
	}#for i

	if(sorder==2)
	{
		dsetwo<-dsetwo[ order(dsetwo[,sorder]), , drop=FALSE];
	}
	else
	{
		dsetwo<-dsetwo[ order(-dsetwo[,sorder]), , drop=FALSE];
	}
	if(ntop>dimv1[1]-nacount)
	{
		ntop <- dimv1[1]-nacount
	}
	return (dsetwo[c(0:ntop), ,drop=FALSE]);
}#diffseall
