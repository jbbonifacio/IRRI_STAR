# ----------------------------
# FROM agricolae package
# with some modification in presenting the output
# ----------------------------

# AAGulles added the parameter order

`SNK.test` <- function (y, trt, DFerror, MSerror, alpha=0.05, group=TRUE,main = NULL, pwOrder = c("trmt", "means")) {
	pwOrder <- match.arg(pwOrder)	# added by AAGulles
	name.y <- paste(deparse(substitute(y)))
    	# name.t <- paste(deparse(substitute(trt)))
	if (is.data.frame(trt)) { 
		name.t <- names(trt) 
	} else { 
		name.t <- paste(deparse(substitute(trt))) 
		index <- which(strsplit(noquote(name.t), split = "")[[1]] == "\"")
		if (length(index) == 2) {
			name.t <- substr(noquote(name.t), start = index[1], stop = index[2])
		}
	} # modified by AAGulles

    	clase<-c("aov","lm")
    	if("aov"%in%class(y) | "lm"%in%class(y)){
    		A<-y$model
    		DFerror<-df.residual(y)
    		MSerror<-deviance(y)/DFerror
    		y<-A[,1]
    		ipch<-pmatch(trt,names(A))
    		if( is.na(ipch)) return(cat("Name: ",trt,"\n",names(A)[-1],"\n"))
    		name.t <-names(A)[ipch]
    		trt<-A[,ipch]
    		name.y <- names(A)[1]
    	}

	if (group && length(unique(trt)) > 26 ) { group <- FALSE }
    	junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
   	means <- STAR::tapply.stat(junto[,1],junto[,2],stat="mean") # change
    	sds <-   STAR::tapply.stat(junto[,1],junto[,2],stat="sd") #change
    	nn <-   STAR::tapply.stat(junto[,1],junto[,2],stat="length") # change
    	means<-data.frame(means,std.err=sds[,2]/sqrt(nn[,2]),replication=nn[,2])
    	names(means)[1:2]<-c(name.t,name.y)
	# row.names(means)<-means[,1]
    	ntr<-nrow(means)
    	Tprob <- qtukey(1-alpha,2:ntr, DFerror)
    	nr <- unique(nn[,2])
	#"Critical Value of Studentized Range")

	# AAGulles suppress the printing of the following: 
	# cat("\nStudy:", main)
    	# cat("\n\nStudent Newman Keuls Test\nfor",name.y,"\n")
    	# cat("\nMean Square Error: ",MSerror,"\n\n")
    	# cat(paste(name.t,",",sep="")," means\n\n")
    	# print(data.frame(row.names = means[,1], means[,-1]))
	if(length(nr) == 1 ) { sdtdif <- sqrt(MSerror/nr)
	} else {
		nr1 <-  1/mean(1/nn[,2])
		sdtdif <- sqrt(MSerror/nr1)
	}
	SNK <- data.frame(rbind(Tprob, Tprob * sdtdif)) 			# modified by AAGulles
	names(SNK)<-2:ntr
	maxWidthTable <- max(nchar(max(round(SNK, 0)))) + 7	# added by AAGulles
	SNK <- cbind(c("Critical Value", "Test Statistics") , SNK)
	rownames(SNK) <- 1:nrow(SNK)
	colnames(SNK)[1] <- "Number of Means"
	maxLabelTable <- max(nchar(colnames(SNK)[1]), max(nchar(SNK[1]))) + 2

	# AAGulles suppress the printing of the following:
	# cat("\nalpha:",alpha,"; Df Error:",DFerror,"\n")
	# cat("\nCritical Range\n")
	# print(SNK)

	# AAGulles added the following:
	# if (!is.null(main)) { cat("\n", main,"\n", sep = "") }
    	# cat("\nStudent Newman Keuls (SNK) Test\n\n", sep = "") # added by AAGulles
   
    	# maxWidth <- max(nchar(round(MSerror)), nchar(DFerror)) + 7
    	# labelWidth <- 29
	# cat(formatC("Alpha", format = "s", width = labelWidth, flag = "-"), formatC(alpha, format = "f", digits = 2, width =  maxWidth, flag = "#"), "\n",sep = "")  
    	# cat(formatC("Error Degrees of Freedom", format = "s", width = labelWidth, flag = "-"), formatC(DFerror, format = "d", width =  maxWidth, flag = "#"), "\n",sep = "")  	
    	# cat(formatC("Error Mean Square", format = "s", width = labelWidth, flag = "-"), formatC(MSerror, format = "f", digits = 4, width =  maxWidth, flag = "#"), "\n",sep = "")  
		
	# if (length(nr) > 1) {
		# cat(formatC("Harmonic Mean of Cell Sizes", format = "s", width = labelWidth, flag = "-"), formatC(nr1, format = "f", digits = 4, width =  maxWidth, flag = "#"), "\n",sep = "")  
		# cat("\nHarmonic Mean of Cell Sizes ", nr1 )
		# cat("\n\nDifferent value for each comparison")
	# }
	# cat("\n")
	#cat(formatC(colnames(SNK)[1], format = "s", width = maxLabelTable, flag = "-"), formatC(colnames(SNK)[2:ncol(SNK)], format = "d", width = maxWidthTable, flag = "#"), "\n", sep = "")
	#for (i in (1:2)) {
	#	cat(formatC(as.character(SNK[i, 1]), format = "s", width = maxLabelTable, flag = "-"), sep = "")
	#	for (j in (2:ncol(SNK))) {
	#		cat(formatC(SNK[i,j], format = "f", digits = 4, width = maxWidthTable, flag = "#"), sep = "")
	#	}
	#	cat("\n")
	#}
     
	# STAR::printDataFrame(SNK, digits = 4)
	# cat("\n")

	if (group) {
		# cat("\nMeans with the same letter are not significantly different.")
		# cat("\n\nGroups, Treatments and means\n")
		output <- order.group(means[,1], means[,2], means[,4], MSerror, 1 ,means[,3], parameter=0.5,snk=1,DFerror,alpha,sdtdif)
		colnames(output)[1] <- name.t  # added by AAGulles
	  	# AAGulles added the following if-else stmt:
	 	if (pwOrder == "trmt") { 
	 	     if (suppressWarnings(all(!is.na(as.numeric(as.character(factor(trimws(output[,1])))))))) {
	 	          output[,1] <- factor(as.numeric(as.character(output[,1])))
	 	     }
			output <- output[order(output[,1]),]
			rownames(output) <- 1:nrow(output)
        	}
		# cat("Summary of the Result:\n", sep = "") # added by AAGulles
		# STAR::printDataFrame(output[,c(1,2,3,5)])	  # added by AAGulles
	  	# cat("Means with the same letter are not significantly different.\n\n") # added by AAGulles
	}

	if (!group) {
		Omeans<-order(means[,2],decreasing = TRUE)
		Ordindex<-order(Omeans)
		comb <-combn(ntr,2)
		nn<-ncol(comb)
		dif<-rep(0,nn)
		LCL<-dif
		UCL<-dif
		sig<-NULL
		pvalue<-dif
		odif<-dif
		for (k in 1:nn) {
			i<-comb[1,k]
			j<-comb[2,k]
			if (means[i, 2] < means[j, 2]){
				comb[1, k]<-j
				comb[2, k]<-i
			}
			dif[k]<-abs(means[i,2]-means[j,2])
		nx<-abs(i-j)+1
		odif[k] <- abs(Ordindex[i]- Ordindex[j])+1
		pvalue[k]<- round(1-ptukey(dif[k]/sdtdif,odif[k],DFerror),6)
		LCL[k] <- dif[k] - SNK[(odif[k]-1) + 1]
		UCL[k] <- dif[k] + SNK[(odif[k]-1) + 1]
			# sig[k]<-" "
			# AAGulles suppress the following if-else stmt:
			# if (pvalue[k] <= 0.001) sig[k]<-"***"
			# else  if (pvalue[k] <= 0.01) sig[k]<-"**"
			# else  if (pvalue[k] <= 0.05) sig[k]<-"*"
			# else  if (pvalue[k] <= 0.1) sig[k]<-"."

			# added by AAGulles
			remark <- "  "					
			if (alpha <= 0.001)  { remark <-"***"	
        		} else {
				if (alpha <= 0.01) { remark <-"**"	
				} else {
					if (alpha <= 0.05) { remark<-"*"	
					} else {
						if (alpha <= 0.1) remark <-"."	
					}
				}
			}  
			if (pvalue[k] <= alpha) { sig[k] <- remark } else { sig[k] <- "   " }
		}
		tr.i <- means[comb[1, ],1]
		tr.j <- means[comb[2, ],1]
		output<-data.frame("Difference" = dif, pvalue=pvalue,sig,LCL,UCL)
		rownames(output)<-paste(tr.i,tr.j,sep="-")
		output <- output[,c("Difference", "pvalue", "sig")]
		#cat("Comparison between treatments means\n")
		#printDataFrame(cbind("Treatment" = rownames(output),output[,1:3]))	# change by AAGulles
		sigResult <- output[output[,"sig"] == remark, c("Difference", "pvalue")]
		colnames(output) <- c("MeanDiff", "Prob", "Sig")
		if (nrow(sigResult) != 0) {
			colnames(sigResult) <- c("Mean Diff", "Prob")
			# cat("Significant Pairwise Mean Comparison at alpha = ", alpha,"\n")
			# STAR::printDataFrame(cbind("Treatment" = rownames(sigResult),sigResult))	# change by AAGulles
		}
		
		# suppress the following stmt by AAGulles
		# print(output)
		# output<-data.frame(trt= means[,1],means= means[,2],M="",N=means[,4],std.err=means[,3])
	}
	invisible(list(method = "Student Newman Keuls (SNK) Test", testStat = SNK, summary = output))
}

