# ----------------------------------------------------------------------
# pairwiseWithin: Function for displaying the mean comparison
# Created by: Alaine A. Gulles for International Rice Research Institute
# Modified by: Alaine A. Gulles
# ----------------------------------------------------------------------

pairwiseWithin <- function(data, respvar, environment, typeTest, nobs1, nobs2, dfError, MSError, f1, f2, f3 = NULL, siglevel, rdata, analysisId ){ #UseMethod("pairwiseWithin")

# pairwiseWithin.default <- function(data, respvar, typeTest, nobs1, nobs2, dfError, MSError, f1, f2, f3 = NULL, siglevel) {
      comparison <- NULL
	if (nlevels(data[,f1]) > 26) {
		for (j in (1:nobs2)) {
			command <- paste("STAR::", typeTest, ".test(data['",respvar,"'][(data['",f2,"']) == levels(data[,'",f2,"'])[[",j,"]]",f3,"], data['",f1,"'][(data['",f2,"']) == levels(data[,'",f2,"'])[[",j,"]]",f3,"], dfError, MSError, alpha = ",siglevel,", group = FALSE, pwOrder = 'trmt')", sep = "")
			cat("\n-----", f2, " = ", levels(data[,f2])[[j]], "-----\n")
			result.pw <- eval(parse(text = command))

			aovPredictions <- data.frame(module = rep("aov",nrow(result.pw$summary)), analysisId = rep(analysisId,nrow(result.pw$summary)),
			                             trait = rep(respvar,nrow(result.pw$summary)), environment = rep(environment,nrow(result.pw$summary)),
			                             designation = rownames(result.pw$summary),
			                             predictedValue = result.pw$summary$MeanDiff,
			                             reliability = result.pw$summary$Prob,
			                             entryType = result.pw$summary$Sig)

			analysisIrdatad$predictions <- rbind(rdata$predictions, aovPredictions)

			aovMetrics <- data.frame(module = rep("aov",2), analysisId = rep(analysisId,2), trait = rep(respvar,2),
			                         environment = rep(environment,2), parameter = c("Critical Value","Test Statistics"),
			                         method = rep(result.pw$method, 2),
			                         value = c(result.pw$tabValue, result.pw$testStat),
			                         stdError = rep(0,2))

			rdata$metrics <- rbind(rdata$metrics, aovMetrics)


		}
	} else {
	      for (j in (1:nobs2)) {
      		# k = 1 + 2*(j - 1)
			        command <- paste("STAR::",typeTest, ".test(data['",respvar,"'][(data['",f2,"']) == levels(data[,'",f2,"'])[[",j,"]]",f3,"], data['",f1,"'][(data['",f2,"']) == levels(data[,'",f2,"'])[[",j,"]]",f3,"], dfError, MSError, alpha = ",siglevel,", group = TRUE, pwOrder = 'trmt')", sep = "")
			        result.pw <- eval(parse(text = command))

	            addComparison <- result.pw$summary[,c(1,3,2,4,5)]
	            addComparison[,f2] <- levels(data[,f2])[[j]]
	            colnames(addComparison)[1] <- f1

	            comparison <- rbind(comparison, addComparison)

	      } ### end stmt -- for (j in 1:nobs2)

		# options(width = 5000)
		# cat("\n", result.pw$method,"\n\n", sep = "")
		# if (typeTest == "LSD" || typeTest == "HSD" || typeTest == "scheffe") {
		# 	maxWidthEntry <- max(nchar(dfError), nchar(round(MSError, 0)), nchar(round(result.pw$tabValue, 0)), nchar(round(result.pw$testStat, 0))) + 7
		# } else {
		# 	maxWidthEntry <- max(nchar(dfError), nchar(round(MSError, 0))) + 7
		# }

		# cat(formatC("Alpha", format = "s", width = 25, flag = "-"), formatC(siglevel, format = "f", digits = 2, width = maxWidthEntry, flag = "#"), "\n", sep = "")
		# cat(formatC("Error Degrees of Freedom", format = "s", width = 25, flag = "-"), formatC(dfError, format = "d", width = maxWidthEntry, flag = "#"), "\n", sep = "")
		# cat(formatC("Error Mean Square", format = "s", width = 25, flag = "-"), formatC(MSError, format = "f", digits = 4, width = maxWidthEntry, flag = "#"), "\n", sep = "")
		if (typeTest == "LSD" || typeTest == "HSD" || typeTest == "scheffe") {
		  # cat(formatC("Critical Value", format = "s", width = 25, flag = "-"), formatC(result.pw$tabValue, format = "f", digits = 4, width =  maxWidthEntry, flag = "#"), "\n", sep = "")
		  # cat(formatC("Test Statistic", format = "s", width = 25, flag = "-"), formatC(result.pw$testStat, format = "f", digits = 4, width =  maxWidthEntry, flag = "#"), "\n\n", sep = "")

		  aovMetrics <- data.frame(module = rep("aov",2), analysisId = rep(analysisId,2), trait = rep(respvar,2),
		                           environment = rep(environment,2), parameter = c("Critical Value","Test Statistics"),
		                           method = rep(result.pw$method, 2),
		                           value = c(result.pw$tabValue, result.pw$testStat),
		                           stdError = rep(0,2))

		  rdata$metrics <- rbind(rdata$metrics, aovMetrics)

		} else {
			# cat("\n")
			#maxWidthLabel <- max(nchar(colnames(result.pw$testStat)[1]), max(nchar(as.character(result.pw$testStat[,1])))) + 2
			#maxWidthEntry <- nchar(max(round(result.pw$testStat[,2:ncol(result.pw$testStat)],0))) + 7
			#cat(formatC(colnames(result.pw$testStat)[1], format = "s", width = maxWidthLabel, flag = "-"), formatC(colnames(result.pw$testStat)[2:ncol(result.pw$testStat)], format = "d", width = maxWidthLabel, flag = "#"), "\n", sep = "")
			#for (i in (1:2)) {
			#	cat(formatC(result.pw$testStat[i,1], format = "s", width = maxWidthLabel, flag = "-"), sep = "")
			#	for (j in (2:ncol(result.pw$testStat))) {
			#		cat(formatC(result.pw$testStat[i,j], format = "d", width = maxWidthLabel, flag = "#"), sep = "")
			#	}
			#	cat("\n")
			#}
			#cat("\n")
		  # STAR::printDataFrame(result.pw$testStat, digits = 4)
			# cat("\n")

			aovMetrics <- data.frame(module = rep("aov",1), analysisId = rep(analysisId,1), trait = rep(respvar,1),
			                         environment = rep(environment,1), parameter = c("Test Statistics"),
			                         method = rep(result.pw$method, 1),
			                         value = c(result.pw$testStat),
			                         stdError = rep(0,1))

			rdata$metrics <- rbind(rdata$metrics, aovMetrics)
		}
		colnames(result.pw$summary)[1] <- f1
		# cat("Summary:", "\n")
		# STAR::printDataFrame(comparison, digits = 4)
		# cat("Means with the same letter are not significantly different\n\n")

		aovPredictions <- data.frame(module = rep("aov",nrow(comparison)), analysisId = rep(analysisId,nrow(comparison)),
		                             trait = rep(respvar,nrow(comparison)), environment = rep(environment,nrow(comparison)),
		                             designation = paste0(comparison[[1]],comparison[[6]]),
		                             predictedValue = comparison$means,
		                             stdError = comparison$std.err,
		                             entryType = comparison$group)

		rdata$predictions <- rbind(rdata$predictions, aovPredictions)

	}

  return(rdata)

} ## END FUNCTION
