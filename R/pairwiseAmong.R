# ----------------------------------------------------------------------
# pairwiseAmong: Function for displaying the mean comparison
# Created by: Alaine A. Gulles for International Rice Research Institute
# Modified by: Alaine A. Gulles
# ----------------------------------------------------------------------

pairwiseAmong <- function(data, respvar, typeTest, trmt, environment, dfError, MSError, siglevel, rdata, analysisId){ #UseMethod("pairwiseAmong")

# pairwiseAmong.default <- function(data, respvar, typeTest, trmt, dfError, MSError, siglevel) {
     options(width = 5000)
     if (nlevels(data[,trmt]) > 26) {
          command <- paste("STAR::",typeTest, ".test(data['",respvar,"'], data['",trmt,"'], dfError, MSError, alpha = ", siglevel,", group = FALSE, pwOrder = 'trmt')", sep = "")
          result.pw <- eval(parse(text = command))

          aovPredictions <- data.frame(module = rep("aov",nrow(result.pw$summary)), analysisId = rep(analysisId,nrow(result.pw$summary)),
                                       trait = rep(respvar,nrow(result.pw$summary)), environment = rep(environment,nrow(result.pw$summary)),
                                       designation = rownames(result.pw$summary),
                                       predictedValue = result.pw$summary$MeanDiff,
                                       reliability = result.pw$summary$Prob,
                                       entryType = result.pw$summary$Sig)

          rdata$predictions <- rbind(rdata$predictions, aovPredictions)

          aovMetrics <- data.frame(module = rep("aov",2), analysisId = rep(analysisId,2), trait = rep(respvar,2),
                                   environment = rep(environment,2), parameter = c("Critical Value","Test Statistics"),
                                   method = rep(result.pw$method, 2),
                                   value = c(result.pw$tabValue, result.pw$testStat),
                                   stdError = rep(0,2))

          rdata$metrics <- rbind(rdata$metrics, aovMetrics)

     } else {
          if (length(typeTest) == 1) {
               command <- paste("STAR::",typeTest, ".test(data['",respvar,"'], data['",trmt,"'], dfError, MSError, alpha = ", siglevel,", group = TRUE, pwOrder = 'trmt')", sep = "")
               result.pw <- eval(parse(text = command))

               aovPredictions <- data.frame(module = rep("aov",nrow(result.pw$summary)), analysisId = rep(analysisId,nrow(result.pw$summary)),
                                            trait = rep(respvar,nrow(result.pw$summary)), environment = rep(environment,nrow(result.pw$summary)),
                                            designation = as.vector(result.pw$summary[[1]]),
                                            predictedValue = result.pw$summary$means,
                                            stdError = result.pw$summary$std.err,
                                            entryType = result.pw$summary$group)

               rdata$predictions <- rbind(rdata$predictions, aovPredictions)

               aovMetrics <- data.frame(module = rep("aov",2), analysisId = rep(analysisId,2), trait = rep(respvar,2),
                                        environment = rep(environment,2), parameter = c("Critical Value","Test Statistics"),
                                        method = rep(result.pw$method, 2),
                                        value = c(result.pw$tabValue, result.pw$testStat),
                                        stdError = rep(0,2))

               rdata$metrics <- rbind(rdata$metrics, aovMetrics)

          } else {
               singleVal <- na.omit(match(typeTest, c("LSD", "HSD", "scheffe")))
               if (length(singleVal) != 0) {
                    summaryTable <- NULL
                    summaryStat <- NULL
                    for (i in (1:length(singleVal))) {
                         command <- paste("STAR::",c("LSD", "HSD", "scheffe")[singleVal[i]], ".test(data['",respvar,"'], data['",trmt,"'], dfError, MSError, alpha = ", siglevel,", group = TRUE, pwOrder = 'trmt')", sep = "")

                         result.pw <- eval(parse(text = command))

                         aovPredictions <- data.frame(module = rep("aov",nrow(result.pw$summary)), analysisId = rep(analysisId,nrow(result.pw$summary)),
                                                      trait = rep(respvar,nrow(result.pw$summary)), environment = rep(environment,nrow(result.pw$summary)),
                                                      designation = as.vector(result.pw$summary[[1]]),
                                                      predictedValue = result.pw$summary$means,
                                                      stdError = result.pw$summary$std.err,
                                                      entryType = result.pw$summary$group)

                         rdata$predictions <- rbind(rdata$predictions, aovPredictions)

                         aovMetrics <- data.frame(module = rep("aov",2), analysisId = rep(analysisId,2), trait = rep(respvar,2),
                                                  environment = rep(environment,2), parameter = c("Critical Value","Test Statistics"),
                                                  method = rep(result.pw$method, 2),
                                                  value = c(result.pw$tabValue, result.pw$testStat),
                                                  stdError = rep(0,2))

                         rdata$metrics <- rbind(rdata$metrics, aovMetrics)

                    }

                    # if (length(singleVal) > 1) {
                    #      maxWidthEntry <- max(nchar(dfError), nchar(round(MSError,0))) + 7
                    #      cat("\n")
                    #      cat(formatC("Alpha", format = "s", width = 25, flag = "-"), formatC(siglevel, format = "f", digits = 2, width = maxWidthEntry , flag = "#"), "\n", sep = "")
                    #      cat(formatC("Error Degrees of Freedom", format = "s", width = 25, flag = "-"), formatC(dfError, format = "d", width = maxWidthEntry , flag = "#"), "\n", sep = "")
                    #      cat(formatC("Error Mean Square", format = "s", width = 25, flag = "-"), formatC(MSError, format = "f", digits = 4, width = maxWidthEntry , flag = "#"), "\n\n", sep = "")
                    #
                    #      maxWidthEntry <- nchar(max(round(summaryStat,0))) + 7
                    #      maxWidthLabel <- max(nchar(rownames(summaryStat))) + 2
                    #
                    #      cat(formatC(paste(rep("-",maxWidthLabel+(length(singleVal)*maxWidthEntry)),collapse = ""), width = maxWidthLabel+(length(singleVal)*maxWidthEntry)), "\n", sep = "")
                    #      cat(formatC("", format = "s", width = maxWidthLabel, flag = "-"), formatC(colnames(summaryStat), format = "s", width = maxWidthEntry, flag = "#"), "\n", sep = "")
                    #      cat(formatC(paste(rep("-",maxWidthLabel+(length(singleVal)*maxWidthEntry)),collapse = ""), width = maxWidthLabel+(length(singleVal)*maxWidthEntry)), "\n", sep = "")
                    #      for (i in (1:2)) {
                    #           cat(formatC(rownames(summaryStat)[i], format = "s", width = maxWidthLabel, flag = "-"), sep = "")
                    #           for (j in (1:ncol(summaryStat))) {
                    #                cat(formatC(summaryStat[i,j], format = "f", digits = 4, width = maxWidthEntry, flag = "#"),sep = "")
                    #           }
                    #           cat("\n")
                    #      }
                    #      cat(formatC(paste(rep("-",maxWidthLabel+(length(singleVal)*maxWidthEntry)),collapse = ""), width = maxWidthLabel+(length(singleVal)*maxWidthEntry)), "\n", sep = "")
                    #
                    #      cat("\nSummary: \n")
                    #      STAR::printDataFrame(summaryTable, digits = 4)
                    #      cat("Means with the same letter are not significantly different\n\n")
                    #
                    #
                    # }

               }

               rangeVal <- na.omit(match(typeTest, c("duncan", "SNK")))
               if (length(rangeVal) != 0) {
                    for (i in (1:length(rangeVal))) {
                         command <- paste("STAR::",c("duncan", "SNK")[rangeVal[i]], ".test(data['",respvar,"'], data['",trmt,"'], dfError, MSError, alpha = ", siglevel,", group = TRUE, pwOrder = 'trmt')", sep = "")
                         result.pw <- eval(parse(text = command))

                         aovPredictions <- data.frame(module = rep("aov",nrow(result.pw$summary)), analysisId = rep(analysisId,nrow(result.pw$summary)),
                                                      trait = rep(respvar,nrow(result.pw$summary)), environment = rep(environment,nrow(result.pw$summary)),
                                                      designation = as.vector(result.pw$summary[[1]]),
                                                      predictedValue = result.pw$summary$means,
                                                      stdError = result.pw$summary$std.err,
                                                      entryType = result.pw$summary$group)

                         rdata$predictions <- rbind(rdata$predictions, aovPredictions)

                         for (j in 2:ncol(result.pw$testStat)){
                           aovMetrics <- data.frame(module = rep("aov",2), analysisId = rep(analysisId,2), trait = rep(respvar,2),
                                                    environment = rep(environment,2), parameter = c(paste0("Critical Value_Number of Means ",colnames(result.pw$testStat)[j]),paste0("Test Statistics_Number of Means ",colnames(result.pw$testStat)[j])),
                                                    method = rep(result.pw$method, 2),
                                                    value = c(result.pw$testStat[1,j], result.pw$testStat[2,j]),
                                                    stdError = rep(0,2))

                           rdata$metrics <- rbind(rdata$metrics, aovMetrics)

                         }
                    }
               }
          }
     }

  return(rdata)

} ## END FUNCTION
