# -------------------------------------------------------------------------------
# STAR
# -------------------------------------------------------------------------------
# ClassInformation
# Created by: Alaine A. Gulles for International Rice Research Institute
# Modified by: Alaine A. Gulles and Rose Imee Zhella A. Morantte 04.16.2013
# Note: Include multiple response variable
# -------------------------------------------------------------------------------

ClassInformation <- function(data, respvar = NULL){ #UseMethod("ClassInformation")

# ClassInformation.default <- function(data, respvar = NULL) {
	result <- STAR::DataAttribute(data)
	result <- result[result[,"TYPE"] == "factor",c(1,3,4)]
	result[,2] <- as.numeric(result[,2])
	result[,2] <- factor(result[,2])
	colnames(result) <- c("FACTOR", "NO. OF LEVELS", "LEVELS")
	#cat("Class Level Information\n")
	# cat("Summary Information\n")
	# STAR::printDataFrame(result)
	#cat("\n")
	if (!is.null(respvar) && length(respvar) == 1) {
		if (length(na.omit(data[,respvar])) == nrow(data)) {
		  read <- used <- nrow(data)
			# cat("Number of Observations Read and Used:",nrow(data),"\n")
		} else {
		  read <- nrow(data)
		  used <- length(na.omit(data[,respvar]))
			# cat("Number of Observations Read:",nrow(data),"\n")
			# cat("Number of Observations Used:",length(na.omit(data[,respvar])),"\n\n")
		}
	}
	if (!is.null(respvar) && length(respvar) > 1) {
		numObsNM = nrow((na.omit(data[,respvar])))
		numObs = nrow((data[,respvar]))

		if (numObs == numObsNM) {
		  read <- used <- numObs
		  # cat("Number of Observations:", numObs,"\n")
		} else {
		  read <- numObs
		  used <- numObsNM
			# cat("Number of Observations:", numObs,"\n")
			# cat("Number of Observations Used:", numObsNM,"\n\n")
		}
	}

	return(list("result" = result, "read" = read, "used" = used))
}

