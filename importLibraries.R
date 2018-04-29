library(logging)
library(jsonlite)
library(lubridate)
library(data.table)
library(stringr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(igraph)
library(scales)
library(Matrix)

# Warnings are turned into errors to avoid weird cases
options(warn = 2)


# Setting up the logging level for debugging
basicConfig()

setLevel('FINEST')
setLevel('INFO', getHandler('basic.stdout'))

# Some constants for the analysis
LAST_DATE = as.Date('2017-10-01')
DATA_FILE = 'reports.RData'
KEY_FILE = 'apiKey.txt'

TOP_COUNTRIES = 12 
TOP_REACTIONS = 6
MIN_INTERACTIONS = 10 
MIN_ASSOCIATIONS = 100 

LIMIT_QUERY = 100

# We save our personal API key for openFDA in a file that is not part of the code for security reason (github)
if (!file.exists(KEY_FILE)) {
	logwarn("API file key '%s' is missing: you will not be able to load the whole data set in one pass", KEY_FILE)
	API_KEY = ''
} else {
	loginfo("Loaded API key from '%s' file", KEY_FILE)
	API_KEY = as.character(read.table(KEY_FILE, stringsAsFactors = FALSE))
}

# Shortcut to get some default palette colours
getColours = function() {
	return(brewer.pal(12,"Set3"))
}

# Remove hypens from date for the REST API queries to openFDA
stripHyphens = function(aDate) {
	return(gsub('-', '', aDate))
}

# This is a recursive function (and complicated) that transforms a JSON as given by jsonlite to a list of data frames
# Report IDs are created within the process to be able to join the nested data frame to their original rows
# We limit ourselves to the first layer of nested data frames which will populate the list of data frames 
# We avoid deliberately the openfda variables as we do not use them later for the sake of simplicity
# The returned list is made of a 'results' data frame, 'reactions' data frame and 'drug' data frame 
stripDataFrame = function(x, aName = 'results', aReportId = NA, aLevel = 1, aMaxDepth = 5) {
	theIndent = paste0(rep('', aLevel), collapse = '____')

	logdebug("%s%02d class(x)=%s, aName=%s", theIndent, aLevel, class(x), aName)

	if (aLevel >= aMaxDepth) {
		return(NULL)
	}

	if (is.data.frame(x)) {
		stopifnot(nrow(x) > 0)

		z = list()

		for (i in colnames(x)) {
			logdebug("%s%02d x is a data.frame class(x[[%s]])=%s, aName=%s", theIndent, aLevel, i, class(x[[i]]), aName)

# avoid processing columns that starts with openfda
			if (startsWith(i, 'openfda')) {
				logdebug("Skipping column %s...", i)
				x[[i]] = NULL
			}
			else {
# Recursive call to stripDataFrame
				y = stripDataFrame(x[[i]], i, aReportId, aLevel + 1)

			 	if (is.null(y)) {
					x[[i]] = NULL
				} else {
					if (is.data.frame(y)) {
						stopifnot(nrow(y) > 0)

						if (length(z) == 0) { z = list(y) } else { logdebug("%s%02d old names=%s", theIndent, aLevel, paste0(names(z), collapse = ', ')); z = c(z, list(y)) }
						x[[i]] = NULL
						names(z)[length(z)] = i
						logdebug("%s%02d last names(z)=%s", theIndent, aLevel, paste0(names(z), collapse = ', '))
					}
				}
			}
		}

# Add a reportid column to the dataframe so that data frames within nested list can be joined by to their original report
		if (!is.na(aReportId)) { x$reportid = aReportId } else { x$reportid = seq.int(1, nrow(x)) }

		if (length(z) == 0) { z = list(x) } else { z = c(list(x), z) }

		names(z)[1] = aName
		logdebug("%s%02d first names(z)=%s", theIndent, aLevel, paste0(names(z), collapse = ', '))

		return(z)
	} else if (is.list(x)) {
		logdebug("%s%02d x is a list of length %d, aName=%s", theIndent, aLevel, length(x), aName)

# Recursive call to stripDataFrame over all the element of the nested list
		y = sapply(seq.int(1, length(x)), function(i, a) { stripDataFrame(a[[i]], aName, i, aLevel + 1) }, a = x)
		z = rbindlist(y, fill = TRUE)

		if (nrow(z) == 0) { return (NULL) } else { return(z) }

	} else {
# The processed columns is neither a nested list or data frame, do not do anything special just return the column to preserve it
		logdebug("%s%02d class(x)=%s length(x)=%d x=%s , aName=%s", theIndent, aLevel, class(x), length(x), paste0(x, collapse = ','), aName)

		return(x)
	}

	return(x)

}
