# This function triggers a REST API call to the openFDA website with a time window passed as arguments
fetchReports = function(aFirstDate, aLastDate, aLimit = LIMIT_QUERY) {
	loginfo("Fetching reports from FDA site between %s and %s...", aFirstDate, aLastDate)

	theQuery = paste0('https://api.fda.gov/drug/event.json?search=receivedate:[', stripHyphens(aFirstDate), '+TO+', stripHyphens(aLastDate), ']', ifelse(aLimit > 0, paste0('&limit=', aLimit), ''), ifelse(nzchar(API_KEY), paste0('&api_key=', API_KEY), ''))

	logdebug("theQuery=%s", theQuery)

# We first receive the JSON from the website
	theReports.json = fromJSON(theQuery, flatten = TRUE)
# We transform it into a list of 3 data tables by calling stripDataFrame
	theResults.l = stripDataFrame(theReports.json$results)

	return(theResults.l)
}

# We create a sequence of consecute periods to loop over
thePeriods.D = seq.Date(LAST_DATE - years(1), LAST_DATE, by = 1)

logdebug("thePeriods=%s", paste0(thePeriods.D, collapse = ', '))

# We test if the variable holding the complete data set does not exist already
if (!exists('theReports.l')) {
# We test if there is a file containing the results available to save downloading time
	if (file.exists(DATA_FILE)) {
		loginfo("Loading data file '%s'", DATA_FILE)
		load(DATA_FILE)
	}
	else {
# If there is no file with an archive of the results then we initialise the containing variable with NULL for a complete download
		loginfo("No previous data file '%s' found: downloading...", DATA_FILE)
		theReports.l = list()
	}
}

# We prepare one data table for each of the three types of data tables that are returned by the REST API call
theResults.dt = data.table(NULL)
theReactions.dt = data.table(NULL)
theDrugs.dt = data.table(NULL)

# We create a flag to test if the final downloaded data needs to be saved at the end because it has been modified
doSave = FALSE

# We iterate over the periods of days 
for (i in seq.int(2, length(thePeriods.D))) {

# We test here if the i - 1 period has already been downloaded in the final list in case we reached a limit with the REST API calls previously
	if (length(theReports.l) < i - 1) {

# We download the list of data frames for the given period
		theReports.l[[i - 1]] = fetchReports(thePeriods.D[i - 1], thePeriods.D[i] - 1)

# We adjust the reportid additional column based of the total of reports loaded so far
		theReports.l[[i - 1]]$results$reportid = theReports.l[[i - 1]]$results$reportid + nrow(theResults.dt)
		theReports.l[[i - 1]]$patient.reaction$reportid = theReports.l[[i - 1]]$patient.reaction$reportid + nrow(theResults.dt)
		theReports.l[[i - 1]]$patient.drug$reportid = theReports.l[[i - 1]]$patient.drug$reportid + nrow(theResults.dt)

# We have to save the data later since we have made some modifications
		doSave = TRUE
	}

# We concatenate the data frames from the current downloaded period to the global data tables 
	theResults.dt = rbind(theResults.dt, theReports.l[[i - 1]]$results, fill = TRUE)
	theReactions.dt = rbind(theReactions.dt, theReports.l[[i - 1]]$patient.reaction, fill = TRUE)
	theDrugs.dt = rbind(theDrugs.dt, theReports.l[[i - 1]]$patient.drug, fill = TRUE)

# We inform the user once every 50 slice loaded
	if ((i - 1) %% 50 == 0) { loginfo("%d results up to slice %d", nrow(theResults.dt), i - 1) }
}

loginfo("Downloaded %d reports...", nrow(theResults.dt))

# We test if the data has been modified within the previous loop
if (doSave) {
	loginfo("Saving report in '%s' file", DATA_FILE)
	save(theReports.l, file = DATA_FILE)
}


# We change the type of some columns in the data table in order to process the data correctly 
theResults.dt$receivedate = as.Date(strptime(theResults.dt$receivedate, "%Y%m%d"))
theResults.dt$receiptdate = as.Date(strptime(theResults.dt$receiptdate, "%Y%m%d"))
theResults.dt$reportduplicate.duplicatesource = as.factor(theResults.dt$reportduplicate.duplicatesource)

# We aggregate the results data table by country to order the countries by number of reports
theTopCountries.dt = theResults.dt %>% group_by(occurcountry) %>% summarise(n = length(reportid)) %>% arrange(-n)

# We aggregate the reaction data table by reaction type to order the reactions by number of reports
theTopReactions.dt = theReactions.dt %>% group_by(reactionmeddrapt) %>% summarise(n = length(reportid)) %>% arrange(-n)

# We set up some keys to our data tables to be able to join and search them efficiently 
setkey(theResults.dt, reportid)
setkey(theReactions.dt, reportid)
setkey(theDrugs.dt, reportid)

# We use the data table merge syntax to combine the results data table and the reaction data table by their reportid
theFullReactions.dt = theReactions.dt[theResults.dt]

# We use the data table merge syntax to combine the reaction data table and the drug data table by their reportid
theFullDrugs.dt = theDrugs.dt[theReactions.dt, allow.cartesian = TRUE]

