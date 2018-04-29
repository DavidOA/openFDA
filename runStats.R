# We create a subset of the full reaction data table limited to the top countries and top reactions for the sake of simplicity
theSubset.dt = as.data.table(theFullReactions.dt %>% filter(occurcountry %in% theTopCountries.dt$occurcountry[1 : TOP_COUNTRIES], reactionmeddrapt %in% theTopReactions.dt$reactionmeddrapt[1 : TOP_REACTIONS]))

# We prepare a data table to hold the pairwise comparisons of the top countries
thePairWises.dt = as.data.table(NULL)

# We iterate over the paiwise combination of countries and run a chi square test to compare their proportions
for (i in combn(theTopCountries.dt$occurcountry[1 : TOP_COUNTRIES], 2, simplify = FALSE)) {
	theTable = xtabs(~ occurcountry + reactionmeddrapt, data = subset(theSubset.dt, occurcountry %in% i))
	theChisq = suppressWarnings(chisq.test(theTable))

# We concatenate all the results of the pairwise comparisons into a single data table
	thePairWises.dt = rbind(thePairWises.dt, data.table('country_one' = i[1], 'country_two' = i[2], pvalue = theChisq$p.value))
} 

# We order the data table by the pvalue column
setkey(thePairWises.dt, 'pvalue')

# We add some significance markers to make the results clearer to the users 
thePairWises.dt$stars = ifelse(thePairWises.dt$pvalue < 0.001, '***', ifelse(thePairWises.dt$pvalue < 0.01, '**', ifelse(thePairWises.dt$pvalue < 0.05, '*', '')))

print(thePairWises.dt)

# We pick a separator character that is not present in the medicinalproduct column
SEP = "`"

# We actually test that the chosen separator is not found the medicinalproduct column
stopifnot(!any(grepl(SEP, theDrugs.dt$medicinalproduct)))


# We return a list containing the number of unique products within the vector of medicinal product x passed to the function, a string that is the concatenations of all the products, the original vector x, and their respective alphabetical order
# The function returns NULL if there are less than 2 products in the vector x since we are interested in interactions here
returnDrugs = function(x) {
	x = unique(x)
	n = length(x)
	p = paste0(x, collapse = SEP)

	if (n >= 2) {
		l = list(rep(n, n), rep(p, n), x, order(x))
		names(l) = c('numproducts', 'listproducts', 'medicinalproduct', 'seqproduct')
		return(l)
	} else {
		return(NULL)
	}
}

# We use the data table syntax to call the returnDrugs function to every group of medicinalproduct by reportid
# We discard all the medicinalproduct that do not have an authorisation number
theInteractions.dt = as.data.table(subset(theDrugs.dt, !is.na(drugauthorizationnumb) & !is.na(medicinalproduct))[, returnDrugs(medicinalproduct), by = c('reportid')])

# We create a data table from all the medicinalproducts that are present in the interactions data table
theUniques.dt = data.table(medicinalproduct = unique(theInteractions.dt$medicinalproduct), key = c('medicinalproduct'))

# We add a second column to create a fast index 
theUniques.dt$index = seq.int(1, nrow(theUniques.dt))

# We merge this two column data table to the interactions data table so that we have now the indices in that table
theInteractions.dt = merge(theInteractions.dt, theUniques.dt, by = c('medicinalproduct'))

loginfo("%d drug-drug interactions detected", nrow(theInteractions.dt))

# We return all the combinations for the list of indices passed in the x argument as a 2 column matrix
findDrugs = function(x) {
	r = t(combn(x[order(x)], 2, simplify = TRUE))
	r = cbind(r, as.integer(r[, 1] == r[, 2]))

	stopifnot(all(r[, 3] == 0))

	l = list(r[, 1], r[, 2])

	names(l) = c('j', 'i')

	return(l)
}

# We call the findDrugs function for every report id in order to generate all the pairwise combinations of drugs within each report
theCells.dt = theInteractions.dt[,  findDrugs(index), by = c('reportid')]

# We check that none of the interactions are with themselves
stopifnot(all(theCells.dt[['i']] != theCells.dt[['j']]))

# We aggregate all the possible interactions to count their occurences in the reports
theHeatMap.dt = theCells.dt[ ,.N, by = c('i', 'j')] 

# We collect the indices of the drugs having at least MIN_INTERACTIONS reports
theNewIndices.dt = data.table(old_i = unique(c(subset(theHeatMap.dt, N >= MIN_INTERACTIONS)[['i']], subset(theHeatMap.dt, N >= MIN_INTERACTIONS)[['j']])), key = c('old_i'))
theNewIndices.dt$new_i = seq.int(1, nrow(theNewIndices.dt))
theNewIndices.dt = merge(theNewIndices.dt, theUniques.dt, by.x = c('old_i'), by.y = c('index'))

# We remap the old indices from all the drugs to the indices of the subset of the drugs from above
theHeatMap.dt = as.data.table(merge(theHeatMap.dt, theNewIndices.dt, by.x = c('i'), by.y = c('old_i'), all.x = TRUE) %>% filter(!is.na(new_i)))
theHeatMap.dt = as.data.table(merge(theHeatMap.dt, theNewIndices.dt, by.x = c('j'), by.y = c('old_i'), all.x = TRUE) %>% filter(!is.na(new_i.y)))
setkeyv(theHeatMap.dt, c('medicinalproduct.x', 'medicinalproduct.y'))

# We rename the new column for the sake of clarity
setnames(theHeatMap.dt, c('new_i.x', 'new_i.y'), c('new_i', 'new_j'))

# We check that none of the interactions are with themselves
stopifnot(all(theHeatMap.dt[['new_i']] != theHeatMap.dt[['new_j']]))

loginfo("%d drug-drug interactions retained above %d threshold", nrow(theHeatMap.dt), MIN_INTERACTIONS)

# We create a data table where a reaction is associated to a medicinal product in at least MIN_ASSOCIATION reports. This data table contains all the reaction-medicinalproduct pairs of interest, which will be used as edges in the graph that we create below
theAssociations.dt = as.data.table(theFullDrugs.dt %>% filter(!is.na(drugauthorizationnumb) & !is.na(medicinalproduct)) %>% group_by(medicinalproduct, reactionmeddrapt) %>% summarise(N = length(reportid)) %>% filter(N >= MIN_ASSOCIATIONS))

loginfo("%d drug-reaction associations retained above %d threshold", nrow(theAssociations.dt), MIN_ASSOCIATIONS)

# We extract the vertices from the edges created above
theVertices.dt = rbind(data.table('name' = unique(theAssociations.dt$medicinalproduct), 'type' = TRUE), data.table('name' = unique(theAssociations.dt$reactionmeddrapt), 'type' = FALSE))

setkeyv(theVertices.dt, c('name'))

# We create a graph from the edges and vertices information found in the data tables 
theGraph = graph_from_data_frame(theAssociations.dt, vertices = theVertices.dt, directed = FALSE)
