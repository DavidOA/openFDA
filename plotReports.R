# We create the first graph that is a weekly time series of the numbers of reports coloured by country and facetted by the top reactions

thePlot = ggplot(theSubset.dt)

thePlot = thePlot + geom_histogram(aes(x = receiptdate, fill = occurcountry), binwidth = 7)
thePlot = thePlot + facet_grid(reactionmeddrapt ~ .)
thePlot = thePlot + ggtitle(sprintf("Distributions of weekly reports for top %d countries and top %d reactions", TOP_COUNTRIES, TOP_REACTIONS))
thePlot = thePlot + scale_fill_manual(values = getColours()[1 : TOP_COUNTRIES])

print(thePlot)
ggsave('weekly.png', width = 10, height  = 6)



# The second map is a heat map that is coloured based on the number of reports containing the pair of drugs found at the intersections

thePlot = ggplot(theHeatMap.dt)

thePlot = thePlot + geom_tile(aes(x = new_i, y = new_j, fill = N))
thePlot = thePlot + scale_x_continuous(breaks = seq.int(1, nrow(theNewIndices.dt)), labels = sapply(str_split(theNewIndices.dt$medicinalproduct, ' '), function(x) x[1]))
thePlot = thePlot + scale_y_continuous(breaks = seq.int(1, nrow(theNewIndices.dt)), labels = theNewIndices.dt$medicinalproduct)
thePlot = thePlot + theme(axis.text.x = element_text(angle = -90, hjust = 0.0, vjust = 0.5), axis.text = element_text(size = 6), panel.grid.minor = element_blank())

# We also show some evidence that one of the high colour drug interaction has been reported in the literature demonstrating the ability of this plot to detect abnormalities 

thePlot = thePlot + annotate('text', x = as.integer(theHeatMap.dt[J('VANCOMYCIN INJECTION, USP', 'CEFEPIME'), 'new_i']), y = as.integer(theHeatMap.dt[J('VANCOMYCIN INJECTION, USP', 'CEFEPIME'), 'new_j']), label = 'http://aac.asm.org/content/61/2/e02089-16.full', size = 2, colour = 'darkgrey')
thePlot = thePlot + scale_fill_gradient2(midpoint = (max(theHeatMap.dt$N) + MIN_INTERACTIONS) / 2, low = muted('blue'), mid = getColours()[12], high = muted('red'))
thePlot = thePlot + theme(axis.title = element_blank())
thePlot = thePlot + ggtitle(sprintf("Drug-drug interactions reported more than %d times", MIN_INTERACTIONS))

print(thePlot)
ggsave('interactions.png', width = 10, height  = 6)

theLayout = layout_nicely(theGraph)

theVertices.dt$x = theLayout[, 1]
theVertices.dt$y = theLayout[, 2]

theEdges.dt = as.data.table(merge(merge(theAssociations.dt, theVertices.dt, by.x = c('medicinalproduct'), by.y = 'name'), theVertices.dt, by.x = ('reactionmeddrapt'), by.y = ('name')))
setnames(theEdges.dt, c('type.x', 'x.x', 'y.x', 'type.y', 'x.y', 'y.y'), c('from_type', 'from_x', 'from_y', 'to_type', 'to_x', 'to_y'))
setkeyv(theEdges.dt, c('medicinalproduct', 'reactionmeddrapt'))


# We annotate some of the edge from evidence found in the literature. Again the potential of this graph is demonstrated

theEdges.dt$label = as.character(NA)
theEdges.dt[J('ARANESP', 'Death'), 'label'] = 'https://www.reuters.com/article/us-amgen-aranesp/amgen-details-higher-death-risk-in-aranesp-trial-idUSN1323853720070416'
theEdges.dt[J('ENBREL', 'Off label use'), 'label'] = 'https://www.fiercepharma.com/sales-and-marketing/amgen-to-shell-out-71m-to-states-aranesp-enbrel-off-label-settlement'

# We create a proper graph with vertices and segments from the data calculated above

thePlot = ggplot(theEdges.dt) 
thePlot = thePlot + geom_segment(aes(x = from_x, y = from_y, xend = to_x, yend = to_y, colour = N))
thePlot = thePlot + scale_colour_gradient2(midpoint = (max(theEdges.dt$N) + MIN_ASSOCIATIONS) / 2, low = muted('blue'), mid = getColours()[12], high = muted('red'))
thePlot = thePlot + geom_text(data = theVertices.dt, aes(x = x, y = y, label = name), size = 2)
thePlot = thePlot + geom_point(data = subset(theEdges.dt, !is.na(label)), aes(x = (from_x + to_x) / 2, y = (from_y + to_y) / 2), size = 2, color = 'white')
thePlot = thePlot + geom_text(data = subset(theEdges.dt, !is.na(label)), aes(x = (from_x + to_x) / 2, y = (from_y + to_y) / 2, label = label ), size = 2, color = 'darkgrey')
thePlot = thePlot + coord_flip()
thePlot = thePlot + theme(axis.ticks = element_blank(), axis.title = element_blank(), axis.text = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
thePlot = thePlot + ggtitle(sprintf("Drug-reaction associations reported more than %d times", MIN_ASSOCIATIONS))

print(thePlot)
ggsave('associations.png', width = 10, height  = 6)


