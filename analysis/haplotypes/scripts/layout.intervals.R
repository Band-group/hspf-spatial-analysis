#' Layout intervals
#'
#' This function takes some intervals expressed as start and end coordinates
#' and computes a level for each interval, such that intervals on the same level
#' do not overlap, taking into account an additional spacer if desired.
#' @param intervals a data frame or similar object with start end end columns
#' @param spacer a named vector of numbers to add to the start and end of each interval
#' @param columns the names of the columns in intervals denoting start and end coordinates.
layout.intervals <- function(
	intervals,
	spacer = c( start = 0, end = 0 ),
	columns = c( "start" = "start", "end" = "end" )
) {
	if( nrow( intervals ) == 0 ) {
		return( integer(0) )
	}
	intervals = intervals[ order( intervals[[ columns['start'] ]] ), ]
	level = rep( NA, nrow( intervals ))
	level[1] = 1
	level.endpoints = ( intervals[1, columns['end']] + spacer['end'] )
	if( nrow( intervals ) > 1 ) {
		for( i in 2:nrow(intervals)) {
			region = intervals[i,]
			end_with_spacer = region[[ columns['end'] ]] + spacer['end'] ;
			# Try to put gene in an existing level
			for( l in 1:length(level.endpoints) ) {
				if( (region[[ columns['start'] ]] - spacer['start']) > level.endpoints[l] ) {
					level[i] = l ;
					level.endpoints[l] = end_with_spacer
					break ;
				}
			}
		
			# Otherwise add a new level
			if( is.na( level[i] )) {
				level.endpoints = c( level.endpoints, end_with_spacer )
				level[i] = length( level.endpoints ) ;
			}
		}
	}
	return( level ) ;
}
