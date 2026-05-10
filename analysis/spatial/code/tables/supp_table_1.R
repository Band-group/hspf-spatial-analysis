library( dplyr )

clean = readr::read_csv( "input/cleanHbSdata.csv", guess_max = 100000 )
colnames(clean)[colnames(clean) == 'source...17'] = "source_type"

piel = readr::read_csv( "input/HbS_survey.csv" )
extended = readr::read_csv( "input/HbSgooglesheet.csv" )

table( clean$Dataset )
dim(clean)
original = (
	clean
	%>% filter( Dataset == 'original' )
	%>% left_join(
		piel[, c("id", "country", "citation", "area_type" )],
		by = c( ID_Piel_OR_PUBMED = "id" )
	)
	%>% mutate(
		spatial_accuracy = `area_type`
	)
)
dim(original)

extended1 = (
	clean
	%>% filter( Dataset == 'extended' & `source_type` == 'genotyping' )
	%>% left_join(
		extended %>% filter( !is.na( hbaa )) %>% select( PMID, hbaa, hbas, country = `ADM-0 (country)`, spatial_accuracy = `Spatial accuracy` ),
		by = c( ID_Piel_OR_PUBMED = "PMID", "hbaa", "hbas" )
	)
	%>% mutate(
		citation = DOI
	)
)
dim(extended1)

extended2 = (
	clean
	%>% filter( Dataset == 'extended' & `source_type` == 'blood_typing' )
	%>% left_join(
		extended %>% filter( !is.na( HbFA )) %>% select( PMID, HbFA, HbFAS, HbFS, country = `ADM-0 (country)`, spatial_accuracy = `Spatial accuracy` ),
		by = c( ID_Piel_OR_PUBMED = "PMID", "HbFA", "HbFAS", "HbFS" )
	)
	%>% mutate( citation = DOI )
)
dim(extended2)
{
	result = (
		dplyr::bind_rows(
			original,
			extended1,
			extended2
		)
	)
	result$source = factor( result$Dataset, levels = c( "original", "extended" ))
	levels( result$source ) = c(
		"Piel et al doi:10.1016/S0140-6736(12)61229-X",
		"Literature search 2011- 2024"
	)

	result = (
		result
		%>% transmute(
			source,
			country,
			`Reported longitude` = original_longitude,
			`Reported latitude` = original_latitude,
			`Adjusted longitude` = longitude,
			`Adjusted latitude` = latitude,
			`Location adjusted?` = "no",
			`Spatial accuracy` = spatial_accuracy,
			`Data type` = gsub( "_", " ", `source_type`, fixed = T ),
			blank1 = '',
			`HbAA` = hbaa,
			`HbAS` = hbas,
			`HbSS` = hbss,
			blank2 = '',
			`HbFA` = HbFA,
			`HbFAS` = HbFAS,
			`HbFS` = HbFS,
			blank3 = '',
			A, S, N,
			blank4 = '',
			`ID (Piel et al, or pubmed)` = ID_Piel_OR_PUBMED,
			`citation or doi` = citation
		)
	)
	result$`Location adjusted?`[
		(result$`Adjusted longitude` != result$`Reported longitude`)
		|
		(result$`Adjusted latitude` != result$`Reported latitude`)
	] = "yes"
}
readr::write_csv( result, "output/tables/table_S1_HbS_survey_points.csv" )
