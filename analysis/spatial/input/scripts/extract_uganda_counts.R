library( dplyr )
paths = list(
	data = "data/uganda/pfsa_data_uganda_wgs.tsv"
)

X = (
	readr::read_tsv( paths$data )
)

X$Pf3D7_02_v3_631190_T_A = c(
	'T|T' = 0,
	'T/T' = 0,
	'T/A' = NA,
	'T|A' = NA,
	'A/A' = 1,
	'A|A' = 1
)[ X$Pf3D7_02_v3_631190_T_A ]

X$Pf3D7_02_v3_814288_C_T = c(
	'C/C' = 0,
	'C|C' = 0,
	'C/T' = NA,
	'C|T' = NA,
	'T/T' = 1,
	'T|T' = 1
)[ X$Pf3D7_02_v3_814288_C_T ]

X$Pf3D7_04_v3_1121472_T_A = c(
	'T|T' = 0,
	'T/T' = 0,
	'T/A' = NA,
	'T|A' = NA,
	'A/A' = 1,
	'A|A' = 1
)[ X$Pf3D7_04_v3_1121472_T_A ]

X$Pf3D7_11_v3_1058035_T_A = c(
	'T|T' = 0,
	'T/T' = 0,
	'T/A' = NA,
	'T|A' = NA,
	'A/A' = 1,
	'A|A' = 1
)[ X$Pf3D7_11_v3_1058035_T_A ]

by_sample = (
	X
	%>% mutate(
		source = "Greenwood_Uganda_2017-2022",
		study = "Greenwood_Uganda_2017-2022",
		country = "Uganda",
		site,
		latitude,
		longitude,
		ID = sample_name, 
		N = 1,
		`Pfsa1:ref`    = 1 - Pf3D7_02_v3_631190_T_A,
		`Pfsa1:nonref` = Pf3D7_02_v3_631190_T_A,
		`Pfsa2:ref`    = 1 - Pf3D7_02_v3_814288_C_T,
		`Pfsa2:nonref` = Pf3D7_02_v3_814288_C_T,
		`Pfsa3:ref` = 1 - Pf3D7_11_v3_1058035_T_A,
		`Pfsa3:nonref` = Pf3D7_11_v3_1058035_T_A,
		`Pfsa4:ref` = 1 - Pf3D7_04_v3_1121472_T_A,
		`Pfsa4:nonref` = Pf3D7_04_v3_1121472_T_A,
		exclude = "no"
	)
	%>% select(
		source, study, country, site, latitude, longitude,
		ID, N,
		`Pfsa1:ref`, `Pfsa1:nonref`,
		`Pfsa2:ref`, `Pfsa2:nonref`,
		`Pfsa3:ref`, `Pfsa3:nonref`,
		`Pfsa4:ref`, `Pfsa4:nonref`,
		`exclude`
	)
)

by_site = (
	by_sample
	%>% group_by(
		source, study, country, site, latitude, longitude
	) %>% summarise(
		N = n(),
		'Pfsa1:ref' = sum(`Pfsa1:ref`, na.rm = T ), `Pfsa1:nonref` = sum( `Pfsa1:nonref`, na.rm = T ),
		'Pfsa2:ref' = sum(`Pfsa2:ref`, na.rm = T ), `Pfsa2:nonref` = sum( `Pfsa2:nonref`, na.rm = T ),
		'Pfsa3:ref' = sum(`Pfsa3:ref`, na.rm = T ), `Pfsa3:nonref` = sum( `Pfsa3:nonref`, na.rm = T ),
		'Pfsa4:ref' = sum(`Pfsa4:ref`, na.rm = T ), `Pfsa4:nonref` = sum( `Pfsa4:nonref`, na.rm = T ),
		exclude = max( exclude )
	)
)

db = DBI::dbConnect( RSQLite::SQLite(), "results/genotypes/hbs-pf-v2.sqlite" )
DBI::dbExecute( db, "DELETE FROM by_site WHERE source == 'Greenwood_Uganda_2017-2022' ")
DBI::dbWriteTable( db, "by_site", by_site, overwrite = FALSE, append = TRUE )
DBI::dbExecute( db, "DELETE FROM by_sample WHERE source == source == 'Greenwood_Uganda_2017-2022' ")
DBI::dbWriteTable( db, "by_sample", by_sample, overwrite = FALSE, append = TRUE )
DBI::dbDisconnect( db )
