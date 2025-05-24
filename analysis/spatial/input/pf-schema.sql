CREATE TABLE `by_sample` (
  `source` TEXT,
  `study` TEXT,
  `datatype` TEXT,
  `country` TEXT,
  `year` TEXT,
  `site` TEXT,
  `latitude` REAL,
  `longitude` REAL,
  `ID` TEXT,
  `exclude` TEXT DEFAULT 'no',
  `locus` TEXT,
  `ref` INTEGER,
  `mixed` INTEGER,
  `nonref` INTEGER
) ;