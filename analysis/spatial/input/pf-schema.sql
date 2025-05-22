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
  `N` INTEGER,
  `Pfsa1:ref` REAL,
  `Pfsa1:nonref` REAL,
  `Pfsa2:ref` REAL,
  `Pfsa2:nonref` REAL,
  `Pfsa3:ref` REAL,
  `Pfsa3:nonref` REAL,
  `Pfsa4:ref` REAL,
  `Pfsa4:nonref` REAL
, exclude DEFAULT 'no');

CREATE TABLE `by_site` (
  `source` TEXT,
  `study` TEXT,
  `datatype` TEXT,
  `country` TEXT,
  `year` TEXT,
  `site` TEXT,
  `latitude` REAL,
  `longitude` REAL,
  `N` INTEGER,
  `Pfsa1:ref` REAL,
  `Pfsa1:nonref` REAL,
  `Pfsa2:ref` REAL,
  `Pfsa2:nonref` REAL,
  `Pfsa3:ref` REAL,
  `Pfsa3:nonref` REAL,
  `Pfsa4:ref` REAL,
  `Pfsa4:nonref` REAL
, exclude DEFAULT 'no');

