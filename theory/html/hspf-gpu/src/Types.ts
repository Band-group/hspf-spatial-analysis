export interface LatLong {
	latitude: number,
	longitude: number
} ;

export interface PixelCoords {
	x: number,
	y: number
} ;

export interface LatLongCoords {
	latlong: LatLong,
	xy: PixelCoords
} ;

export type PfsaDataKey =
	'pfsa1p' | 'pfsa1m' | 'pfsa1N' |
	'pfsa2p' | 'pfsa2m' | 'pfsa2N' |
	'pfsa3p' | 'pfsa3m' | 'pfsa3N' |
	'pfsa4p' | 'pfsa4m' | 'pfsa4N' |
	'pfsa13mm' | 'pfsa13mp' | 'pfsa13pm' | 'pfsa13pp' | 'pfsa13N';

export interface PfsaCounts {
	country: Country,
	admin1: string,
	latlong: LatLong,
	xy: PixelCoords,
	pfsa1p: number,
	pfsa1m: number,
	pfsa1N: number,
	pfsa2p: number,
	pfsa2m: number,
	pfsa2N: number,
	pfsa3p: number,
	pfsa3m: number,
	pfsa3N: number,
	pfsa4p: number,
	pfsa4m: number,
	pfsa4N: number,
	pfsa13mm: number,
	pfsa13mp: number,
	pfsa13pm: number,
	pfsa13pp: number,
	pfsa13N: number
} ;

export const colours = {
	"Morocco": "#292933",
	"Mauritania": "#090953",
	"Gambia": "#0c0c83",
	'Senegal': "#2323f6",
	'Guinea-Bissau': "#0000CD",
	'Guinea': "#3a3a9f",
	"Mali": "#42426F",
	"Sierra Leone": "#42628F",
	"Liberia": "#42628F",
	"Burkina_Faso": "#377EB8",
	"Burkina Faso": "#377EB8",
	"IvoryCoast": "#2ecdab",
	"Ivory Coast": "#2ecdab",
	"Cote_dIvoire": "#2ecdab",
	"Cote d'Ivoire": "#2ecdab",
	"Ghana": "#03B4CC",
	"Benin": "#03cc53",
	"Nigeria": "#a57d0f",
	"Niger": "#c57d0f",
	"Chad": "#fecb00",
	"Cameroon": "#007a5e",
	"Gabon": "#009E60",
	"Republic of the Congo": "#dc241f",
	"Democratic_Republic_of_the_Congo": "#ef3340",
	"Democratic Republic of the Congo": "#ef3340",
	"Congo": "#ef3340",
	"Rwanda": "#e5be01",
	"Zambia": "#A4081C",
	"Sudan": "#c59d0f",
	"Uganda": "#fcdc04",
	"Malawi": "#A65628",
	"Tanzania": "#EE5C42",
	"United Republic of Tanzania": "#EE5C42",
	"Mozambique": "#EE5C42",
	"Kenya": "#FF7F00",
	"Ethiopia": "#d1cd0c",
	"Madagascar": "#c800ff",
	'Bangladesh': "#444444",
	'Myanmar': "#444444",
	'Laos': "#444444",
	'Thailand': "#444444",
	'Cambodia': "#444444",
	'Vietnam': "#444444",
	'Indonesia': "#444444",
	'PNG': "#444444",
	'South Africa': "#23f623",
	'eSwatini': "#23f623",
	"other": "#AAAAAA"
} as const;

export type Country = keyof typeof colours;

