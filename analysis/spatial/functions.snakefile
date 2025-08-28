def srcdir(x):
	return x

# dict_product from StackOverflow:
# https://stackoverflow.com/questions/5228158/cartesian-product-of-a-dictionary-of-lists/40623158#40623158
import itertools
def dict_product(dicts):
	"""
	>>> list(dict_product(dict(number=[1,2], character='ab')))
	[{'character': 'a', 'number': 1},
	 {'character': 'a', 'number': 2},
	 {'character': 'b', 'number': 1},
	 {'character': 'b', 'number': 2}]
	"""
	return list(dict(zip(dicts, x)) for x in itertools.product(*dicts.values()))


def get_area_definitions( names ):
	defs = {
		'global': None,
		'africa': [
			'Gambia', 'Senegal', 'Mali', 'Benin', 'Burkina Faso', 'Ivory Coast', 'Ghana', 'Guinea', 'Mauritania', 'Nigeria', 'Senegal', 'Togo',
			'Central African Republic', 'Angola', 'Cameroon', 'Gabon', 'Republic of the Congo', 'Democratic Republic of the Congo',
			'Ethiopia', 'Kenya', 'Madagascar', 'Malawi', 'Mozambique', 'Rwanda', 'Uganda', 'United Republic of Tanzania', 'Zambia'
		],
		# Western africa region matching Pfsa2/4 distribution split
		'waf': [ 'Mauritania', 'Senegal', 'Gambia', 'Guinea', 'Mali', 'Burkina Faso', 'Ivory Coast', 'Ghana', 'Togo', 'Benin', 'Nigeria', 'Cameroon', 'Gabon'  ],
		'wwaf': [ 'Mauritania', 'Senegal', 'Gambia', 'Guinea', 'Mali' ],
		'ewaf': [ 'Burkina Faso', 'Ivory Coast', 'Ghana', 'Togo', 'Benin', 'Nigeria', 'Cameroon', 'Gabon'  ],
		# Central africa region.  Not needed.
		'caf': [ 'Democratic Republic of the Congo', 'Zambia', 'Gabon' ],
		# Eestern africa region matching Pfsa2/4 distribution split
		'DRC+eaf': [ 'Democratic Republic of the Congo', 'Sudan', 'Ethiopia', 'Kenya', 'Rwanda', 'Uganda', 'Malawi', 'Zambia', 'Mozambique', 'United Republic of Tanzania',  'Madagascar' ],
		'eaf': [ 'Kenya', 'Rwanda', 'Uganda', 'Malawi', 'Zambia', 'Mozambique', 'United Republic of Tanzania' ],
		#
		'gambia+senegal': [ 'Gambia', 'Senegal' ],
		'mali': [ 'Mali' ],
		'ghana': [ 'Ghana' ],
		'ghana+burkina+togo': [ 'Ghana', 'Burkina Faso', 'Togo' ],
		'ghana+burkina+togo+benin+ivorycoast': [ 'Ghana', 'Burkina Faso', 'Togo', 'Ivory Coast', 'Benin' ],
		'uganda': [ 'Uganda' ],
		'tanzania': [ 'United Republic of Tanzania' ],
		'tanzania+kenya+uganda+rwanda': [ 'United Republic of Tanzania', 'Kenya', 'Uganda', 'Rwanda' ],
		'DRC': [ 'Democratic Republic of the Congo' ]
	}
	result = {}
	for name in names:
		result[name] = defs[name]
	return result

def remove_keys( dictionary, keys_to_remove ):
	result = dictionary.copy()
	for key in keys_to_remove:
		result.pop( key )
	return result

