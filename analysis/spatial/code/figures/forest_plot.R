library( ggplot2 )
X = readr::read_tsv( "output/hspf/fixed-r0=25.0-sigma0=0.6-fc=none/all_hspf_analyses_summary.tsv" )
X$area = factor( X$area, levels = rev(c( "global", "africa", "waf", "eaf", "gambia", "gambia+senegal", "mali", "ghana+burkina+togo", "DRC", "tanzania" )))

country.names = c(
	"global" = "Global",
	"africa" = "Africa",
	"waf" = "west Africa",
	"eaf" = "east Africa",
	"gambia" = "Gambia",
	"gambia+senegal" = "Gambia + Senegal",
	"mali" = "Mali",
	"ghana" = "Ghana",
	"ghana+burkina+togo" = "Ghana, Burkina Faso and Togo",
	"DRC" = "DRC",
	"tanzania" = "Tanzania"
)
X$country = factor( country.names[ as.character(X$area) ], levels = rev(country.names) )
p = (
	ggplot(  data = X %>% filter( allele == 'Pfsa1_+' & model == 'bym2' ))
	+ geom_segment( aes( x = pmax(beta.q2.5,-10), xend = pmin(beta.q97.5,30), y = country, yend = country ), linewidth = 4, col = 'grey' )
#	+ facet_grid( allele ~ model )
	+ xlim( -10, 30 )
	+ geom_point(aes( x =  beta.mean, y = country), size = 4) + theme_minimal(16)
	+ geom_vline( xintercept = 0, linetype = 2 )
	+ xlab( "Slope estimate (log-odds scale)" )
	+ theme(axis.text.y = element_text( size = 14, colour = 'black' ))
)
print(p)

ggsave( p, file = "output/figures/for_slides/forest_plot.pdf", width = 7, height = 5 )
