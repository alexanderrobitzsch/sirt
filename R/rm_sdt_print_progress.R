## File Name: rm_sdt_print_progress.R
## File Version: 0.02
## File Last Change: 2017-10-02 18:40:40

rm_sdt_print_progress <- function( dev, dev0, c.rater, c.rater0, d.rater, d.rater0,
	tau.item, tau.item0, a.item, a.item0, mu, sigma, iter, digits_parm=6, digits_trait=3, digits_deviance=4 )
{
	#--- print deviance
	rm_facets_print_progress_deviance( dev=dev, dev0=dev0, digits_deviance=digits_deviance, iter=iter ) 
	#--- item and rater parameters
	rm_facets_print_progress_parameter( parm=c.rater, parm0=c.rater0, parmlabel="c.rater", digits_parm=digits_parm )
	rm_facets_print_progress_parameter( parm=d.rater, parm0=d.rater0, parmlabel="d.rater", digits_parm=digits_parm )
	rm_facets_print_progress_parameter( parm=tau.item, parm0=tau.item0, parmlabel="tau.item", digits_parm=digits_parm )									
	rm_facets_print_progress_parameter( parm=a.item, parm0=a.item0, parmlabel="a.item", digits_parm=digits_parm )									
	#--- trait distribution
	rm_facets_print_progress_trait_distribution( parm=mu, parmlabel="Trait M  ", digits_trait=digits_trait )	
	rm_facets_print_progress_trait_distribution( parm=sigma, parmlabel="Trait SD ", digits_trait=digits_trait )	
	utils::flush.console()			
}
		
