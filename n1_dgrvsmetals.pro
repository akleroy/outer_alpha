pro n1_dgrvsmetals,$
	targetlist=targetlist,$
	mwdgr=mwdgr,$
	mwmetal=mwmetal,$
	outfile=outfile

;+
;
; First hack at code to use linear DGR(Z) for estimating a_co.
; Code takes galaxies with gradients from Moustakas+ 2010, Tb 8
; scales the MW DGR and Z and then turns dust mass surface densities
; into gas mass surface densities.  Then the HI is subtracted and 
; we use the remainder to calculate a_co.
;
; Since we're doing this with radial metallicity gradients, we will
; do the whole thing with radial profiles.
;
; Written by Karin - 8/14/12
;
;-

	datadir = 'data/'
	readcol,targetlist,galname,format='A'
	ntarg = n_elements(galname)

	for i=0,ntarg-1 do BEGIN

		; restore the sampled structure
		restore,datadir+galname[i]+'_samp.sav'

	endfor

end
