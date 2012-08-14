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

		if gstruct.metal_source eq 'Not in M10 Table 8' then goto,skip

		; convert the hi into mass surface density
		fac = 1.36*mh*pc*pc/ms ; accounts for He
		hi = gstruct.hi*fac
		hiunc = gstruct.hi_unc*fac

		; scale MW DGR with metallicity
		mw_oh = 10.^(mwmetal-12d)
		oh = 10.^(gstruct.metal)
		scldgr = mwdgr*oh/mw_oh

		; convert dust mass surface density into gas with DGR
		


		stop

		skip:
	endfor

end
