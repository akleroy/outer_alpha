pro n1_dgrvsmetals,$
	targetlist=targetlist,$
	just=just,$
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

	@code/constants.bat	

	datadir = 'data/'
	
	if keyword_set(just) eq 0 then BEGIN
		readcol,targetlist,galname,format='A'
	endif else BEGIN
		galname = [just]
	endelse

	ntarg = n_elements(galname)

;	plot,findgen(10),findgen(10),/xlog,/ylog,$
;		xr=[3d-3,3d-2],yr=[1d-1,1d3],/xs,/ys,/nodata
	plot,findgen(10),findgen(10),/ylog,$
		xr=[8.0,9.0],yr=[1d-1,1d3],/ys,/xs,/nodata

	for i=0,ntarg-1 do BEGIN

		; restore the sampled structure
		restore,datadir+galname[i]+'_samp.sav'

		if gstruct.metal_source eq 'Not in M10 Table 8' then BEGIN
			print,'Skipping '+galname[i]
			goto,skip
		endif

		; make radial profiles
		bin_prof,gstruct.r25,gstruct.co_hera,binsize=0.1,$
			medprof=medcoprof,meanprof=meancoprof,xmid_bin=xout

			stop


		; convert the hi into mass surface density
		fac = 1.36*mh*pc*pc/ms ; accounts for He
		hi = gstruct.hi*fac
		hiunc = gstruct.hi_unc*fac

		; scale MW DGR with metallicity
		mw_oh = 10.^(mwmetal-12d)
		oh = 10.^(gstruct.metal-12d)
		scldgr = mwdgr*oh/mw_oh

		; convert dust mass surface density into gas with DGR
		totgas = gstruct.sigdust/scldgr	
		
		h2 = totgas-hi

		; conversion factor!
		aco = h2/(gstruct.co_hera/0.7d)

;		oplot,scldgr,aco,ps=3
		oplot,gstruct.metal,aco,ps=3

		stop

		skip:
	endfor

end
