pro n1_dgrvsmetals,$
	targetlist=targetlist,$
	just=just,$
	mwdgr=mwdgr,$
	mwmetal=mwmetal,$
	r25bin=r25bin,$
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
	if ntarg gt 1 then BEGIN
		clr = scale_vector(findgen(ntarg),0,255)
	endif else BEGIN
		clr = 255
	endelse

;	plot,findgen(10),findgen(10),/xlog,/ylog,$
;		xr=[3d-3,3d-2],yr=[1d-1,1d3],/xs,/ys,/nodata
	plot,findgen(10),findgen(10),/ylog,$
		xr=[7.7,9.0],yr=[1d-1,1d3],/ys,/xs,/nodata

	j=0
	loadct,4
	for i=0,ntarg-1 do BEGIN

		; restore the sampled structure
		restore,datadir+galname[i]+'_samp.sav'

		if gstruct.metal_source eq 'Not in M10 Table 8' then BEGIN
			print,'Skipping '+galname[i]
			goto,skip
		endif

		r25 = gstruct.r25
		co = gstruct.co_hera
		hi = gstruct.hi
		dust = gstruct.sigdust
		metal = gstruct.metal

		ok = where(finite(co) and finite(hi) and finite(dust))

		; make radial profiles
		bin_prof,r25[ok],co[ok],binsize=r25bin,$
			medprof=medcoprof,meanprof=meancoprof,xmid_bin=xout,$
			stdprof=stdcoprof

		bin_prof,r25[ok],hi[ok],binsize=r25bin,$
			medprof=medhiprof,meanprof=meanhiprof,xmid_bin=xout2,$
			stdprof=stdhiprof

		bin_prof,r25[ok],dust[ok],binsize=r25bin,$
			medprof=meddustprof,meanprof=meandustprof,xmid_bin=xout3,$
			stdprof=stddustprof

		bin_prof,r25[ok],metal[ok],binsize=r25bin,$
			medprof=medmetalprof,meanprof=meanmetalprof,xmid_bin=xout4
			
		; convert the hi into mass surface density
		fac = 1.36*mh*pc*pc/ms ; accounts for He
		hi = meanhiprof*fac
		hi_std = stdhiprof*fac

		; scale MW DGR with metallicity
		mw_oh = 10.^(mwmetal-12d)
		oh = 10.^(meanmetalprof-12d)
		scldgr = mwdgr*oh/mw_oh
		scldgr_x2 = (mwdgr*2d)*oh/mw_oh
		scldgr_d2 = (mwdgr/2d)*oh/mw_oh

		; convert dust mass surface density into gas with DGR
		totgas = meandustprof/scldgr	
		totgas_std = stddustprof/scldgr
		totgas_dgrx2 = meandustprof/scldgr_x2
		totgas_dgrd2 = meandustprof/scldgr_d2

		h2 = totgas-hi
		h2_dgrx2 = totgas_dgrx2 - hi
		h2_dgrd2 = totgas_dgrd2 - hi

		; conversion factor!
		aco = h2/(meancoprof/0.7d)
		aco_dgrx2 = h2_dgrx2/(meancoprof/0.7d)
		aco_dgrd2 = h2_dgrd2/(meancoprof/0.7d)

;		oplot,scldgr,aco,ps=3
		oplot,meanmetalprof,aco,color=clr[i],ps=-4
		oplot,meanmetalprof,aco_dgrx2,color=clr[i],linestyle=2
		oplot,meanmetalprof,aco_dgrd2,color=clr[i],linestyle=2

		xyouts,[0.8],[0.8]-(j*0.05),galname[i],color=clr[i],/normal
		j += 1
		skip:
	endfor
stop
end
