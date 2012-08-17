pro paper1_dgrvsmetal,$
	saved=saved,$
	targetlist=targetlist,$
	just=just,$
	mwmetal=mwmetal,$
	mwdgr=mwdgr,$
	goodsol=goodsol,$
	thetacut=thetacut

	; restore solutions
	restore,saved

	; if you just want to look at one galaxy
	if keyword_set(just) then BEGIN
		galname = [just]
	endif else BEGIN
		readcol,targetlist,galname,format='A'
	endelse

	; color code by galaxy
	ntarg = n_elements(galname)
    if ntarg gt 1 then BEGIN
       clr = scale_vector(findgen(ntarg),0,255)
    endif else BEGIN
       clr = 255
    endelse

	; generate plot stuff
	if keyword_set(mwmetal) eq 0 then mwmetal = 8.5
	if keyword_set(mwdgr) eq 0 then mwdgr = 1d-2
	plot_oh = findgen(20)/10 + 7.5
	oh = 10^(plot_oh-12d)
	mw_oh = 10.^(mwmetal-12d)
	scldgr = mwdgr*oh/mw_oh
	scldgr_d2 = 0.5*mwdgr*oh/mw_oh
	scldgr_x2 = 2.*mwdgr*oh/mw_oh

	; where to find sampled galaxies
	datadir = 'data/'

	; can skip all the stuff if you just want good solutions
	if keyword_set(goodsol) then BEGIN
		restore,'allgoodsol.sav'
		goto,goodsol
	endif

	; initial plot
	plot,findgen(10),findgen(10),xr=[8.0,9.0],yr=[1d-3,1d-1],/ylog,$
		/xs,/ys,/nodata
	
	oplot,plot_oh,scldgr
	oplot,plot_oh,scldgr_d2,linestyle=2
	oplot,plot_oh,scldgr_x2,linestyle=2

	; set up vector to hold indices of galaxies with M10 gradients
	; and good solutions from paper 1 and only their major axes if
	; inclination is high
	allgood = [0]

	; loop through galaxies
	j = 0
	loadct,4
	for i=0,ntarg-1 do BEGIN

		; get galaxy orientation parameters
		s = kingfish_galaxies(galname[i])

		restore,datadir+galname[i]+'_samp.sav'

		if gstruct.metal_source eq 'Not in M10 Table 8' then BEGIN
			print,'Skipping '+galname[i]
			goto,skip
		endif

		ok = where(allsol.gal eq galname[i],ct)
		if ct eq 0 then goto,skip
		deproject,allsol.ra_cen,allsol.dec_cen,gal=s,rgrid=rgrid,tgrid=tgrid,/vector

		; pick out the good solutions
		if s.incl_deg gt 60 and keyword_set(thetacut) then BEGIN
			good = where(allsol.gal eq galname[i] and allsol.aco_unc lt 0.2 and $
				abs(cos(tgrid)) gt 0.5,gct)
		endif else BEGIN
			good = where(allsol.gal eq galname[i] and allsol.aco_unc lt 0.2,gct)
		endelse

		oplot,allsol[ok].metal,allsol[ok].dgr,ps=3,color=clr[i] ;,/ylog,yr=[1d-4,1d0]
		if gct gt 0 then BEGIN
			oplot,allsol[good].metal,allsol[good].dgr,color=clr[i],ps=2
			allgood = [allgood,good]
		endif
		xyouts,[0.8],[0.8]-(j*0.05),galname[i],/normal,color=clr[i]
		j += 1

		skip:
	endfor

	; remove point 1 which was a placeholder
	allgood = allgood[1:n_elements(allgood)-1]

	; save the good solutions in a separate file
	allgoodsol = allsol[allgood]
	save,file='allgoodsol.sav',allgoodsol

	goodsol:
	; plot
	psopen,'fig1.eps',/encapsulated,/portrait,/color,$
		xsize=8,ysize=5,/inches

	plot,allgoodsol.metal,alog10(allgoodsol.dgr),ps=1,xr=[8.0,9.0],yr=[-3,-1],$
		tit='!6',xtit='12+log(O/H) (PT05)',$
		ytit='Log(DGR)',$
		xthick=5,ythick=5,thick=3,charthick=5,charsize=1.3

	blah = linfit(allgoodsol.metal,alog10(allgoodsol.dgr))
	x = findgen(20)/10 + 7.5
	oplot,x,blah[0]+blah[1]*x,linestyle=2,thick=5

	meanoff = (mean(alog10(allgoodsol.dgr)-allgoodsol.metal))

	oplot,x,x+meanoff,thick=5
	oplot,x,x+meanoff-alog10(2d),linestyle=1,thick=5
	oplot,x,x+meanoff+alog10(2d),linestyle=1,thick=5

	at85 = 8.5+meanoff

	print,'Rank Correlation: '+string((r_correlate(allgoodsol.metal,alog10(allgoodsol.dgr)))[0])
	print,'log10(DGR) at 12+log(O/H) = 8.5, for linear slope: '+string(at85)
	print,'Best fit slope: '+string(blah[1])
	print,'Best fit intercept: '+string(blah[0])

	restore,'data/prof_struct.idl'
	plotsym,0,/fill
	hidom = where(big_struct.hi gt big_struct.co21*6d)
	oplot,big_struct[hidom].metals,alog10(big_struct[hidom].local_dgr),ps=8,$
		thick=5,color=getcolor('forest')

	p1resid = alog10(allgoodsol.dgr)-(allgoodsol.metal+meanoff)
	locresid = alog10(big_struct[hidom].local_dgr) - (big_struct[hidom].metals+meanoff)

	plots,[0.17],[0.88],/normal,ps=1,thick=3
	plots,[0.17],[0.83],/normal,ps=8,color=getcolor('forest')
	xyouts,[0.185],[0.87],'Paper 1',charthick=5,/normal,charsize=1.3
	xyouts,[0.185],[0.82],'HI Dominated Regions',charthick=5,/normal,charsize=1.3

	plots,[0.17,0.22],[0.28,0.28],/normal,thick=5
	plots,[0.17,0.22],[0.23,0.23],/normal,thick=5,linestyle=1
	plots,[0.17,0.22],[0.18,0.18],/normal,thick=5,linestyle=2
	xyouts,[0.23],[0.27],'Linear DGR(Z)',charthick=5,/normal,charsize=1.3
	xyouts,[0.23],[0.22],'Linear DGR(Z) !MX!X/!M/!X 2',charthick=5,/normal,charsize=1.3
	xyouts,[0.23],[0.17],'Best Fit',charthick=5,/normal,charsize=1.3

	plothist,p1resid,bin=0.1,pos=[0.65,0.26,0.9,0.5],/noerase,charthick=5,$
		thick=5,xthick=5,ythick=5,xr=[-0.5,0.5],/xs,$
		xtit='DGR Residual (dex)',/ylog,yr=[1d0,1d2]
	plothist,locresid,bin=0.1,pos=[0.65,0.26,0.9,0.5],/overplot,thick=5,$
		color=getcolor('forest'),/ylog
	oplot,alog10([0.5,0.5]),[0.1,1d3],linestyle=1,thick=5	
	oplot,alog10([2.,2.]),[0.1,1d3],linestyle=1,thick=5

	psclose

stop
end
	
