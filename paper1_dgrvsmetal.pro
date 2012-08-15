pro paper1_dgrvsmetal,$
	saved=saved,$
	targetlist=targetlist,$
	just=just,$
	mwmetal=mwmetal,$
	mwdgr=mwdgr,$
	goodsol=goodsol

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
	plot_oh = findgen(20)/10 + 7.5
	oh = 10^(plot_oh-12d)
	mw_oh = 10.^(mwmetal-12d)
	scldgr = mwdgr*oh/mw_oh
	scldgr_d2 = 0.5*mwdgr*oh/mw_oh
	scldgr_x2 = 2.*mwdgr*oh/mw_oh

	; where to find sampled galaxies
	datadir = 'data/'

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
	; and good solutions from paper 1
	allgood = [0]

	; loop through galaxies
	j = 0
	loadct,4
	for i=0,ntarg-1 do BEGIN

		restore,datadir+galname[i]+'_samp.sav'
		if gstruct.metal_source eq 'Not in M10 Table 8' then BEGIN
			print,'Skipping '+galname[i]
			goto,skip
		endif

		ok = where(allsol.gal eq galname[i],ct)
		if ct eq 0 then goto,skip

		; pick out the good solutions
		good = where(allsol.gal eq galname[i] and allsol.aco_unc lt 0.2,gct)

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
	plot,allgoodsol.metal,alog10(allgoodsol.dgr),ps=1,xr=[8.0,9.0],yr=[-3,-1]
	blah = linfit(allgoodsol.metal,alog10(allgoodsol.dgr))
	x = findgen(20)/10 + 7.5
	oplot,x,blah[0]+blah[1]*x,linestyle=2

	oplot,x,x+(mean(alog10(allgoodsol.dgr)-allgoodsol.metal))
	oplot,x,x+(mean(alog10(allgoodsol.dgr)-allgoodsol.metal))-alog10(2d),linestyle=1
	oplot,x,x+(mean(alog10(allgoodsol.dgr)-allgoodsol.metal))+alog10(2d),linestyle=1

	at85 = 8.5+(mean(alog10(allgoodsol.dgr)-allgoodsol.metal))

	print,'Rank Correlation: '+string((r_correlate(allgoodsol.metal,alog10(allgoodsol.dgr)))[0])
	print,'log10(DGR) at 12+log(O/H) = 8.5, for linear slope: '+string(at85)
	print,'Best fit slope: '+string(blah[1])
	print,'Best fit intercept: '+string(blah[0])

stop
end
	
