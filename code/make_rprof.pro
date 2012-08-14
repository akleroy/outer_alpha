pro make_rprof, just=just

  readcol, "targets.list", gals, format='A'
  ngals = n_elements(gals)
  working_dir = "data/"
  working_ext = "_samp.sav"

  for i = 0, ngals-1 do begin
     g = gals[i]
     if n_elements(just) gt 0 then $
        if (just ne g) then continue
     restore, working_dir+g+working_ext, /v

;    &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
;    SPECIFY PROFILE
;    &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     xmin=0.0
     xmax=1.0
     binsize=0.05

;    &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
;    MAKE PROFILE
;    &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     x = gstruct.r25
     y1 = gstruct.co_hera
     y2 = gstruct.hi
     y3 = gstruct.sigdust

     use_ind = where(finite(x) and $
                     finite(y1) and $
                     finite(y2) and $
                     finite(y3))
     
     x = x[use_ind]
     y1 = y1[use_ind]
     y2 = y2[use_ind]
     y3 = y3[use_ind]

     bin_prof, x, y1 $
               , xmin=xmin $
               , xmax=xmax $
               , xmid=xmid $
               , binsize=binsize $
               , meanprof=co_mean $
               , madlogprof=co_madlog

     bin_prof, x, y2 $
               , xmin=xmin $
               , xmax=xmax $
               , xmid=xmid $
               , binsize=binsize $
               , meanprof=hi_mean $
               , madlogprof=hi_madlog

     bin_prof, x, y3 $
               , xmin=xmin $
               , xmax=xmax $
               , xmid=xmid $
               , binsize=binsize $
               , meanprof=dust_mean $
               , madlogprof=dust_madlog

;    &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
;    PLOT
;    &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     circle
     plot, xmid, co_mean/median(co_mean), /ylo, ps=-8 $
           , yrange=[1e-2, 1e2], title=g
     oplot, xmid, hi_mean/median(hi_mean), ps=-8, color=getcolor('red')
     oplot, xmid, dust_mean/median(dust_mean), ps=-8, color=getcolor('blue')

     ch = get_kbrd(1)
;     stop
  endfor

end
