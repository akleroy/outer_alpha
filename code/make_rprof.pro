pro make_rprof, just=just

  @constants.bat

  readcol, "../targets.list", gals, format='A'
  ngals = n_elements(gals)
  working_dir = "../data/"
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
     binsize=0.10

     nan = !values.f_nan
     empty_struct = { $
                    gal:"", $
                    rmid:nan, $
                    rmin:nan, $
                    rmax:nan, $
                    hi:nan, $
                    co:nan, $
                    dust:nan, $
                    metals:nan $
                    }

;    &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
;    MAKE PROFILE
;    &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     x = gstruct.r25
     y1 = gstruct.co_hera
     y2 = gstruct.hi
     y3 = gstruct.sigdust
     y4 = gstruct.metal

     use_ind = where(finite(x) and $
                     finite(y1) and $
                     finite(y2) and $
                     finite(y3) and $
                     finite(y4), use_ct)

     if use_ct eq 0 then continue
     
     x = x[use_ind]
     y1 = y1[use_ind]
     y2 = y2[use_ind]
     y3 = y3[use_ind]
     y4 = y4[use_ind]

     bin_prof, x, y1 $
               , xmin=xmin $
               , xmax=xmax $
               , xmid=xmid $
               , countprof=ctprof $
               , binsize=binsize $
               , meanprof=co_mean $
               , madprof=co_mad

     bin_prof, x, y2 $
               , xmin=xmin $
               , xmax=xmax $
               , xmid=xmid $
               , binsize=binsize $
               , meanprof=hi_mean $
               , madprof=hi_mad

     bin_prof, x, y3 $
               , xmin=xmin $
               , xmax=xmax $
               , xmid=xmid $
               , binsize=binsize $
               , meanprof=dust_mean $
               , madprof=dust_mad

     bin_prof, x, y4 $
               , xmin=xmin $
               , xmax=xmax $
               , xmid=xmid $
               , binsize=binsize $
               , meanprof=metal_mean $
               , madprof=metal_mad

     nprof = n_elements(xmid)
     for j = 0, nprof-1 do begin
        this_struct = (replicate(empty_struct,1))[0]
        this_struct.gal = g
        this_struct.rmin = xmid[j]-binsize*0.5
        this_struct.rmax = xmid[j]+binsize*0.5
        this_struct.rmid = xmid[j]
        this_struct.co = co_mean[j]
        this_struct.hi = hi_mean[j]*mh*1.36/ms*pc*pc
        this_struct.dust = dust_mean[j]
        this_struct.metals = metal_mean[j]
        if n_elements(big_struct) eq 0 then begin
           big_struct = [this_struct]
        endif else begin
           big_struct = [big_struct, this_struct]
        endelse
     endfor

;    &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
;    CALCULATE XCO
;    &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

;  dgr = 1e-2*10.^(big_struct.metals-8.5)
  dgr = 1.5e-2*10.^(big_struct.metals-8.5)
  h2 = big_struct.dust/dgr - big_struct.hi
  alpha = h2 / big_struct.co
  plot, big_struct.metals, alpha, ps=8, /ylo, yrange=[1,1e3]

;    &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
;    PLOT
;    &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  psfile = "alpha_co
  ps, 
  circle, /fill
  plot, xmid, (co_mean/median(co_mean)) $
        , /ylo, ps=-8, yrange=[1e-2, 1e2], title=g
  oplot, xmid, hi_mean/median(hi_mean) $
         , ps=-8, color=getcolor('red')
     oplot, xmid, dust_mean/median(dust_mean) $
            , ps=-8, color=getcolor('blue')
     
     ;ch = get_kbrd(1)
  endfor




  stop

end
