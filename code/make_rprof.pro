pro make_rprof, just=just

  @constants.bat

  readcol, "../targets.list", gals, format='A'
  ngals = n_elements(gals)
  working_dir = "../data/"
  working_ext = "_samp.sav"

  for i = 0, ngals-1 do begin
     g = gals[i]
     s = things_galaxies(g)
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
                    metals:nan, $
                    local_dgr:nan $
                    }

;    &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
;    MAKE PROFILE
;    &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%     

     x = gstruct.r25

     ra = gstruct.ra
     dec = gstruct.dec
     deproject, ra, dec, gal=s $
                , RGRID = rad, TGRID = theta $
                , /vector

     y1 = gstruct.co_hera
     y2 = gstruct.hi*mh*1.36/ms*pc*pc
     y3 = gstruct.sigdust
     y4 = gstruct.metal

     use_ind = where(finite(x) and $
                     finite(y1) and $
                     finite(y2) and $
                     finite(y3) and $
                     finite(y4) and $
                     (abs(cos(theta)) gt 0.5 or $
                      s.incl_deg lt 60.) $
                     , use_ct)

     if use_ct eq 0 then continue
     
     x = x[use_ind]
     y1 = y1[use_ind]
     y2 = y2[use_ind]
     y3 = y3[use_ind]
     y4 = y4[use_ind]

     dgr = y3/(y2+y1*6.3)     
     
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

     bin_prof, x, dgr $
               , xmin=xmin $
               , xmax=xmax $
               , xmid=xmid $
               , binsize=binsize $
               , medprof=dgr_med $
               , madprof=dgr_mad

     nprof = n_elements(xmid)
     for j = 0, nprof-1 do begin
        this_struct = (replicate(empty_struct,1))[0]
        this_struct.gal = g
        this_struct.rmin = xmid[j]-binsize*0.5
        this_struct.rmax = xmid[j]+binsize*0.5
        this_struct.rmid = xmid[j]
        this_struct.co = co_mean[j]
        this_struct.hi = hi_mean[j]
        this_struct.dust = dust_mean[j]
        this_struct.metals = metal_mean[j]
        this_struct.local_dgr = dgr_med[j]
        if n_elements(big_struct) eq 0 then begin
           big_struct = [this_struct]
        endif else begin
           big_struct = [big_struct, this_struct]
        endelse
     endfor

  endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CALCULATE XCO
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

;  dgr = 1e-2*10.^(big_struct.metals-8.5)
  solar_dgr = 1.7e-2
  dgr_pred = solar_dgr*10.^(big_struct.metals-8.5)
  h2 = big_struct.dust/dgr_pred - big_struct.hi
  alpha = h2 / big_struct.co

  fid_x = findgen(101)/100. + 8.0
  fid_dgr = solar_dgr*10.^(fid_x-8.5)
  fid_y = calc_alpha(fid_dgr/solar_dgr, alpha_0=4.4 $
                     , av_0 = 2.3)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DUST-TO-GAS RATIO CHECK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  dgr_fid = big_struct.dust/ $
            (big_struct.co*6.3 + big_struct.hi)

  dgr_hi = big_struct.dust / big_struct.hi
  
  plot, big_struct.metals, dgr_fid $
        , ps=1, /ylo, xrange=[8,9], yrange=[1e-4, 1e0]
  circle, /fill
  for k = 0, ngals-1 do begin
     ind = where(big_struct.gal eq gals[k] and $
                 big_struct.co*6.3*3. lt big_struct.hi, ct)
     if ct eq 0 then continue
     plot, [big_struct[ind].metals], [dgr_fid[ind]], ps=8 $
           , xrange=[8,9], yrange=[1e-3, 1e0], /ylo
     print, gals[k]
     xfid = findgen(101)/100.+8.0
     yfid = 1.7e-2*10.^(xfid-8.5)
     oplot, [xfid], [yfid], thick=3, lines=1
                                ;ch = get_kbrd(1) & $     
  endfor

  save, file="../data/prof_struct.idl" $
        , big_struct
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PLOT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  psfile = "../plots/alpha_co.eps"
  ps, file=psfile, /def, xsize=12, ysize=8, /color, /encapsulated, /ps
  circle, /fill
  ind = where(big_struct.gal ne "ngc2841")
  plot, big_struct[ind].metals, alpha[ind] $
        , psym=8 $
        , ytitle="!7a!6!dCO!n [M!d!9n!6!n pc!u-2!n (K km s!u-1!n)!u-1!n]" $
        , xtitle="!6Metallicity [12 + log O/H]" $
        , charthick=3, charsize=2, /ylo, yrange=[1e-1, 1e3]
  oplot, [8.0, 9.0], [4.4, 4.4], lines=2, thick=3  
  oplot, fid_x, fid_y, thick=4, lines=1
  ps, /x
  spawn, 'gv '+psfile+' &'




  stop

end
