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
                    hi_unc:nan, $
                    co21:nan, $
                    co21_unc:nan, $
                    co10:nan, $
                    co10_unc:nan, $
                    dust:nan, $
                    dust_unc:nan, $
                    metals:nan, $
                    metals_unc:nan, $
                    local_dgr:nan, $
                    local_dgr_mad:nan $
                    }

     oversamp = 6

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
     e1 = gstruct.co_hera_unc

     y2 = gstruct.hi*mh*1.36/ms*pc*pc
     e2 = gstruct.hi_unc*mh*1.36/ms*pc*pc

     y3 = gstruct.sigdust
     e3 = gstruct.sigdust_unc

     y4 = gstruct.metal
     e4 = gstruct.metal_unc

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
     e1 = e1[use_ind]

     y2 = y2[use_ind]
     e2 = e2[use_ind]

     y3 = y3[use_ind]
     e3 = e3[use_ind]

     y4 = y4[use_ind]
     e4 = e4[use_ind]

     dgr = y3/(y2+y1*6.3)
     
     bin_prof, x, y1 $
               , unc=e1 $
               , xmin=xmin $
               , xmax=xmax $
               , xmid=xmid $
               , countprof=co21_ctprof $
               , binsize=binsize $
               , meanprof=co21_mean $
               , madprof=co21_mad $
               , oversamp=oversamp $
               , errprof=co21_errprof

     bin_prof, x, y2 $
               , unc=e2 $
               , xmin=xmin $
               , xmax=xmax $
               , xmid=xmid $
               , countprof=hi_ctprof $
               , binsize=binsize $
               , meanprof=hi_mean $
               , madprof=hi_mad $
               , oversamp=oversamp $
               , errprof=hi_errprof

     bin_prof, x, y3 $
               , unc=e3 $
               , xmin=xmin $
               , xmax=xmax $
               , xmid=xmid $
               , countprof=dust_ctprof $
               , binsize=binsize $
               , meanprof=dust_mean $
               , madprof=dust_mad $
               , oversamp=oversamp $
               , errprof=dust_errprof

     bin_prof, x, y4 $
               , unc=e4 $
               , xmin=xmin $
               , xmax=xmax $
               , xmid=xmid $
               , countprof=metal_ctprof $
               , binsize=binsize $
               , meanprof=metal_mean $
               , madprof=metal_mad $
               , oversamp=oversamp $
               , errprof=metal_errprof

     bin_prof, x, dgr $
               , xmin=xmin $
               , xmax=xmax $
               , xmid=xmid $
               , countprof=dgr_ctprof $
               , binsize=binsize $
               , medprof=dgr_med $
               , madprof=dgr_mad $
               , oversamp=oversamp

     nprof = n_elements(xmid)
     for j = 0, nprof-1 do begin
        this_struct = (replicate(empty_struct,1))[0]
        this_struct.gal = g
        this_struct.rmin = xmid[j]-binsize*0.5
        this_struct.rmax = xmid[j]+binsize*0.5
        this_struct.rmid = xmid[j]

        this_struct.co21 = co21_mean[j]        
        this_struct.co21_unc = co21_errprof[j]

        this_struct.hi = hi_mean[j]
        this_struct.hi_unc = hi_errprof[j]

        this_struct.dust = dust_mean[j]
        this_struct.dust_unc = dust_errprof[j]

        this_struct.metals = metal_mean[j]
        this_struct.metals_unc = nan

        this_struct.local_dgr = dgr_med[j]
        this_struct.local_dgr_mad = dgr_mad[j]

        if n_elements(big_struct) eq 0 then begin
           big_struct = [this_struct]
        endif else begin
           big_struct = [big_struct, this_struct]
        endelse
     endfor

  endfor

  save, file="../data/prof_struct.idl" $
        , big_struct

  stop

end
