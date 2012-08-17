pro bin_prof, x, y $
              , unc=unc $
              , xmin = xmin_in $
              , xmax = xmax_in $
              , binsize = binsize_in $
              , irregular = irregular $
              , oversamp=oversamp_in $
              , medprof=medprof $
              , meanprof=meanprof $
              , madprof=madprof $
              , madlogprof=madlogprof $
              , stdprof=stdprof $
              , stdmeanprof=stdmeanprof $
              , errprof=errprof $
              , maxprof=maxprof $
              , minprof=minprof $
              , countprof=countprof $
              , outfile=outfile $
              , xmid_bin=xmid_bin $
              , xgeo_bin=xgeo_bin $
              , xmeanprof = xmeanprof $
              , percentile=percentile $
              , lopercprof=lopercprof $
              , hipercprof=hipercprof $
              , writegeo=writegeo $
              , writexmean=writexmean

;+
; NAME:
;
;  bin_prof
;
; PURPOSE:
;
;  Eventually, a generic routine for extracting binned profiles. Pretty damn
;  close right now.
;
; CATEGORY:
;
;  Science tool. Documentation to follow.
;
; CALLING SEQUENCE:
;
;  bin_prof, x, y $
;              , xmin = xmin_in, xmax = xmax_in $
;              , binsize=binsize_in $
;              , medprof=medprof, meanprof=meanprof $
;              , madprof=madprof, madlogprof=madlogprof $
;              , stdprof=stdprof, errprof=errprof $
;              , maxprof=maxprof, minprof=minprof $
;              , countprof=countprof, oversamp=overamp $
;              , percentile=percentile, lopercprof=lopercprof $
;              , hipercprof=hipercprof
; INPUTS:
;
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OPTIONAL OUTPUTS:
;
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;
; Added MIN/MAX measurement
;
;-

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ERROR CHECKING
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; IF THERE IS NO DATA OR NOT-MATCHING DATA, SCREAM AND QUIT
  if (n_elements(y) ne n_elements(x)) or (n_elements(x) eq 0) then begin
      message, "Mismatched or missing data. Returning.", /informational
      return
  endif


; OVERSAMPLING FACTOR - DEFAULTS TO ONE
  if n_elements(oversamp_in) eq 0 then begin
      message, 'No oversampling supplied. Assuming none.', /informational
      oversamp = 1.
  endif else begin
      oversamp = oversamp_in
  endelse

; JUST A CHECK
  if n_elements(percentile) gt 0 then begin
      message, 'Returning upper/lower brackets for range '+ $
        string(percentile), /informational
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CONSTRUCT THE BINS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(xmax_in) eq 0 then begin
      if keyword_set(irregular) then begin
          message, 'Bin maxima required with /IRREGULAR flag. Returning.', /info
          return
      endif
      message, 'No maximum specified. Defaulting to data.', /informational
      xmax = max(x, /nan)
  endif else begin
      xmax = xmax_in
  endelse

  if n_elements(xmin_in) eq 0 then begin
      if keyword_set(irregular) then begin
          message, 'Bin minima required with /IRREGULAR flag. Returning.', /info
          return
      endif
      message, 'No minimum specified. Defaulting to data.', /informational
      xmin = min(x, /nan)
  endif else begin
      xmin = xmin_in
  endelse

  if keyword_set(irregular) eq 0 then begin
      deltax = (xmax - xmin)
      
      if n_elements(binsize_in) eq 0 then begin
          message, "No binsize specified. Making 10 bins.", /info
          binsize = deltax / 10.
      endif else begin
          binsize = binsize_in
      endelse

;     MAKE THE BINS
      nbins = ceil(deltax / binsize) 
      xmin_bin = findgen(nbins)*binsize+xmin 
      xmax_bin = xmin_bin+binsize < xmax 
      xmid_bin = (xmin_bin+xmax_bin)*0.5      

  endif else begin
      
      nbins = n_elements(xmax_in)
      xmin_bin = xmin
      xmax_bin = xmax
      xmid_bin = (xmin_bin+xmax_bin)*0.5
      xgeo_bin = sqrt(xmin_bin*xmax_bin)

  endelse

; TO DO: CONSISTENCY CHECK XMAX, XMIN, BINSIZE, NBINS

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; INITIALIZE PROFILES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; INITIALIZE THE OUTPUT AND FILL WITH NANS
  xmeanprof = dblarr(nbins)*!values.f_nan
  medprof = dblarr(nbins)*!values.f_nan
  meanprof= dblarr(nbins)*!values.f_nan
  madprof = dblarr(nbins)*!values.f_nan
  madlogprof = dblarr(nbins)*!values.f_nan
  stdprof = dblarr(nbins)*!values.f_nan
  stdmeanprof = dblarr(nbins)*!values.f_nan
  errprof = dblarr(nbins)*!values.f_nan
  maxprof = dblarr(nbins)*!values.f_nan
  minprof = dblarr(nbins)*!values.f_nan
  countprof = lonarr(nbins)*!values.f_nan
  if n_elements(percentile) gt 0 then begin
      hipercprof = lonarr(nbins)*!values.f_nan
      lopercprof = lonarr(nbins)*!values.f_nan
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LET'S MAKE US SOME DAMN PROFILES!
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; LOOP AND BIN
  for ii = 0, nbins-1 do begin

;     FIND THE PIXELS WITHIN THE PRESENT RING
      binind = where(x gt xmin_bin[ii] AND x lt xmax_bin[ii], binct)
      
;     NOTE THE COUNTS IN THE COUNTS PROFILE
      countprof[ii] = binct

;     IF WE HAVE SOME PIXELS, BIN... AND PROFIT!
      if (binct gt 1) then begin   
          indep = x[binind]
          data = y[binind]         

;         ... MEAN AND MEDIAN VALUES
          xmeanprof[ii] = mean(indep,/nan)
          medprof[ii]  = median(data,/even)
          meanprof[ii] = mean(data, /nan)
          
;        ... MAX AND MIN VALUES
          minprof[ii] = min(data, /nan)
          maxprof[ii] = max(data, /nan)         
          
;         ... UNCERTAINTY ESTIMATES
          if (binct gt 5) then begin
              madprof[ii] = mad(data)
              madlogprof[ii] = mad(alog10(data))
              stdprof[ii] = stddev(data)
              stdmeanprof[ii] = stdprof[ii] / sqrt(countprof[ii]/oversamp)
              if n_elements(unc) gt 0 then begin
                 if n_elements(unc) eq 1 then begin
                    errprof[ii] = unc / sqrt(countprof[ii]/oversamp)
                 endif else begin
                    errprof[ii] = $
                       sqrt(total(unc[binind]^2) / oversamp) $
                       / countprof[ii]
                 endelse
              endif else begin
                 errprof[ii] = !values.f_nan
              endelse
          endif

;        ... IF REQUESTED, THEN SORT AND RETURN PERCENTILE VALUES
;        (THIS IS COMPUTATIONALLY A BIT INTENSIVE, SO NORMALLY WE AVOID IT)
          if n_elements(percentile) gt 0 then begin             
              data = data[sort(data)]
              loind = ceil(binct*(0.5-percentile/2.))-1
              hiind = ceil(binct*(0.5+percentile/2.))-1
              lopercprof[ii] = data[loind]
              hipercprof[ii] = data[hiind]
          endif
      endif
  endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE ASCII OUTPUT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%     

  xstattag = 'xmid'
  if keyword_set(writegeo) then xstattag = 'xgeo'
  if keyword_set(writexmean) then xstattag = 'xmean'  

  if n_elements(outfile) ne 0 then begin
      openw,1,outfile
      if n_elements(percentile) eq 0 then begin
          printf,1, $
            'xmin|'+xstattag+'|xmax|'+$
            'counts|mean|median|'+$
            'stddev|mad|err|madlog|'+$
            'min|max'
      endif else begin
          printf,1, $
            'xmin|'+xstattag+'|xmax|'+$
            'counts|mean|median|'+$
            'stddev|mad|err|madlog|'+$
            'perc|loperc|hiperc|'+$
            'min|max'
      endelse

      for kk = 0, nbins - 1 do begin
;         MEASUREMENTS OF THE INDEPENDENT VARIABLE
          rstring = string(xmin_bin[kk])+'|'
          if keyword_set(writegeo) then $
            rstring += string(xgeo_bin[kk]) $
          else if keyword_set(writexmean) then $
            rstring += string(xmeanprof[kk]) $
          else $
            rstring += string(xmid_bin[kk])
          rstring += '|'+string(xmax_bin[kk])+'|'

;         PROFILE MEASUREMENTS
          pstring = $
            string(countprof[kk])+'|'+string(meanprof[kk])+'|'+$
            string(medprof[kk])+'|'
;         ERROR MEASUREMENTS
          estring = $
            string(stdprof[kk])+'|'+ string(madprof[kk])+'|'+$
            string(errprof[kk])+'|'+string(madlogprof[kk])+'|'
;         MIN/MAX MEASUREMENTS
          mmstring = $
            string(minprof[kk])+'|'+string(maxprof[kk])
;         PERCENTILE MEASUREMENTS (IF REQUESTED)
          if n_elements(percentile) eq 0 then begin
              printf,1,rstring+pstring+estring+mmstring
          endif else begin
              percstring = $
                '|'+string(percentile)+'|'+ string(lopercprof[kk])+'|'+$
                string(hipercprof[kk])+'|'
              printf,1,rstring+pstring+estring+percstring+mmstring
          endelse
      endfor
      close,1
  endif
  
end                             ; OF RPROF
