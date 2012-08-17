function calc_alpha, dgr_prime $
                     , alpha_0 = alpha_0 $
                     , av_0 = av_0 $
                     , ico = ico $
                     , w10_alpha = w10_alpha $
                     , f12_alpha = f12_alpha $
                     , n12_alpha = n12_alpha $
                     , cap_alpha = cap_alpha

; ALPHA_CO AT SOLAR
  if n_elements(alpha_0) eq 0 then $
     alpha_0 = 6.3

; DON'T LET ALPHA GET ABOVE THIS VALUE
  if n_elements(cap_alpha) eq 0 then $
     cap_alpha = 100.
          
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; WOLFIRE ET AL. (2010)
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
  
; SHIELDING LAYER
  g0 = 1
  nh = 1
  delta_av = (0.53 - 0.045 * alog(g0/nh) - 0.097 * alog(dgr_prime))
  delta_av_0 = (0.53 - 0.045 * alog(g0/nh) - 0.097 * alog(1.0))

; SHIELDING IN THIS CLOUD
  if n_elements(av_0) eq 0 then $
     av_0 = 5.26 * (7.5d21/2.0/1d22) ; WOLFIRE + KRUMHOLZ

; AV IN THIS CLOUD
  av = av_0*dgr_prime

  x = -4.0 * delta_av / av
  fco_to_fh2 = exp(x)

  x_0 = -4.0 * delta_av_0 / av_0
  fco_to_fh2_0 = exp(x_0)

  w10_alpha = alpha_0 * (1./fco_to_fh2) * fco_to_fh2_0

  w10_alpha = w10_alpha < cap_alpha

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; FELDMAN (2012)
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

; COULD ALSO SERVE PRETTY WELL AS THE WILSON ET AL. (1996) VALUE

  f12_alpha = alpha_0/dgr_prime^(0.7)

  f12_alpha = f12_alpha < cap_alpha

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; NARAYANAN (2012)
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

  if n_elements(ico) gt 0.0 then begin
     n12_alpha = ((10.7*ico^(-0.32)) < 6.3)/dgr_prime^(0.65)
  endif else begin
     n12_alpha = !values.f_nan
  endelse

  n12_alpha = n12_alpha < cap_alpha

; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
; CLEAN UP SPOTS WHERE WE DON'T HAVE A DUST-TO-GAS RATIO
; %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

  nodust = where(dgr_prime eq 0 or finite(dgr_prime) eq 0, nodust_ct)

  if nodust_ct gt 0 then begin
     w10_alpha[nodust] = !values.f_nan
     f12_alpha[nodust] = !values.f_nan
     n12_alpha[nodust] = !values.f_nan
  endif

  return, w10_alpha

end
