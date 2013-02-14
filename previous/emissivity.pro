; this relies on throughput.pro, and like that is ONLY meant for 1 um and redder
; assumes everyting in <names> is at thesame temperature
; call this function multiple times and add the results for different temps
; input lambda is in nm
; temperature in C
; delta is subtracted from the reflectivity of the coating before it
; is returned by the throughput function
; dustfac is a fraction of the surface with emissivity=1, possibly
; from dust or other contamination
function emissivity, lambdavec, names, temperature, delta=delta, dustfac=dustfac, fileroot=fileroot, throughput_eff=throughput_eff

  if not(keyword_set(dustfac)) then $
     dustfac = 0.0D

  nlam = n_elements(lambdavec)
  nsurface = n_elements(names)
  output = fltarr(nlam)
  throughput_eff=fltarr(nlam)
  
  tempK = temperature+273.

  ; convert to cm from nm
  lambda = lambdavec*100./1e9
  ; 2*h*c*c units W/(m^2 ster cm^(-4))
  c1 = 1.191044e-8 
  ; hc/k units (K cm)
  c2 = 1.438769
  ; hc units W cm s
  c3 = 1.98645e-23
  intensity = (c1/(lambda^5))*(1./(exp(c2/(lambda*tempK)) - 1.0))
  ; this is in photons s^-1 m^-2 ster^-2 cm^-1
  nphotons = intensity/(c3/lambda)            
  ; convert to photons s^-2 m^-2 arcsec^-2 nm^-1
  nphotons=nphotons*100./(206265.*206265*1e9)

  for i=0, nsurface-1 do begin
     thvec = throughput(lambdavec, names[i], fileroot=fileroot, delta=delta)
     epsilon = 1.0-thvec
;     ifix = where(epsilon lt 0.2, nfix)
;     epsilon[ifix] = 0.2
; account for particle contamination, assume particles have emiss=0.5
     epsilon_eff = (1.-dustfac)*epsilon + 1.0*dustfac
     onesurf=nphotons*epsilon_eff
     output = output+onesurf
     if (i eq 0) then begin
        throughput_eff = 1.-epsilon_eff
     endif else begin
        throughput_eff = throughput_eff*(1.0-epsilon_eff)
     endelse
  endfor
  return, output
end
     

  
