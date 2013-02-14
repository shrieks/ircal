; input lambda is in nm
; temperature in C
; delta is subtracted from the reflectivity of the coating before it
; is returned by the throughput function
; dustfrac is a fraction of the surface with emissivity=1, possibly
; from dust or other contamination
; set allemis to 1 to add background radiation of emissivity 1.0 from each
; surface without seting the throughput to zero. This is to account for, e.g.,
; the oversized cold stop in ircal
function background, lambdavec, namelist, templist, skyemis, deltalist, dustfrac=dustfrac, edust=edust, fileroot=fileroot, total_trans=total_trans, allemis=allemis

  if not(keyword_set(dustfrac)) then $
     dustfrac = 0.0D
  if not(keyword_set(edust)) then $
     edust = 1.0

  nlam = n_elements(lambdavec)
  nsurface = n_elements(namelist)
  total_emis = fltarr(nlam)
  total_trans = fltarr(nlam)
  throughput_eff=fltarr(nlam)

  ; convert to cm from nm
  lambda = lambdavec*100./1e9
  ; 2*h*c*c units W/(m^2 ster cm^(-4))
  c1 = 1.191044e-8 
  ; hc/k units (K cm)
  c2 = 1.438769
  ; hc units W cm s
  c3 = 1.98645e-23

  for i=0, nsurface-1 do begin
;     if (i eq 11) then print, foo
     thvec = throughput(lambdavec, namelist[i], fileroot=fileroot, delta=deltalist[i])
     epsilon = 1.0-thvec
     ; account for particle contaminatio 
     epsilon_eff = (1.-dustfrac)*epsilon + edust*dustfrac
     th_eff = 1.0-epsilon_eff
     if (keyword_set(allemis)) then $
        epsilon_eff = epsilon_eff*0.0+1.0
;     th_eff = thvec-dustfrac
;     epsilon_eff = 1.0-th_eff

     ; now calculate n photons
     tempK = templist[i]+273.
     intensity = (c1/(lambda^5))*(1./(exp(c2/(lambda*tempK)) - 1.0))
     ; this is in photons s^-1 m^-2 ster^-2 cm^-1 from a blackbody
     nphotons = intensity/(c3/lambda)            
     ; convert to photons s^-2 m^-2 arcsec^-2 nm^-1
     nphotons=nphotons*100./(206265.*206265*1e9)
     ; correct emitted photons from a BB to effective emissivity
     onesurf=nphotons*epsilon_eff
     
     if (i eq 0) then begin
        total_trans = th_eff
        total_emis = skyemis + onesurf
     endif else begin 
        total_emis = total_emis*th_eff+onesurf
        total_trans = total_trans*th_eff
     endelse
  endfor
  return, total_emis
end
