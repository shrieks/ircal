; add <scalefac>*(1-transmissionvec) times a Planck function made for <temp> 
; to existing <skyvec>. wavevec in nm. assumes transmissionvec
; and skyvec use the same wavelength vector
function makebbsky, transmissionvec, skyvec, wavevecin, temp, floor, nadd=nadd, smooth=smooth

  if (keyword_set(smooth)) then begin
     trans = convgauss(transmissionvec, smooth)
  endif else begin
     trans = transmissionvec
  endelse
  emis = 1.0-trans
  iboost = where(emis lt floor)
  emis[iboost] = floor

  ; convert to cm
  lambda = wavevecin*100/1e9 
  ; 2*h*c*c units W/(m^2 ster cm^(-4))
  c1 = 1.191044e-8 
  ; hc/k units (K cm)
  c2 = 1.438769
  ; hc units W cm s
  c3 = 1.98645e-23
  intensity = (c1/(lambda^5))*(1./(exp(c2/(lambda*temp)) - 1.0))
  ; this is in photons s^-1 m^-2 ster^-2 cm^-1
  nphotons = intensity/(c3/lambda)            
  ; convert to photons s^-2 m^-2 arcsec^-2 nm^-1
  nphotons=nphotons*100./(206265.*206265*1e9)   
  
  nadd = nphotons*emis
  out = skyvec+nadd
  return, out
end
