; this relies on throughput.pro, and like that is ONLY meant for 1 um and redder
; assumes everyting in <names> is at thesame temperature
; call this function multiple times and add the results for different temps
; input lambda is in nm
; temperature in C
; delta is subtracted from the reflectivity of the coating before it
; is returned by the throughput function
; dustfrac is a fraction of the surface with emissivity=1, possibly
; from dust or other contamination
function emissivity2, lambdavec, names, temp1=temp1, temp2=temp2, n_at_temp1=n_at_temp1, $
  delta=delta, dustfrac=dustfrac, dustemis = dustemis, $
  fileroot=fileroot, throughput_eff=throughput_eff

  if N_params() LT 1 THEN BEGIN
    print,'syntax - xxx = EMISSIVITY( lambdavec, names, temperature, delta=delta, dustfrac=dustfrac, dustemis=dustemis, fileroot=fileroot, throughput_eff=throughput_eff)'
    print,'  -this relies on throughput.pro, and is ONLY meant for 1 um and redder'
    print,'  -returns the total emission from all optics specified'
    print,'  -this assumes everyting in <names> is at the same temperature'
    print,'  -call this function multiple times and add the results for different temps'
    print,'  -input <lambdavec> is in nm'
    print,'  -input <names> is a list of string surface names'
    print,'  -input <temp1> is in degree C'
    print,'  -optional <temp2> is in degree C; see <n_at_temp1> for more details'
    print,'  -optional <n_at_temp1> is the number of optics at <temp1> '
    print,'    -if <n_at_temp1> unspecified and <temp2> specified, 1st 3 optics are assumed to be '
    print,'       at <temp1> (= Keck system);'
    print,'    -if both <n_at_temp1> and <temp2> unspecified, all optics are assumed to be at <temp1>'
    print,'  -optional <delta> is subtracted from the reflectivity of the coating before it is '
    print,'     returned by the throughput function; currently only applied to Aluminum'
    print,'  -optional <dustfrac> is the fraction of the surface covered with dust or other '
    print,'     contaminate, with a unique emissivity (specified by <dustemis>)'
    print,'  -optional <dustemis> is the emissivity of the dust/contaminate on the optics; equal '
    print,'     to 1.0 if unspecified'
    print,'  -optional <fileroot> gives location of data files (reflectivity curves, etc)'
    print,'  -optional <throughput_eff> returns the total effective throughput of system'
    stop
  ENDIF

  if not(keyword_set(dustfrac)) then $
     dustfrac = 0.0D
  IF not(keyword_set(dustemis)) then $
     dustemis = 1.0D
     
  nlam = n_elements(lambdavec)
  nsurface = n_elements(names)
  output = fltarr(nlam)
  throughput_eff=fltarr(nlam)
 
  ; IF n_at_temp is not set OR if it is set to an unreal value, set to all surfaces and 
  IF not(keyword_set(n_at_temp1)) THEN BEGIN
     IF keyword_set(temp2) THEN BEGIN
        n_at_temp1 = 3
        print,' number of surfaces at temp1 not specified - assuming Keck system ' +$
          'of 3 surfaces at temp1'
     ENDIF ELSE BEGIN
        n_at_temp1 = nsurface 
        print,' applying temperature 1 to all surfaces '
     ENDELSE
  ENDIF ELSE BEGIN
     IF n_at_temp1 GT nsurface OR n_at_temp1 LT 1 THEN BEGIN
        n_at_temp1 = nsurface
        print,' number of surfaces supplied to be at temp1 is not realistic '+$
               '- setting to all surfaces'
     ENDIF
  ENDELSE

  tempK1 = temp1+273.
  IF keyword_set(temp2) THEN tempK2 = temp2+273.

  ; convert to cm from nm
  lambda = lambdavec*100./1e9
  ; 2*h*c*c units W/(m^2 ster cm^(-4))
  c1 = 1.191044e-8 
  ; hc/k units (K cm)
  c2 = 1.438769
  ; hc units W cm s
  c3 = 1.98645e-23
  intensity1 = (c1/(lambda^5))*(1./(exp(c2/(lambda*tempK1)) - 1.0))
  ; this is in photons s^-1 m^-2 ster^-2 cm^-1
  nphotons1 = intensity1/(c3/lambda)            
  ; convert to photons s^-2 m^-2 arcsec^-2 nm^-1
  nphotons1=nphotons1*100./(206265.*206265.*1.e9)
  IF keyword_set(temp2) THEN BEGIN
     intensity2 = (c1/(lambda^5))*(1./(exp(c2/(lambda*tempK2)) - 1.0))
     nphotons2 = intensity2/(c3/lambda)
     nphotons2=nphotons2*100./(206265.*206265.*1.e9)
  ENDIF

  for i=0, nsurface-1 do begin

     IF (strlowcase(strtrim(names[i],2)) EQ 'aluminum') AND (keyword_set(delta)) THEN $
        delta_i = delta ELSE $
        delta_i = 0.d

     thvec = throughput(lambdavec, names[i], fileroot=fileroot, delta=delta_i)
     epsilon = 1.0-thvec
;     ifix = where(epsilon lt 0.2, nfix)
;     epsilon[ifix] = 0.2
; account for particle contamination, assume particles have emiss=dustemis ;emiss=0.5
     epsilon_eff = (1.-dustfrac)*epsilon + dustfrac*dustemis
     tau_eff = 1.-epsilon_eff ; new effective throughput/transmissivity

     IF i GE n_at_temp1 AND keyword_set(temp2) THEN $
        nphotons = nphotons2 ELSE $ ; surfaces at temp 2
        nphotons = nphotons1        ; surfaces at temp 1

     onesurf=nphotons*epsilon_eff ; apply the effective emissivity
     ; output of optic[i] = (output of previous optics) times (throughput of optic[i]) plus emissivity of optic[i]
     output = output*tau_eff + onesurf
     if (i eq 0) then begin
        throughput_eff = tau_eff
     endif else begin
        throughput_eff = throughput_eff * tau_eff
     endelse
  endfor
  return, output
end
     

  
