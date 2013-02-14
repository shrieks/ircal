; uses getcoatingdata, and is ONLY meant to work at wavelengths beyond 1 um
; names is a list of coatings to be multipled to give a throughput 
; lambda in nm
; delta is subtracted from the reflectivity 
function throughput, lambda, namelist, delta=delta, fileroot=fileroot
  
  if not(keyword_set(delta)) then $
     delta = 0.0D

  getcoatingdata, fileroot=fileroot, namelist, wavearr, coatarr, splinearr
  
  nvals = n_elements(lambda)
  outvec = dblarr(nvals)+1.0

  nnames = n_elements(namelist)
  
  for j = 0, nnames-1 do begin
     name = namelist[j]
 ;    stridx = stregex(name, allnamelist)
     ic = where(name eq namelist, nmatch)
;     print, j, ic[0]
     wvec = *wavearr[ic[0]]
     cvec = *coatarr[ic[0]]
     spline = *splinearr[ic[0]]
     if (max(spline) gt 0) then begin
        cval = spl_interp(*wavearr[ic[0]], *coatarr[ic[0]], *splinearr[ic[0]], lambda)
     endif else begin
        ; if the spline is null the input vector was a constant
        print, "Null spline,", namelist[j]
        cval = outvec*0.0+min(cvec)
     endelse
     print, name, max(cval), min(wvec), max(wvec)
     cval = cval-delta
     outvec = outvec*cval
  endfor
  return, outvec
end
