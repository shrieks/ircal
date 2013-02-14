; read in coating data. namelist is the list of coatings, wavearr and
; coatarr are arrays fo pointers to vectors that contain the
; wavelength vectors in nm and the coating data, in the order
; given in <namelist>. The splines are returned in splinearr 

pro getcoatingdata, namelist, wavearr, coatarr, splinearr, fileroot=fileroot, irdust=irdust

  if not(keyword_set(fileroot)) then $
     fileroot = '/b/ircal/throughput'
  
  if not(keyword_set(irdust)) then $
     irdust = 0.8 ; note this is transmission, not as in spreadsheet

; namelist = ['MEMS_window', 'PICNIC_QE', 'IR_transmission', 'Kprime', 'Hband', 'Jband', 'ALPAO', 'Aluminum', 'AluminumPlusDust', 'ARIR', 'CaF', 'FSG98', 'H2RG_QE', 'LGSWFS_Dich', 'NaHG', 'NaR_Splitter', 'TT_Dichroic_Refl','X1_Silver_Extrap', 'deimos_qe_data', 'DEIMOS_Camera', 'Grating_900_550']

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; no, make this more general for other throughput and put these in files like
; any other coating data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; other stuff
;  oldao_1turn = 0.97
;  oldao_tiptilt = 0.96
;  oldao_oap1 = 0.96
;  oldao_dm = 0.93
;  oldao_oap2 = 0.94
;  othernamelist = ['oldao_1turn', 'oldao_tiptilt', 'oldao_oap1', 'oldao_dm', 'oldao_oap2']
;  otherlist = [oldao_1turn, oldao_tiptilt, oldao_oap1, oldao_dm, oldao_oap2]

;  alllist = [namelist, othernamelist]
  nlist = n_elements(namelist)
;  nall = n_elements(alllist)
;  nother = n_elements(othernamelist)

  wavearr = ptrarr(nlist)
  coatarr = ptrarr(nlist)
  splinearr = ptrarr(nlist)

  for i=0,nlist-1 do begin
     ; 
     if (stregex(namelist[i], 'PlusDust') gt -1) then begin 
        foolist = stregex(namelist[i], '(^[a-zA-Z0-9]*)(PlusDust)', /subexpr, /extract)
        filename = fileroot + '/' + foolist[1] + '.txt'
        scalefac = irdust
     endif else begin 
        ; just the file, no modifiers
        filename = fileroot + '/' + namelist[i] + '.txt'
        scalefac = 1.0D
     endelse

     cmdstr = 'ls ' + filename
     spawn, cmdstr, flist
     if (strlen(flist) lt strlen(namelist[i])) then begin
        print,'file for ' + namelist[i] + ' not found'
     endif
     readcol, filename, /silent ,format='f,f', wavevec, coatvec
     ; multiply by any scale factor
     coatvec = coatvec*scalefac
     
     ; H2RG_QE is in nm, PICNIC_QE is also in nm
     ; al others in microns, need to convert
     if (not(stregex(namelist[i], 'H2RG_QE') eq 0 OR stregex(namelist[i], 'PICNIC_QE') eq 0)) then $
        wavevec = wavevec*1000.

     ; everything given in percent except IR_transmission, want fractions always
     if (stregex(namelist[i], 'IR_transmission') lt 0) then $
        coatvec = coatvec/100.

     spline = spl_init(wavevec, coatvec, /double)
     wavearr[i] = ptr_new(wavevec)
     coatarr[i] = ptr_new(coatvec)
     splinearr[i] = ptr_new(spline)
  endfor

  ; NOTE this is currently the MEMS_window wavelength vec, make sure 
  ; if you change the order in namelist
;  wvec = *wavearr[0]
  ; now add these one at a time
;  for i=0, nother-1 do begin
;     cvec = wvec*0.0+otherlist[i]
;     wavearr[nlist+i] = ptr_new(wvec)
;     coatarr[nlist+i] = ptr_new(cvec)
;     ospline = spl_init(wvec, cvec, /double)
;     splinearr[nlist+i] = ptr_new(ospline)
;  endfor
end    
