pro coatout, namelist

fileroot='./emissdatafiles'
getcoatingdata, namelist, warr, carr, sarr, fileroot=fileroot

nnames = n_elements(namelist)

for i = 1, nnames-1 do begin
   forprint, *warr[i], *carr[i], *sarr[i], format='(f9.3,1x,f11.8,1x,e14.4)', $
             textout='./emissdatafiles/'+namelist[i]+'.out', $
             comment = "# lambda(nm)  coeff spline
end


end
@getcoatingdata
