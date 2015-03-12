;Procedure to download dayly data from HMI instrument

for i=0, 28 do begin

if i eq 0 then begin

T0= '2011/02/'+ strtrim(i,1) +'T0:0:0-2011/02/'+ strtrim(i,1) +'T0:1:0'

endif else begin

T0= [T0, '2011/02/'+ strtrim(i,1) +'T0:0:0-2011/02/'+ strtrim(i,1) +'T0:1:0']

endelse

endfor

for i=0,28 do begin 
bar,i,28
print,t0[i], t1[i]
downr = vso_search(date=T0[i], physobs='intensity', inst='hmi')
status=vso_get(downr,out_dir='../data/')
endfor
end