;data path
dpath = '../data/'

;find data
files = file_search(dpath+'*.fits')

ff=0

;reading data and pre analysis
read_sdo,files[ff],hdr,data
aia_prep,hdr,data,hdr2,data2,/use_hdr,scale_ref=.504

;image size
sd = size(data)

;center and sun radiii
xcen = hdr2.CRPIX1
ycen = hdr2.CRPIX2
rsun = hdr2.RSUN_OBS
scale= hdr2.CDELT1

;correct limb darkening
darklimb_correct,data2,data,lambda=6173,$
limbxyr=[xcen, ycen, rsun/scale] 

;make a circle to eliminate limb pixels
dist_circle,dist_grid, sd[1], xcen,ycen                               
dist_grid=dist_grid/1850.         
outside=where(dist_grid gt 1., complement=inside)                            
dist_grid(outside)=0.
dist_grid(inside)=1.  
data = dist_grid*data


;mark outside pixels
mk = where (data lt 10000.,complement= inmk)
data(mk) = -0.

;Converting to map
index2map, hdr2, data, map

;mean ans dev standard
S_mean = mean(data(inmk))
s_sdev = stddev(data(inmk))

;reversal the image
ndata = 1.- data +2.*s_mean
ndata(mk) = -0

;mask image
mask = intarr(sd[1:2])
armk = where(ndata gt s_mean+ 3.*s_sdev)
mask(armk) = 1

;Label each region
data_lab = label_region(mask, /all) 

goto,jump
;eliminate small regions with less than lmpx number of pixels
lmpx=100
for i=1,max(data_lab) do begin
  tmp = where(data_lab eq i,ntmp)

if (ntmp gt 0) then begin
  if (ntmp lt lmpx) then mask(tmp) = 0
endif
endfor ;for i

;label final regions
mask2 = label_region(mask, /all) 

jump:

;find center of mass for each region
pos_arcsec=fltarr(2, max(mask2)) 
x_ar = fltarr(2, max(mask2)) 
y_ar = fltarr(2, max(mask2)) 

for i=0,n_elements(pos_arcsec[0,*])-1 do begin

tmp = where_xyz(mask2 eq i+1, xind=xind, yind=yind) 

;box AR
x_ar[*,i] = ([min(xind), max(xind)] - xcen) * scale
y_ar[*,i] = ([min(yind), max(yind)] - ycen) * scale

rmass_x = 0 & rmass_y = 0
for j=0, n_elements(tmp)-1 do begin

rmass_x += ndata(tmp[j])*(xind[j] - xcen)  
rmass_y += ndata(tmp[j])*(yind[j] - ycen)  
 
endfor

masstot = total(ndata(tmp)) 
pos_arcsec[*,i] = [(rmass_x/masstot)*scale, (rmass_y/masstot)*scale]

endfor;for i


;arcsec to angle
pos_angle = arcmin2hel(pos_arcsec[0,*]/60., pos_arcsec[1,*]/60., radio=rsun)


;-----plotting subregions
file_mkdir,'../results/'
tmp = STRJOIN(STRSPLIT(map.time, '.', /EXTRACT), '_')
tmp = STRJOIN(STRSPLIT(tmp, ' ', /EXTRACT), '_')
name = STRJOIN(STRSPLIT(tmp, ':', /EXTRACT), '_')
fname='../results/'+ name +'.eps' 
PS_start,filename=fname,/nomatch,/QUIET 
device,xsize=16, ysize=20, xoffset=0., yoffset=0. ;SIZE OUTPUT FILE
nc=100                   
loadct,3, ncolors=nc 
!p.font=7

!x.ticklen = -.02
!y.ticklen = -.02
!p.charsize=0.7
plot_map, map, ncolors=nc, grid= 20, title = map.time
!p.charsize=0.9
for i=0, n_elements(pos_arcsec[0,*])- 1 do begin

pl=20
sub_map,map,smap, xrange=[x_ar[0,i]-pl, x_ar[1,i]+ pl], yrange= [y_ar[0,i]- pl, y_ar[1,i]+ pl], /noplot

tit='R'+ strtrim(i,1) + ' (' +  strtrim(pos_arcsec[0,i] ,1) +', '+ strtrim(pos_arcsec[1,i] ,1) +')'
plot_map,smap,ncolors=nc, title = tit, /grid, multi=[2,3]


endfor ;for i
!p.multi=0
;CLOSE PS	
PS_end,/nomessage 
cgps2pdf,fname,/delete_ps,unix_convert_cmd='epstopdf',/silent ;CONVERT TO PDF



;write table
head = ['Time', 'Index AR', 'Xpos arcsec', 'Ypos arcsec', 'Latitud','Longitud']
ntable = '../results/'+ name+'.csv'
tim = replicate(map.time, n_elements( pos_arcsec[0,*]))
index = indgen(n_elements( pos_arcsec[0,*]))
write_csv, ntable, tim, index, pos_arcsec[0,*], pos_arcsec[1,*], pos_angle[0,*], pos_angle[1,*], header=head



END