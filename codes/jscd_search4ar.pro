;data path
dpath = '../data/'

;find data
files = file_search(dpath+'*.fits')

ff=0

;reading data and pre analysis
read_sdo,files[ff],hdr,data
aia_prep,hdr,data,hdr2,data2,/use_hdr,scale_ref=.504
index2map,hdr2,data2,map



;CRPIX1 CRPIX2 Sun center on CCD
plot,(map.data)[2444,*]

;Falta quitar el limb darkening

;WATERSHED(map.data)     



END