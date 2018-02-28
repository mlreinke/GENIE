pro read_fiesta_eqdsk,fname,efit,skip_end=skip_end

openr,lun,fname,/get_lun
;str11 = '(a50,3i4)'  ; modifed to cope with 129 x 129 files
str11 = '(a48,3i4)'
str1 = '(2x,a5,4x,a5,a10,a6,2x,a4,a2,8x,3i4)'
str3 = '(2i5)'
str4 = '(5e17.9)';changed 5e16.9 to 5e17.9
str5 = '(4i5)'
str6 = '(2i5,i6,i5)'
str7 = '(i5,e16.9,i5)'
str8 = '(1x,a42,1x,a3)'
str9 = '(a10,5e16.9)'
str10 = '(i4)'

strtitle = ''
readF, LUN,strtitle,idum,nw,nh,FORMAT = str11
nw = fix(nw) & nh=fix(nh)
print,'nw,nh',nw,nh

readF, LUN,xdim,zdim,rzero,rgrid,zmid,FORMAT = str4

readF, LUN,rmaxis,zmaxis,ssim,ssib,bcentr ,FORMAT = str4
readF, LUN,cpasma,ssim,xdum,rmaxis,xdum ,FORMAT = str4
readF, LUN, zmaxis,xdum,ssib,xdum,xdum ,FORMAT = str4

fpol = fltarr(nw)
pres = fltarr(nw)
ffprim = fltarr(nw)
pprime = fltarr(nw)
psirz = fltarr(nw,nh)
psirzt = fltarr(nh)
qpsi = fltarr(nw)

readF, LUN,fpol ,FORMAT = str4
readF, LUN,pres ,FORMAT = str4
readF, LUN,ffprim ,FORMAT = str4
readF, LUN,pprime ,FORMAT = str4
for i=0,nh-1 do begin
  readF, LUN,psirzt,FORMAT = str4
  psirz(*,i) = psirzt
endfor
readF, LUN,qpsi ,FORMAT = str4
readF, LUN, nbbbs,limitr ,FORMAT = str3

nbbbs = fix(nbbbs)
limitr = fix(limitr)
ii =  2*indgen(nbbbs)
rbbbs = fltarr(nbbbs)
zbbbs = fltarr(nbbbs)
for i=0,nbbbs-1 do begin
   readF, LUN, rb1,zb1,FORMAT = str4
   rbbbs(i) = rb1
   zbbbs(i) = zb1
endfor

xlim = fltarr(limitr)
ylim = fltarr(limitr)
for i=0,limitr-1 do begin
   readF, LUN, xb1,yb1,FORMAT = str4
   xlim(i) = xb1
   ylim(i) = yb1
endfor

Close,lun
FREE_LUN, LUN

efit = create_struct('nw',nw,'nh',nh,'nprof',nw, $
       'xdim',xdim,'zdim',zdim,'rzero',rzero,'rmin',rgrid,'zmid',zmid,$
 'rmaxis',rmaxis,'zmaxis',zmaxis,'ssimag',ssim,'ssibry',ssib,'bcentr',bcentr,$
'cpasma',cpasma,'fpol',fpol,'pres',pres,'ffprim',ffprim,'pprime',pprime,$
'psirz',psirz,'qpsi',qpsi,'nbbbs',nbbbs,'limitr',limitr,$
'rbbbs',rbbbs,'zbbbs',zbbbs,'xlim',xlim,'ylim',ylim)


end
