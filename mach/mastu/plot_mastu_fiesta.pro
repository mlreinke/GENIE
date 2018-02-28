; Purpose: Plot MAST-U equilibrium from the FIESTA code
; Author: Devon Battaglia

pro plot_mastu_fiesta,fiesta_filename,ps_filename=ps_filename, $
                      rrange=rrange, zrange=zrange

  if not keyword_set(ps_filename) then ps_filename = 'mastu_fiesta_eq.ps'
  if not keyword_set(rrange) then rrange = [0,2.2]
  if not keyword_set(zrange) then zrange = [-2.5,2.5]

; Load MAST-U limiter, coil and wall structure
  
  mastu_limiter_surface = '/u/dbattagl/LRDFIT/lrdfit3/MAST/2017/verStandard/device_definition/limiter_surface'
  a = read_ascii(mastu_limiter_surface)
  xlim_mastu = (a.field1)(0,*)
  ylim_mastu = (a.field1)(1,*)
  
  restore,'/u/dbattagl/LRDFIT/lrdfit3/MAST/diagsys/diagsys_reduced.MAST.2017.Standard'
  s_mastu = S
  ;help,S.full.cedata,/str
 
; Load FIESTA equilibrium

  read_fiesta_eqdsk,fiesta_filename,g
  
  psirz = g.psirz
  psi_axis = g.ssimag
  psi_bndry = g.ssibry
  xlim = g.xlim
  ylim = g.ylim
  xdim = g.xdim
  zdim = g.zdim
  nw = g.nw
  nh = g.nh
  rmin = g.rmin
  rmaxis = g.rmaxis
  zmaxis = g.zmaxis
  rbbs = g.rbbbs
  zbbs = g.zbbbs 
  ;help,g_mastu,/str

; Flux surface plot computations
  
  ;Normalize psi 
  psin = (psirz - psi_axis)/(psi_bndry - psi_axis)
  
  ; Define x,y axis of contour plots
  xaxis = (xdim) * findgen(nw)/(nw-1) + rmin
  yaxis = (zdim) * findgen(nh)/(nh-1) - (zdim)/2 ;assume zmid = 0
  
  ;Find x-points (assume at max/min z points on boundary) 
  zxpnt_lower = min(zbbs,ind)
  rxpnt_lower = rbbs(ind)
  zxpnt_upper = max(zbbs,ind)
  rxpnt_upper = rbbs(ind)
  
  ;Psin for plotting surfaces inside separatrix (outside separatrix = 1)
  psin_inside = psin
  psin_inside(*,where( (yaxis lt zxpnt_lower) or (yaxis gt zxpnt_upper))) = 1.0
  psin_inside(where(psin_inside gt 1)) = 1.0
 
  ;Psin for SOL and PFR (inside separatrix = 1)
  psin_outside = psin
  psin_outside(where(psin_inside lt 1)) = 1.0
  
  ;Psin for SOL and PFR where outside limiter = 10
  xlimind = interpol(findgen(nw),xaxis,xlim_mastu)
  ylimind = interpol(findgen(nh),yaxis,ylim_mastu)
  inside_lim = polyfillv(xlimind,ylimind,nw,nh)
  psin_sol_pfr = 0.0*psin+10.0
  psin_sol_pfr(inside_lim) = psin_outside(inside_lim)
  

; Set plotting attributes  
 
    !p.multi=0
    set_plot, 'ps'
    device,/portrait,/inches,xoffset=0,xsize=8.5,yoffset=0,ysize=11
    device,/color,filename=ps_filename
    device,/helvetica,font_index=5,font_size=10
 
    pthick=!p.thick
    pcharthick=!p.charthick
    xthick=!x.thick
    ythick=!y.thick
    pcharsize = !p.charsize
      
    !p.thick=2
    !p.charthick=4
    !x.thick=6
    !y.thick=6
    !p.charsize=2. 

; Plot x and y axis
  ;'right_aspect' found in:
  ; /p/nstxops/cvs/idl_cvs
  ; /u/bdavis/GArelease/general
  ; /p/nstxops/cvs/ga/efitviewer_nstx/general
  
  loadct,0,/silent
  right_aspect,rrange,zrange,position=position
  plot,rrange,zrange,position=position,/nodata,/xstyle,/ystyle,ytitle = 'Z (m)', xtitle = 'R (m)'

; Plot flux surfaces

  ;Separatrix and strikepoints
  contour,psin,xaxis,yaxis,/overplot,xrange=rrange,yrange=yrange,/xstyle,/ystyle,$
    level = [1.0],thick=6,color=0

  ;Inner flux surfaces
   nsurfaces = 9
   delta_surface = 1.0/(nsurfaces + 1.0)
   contour,psin_inside,xaxis,yaxis,/overplot,xrange=rrange,yrange=yrange,/xstyle,/ystyle,$
    level = delta_surface*(1.0+findgen(nsurfaces)),thick=4,color=30
  
  ;SOL
   nsurfaces = 5
   delta_surface = 0.02 ;psin
   contour,psin_sol_pfr,xaxis,yaxis,/overplot,xrange=rrange,yrange=yrange,/xstyle,/ystyle,$
    level = 1.0+delta_surface*(1.0+findgen(nsurfaces)),thick=4,color=90
  
  ;PFR
   nsurfaces = 5
   delta_surface = 0.02 ;psin
   contour,psin_sol_pfr,xaxis,yaxis,/overplot,xrange=rrange,yrange=yrange,/xstyle,/ystyle,$
    level = (1.0-delta_surface*(1.0+nsurfaces))+0.02*(1.0+findgen(nsurfaces)),thick=4,color=90
 
  ;Magnetic axis
   oplot,[rmaxis],[zmaxis],psym=1,color=0

; Plot MAST-U wall elements
;  Note: at this point, only contains toroidally continuous elements
  
  s = s_mastu
  passive_number = s.full.cedata.passive_numbers
  
  ind = where(passive_number ge 1.0)
  
  r	= (s.full.cedata.r)(ind)
  z	= (s.full.cedata.z)(ind)
  w	= ((s.full.cedata.w)(ind))/2.
  h	= ((s.full.cedata.h)(ind))/2.
  rv	= [[r],[r+w],[r-w],[r-w],[r+w]]
  zv	= [[z],[z+h],[z+h],[z-h],[z-h]]
   
   loadct, 0, /silent
   
   FOR i=0, N_ELEMENTS(r)-1 DO BEGIN
     rvec	= TRANSPOSE([[rv(i,1:4)],[rv(i,1)]])
     zvec	= TRANSPOSE([[zv(i,1:4)],[zv(i,1)]])
     POLYFILL,rvec,zvec,color=120,NOCLIP=0
     OPLOT,rvec,zvec,color=120
   ENDFOR

; Plot MAST-U coils
   
   s = s_mastu
   coil_name = s.full.cedata.coil_names

   ind = where(strmatch(coil_name,'P4*',/fold_case) or strmatch(coil_name,'D*',/fold_case) or $
               strmatch(coil_name,'P5*',/fold_case) or strmatch(coil_name,'P1*',/fold_case)or $
	       strmatch(coil_name,'PX*',/fold_case) or strmatch(coil_name,'P6*',/fold_case)or $
	       strmatch(coil_name,'PC*',/fold_case))
	       
   r	= (s.full.cedata.r)(ind)
   z	= (s.full.cedata.z)(ind)
   w	= ((s.full.cedata.w)(ind))/2.
   h	= ((s.full.cedata.h)(ind))/2.
   rv	= [[r],[r+w],[r-w],[r-w],[r+w]]
   zv	= [[z],[z+h],[z+h],[z-h],[z-h]]
  
   loadct,0,/silent

   FOR i=0, N_ELEMENTS(r)-1 DO BEGIN
     rvec	= TRANSPOSE([[rv(i,1:4)],[rv(i,1)]])
     zvec	= TRANSPOSE([[zv(i,1:4)],[zv(i,1)]])
     POLYFILL,rvec,zvec,color=120,NOCLIP=0
     OPLOT,rvec,zvec,color=0
   ENDFOR

; Plot limiter

  loadct,0,/silent
  oplot,xlim_mastu,ylim_mastu,color=0,thick=4,linestyle=2
 
 ; Reset plotting attributes
    
    device,/close
    set_plot,'x'
    !p.background=255
    !p.color=0
    !p.thick=pthick
    !p.charthick=pcharthick
    !x.thick=xthick
    !y.thick=ythick

end
