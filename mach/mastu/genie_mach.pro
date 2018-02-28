;+
;NAME:
;	VESSEL_PLOT
;
;PURPOSE:
;	This procedures sends a properly scaled plot of a poloidal  vacuum vessel/tile cross-section 
;	to either the postscript or x-windows device.
;
;CALLING SEQUENCE:
;	VESSEL_PLOT
;
;OPTIONAL INPUTS:
;	n:		INT:	The window number to plot to
;	color:		INT:	Color of the vessel plot lines (not color of plot) DEFAULT: 0
;	title:		STR:	String of the plot title DEFAULT: "MAST-U"
;	thick:		INT:	Thickness of the vessel plot lines DEFAULT: 0
;	shot:		LON:	Shot number to determine which vessel cross-section to use DEFAULT: 135440
;
;KEYWORD PARAMETERS:
;	force:		/force will force window number n to reopen and resize it accordingly
;	oplot:		/oplot will oplot the vessel structure on the currently selected plotting device
;	div:		/div will create a properly scaled zoomed plot of the lower divertor
;	edge:		/edge will create a properly scaled zoomed plot of the outer edge
;
;OUTPUTS:
;	All outputs are sent to the currently selected plotting device.  VESSEL_PLOTS will do nothing if a device other than
;	the X or PS is selected.
;
;PROCEDURE:
;	The vessel cross-section plots are taken from:
;	/u/dbattagl/LRDFIT/lrdfit3/MAST/diagsys/diagsys_reduced.MAST.2017.Standard
;	rlim=s.field.lds.rlimh
;	zlim=s.field.lds.zlimh	
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke, Fall 2006
;	5-10-07:	ML Reinke - added the new cyropump vessel cross-section and made it the default
;	7-02-15:	ML Reinke - modified for use with NSTX
;	11-02-17:	ML Reinke - modified for use with MAST-U
;	2-28-18:	ML Reinke - modified to use limiter/efit from UDA
;-

PRO vessel_plot,ps=ps,n=n, force=force,title=title,x_offset=x_offset,y_offset=y_offset,color=color,thick=thick,oplot=oplot,$
		shot=shot,d_old=d_old,div=div,edge=edge,pedestal=pedestal,wall=wall
	IF NOT keyword_set(shot) THEN shot=50000
	IF NOT keyword_set(title) THEN BEGIN
		title='MAST-U'
	ENDIF
	IF NOT keyword_set(x_offset) THEN x_offset=0.0
	IF NOT keyword_set(y_offset) THEN y_offset=0.0
	
	CASE !d.name OF
		'PS': 	ps=1
		'X':	ps=0
		ELSE:	RETURN
	ENDCASE


	IF NOT keyword_set(div) AND NOT keyword_set(edge) AND NOT keyword_set(pedestal) AND NOT keyword_set(wall)  THEN BEGIN
		yrange=[-2.2,2.2]
		xrange=[0.0,1.9]
		IF keyword_set(ps) THEN BEGIN
			xsize=4.75
			ysize=8.98
		ENDIF ELSE BEGIN
			xsize=451.0
			ysize=890.0
		ENDELSE
	ENDIF ELSE BEGIN
		IF keyword_set(div) THEN BEGIN
			yrange=[-2.1,-0.6]
			xrange=[0.2,1.9]
			IF keyword_set(ps) THEN BEGIN
				xsize=7.0
				ysize=6.3
			ENDIF ELSE BEGIN
				xsize=672.0
				ysize=650.0
			ENDELSE
		ENDIF
; 		IF keyword_set(edge) THEN BEGIN
; 			yrange=[-0.35,0.35]
; 			xrange=[0.574,0.95]
; 			IF keyword_set(ps) THEN BEGIN
; 				xsize=6.0
; 				ysize=6.0*810.0/500.0
; 			ENDIF ELSE BEGIN
; 				xsize=500
; 				ysize=818
; 			ENDELSE	
; 		ENDIF
; 		IF keyword_set(wall) THEN BEGIN
; 			yrange=[-0.35,0.35]
; 			xrange=[0.574,0.95]-0.174
; 			IF keyword_set(ps) THEN BEGIN
; 				xsize=6.0
; 				ysize=6.0*810.0/500.0
; 			ENDIF ELSE BEGIN
; 				xsize=500
; 				ysize=818
; 			ENDELSE	
; 		ENDIF
; 		IF keyword_set(pedestal) THEN BEGIN
; 			yrange=[-0.1,0.1]
; 			xrange=[0.822,0.93]
; 			IF keyword_set(ps) THEN BEGIN
; 				xsize=6.0
; 				ysize=6.0*810.0/500.0
; 			ENDIF ELSE BEGIN
; 				xsize=500
; 				ysize=818
; 			ENDELSE	
; 		ENDIF
	ENDELSE

	IF NOT keyword_set(oplot) THEN BEGIN
		IF NOT keyword_set(ps) THEN BEGIN
			IF NOT keyword_set(n) THEN n=19
			device, window_state=var
			IF var[n] EQ 0 OR keyword_set(force) EQ 1 THEN window,n,xsize=xsize,ysize=ysize,xpos=0,ypos=670,title='vessel cx,'+num2str(n) $
				ELSE wset,n
		ENDIF ELSE BEGIN
			d_old=!d
			device, xsize=xsize, ysize=ysize, /inches
		ENDELSE
		plot, [0],[0],title=title,chars=1.3,xrange=xrange,yrange=yrange,xtit='R (m)',ytit='Z (m)',/xsty, /ysty,/iso
	ENDIF

	limiter=getdata('GEOM::/limiter/efit',shot)    
	rlim=limiter.data.r
	zlim=limiter.data.z
	oplot,rlim,zlim

	;IF NOT keyword_set(oplot) THEN IF keyword_set(ps) THEN device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm
END

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
;print,'nw,nh',nw,nh

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
	



;+
;NAME:
;	FS_PLOT
;
;PURPOSE:
;	This procedure OPLOTS flux surfaces, presumabley onto a window where VESSEL_PLOT has been called.
;
;CALLING SEQUENCE
;	FS_PLOT,shot,time
;
;INPUTS:
;	shot:	LONG	shot number
;	time:	FLT	time point [sec]
;
;MODIFICATION HISTORY:
;	6/18/11		M.L. Reinke - modified to allow non-ANALYSIS EFIT's to be used
;	7/2/15		M.L. Reinke - modified to work with NSTX EFIT system
;-

PRO fs_plot,shot,time,num=num,color=color,style=style,thick=thick,sol=sol,nsol=nsol,rad_sol=rad_sol,invssel=invessel,tree=tree
	restore,'/u/mreinke/GENIE/mach/mastu/mastu_limiter.dat'	;replace with something that is more official at some point
	xlim_mastu=rlim
	ylim_mastu=zlim

	fiesta_filename='/home/mreinke/GENIE/mach/mastu/g20171023'
	
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
  



  	; Plot flux surfaces

  	;Separatrix and strikepoints
  	contour,psin,xaxis,yaxis,/overplot,xrange=rrange,yrange=yrange,/xstyle,/ystyle,$
    		level = [1.0],thick=1,color=200

  	;Inner flux surfaces
   	nsurfaces = 9
   	delta_surface = 1.0/(nsurfaces + 1.0)
   	contour,psin_inside,xaxis,yaxis,/overplot,xrange=rrange,yrange=yrange,/xstyle,/ystyle,$
    		level = delta_surface*(1.0+findgen(nsurfaces)),thick=1,color=100
  
  	;SOL
   	nsurfaces = 5
   	delta_surface = 0.02 ;psin
   	contour,psin_sol_pfr,xaxis,yaxis,/overplot,xrange=rrange,yrange=yrange,/xstyle,/ystyle,$
    		level = 1.0+delta_surface*(1.0+findgen(nsurfaces)),thick=1,color=30
  
  	;PFR
   	nsurfaces = 5
   	delta_surface = 0.02 ;psin
   	contour,psin_sol_pfr,xaxis,yaxis,/overplot,xrange=rrange,yrange=yrange,/xstyle,/ystyle,$
    		level = (1.0-delta_surface*(1.0+nsurfaces))+0.02*(1.0+findgen(nsurfaces)),thick=1,color=100
 
  	;Magnetic axis
  	oplot,[rmaxis],[zmaxis],psym=1,color=0
END

;+
;NAME:
;	LINE_PATH_PLOTS
;	
;PURPOSE:
;	This procedure plots the line of sight of an arbitrary number of POS vectors
;	in both poloidal and toroidal projections.
;
;CALLING SEQUENCE:
;	LINE_PATH_PLOTS,pos
;
;INPUTS:
;	pos: 	FLTARR 2D (m,4) defines the diagnostic view(s)
;
;OPTIONAL INPUTS:
;	shot:	LONINT of the shot number to display flux surface information DEFAULT=1050426022 (everyone's favorite)
;	tpt:	FLT of the time point to display flux surface information DEFAULT=1.0
;	thick:	FLT of the line thickness DEFAULT=2.0
;	fovpos:	FLTARR of pos vectors for the field of view error bars (see GENPOS)
;	rzbnd:	FLTARR boundary [r_in, r_out, abs(z)] for line of sight [DEFAULT = [0.44,1.0,0.6]]
;
;KEYWORD PARAMETERS:
;	div:	/div will make the poloidal view zoomed into the divertor region
;	tang:	/tang will display tangency points views that are constrained to a single Z-plane
;	verb:	/verb display line of sight info (same as LINE_BR)
;	debug:	/debug stops at various locations
;	ps:	/ps supresses WSET command.  Postscript device must be setup prior to call
;
;OUTPUTS:
;	Outputs to the X-windows plotting device will be made on windows 19 and 20.  If currently being
;	used, close so that LINE_PATH_PLOTS resizes them appropriately.  Note that switching between divertor
;	and full view requires the window to be closed for resizing.
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - Spring 2006
;	8-24-06:	ML Reinke - added the rzbnd optional input 
;	7-13-07:	ML Reinke - made this program use the vessel cross-section from VESSEL_PLOT in MLR_FUNCTIONS
;-


PRO line_path_plots,pos,shot=shot,tpt=tpt,ps=ps,verb=verb,div=div,debug=debug,thick=thick,tang=tang,fovpos=fovpos,rzbnd=rzbnd,nsol=nsol,$
		sol=sol,tor=tor,edge=edge,bspline=bspline,n_s=n_s,inv=inv,col=col,tree=tree,win=win,oplot=oplot

	IF NOT keyword_set(tor) THEN notor=1 ELSE notor=0
	IF NOT keyword_set(col) THEN col=200
	IF NOT keyword_set(shot) THEN shot=135440	;shot for flux surfraces
	IF NOT keyword_set(tpt) THEN tpt =1.0		;time point for flux surfaces
	IF NOT keyword_set(thick) THEN thick=2.0
	IF NOT keyword_set(n_s) THEN n_s=100	;set number of segments to define line of sight
	IF NOT keyword_set(rzbnd) THEN rzbnd=[0.26,1.8,2.06]	;set path's [r_min,r_max,abs(z_max)] limit
	IF NOT keyword_set(nsol) THEN nsol=4
	IF NOT keyword_set(tree) THEN tree='efit01'
	
	;determine the number of views
	xpos=size(pos)	
	IF xpos[0] EQ 1 THEN n_view=1 ELSE n_view=xpos[2]

	IF NOT keyword_set(oplot) THEN BEGIN
		IF shot GT 0 THEN BEGIN
			vessel_plot,shot=shot,div=div,n=win
			fs_plot,shot,tpt,tree=tree
	        ENDIF ELSE BEGIN
			vessel_plot,div=div,/isnstxu,n=win
                    ENDELSE	
	ENDIF
	;debug information
	IF keyword_set(debug) THEN STOP

	;start plotting each view
	smin=fltarr(n_view)
	FOR i=0L,n_view-1 DO BEGIN
		IF keyword_set(verb) THEN print, ' i = '+num2str(i,1)

		;determine where the line of sight terminates in s due to radial boundries
		IF pos[2,i] LT rzbnd[0] THEN BEGIN	;if r_tang < r_min
			svec=line_s(pos[*,i],r=rzbnd[0])
			smin_r=svec[1]	 ;choose the smallest s where r=rmin
		ENDIF ELSE BEGIN			;if r_tang > r_min
			svec=line_s(pos[*,i],r=rzbnd[1])
			smin_r=svec[0]	 ;choose the largest (and non-negative) s where r=rmax
 		ENDELSE

		;determine where the line of sight terminates in s due to z boundries
		IF pos[3,i] NE 0 THEN BEGIN
			IF line_s(pos[*,i],z=rzbnd[2]) LT 0 THEN $
				smin_z=line_s(pos[*,i],z=-1.0*rzbnd[2]) ELSE $
				smin_z=line_s(pos[*,i],z=rzbnd[2])
		ENDIF ELSE smin_z=2.0*smin_r	;if no inclination, smin_z set > smin_r	
		IF smin_r GT smin_z THEN smin[i] = smin_z ELSE smin[i]=smin_r		;set minimum s
		s=make(0,smin[i],n_s)
		IF keyword_set(verb) THEN print, s		
		;there is now array of points that go from the pinhole to the closest limiting
		;surface and is divided into n_s parts.  

		;plot the line of sight on the polodial plane
		
		IF total(abs(pos[3,*])) NE 0 THEN BEGIN
			rplot=line_r(s,pos[*,i])
			zplot=line_z(s,pos[*,i])
			IF keyword_set(inv) THEN BEGIN
				good=line_invessel(rplot,zplot,shot=shot)
				tmp=where(good EQ 0)
				oplot,rplot[0:tmp[0]],zplot[0:tmp[0]],color=col,thick=thick
			ENDIF ELSE oplot,rplot,zplot,color=col,thick=thick
		ENDIF
		IF keyword_set(tang) AND NOT keyword_set(fovpos) THEN BEGIN
			makesym,10
			oplot, [pos[2,i]],[pos[1,i]-tan(pos[3,i])*sqrt(pos[0,i]^2-pos[2,i]^2)],psym=8,color=col
		ENDIF
		IF keyword_set(fovpos) THEN BEGIN
			out=pos2width(pos,fovpos)
			oploterror,[pos[2,i]],[pos[1,i]-tan(pos[3,i])*sqrt(pos[0,i]^2-pos[2,i]^2)],out.r[i],$
				out.z[i],psym=8,color=col,errcolor=col
		ENDIF

	ENDFOR

	IF NOT keyword_set(notor) THEN BEGIN
		;if plotting to x-windows then setup top-down view
		IF NOT keyword_set(ps) THEN BEGIN
			device, window_state=var
			IF var[20] EQ 0 THEN window,20,xsize=585,ysize=785,xpos=500,ypos=670,title='vessel cx,20' ELSE wset,20
		ENDIF ELSE BEGIN
			device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm
		ENDELSE
			
		alpha=make(0,2*!pi,100)

		;plot inner wall
		x=0.44*cos(alpha)
		y=0.44*sin(alpha)	
		plot,x,y,xr=[-0.1,1.0],yr=[-1.0,0.6],title='SHOT: '+num2str(shot,1)+' TIME: '+num2str(tpt,dp=3),chars=1.3,/xsty,/ysty

		;plot outerlimiter
		rlim=max(r_rf)
		x=rlim*cos(alpha)
		y=rlim*sin(alpha)	
		oplot,x,y	

		;plot inner LCFS
		y=min(rb_plt)*sin(alpha)
		x=min(rb_plt)*cos(alpha)
		oplot,x,y,color=100,thick=2.0

		;plot outer LCFS
		y=max(rb_plt)*sin(alpha)
		x=max(rb_plt)*cos(alpha)
		oplot,x,y,color=100,thick=2.0

		FOR i=0,n_view-1 DO BEGIN
			ipos=reform(pos[*,i])
			s=make(0,smin[i],n_s)
			xpts=line_r(s,ipos)*cos(line_th(s,ipos,/cont)-!pi/2.0)
			ypts=line_r(s,ipos)*sin(line_th(s,ipos,/cont)-!pi/2.0)
					

			;plot line of sight on z=const plane
			oplot,xpts,ypts,color=200,thick=2.0
	
		ENDFOR
	ENDIF 
	IF keyword_set(debug) THEN stop

END
