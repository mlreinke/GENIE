
PRO wspec_load_spred,shot,specbr,lam,time,short=short,status=status
        mdsopen,'passivespec',shot
        t=mdsvalue('\PASSIVESPEC::TOP.SPRED.TIMES',/quiet,status=tstatus)
        d=mdsvalue('_sig=\PASSIVESPEC::TOP.SPRED.RAWDATA.SPECTRA',/quiet,status=dstatus) ; [nwl nch nt]
     	C=mdsvalue('\PASSIVESPEC::TOP.SPRED.POLY_COEFFS',/quiet,status=lstatus)
	mdsclose,'passivespec',shot
	nl=1024
	IF tstatus AND dstatus AND dstatus THEN BEGIN
		time=float(t[1:n(t)-1])	;remove first and last frame
	        IF keyword_set(short) THEN BEGIN
			specbr=float(reform(d[*,0,*])) 
			lam=poly(indgen(nl),C[0:3])
	        ENDIF ELSE BEGIN
			specbr=float(reform(d[*,1,*]))
			lam=poly(indgen(nl),C[4:7])/10.0
                ENDELSE
		specbr=specbr[*,1:*]
		status=1
       ENDIF ELSE status=0
END


PRO wspec_load_vuv,shot,specbr,lam,time,xeus=xeus,loweus=loweus,status=status	
	IF keyword_set(xeus) THEN path='\PSPEC_PC::TOP.XEUS'
	IF keyword_set(loweus) THEN path='\PSPEC_PC::TOP.LOWEUS'	
	mdsopen,'pspec_pc',shot
	specbr=float(mdsvalue('_sig='+path+':IMAGE',/quiet,status=status))
	lam=mdsvalue('dim_of(_sig,0)',/quiet)
	time=float(mdsvalue('dim_of(_sig,2)',/quiet))/1.0e3
	mdsclose,'pspec_pc',shot
END

PRO wspec_load_xeus,shot,specbr,lam,time,status=status
	wspec_load_vuv,shot,ispecbr,lam,time,/xeus,status=status
;	IF NOT status THEN BEGIN
;		load_xeus,shot,cnts
;		nframes=n(cnts[0,*])+1
;		xeus_lam_time,lam
;		specbr=cnts/65535.0
;		time=-0.05+indgen(nframes)*1.977/1.0e3
;	ENDIF
;	status=1
	x=size(ispecbr)
	specbr=fltarr(x[1],x[3])
	FOR i=0,x[1]-1 DO FOR j=0,x[3]-1 DO specbr[i,j]=median(ispecbr[i,*,j])
END

PRO wspec_load_loweus,shot,specbr,lam,time,status=status
	wspec_load_vuv,shot,ispecbr,lam,time,/loweus,status=status
;	IF NOT status THEN BEGIN
;		load_loweus,shot,cnts
;		nframes=n(cnts[0,*])+1
;		lowes_lam_time,lam
;		specbr=cnts/65535.0
;		time=-0.05+indgen(nframes)*1.977/1.0e3
;	ENDIF
;	status=1
	x=size(ispecbr)
	specbr=fltarr(x[1],x[3])
	FOR i=0,x[1]-1 DO FOR j=0,x[3]-1 DO specbr[i,j]=median(ispecbr[i,*,j])
END

PRO wspec_load_dalpha,shot,time,dalpha,status=status
	mdsopen,'passivespec',shot
       	dalpha=mdsvalue('_sig=\PASSIVESPEC::TOP.FILTERED_VIS.BAYC_OPIPE:DALPHA_EIES',/quiet,status=status)
        time=mdsvalue('dim_of(_sig,0)',/quiet)
        mdsclose,'passivespec',shot
END

PRO wspec_load_corextomo,shot,time,br,ch=ch,status=status
	IF NOT keyword_set(ch) THEN ch=5
	IF ch LT 10 THEN ch_str='0' ELSE ch_str=''
	pt='HDOWN.CHORD_'+ch_str+num2str(ch,1)
	mdsopen,'usxr',shot
        br=mdsvalue('_sig=\USXR::TOP.'+pt,/quiet,status=status)
        time=mdsvalue('dim_of(_sig,0)',/quiet)
        mdsclose,'usxr',shot
	

END

PRO wspec_load_dens,shot,time,dens,status=status
        mdsopen,'activespec',shot
        dens=mdsvalue('_sig=\ACTIVESPEC::TOP.MPTS.OUTPUT_DATA.BEST:FIT_NE',/quiet,status=status)
        time=mdsvalue('dim_of(_sig,0)',/quiet)
        mdsclose,'activespec',shot
	dens/=1.0e14	;put in 10^20/m^3

	IF status THEN BEGIN
		x=size(dens,/dimensions)
		nt=x[0]
		nr=x[1]
		dens=total(dens,2)/nr
	ENDIF
END

PRO wspec_load_nbipow,shot,time,pow,status=status
        mdsopen, 'nbi', shot
        pow=mdsvalue('_sig=\NBI::TOP.ANALYSIS:P_INJ',/quiet,status=status) ;[MW]
        time=mdsvalue('dim_of(_sig,0)',/quiet)
        mdsclose
END

PRO wspec_load_hhfwpow,shot,time,pow,status=status
        mdsopen, 'rf', shot
        pow=mdsvalue('_sig=\RF::TOP.HHFW.POWERSIGS:HHFW_POWER',/quiet,status=status)/1.0e6 ; [MW]
        time=mdsvalue('dim_of(_sig,0)',/quiet)
        mdsclose
END

;6/15/12 - moved line list into GENIE
PRO wspec_load_linestr,line,convert=convert
	;csvpath='/usr/local/cmod/codes/spectroscopy/vuv/vuv_line_list.csv'
	;datpath='/usr/local/cmod/codes/spectroscopy/vuv/vuv_line_list.dat'
	
	;csvpath='/usr/local/cmod/idl/GENIE/IMPSPEC/vuv_line_list.csv'
	;datpath='/usr/local/cmod/idl/GENIE/IMPSPEC/vuv_line_list.dat'
	;IF keyword_set(convert) THEN csv_convert,csvpath,datpath
	;restore,datpath

	GENIE_PATH=getenv('GENIE_PATH') & if GENIE_PATH eq '' then GENIE_PATH='/u/mreinke/GENIE/'
	tsvpath=GENIE_PATH+'IMPSPEC/vuv_line_list.tsv'
	IF keyword_set(convert) THEN print,'using tsv_read instead of convert...'
	tsv_read,tsvpath

	list=strtrim(elem[0])
	ielem=strtrim(elem[0])
	FOR i=1,n(elem) DO BEGIN
		IF strtrim(elem[i]) NE ielem THEN list=list+','+strtrim(elem[i])
		ielem=strtrim(elem[i])
	ENDFOR
	line={lam:wave,elem:elem,cs:cs,note:note,iso:iso,source:source,priority:priority,list:list}

END

PRO wspec_max_tplot,u
	xr=[u.plot.x0[0],u.plot.x1[0]]
	IF size(u.dat.dalpha,/type) EQ 8 THEN BEGIN
		tmp=where(u.dat.dalpha.time GE xr[0] AND u.dat.dalpha.time LE xr[1])
		IF tmp[0] NE -1 THEN u.plot.y1[0]=max(u.dat.dalpha.dalpha[tmp])*1.05
        ENDIF
	IF size(u.dat.xtomo,/type) EQ 8 THEN BEGIN
		tmp=where(u.dat.xtomo.time GE xr[0] AND u.dat.xtomo.time LE xr[1])
		IF tmp[0] NE -1 THEN u.plot.y1[1]=max(u.dat.xtomo.br[tmp])*1.05
        ENDIF
	IF size(u.dat.dens,/type) EQ 8 THEN BEGIN
		tmp=where(u.dat.dens.time GE xr[0] AND u.dat.dens.time LE xr[1])
		IF tmp[0] NE -1 THEN u.plot.y1[2]=max(u.dat.dens.dens[tmp])*1.05
        ENDIF	
	IF size(u.dat.nbi,/type) EQ 8 THEN BEGIN
		tmp=where(u.dat.nbi.time GE xr[0] AND u.dat.nbi.time LE xr[1])
		IF tmp[0] NE -1 THEN u.plot.y1[3]=max(u.dat.nbi.pow[tmp])*1.05
        ENDIF
	IF u.plot.y1[3] LT 1.0 THEN u.plot.y1[3]=1.1
	ymax=0.1
	FOR i=0,n(u.line.plot) DO BEGIN
		IF u.line.plot[i] THEN BEGIN
			x=*u.dat.br[i+4]
			y=*u.dat.br[i]
			tmp=where(x GE xr[0] AND x LE xr[1])
			IF tmp[0] NE -1 THEN ymax=ymax > max(y[tmp])
                ENDIF
        ENDFOR
	u.plot.y1[4]=ymax*1.05
END

PRO wspec_calc_br,u,line
	lamc=0.5*(u.line.lam0[line]+u.line.lam1[line])
	spec='none'
	IF lamc GT 10 AND lamc LT 90 THEN spec='xeus'
	IF lamc GT 100 AND lamc LT 300 THEN BEGIN
		IF u.stat.low THEN spec='loweus'
	ENDIF
	IF lamc GT 300 THEN spec='spredl'

	IF spec EQ 'xeus' AND u.stat.xeus NE 1 THEN spec='x'
	IF spec EQ 'loweus' AND u.stat.low NE 1 THEN spec='x'
	IF spec EQ 'spredl' AND u.stat.spredl NE 1 THEN spec='x'


	CASE spec OF 
		'xeus' : BEGIN
			ntime=n(u.dat.xeus.time)+1
			*u.dat.br[4+line]=u.dat.xeus.time
			tmp=where(u.dat.xeus.lam GE u.line.lam0[line] AND u.dat.xeus.lam LE u.line.lam1[line])
			br=sum_array(u.dat.xeus.specbr[tmp,*],/j)-0.5*(u.dat.xeus.specbr[tmp[0],*]+u.dat.xeus.specbr[last(tmp),*])*(n(tmp)+1)
			*u.dat.br[line]=br
                END
		'loweus' : BEGIN
			ntime=n(u.dat.low.time)+1
			*u.dat.br[4+line]=u.dat.low.time
			tmp=where(u.dat.low.lam GE u.line.lam0[line] AND u.dat.low.lam LE u.line.lam1[line])
			br=sum_array(u.dat.low.specbr[tmp,*],/j)-0.5*(u.dat.low.specbr[tmp[0],*]+u.dat.low.specbr[last(tmp),*])*(n(tmp)+1)
			*u.dat.br[line]=br
                END
		'spredl' : BEGIN
			ntime=n(u.dat.spredl.time)+1
			*u.dat.br[4+line]=u.dat.spredl.time
			tmp=where(u.dat.spredl.lam GE u.line.lam0[line] AND u.dat.spredl.lam LE u.line.lam1[line])
			br=sum_array(u.dat.spredl.specbr[tmp,*],/j)-0.5*(u.dat.spredl.specbr[tmp[0],*]+u.dat.spredl.specbr[last(tmp),*])*(n(tmp)+1)
			*u.dat.br[line]=br
                END
		'spreds' : BEGIN
			ntime=n(u.dat.spreds.time)+1
			*u.dat.br[4+line]=u.dat.spreds.time
			tmp=where(u.dat.spreds.lam GE u.line.lam0[line] AND u.dat.spreds.lam LE u.line.lam1[line])
			br=sum_array(u.dat.spreds.specbr[tmp,*],/j)-0.5*(u.dat.spreds.specbr[tmp[0],*]+u.dat.spreds.specbr[last(tmp),*])*(n(tmp)+1)
			*u.dat.br[line]=br
                END
		'helike' : BEGIN
			ntime=n(u.dat.hrxhe.time)+1
			*u.dat.br[4+line]=u.dat.hrxhe.time
			tmp=where(u.dat.hrxhe.lam GE u.line.lam0[line] AND u.dat.hrxhe.lam LE u.line.lam1[line])
			br=sum_array(u.dat.hrxhe.specbr[tmp,*],/j)-0.5*(u.dat.hrxhe.specbr[tmp[0],*]+u.dat.hrxhe.specbr[last(tmp),*])*(n(tmp)+1)
			*u.dat.br[line]=br
                 END
		ELSE : BEGIN
			;print, 'cannot calc_br for case : '+spec
			*u.dat.br[4+line]=make(0,u.stat.tmax,100)	;read in zero array
			*u.dat.br[line]=fltarr(100)
		END
        ENDCASE
	u.line.spec[line]=spec

END

PRO wspec_plot_time,u
	IF u.stat.ps THEN BEGIN
		xsize=6.0
		ysize=9.5
		ls=1.1
		tit=num2str(u.shot,1)+' t='+num2str(u.stat.time,dp=3)
		d_old=!d
		lcol=u.line.pscol
		device,xsize=xsize,ysize=ysize,/in
	ENDIF ELSE BEGIN
		ls=1.2
		widget_control,u.id.draw1,get_value=draw_win
		window,0,xsize=u.plot.asize[0],ysize=u.plot.asize[1],/pixmap
		lcol=u.line.col
		tit=''
	ENDELSE
	xr=[u.plot.x0[0],u.plot.x1[0]]
	tau=make(xr[0],xr[1],u.plot.ntau)
	nt=u.plot.nt
	nm=u.plot.nm
	low=0.07
	del=(0.96-low)/5
	tau0=[u.stat.time,u.stat.time]
	wspec_max_tplot,u

	i=4
	yr=[u.plot.y0[4-i],u.plot.y1[4-i]]
	pos=[0.13,low+i*del,0.9,low+(i+1)*del]
	plot,[0],[0],xr=xr,pos=pos,/xsty,chars=1.0*ls,xticks=xt,xminor=nm,xtickname=replicate(' ',nt+2),yr=yr,ysty=1,ytit='[AU]'
	xyouts,1.05*(xr[1]-xr[0])+xr[0],0.37*(yr[1]-yr[0])+yr[0],'DALPHA',chars=1.0*ls,orient=90
	IF size(u.dat.dalpha,/type) EQ 8 THEN oplot,tau,interpol(u.dat.dalpha.dalpha,u.dat.dalpha.time,tau)
	IF u.stat.tline THEN oplot,tau0,yr,linestyle=1.0,color=200
	IF u.stat.ps THEN xyouts,0.8*(xr[1]-xr[0])+xr[0],1.05*(yr[1]-yr[0])+yr[0],num2str(u.shot,1),chars=0.8*ls
	i=3
	yr=[u.plot.y0[4-i],u.plot.y1[4-i]]
	pos=[0.13,low+i*del,0.9,low+(i+1)*del]
	plot,[0],[0],xr=xr,pos=pos,/xsty,chars=1.0*ls,xticks=xt,xminor=nm,xtickname=replicate(' ',nt+2),yr=yr,/ysty,/noerase,ytit='[AU]'
	xyouts,1.05*(xr[1]-xr[0])+xr[0],0.2*(yr[1]-yr[0])+yr[0],'Core XTOMO',chars=1.0*ls,orient=90
	IF size(u.dat.xtomo,/type) EQ 8 THEN oplot,tau,interpol(u.dat.xtomo.br,u.dat.xtomo.time,tau)
	IF u.stat.tline THEN oplot,tau0,yr,linestyle=1.0,color=200

	i=2
	yr=[u.plot.y0[4-i],u.plot.y1[4-i]]
	pos=[0.13,low+i*del,0.9,low+(i+1)*del]
	plot,[0],[0],xr=xr,pos=pos,/xsty,chars=1.0*ls,xticks=xt,xminor=nm,xtickname=replicate(' ',nt+2),yr=yr,/ysty,/noerase,ytit='[10!u20!n m!u-3!n]'
	xyouts,1.05*(xr[1]-xr[0])+xr[0],0.45*(yr[1]-yr[0])+yr[0],'DENS',chars=1.0*ls,orient=90
	IF size(u.dat.dens,/type) EQ 8 THEN oplot,tau,interpol(u.dat.dens.dens,u.dat.dens.time,tau)
	IF u.stat.tline THEN oplot,tau0,yr,linestyle=1.0,color=200

	i=1
	yr=[u.plot.y0[4-i],u.plot.y1[4-i]]
	pos=[0.13,low+i*del,0.9,low+(i+1)*del]
	plot,[0],[0],xr=xr,pos=pos,/xsty,chars=1.0*ls,xticks=xt,xminor=nm,xtickname=replicate(' ',nt+2),yr=yr,/ysty,/noerase, ytit='[MW]'
	xyouts,1.05*(xr[1]-xr[0])+xr[0],0.25*(yr[1]-yr[0])+yr[0],'RF',chars=1.0*ls,orient=90,color=196
	xyouts,1.05*(xr[1]-xr[0])+xr[0],0.6*(yr[1]-yr[0])+yr[0],'NBI',chars=1.0*ls,orient=90
	IF size(u.dat.rf,/type) EQ 8 THEN BEGIN
		oplot,tau,interpol(u.dat.rf.pow,u.dat.rf.time,tau),color=196
	ENDIF
	IF size(u.dat.nbi,/type) EQ 8 THEN BEGIN
		oplot,tau,interpol(u.dat.nbi.pow,u.dat.nbi.time,tau)
	ENDIF
	IF u.stat.tline THEN oplot,tau0,yr,linestyle=1.0,color=200

	i=0
	yr=[u.plot.y0[4-i],u.plot.y1[4-i]]
	pos=[0.13,low+i*del,0.9,low+(i+1)*del]
	plot,[0],[0],xr=xr,pos=pos,/xsty,chars=1.0*ls,xticks=xt,xminor=nm,xtit='Time [sec]',yr=yr,/ysty,/noerase,ytit='[AU]'
	xyouts,1.05*(xr[1]-xr[0])+xr[0],0.15*(yr[1]-yr[0])+yr[0],'Line Brightness',chars=1.0*ls,orient=90
	IF u.stat.tline THEN oplot,tau0,yr,linestyle=1.0,color=200
	ylab=0.9*(yr[1]-yr[0])+yr[0]
	xlab=0.85*(xr[1]-xr[0])+xr[0]
	IF u.line.plot[0] THEN BEGIN
		oplot,*u.dat.br[4],*u.dat.br[0],color=lcol[0]
		xyouts,xlab,ylab,u.line.label[0],color=lcol[0],chars=0.8*ls
		ylab-=0.1*(yr[1]-yr[0])
	ENDIF
	IF u.line.plot[1] THEN BEGIN
		oplot,*u.dat.br[5],*u.dat.br[1],color=lcol[1]
		xyouts,xlab,ylab,u.line.label[1],color=lcol[1],chars=0.8*ls
		ylab-=0.1*(yr[1]-yr[0])
	ENDIF
	IF u.line.plot[2] THEN BEGIN
		oplot,*u.dat.br[6],*u.dat.br[2],color=lcol[2]
		xyouts,xlab,ylab,u.line.label[2],color=lcol[2],chars=0.8*ls
		ylab-=0.1*(yr[1]-yr[0])
	ENDIF
	IF u.line.plot[3] THEN BEGIN
		oplot,*u.dat.br[7],*u.dat.br[3],color=lcol[3]
		xyouts,xlab,ylab,u.line.label[3],color=lcol[3],chars=0.8*ls
		ylab-=0.1*(yr[1]-yr[0])
	ENDIF

	;dat={dalpha:dalpha,xtomo:xtomo,dens:dens,nbi:nbi,rf:rf}
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.asize[0],u.plot.asize[1],0,0,0]
	ENDIF ELSE device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm

	u.plot.x[0]=!X & u.plot.y[0]=!Y & u.plot.p[0]=!P ; save context
END

FUNCTION wspec_maxsbr,u
	xr=[u.plot.x0[1],u.plot.x1[1]]
	ymax=0.0
	IF u.stat.low AND u.plot.low THEN BEGIN
		index=ipt(u.stat.time,u.dat.low.time)	
		tmp=where(u.dat.low.lam GE xr[0] AND u.dat.low.lam LE xr[1])
		IF tmp[0] NE -1 THEN ymax=ymax > max(u.dat.low.specbr[tmp,index])
        ENDIF
	IF u.stat.xeus AND u.plot.xeus THEN BEGIN
		index=ipt(u.stat.time,u.dat.xeus.time)	
		tmp=where(u.dat.xeus.lam GE xr[0] AND u.dat.xeus.lam LE xr[1])
		IF tmp[0] NE -1 THEN ymax=ymax > max(u.dat.xeus.specbr[tmp,index])
        ENDIF
	IF u.stat.spreds AND u.plot.spreds THEN BEGIN
		index=ipt(u.stat.time,u.dat.spreds.time)	
		tmp=where(u.dat.spreds.lam GE xr[0] AND u.dat.spreds.lam LE xr[1])
		IF tmp[0] NE -1 AND index NE -1 THEN ymax=ymax > max(u.dat.spreds.specbr[tmp,index])
        ENDIF
	IF u.stat.spredl AND u.plot.spredl THEN BEGIN
		index=ipt(u.stat.time,u.dat.spredl.time)	
		tmp=where(u.dat.spredl.lam GE xr[0] AND u.dat.spredl.lam LE xr[1])
		IF tmp[0] NE -1 AND index NE -1 THEN ymax=ymax > max(u.dat.spredl.specbr[tmp,index])
        ENDIF
	RETURN,ymax
END

PRO wspec_plot_spec,u
	IF u.stat.ps THEN BEGIN
		xsize=7.0
		ysize=4.5
		ls=0.9
		col=u.plot.pscol 
		lcol=u.line.pscol
		tit=num2str(u.shot,1)+' t='+num2str(u.stat.time,dp=3)
		d_old=!d
		device,xsize=xsize,ysize=ysize,/in
	ENDIF ELSE BEGIN
		ls=1.0
		col=u.plot.col
		lcol=u.line.col
		widget_control,u.id.draw2,get_value=draw_win
		window,0,xsize=u.plot.bsize[0],ysize=u.plot.bsize[1],/pixmap
		tit=''
        ENDELSE
	xr=[u.plot.x0[1],u.plot.x1[1]]
	IF u.plot.log THEN BEGIN
	 ylog=1
	ENDIF
	IF u.plot.auto THEN BEGIN
		yr=[0,wspec_maxsbr(u)*1.05]
		u.plot.y0[5]=yr[0]
		u.plot.y1[5]=yr[1]
		widget_control, u.id.sbr0,set_value=num2str(u.plot.y0[5],dp=1)
		widget_control, u.id.sbr1,set_value=num2str(u.plot.y1[5],dp=2)
        ENDIF ELSE BEGIN 
		yr=[u.plot.y0[5],u.plot.y1[5]]
	ENDELSE

	plot,[0],[0],xr=xr,yr=yr,ylog=ylog,/xsty,/ysty,xtit='Wavelength [Ang]',ytit='Spectral Brightness [AU]',chars=1.2*ls


	IF u.stat.xeus AND u.plot.xeus THEN BEGIN
		index=ipt(u.stat.time,u.dat.xeus.time)
		IF index NE -1 THEN oplot,u.dat.xeus.lam,u.dat.xeus.specbr[*,index],col=col[0]
		u.dat.xeus.index=index
        ENDIF 

	IF u.stat.low AND u.plot.low THEN BEGIN
		index=ipt(u.stat.time,u.dat.low.time)
		IF index NE -1 THEN oplot,u.dat.low.lam,u.dat.low.specbr[*,index],col=col[0]
		u.dat.low.index=index
        ENDIF
 
	IF u.stat.spreds AND u.plot.spreds THEN BEGIN
		index=ipt(u.stat.time,u.dat.spreds.time)
		IF index NE -1 THEN oplot,u.dat.spreds.lam,u.dat.spreds.specbr[*,index],col=col[0]
		u.dat.spreds.index=index
        ENDIF 

	IF u.stat.spredl AND u.plot.spredl THEN BEGIN
		index=ipt(u.stat.time,u.dat.spredl.time)
		IF index NE -1 THEN oplot,u.dat.spredl.lam,u.dat.spredl.specbr[*,index],col=col[0]
		u.dat.spredl.index=index
        ENDIF 

	ylab=1.02*(yr[1]-yr[0])+yr[0]
	xlab=0.1*(xr[1]-xr[0])+xr[0]	
	FOR i=0,n(u.line.plot) DO BEGIN
		IF u.line.plot[i] THEN BEGIN
			CASE u.line.spec[i] OF
				'xeus' : BEGIN
					IF u.stat.xeus THEN BEGIN
						tmp=where(u.dat.xeus.lam GE u.line.lam0[i] AND u.dat.xeus.lam LE u.line.lam1[i])
						index=ipt(u.stat.time,u.dat.xeus.time)
						IF index NE -1 THEN oplot,u.dat.xeus.lam[tmp],u.dat.xeus.specbr[tmp,index],col=lcol[i]
					ENDIF
	                        END

				'loweus' : BEGIN
					IF u.stat.low THEN BEGIN
						tmp=where(u.dat.low.lam GE u.line.lam0[i] AND u.dat.low.lam LE u.line.lam1[i])
						index=ipt(u.stat.time,u.dat.low.time)
						IF index NE -1 THEN oplot,u.dat.low.lam[tmp],u.dat.low.specbr[tmp,index],col=lcol[i]
					ENDIF
                                END
				'spreds' : BEGIN
					IF u.stat.spreds THEN BEGIN
						tmp=where(u.dat.spreds.lam GE u.line.lam0[i] AND u.dat.spreds.lam LE u.line.lam1[i])
						index=ipt(u.stat.time,u.dat.spreds.time)
						IF index NE -1 THEN oplot,u.dat.spreds.lam[tmp],u.dat.spreds.specbr[tmp,index],col=lcol[i]
					ENDIF
				END
				'spredl' : BEGIN
					IF u.stat.spredl THEN BEGIN
						tmp=where(u.dat.spredl.lam GE u.line.lam0[i] AND u.dat.spredl.lam LE u.line.lam1[i])
						index=ipt(u.stat.time,u.dat.spredl.time)
						IF index NE -1 THEN oplot,u.dat.spredl.lam[tmp],u.dat.spredl.specbr[tmp,index],col=lcol[i]
					ENDIF
				END
				else  : ;print, u.line.spec[i]
                      	ENDCASE
			IF u.stat.ps THEN xyouts,xlab,ylab,u.line.label[i],color=lcol[i]
			xlab+=0.2*(xr[1]-xr[0])

               ENDIF
	ENDFOR
	IF u.plot.label THEN BEGIN
		widget_control,u.id.lelem,get_value=lelem
		eplot=strsplit(lelem,',',/extract)
		FOR order=1,u.plot.order DO BEGIN
			tmp=where( ((order*u.dat.line.lam) GT xr[0]) AND ((order*u.dat.line.lam) LT xr[1]) )
		 	IF tmp[0] NE -1 THEN BEGIN
				FOR i=0,n(tmp) DO BEGIN
					IF total(where(eplot EQ strtrim(u.dat.line.elem[tmp[i]]))) NE -1 AND u.dat.line.priority[tmp[i]] LE u.plot.plabel THEN BEGIN
						oplot, order*u.dat.line.lam[tmp[i]]*[1,1],yr,linestyle=1.0,col=col[order-1]
						xyouts,order*u.dat.line.lam[tmp[i]]-0.005*(xr[1]-xr[0]),yr[0]+(yr[1]-yr[0])*0.99,u.dat.line.elem[tmp[i]]+' '+u.dat.line.cs[tmp[i]],orient=90,align=1,col=col[order-1]
					ENDIF
				ENDFOR
			 ENDIF
		ENDFOR
	ENDIF
	IF u.stat.ps THEN xyouts,1.03*(xr[1]-xr[0])+xr[0],0.1*(yr[1]-yr[0])+yr[0],num2str(u.shot,1)+'     t='+num2str(u.stat.time,dp=3),chars=0.8*ls,orient=90

	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.bsize[0],u.plot.bsize[1],0,0,0]
	ENDIF ELSE device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm
	
	u.plot.x[1]=!X & u.plot.y[1]=!Y & u.plot.p[1]=!P ; save context

END



PRO wspec_plot,u
	IF u.stat.tline THEN wspec_plot_time,u
	wspec_plot_spec,u
END

PRO wspec_load_data,u
	widget_control,/hourglass
	shot=u.shot
	u.stat.tmax=0.0
	u.stat.lmin=1500.0
	u.stat.lmax=0.0

	;load d-alpha
	wspec_load_dalpha,shot,time,dalpha,status=status
	IF status THEN dalpha={time:time,dalpha:dalpha} ELSE dalpha=-1

	;load core xtomo brightness
	wspec_load_corextomo,shot,time,br,status=status
	IF status THEN xtomo={time:time,br:br*1.0e-3} ELSE xtomo=-1

	;load density
	wspec_load_dens,shot,time,dens,status=status
	IF status THEN dens={time:time,dens:dens} ELSE dens=-1

	;load NBI power
	wspec_load_nbipow,shot,time,nbipow,status=status
	IF status THEN nbi={time:time,pow:nbipow} ELSE nbi=-1
	nbit=time

	;load HHFW power
	wspec_load_hhfwpow,shot,time,rfpow,status=status
	IF status THEN rf={time:time,pow:rfpow}  ELSE rf=-1
	widget_control,u.id.message,set_value='SHOT '+num2str(shot,1)+': Time Trace Data Loaded',/append

	;load short wavelength SPRED
	IF u.stat.load.s THEN wspec_load_spred,shot,specbr,lam,time,/short,status=status ELSE status=0
	IF status THEN BEGIN
		spreds={specbr:specbr,lam:lam,time:time,index:ipt(time,u.stat.time)} 
		u.stat.spreds=1
		if max(time) gt u.stat.tmax then u.stat.tmax=max(time)
		IF min(lam) LT u.stat.lmin THEN u.stat.lmin=min(lam)
		IF max(lam) GT u.stat.lmax THEN u.stat.lmax=max(lam)
		widget_control,u.id.message,set_value='SHOT '+num2str(shot,1)+': SPRED (s) Loaded',/append
		;widget_control,u.id.spreds,set_button=1
        ENDIF ELSE BEGIN
		spreds=-1
		u.stat.spreds=0
		if max(time) gt u.stat.tmax then u.stat.tmax=max(time)
		IF u.stat.load.s THEN widget_control,u.id.message,set_value='SHOT '+num2str(shot,1)+': No SPRED (s) Data',/append
		widget_control,u.id.spreds,set_button=0
        ENDELSE

	;load long wavelength SPRED
	IF u.stat.load.s THEN wspec_load_spred,shot,specbr,lam,time,status=status ELSE status=0
	IF status THEN BEGIN
		spredl={specbr:specbr,lam:lam,time:time,index:ipt(time,u.stat.time)} 
		u.stat.spredl=1
		if max(time) gt u.stat.tmax then u.stat.tmax=max(time)
		IF min(lam) LT u.stat.lmin THEN u.stat.lmin=min(lam)
		IF max(lam) GT u.stat.lmax THEN u.stat.lmax=max(lam)
		widget_control,u.id.message,set_value='SHOT '+num2str(shot,1)+': SPRED (l) Loaded',/append
		;widget_control,u.id.spredl,set_button=1
        ENDIF ELSE BEGIN
		spredl=-1
		u.stat.spredl=0
		IF u.stat.load.s THEN widget_control,u.id.message,set_value='SHOT '+num2str(shot,1)+': No SPRED (long) Data',/append
		widget_control,u.id.spredl,set_button=0
        ENDELSE

	;load XEUS spectra
	IF u.stat.load.x THEN wspec_load_xeus,shot,specbr,lam,time,status=status ELSE status=0
	IF status THEN BEGIN
		xeus={specbr:specbr,lam:lam,time:time,index:ipt(time,u.stat.time)} 
		u.stat.xeus=1
		if max(time) gt u.stat.tmax then u.stat.tmax=max(time)
		IF min(lam) LT u.stat.lmin THEN u.stat.lmin=min(lam)
		IF max(lam) GT u.stat.lmax THEN u.stat.lmax=max(lam)
		widget_control,u.id.message,set_value='SHOT '+num2str(shot,1)+': XEUS Loaded',/append
		widget_control,u.id.xeus,set_button=1
        ENDIF ELSE BEGIN
		xeus=-1
		u.stat.xeus=0
		IF u.stat.load.s THEN widget_control,u.id.message,set_value='SHOT '+num2str(shot,1)+': No XEUS Data',/append
		widget_control,u.id.xeus,set_button=0
	ENDELSE

	;load LoWEUS spectra
	IF u.stat.load.l THEN wspec_load_loweus,shot,specbr,lam,time,status=status ELSE status=0
	IF status THEN BEGIN
		low={specbr:specbr,lam:lam,time:time,index:ipt(time,u.stat.time)} 
		u.stat.low=1
		if max(time) gt u.stat.tmax then u.stat.tmax=max(time)
		IF min(lam) LT u.stat.lmin THEN u.stat.lmin=min(lam)
		IF max(lam) GT u.stat.lmax THEN u.stat.lmax=max(lam)
		widget_control,u.id.message,set_value='SHOT '+num2str(shot,1)+': LoWEUS Loaded',/append
		widget_control,u.id.low,set_button=1
        ENDIF ELSE BEGIN
		low=-1
		u.stat.low=0
		IF u.stat.load.l THEN widget_control,u.id.message,set_value='SHOT '+num2str(shot,1)+': No LoWEUS Data',/append
		widget_control,u.id.low,set_button=0
        ENDELSE

	;print,'u.stat.tmax='+num2str(u.stat.tmax)
	widget_control,u.id.t_slider,set_slider_max=u.stat.tmax*1000

	;initialize brightness ptrarr
	nlines=n(u.line.label)+1
	IF u.stat.dat THEN BEGIN
		heap_free,u.dat.br 		;release br
		heap_free,u.dat.sbr		;release sbr
	ENDIF
	br=ptrarr(nlines*2, /allocate_heap)
	sbr=ptrarr(nlines*2,/allocate_heap)

	;load LINE table
	wspec_load_linestr,line
	;widget_control,u.id.lelem,set_value=line.list

	u.stat.dat=1 
	dat={low:low,xeus:xeus,spredl:spredl,spreds:spreds,dalpha:dalpha,xtomo:xtomo,dens:dens,nbi:nbi,rf:rf,line:line,br:br,sbr:sbr}
	u={id:u.id,shot:u.shot,stat:u.stat,plot:u.plot,line:u.line,dat:dat}
	widget_control,u.id.base, set_uvalue=u

	FOR i=0,n(u.line.plot) DO IF u.line.plot[i] THEN wspec_calc_br,u,i
	wspec_plot_time,u
	wspec_plot_spec,u
END

PRO wspec_reset_lzoom,u
	widget_control,u.id.lowz,set_button=0
	widget_control,u.id.midz,set_button=0
	widget_control,u.id.highz,set_button=0
	widget_control,u.id.outz,set_button=0
END

PRO w_spec_event,event
	widget_control,event.top,get_uvalue=u
	id = u.id
	tag = tag_names(event,/st)
	button=' '
	idtags=tag_names(id)
	FOR i=0,n(idtags) DO IF id.(i) EQ event.id THEN ename=idtags[i]
	CASE tag OF
		"WIDGET_BASE" : BEGIN

		END
		"WIDGET_BUTTON": BEGIN
			widget_control,event.id,get_value=button,get_uvalue=uvalue
			IF NOT keyword_set(uvalue) THEN uvalue='none'
			CASE ename OF
				"LOAD" : wspec_load_data,u
				"QUIT" : BEGIN
					IF u.stat.dat THEN BEGIN
						heap_free,u.dat.br
						heap_free,u.dat.sbr
					ENDIF
					widget_control,event.top,/destroy
					!except=1
					heap_gc
				END
				"PRINT" : BEGIN
					psplot
					u.stat.ps=1
					wspec_plot_time,u
					wspec_plot_spec,u
					psc
					u.stat.ps=0
					xwplot
					set_plot,'x'	
				END
				"STOP" : BEGIN
					stop
                                END
				"SPREDS" : BEGIN
					IF u.stat.dat THEN BEGIN
						IF event.select EQ 1 AND u.stat.spreds THEN u.plot.spreds=1 ELSE u.plot.spreds=0
						wspec_plot_spec,u
                                        ENDIF
                                END
				"SPREDL" : BEGIN
					IF u.stat.dat THEN BEGIN
						IF event.select EQ 1 AND u.stat.spredl THEN u.plot.spredl=1 ELSE u.plot.spredl=0
						wspec_plot_spec,u
                                        ENDIF
                                END
				"XEUS" : BEGIN
					IF u.stat.dat THEN BEGIN
						IF event.select EQ 1 AND u.stat.xeus THEN u.plot.xeus=1 ELSE u.plot.xeus=0
						wspec_plot_spec,u
                                        ENDIF
                                END
				"LOW" : BEGIN
					IF u.stat.dat THEN BEGIN
						IF event.select EQ 1 AND u.stat.low THEN u.plot.low=1 ELSE u.plot.low=0
						wspec_plot_spec,u
                                        ENDIF
                                END
        			"SETCURR": BEGIN
            				u.shot = mdsvalue('current_shot("cmod")')
            				widget_control,u.id.shotid,set_value=num2str(u.shot,1)
					wspec_load_data,u
                                END
				'LOWZ' : BEGIN
					wspec_reset_lzoom,u
					u.plot.x0[1]=10.0
					u.plot.x1[1]=65.0
					widget_control, u.id.lam0,set_value=num2str(u.plot.x0[1],dp=1)
					widget_control, u.id.lam1,set_value=num2str(u.plot.x1[1],dp=1)
					wspec_plot_spec,u
                                END
				'MIDZ' : BEGIN
					wspec_reset_lzoom,u
					u.plot.x0[1]=140.0
					u.plot.x1[1]=270.0
					widget_control, u.id.lam0,set_value=num2str(u.plot.x0[1],dp=1)
					widget_control, u.id.lam1,set_value=num2str(u.plot.x1[1],dp=1)
					wspec_plot_spec,u
                                END
				'HIGHZ' : BEGIN
					wspec_reset_lzoom,u
					u.plot.x0[1]=110.0
					u.plot.x1[1]=135.0
					widget_control, u.id.lam0,set_value=num2str(u.plot.x0[1],dp=1)
					widget_control, u.id.lam1,set_value=num2str(u.plot.x1[1],dp=1)
					wspec_plot_spec,u
                                END
				'OUTZ' : BEGIN
					wspec_reset_lzoom,u
					u.plot.x0[1]=u.stat.lmin
					u.plot.x1[1]=u.stat.lmax
					widget_control, u.id.lam0,set_value=num2str(u.plot.x0[1],dp=1)
					widget_control, u.id.lam1,set_value=num2str(u.plot.x1[1],dp=1)
					wspec_plot_spec,u
                                END
				'SPRED' : BEGIN
					wspec_reset_lzoom,u
					u.plot.x0[1]=100
					u.plot.x1[1]=1200
					widget_control, u.id.lam0,set_value=num2str(u.plot.x0[1],dp=1)
					widget_control, u.id.lam1,set_value=num2str(u.plot.x1[1],dp=1)
					wspec_plot_spec,u
                                END
				'LSPRED' : IF event.select EQ 1 THEN u.stat.load.s=1 ELSE u.stat.load.s=0
				'LLOW' : IF event.select EQ 1 THEN u.stat.load.l=1 ELSE u.stat.load.l=0
				'LXEUS' : IF event.select EQ 1 THEN u.stat.load.x=1 ELSE u.stat.load.x=0
				'APLOT' : BEGIN
					IF event.select EQ 1 THEN u.line.plot[0]=1 ELSE u.line.plot[0]=0
					wspec_calc_br,u,0
					wspec_plot_time,u
                                 END
				'BPLOT' : BEGIN
					IF event.select EQ 1 THEN u.line.plot[1]=1 ELSE u.line.plot[1]=0
					wspec_calc_br,u,1
					wspec_plot_time,u
                                 END
				'CPLOT' : BEGIN
					IF event.select EQ 1 THEN u.line.plot[2]=1 ELSE u.line.plot[2]=0
					wspec_calc_br,u,2
					wspec_plot_time,u
                                 END
				'DPLOT' : BEGIN
					IF event.select EQ 1 THEN u.line.plot[3]=1 ELSE u.line.plot[3]=0
					wspec_calc_br,u,3
					wspec_plot_time,u
                                 END
				'AELEM' : BEGIN
					IF event.select EQ 1 THEN BEGIN
						IF u.stat.dat THEN BEGIN
							widget_control,u.id.lelem,set_value=u.dat.line.list
							wspec_plot_spec,u
						ENDIF
						
						widget_control, u.id.aelem,set_button=0
					ENDIF
				END
				'OELEM' : BEGIN
					IF event.select EQ 1 THEN u.plot.order=2 ELSE u.plot.order=1
					wspec_plot_spec,u
				END
				'LOGSBR' : BEGIN
					IF event.select EQ 1 THEN u.plot.log=1 ELSE u.plot.log=0
					wspec_plot_spec,u
				END
				'AUTOSBR' : BEGIN
					IF event.select EQ 1 THEN u.plot.auto=1 ELSE u.plot.auto=0
					wspec_plot_spec,u
				END
				'LABEL' : BEGIN
					IF event.select EQ 1 THEN u.plot.label=1 ELSE u.plot.label=0
					wspec_plot_spec,u
				END
				'TLINE' : BEGIN
					IF event.select EQ 1 THEN u.stat.tline=1 ELSE u.stat.tline=0
					wspec_plot_time,u
                                END
				'LZOOM' : BEGIN
					wspec_plot_spec,u
					cursor,x1,y1
					wait,0.1
					cursor,x2,y2
					u.plot.x0[1]=x1
					u.plot.x1[1]=x2
					widget_control, u.id.lam0,set_value=num2str(u.plot.x0[1],dp=1)
					widget_control, u.id.lam1,set_value=num2str(u.plot.x1[1],dp=1)
					wspec_plot_spec,u
                                 END
				'TZOOM' : BEGIN
					wspec_plot_time,u
					cursor,x1,y1
					wait,0.1
					cursor,x2,y2
					u.plot.x0[0]=x1
					u.plot.x1[0]=x2
					widget_control, u.id.tau0,set_value=num2str(u.plot.x0[0],dp=1)
					widget_control, u.id.tau1,set_value=num2str(u.plot.x1[0],dp=1)
					wspec_plot_time,u
				END
				'AFIND' : BEGIN
					wspec_plot_spec,u
					widget_control,u.id.message,set_value='Select Lower Wavelength',/append
					cursor,x1,y1
					wait,0.15
					widget_control,u.id.message,set_value='Select Upper Wavelength',/append
					cursor,x2,y2
					u.line.lam0[0]=x1
					u.line.lam1[0]=x2
					widget_control, u.id.alam0,set_value=num2str(u.line.lam0[0],dp=2)
					widget_control, u.id.alam1,set_value=num2str(u.line.lam1[0],dp=2)
					wspec_plot_time,u
					wspec_calc_br,u,0
					u.line.plot[0]=1
					widget_control,u.id.aplot,set_button=1
					wspec_plot_time,u
					wspec_plot_spec,u
                                END
				'BFIND' : BEGIN
					wspec_plot_spec,u
					widget_control,u.id.message,set_value='Select Lower Wavelength',/append
					cursor,x1,y1
					wait,0.15
					widget_control,u.id.message,set_value='Select Upper Wavelength',/append
					cursor,x2,y2
					u.line.lam0[1]=x1
					u.line.lam1[1]=x2
					widget_control, u.id.blam0,set_value=num2str(u.line.lam0[1],dp=2)
					widget_control, u.id.blam1,set_value=num2str(u.line.lam1[1],dp=2)
					wspec_plot_time,u
					wspec_calc_br,u,1
					u.line.plot[1]=1
					widget_control,u.id.bplot,set_button=1
					wspec_plot_time,u
					wspec_plot_spec,u
				END
				'CFIND' : BEGIN
					wspec_plot_spec,u
					widget_control,u.id.message,set_value='Select Lower Wavelength',/append
					cursor,x1,y1
					wait,0.15
					widget_control,u.id.message,set_value='Select Upper Wavelength',/append
					cursor,x2,y2
					u.line.lam0[2]=x1
					u.line.lam1[2]=x2
					widget_control, u.id.clam0,set_value=num2str(u.line.lam0[2],dp=2)
					widget_control, u.id.clam1,set_value=num2str(u.line.lam1[2],dp=2)
					wspec_plot_time,u
					wspec_calc_br,u,2
					u.line.plot[2]=1
					widget_control,u.id.cplot,set_button=1
					wspec_plot_time,u
					wspec_plot_spec,u
				END
				'DFIND' : BEGIN
					wspec_plot_spec,u
					widget_control,u.id.message,set_value='Select Lower Wavelength',/append
					cursor,x1,y1
					wait,0.15
					widget_control,u.id.message,set_value='Select Upper Wavelength',/append
					cursor,x2,y2
					u.line.lam0[3]=x1
					u.line.lam1[3]=x2
					widget_control, u.id.dlam0,set_value=num2str(u.line.lam0[3],dp=2)
					widget_control, u.id.dlam1,set_value=num2str(u.line.lam1[3],dp=2)
					wspec_plot_time,u
					wspec_calc_br,u,3
					u.line.plot[3]=1
					widget_control,u.id.dplot,set_button=1
					wspec_plot_time,u
					wspec_plot_spec,u
				END
				ELSE:
			ENDCASE
		END
		"WIDGET_SLIDER": BEGIN
			widget_control,event.id,get_value=slider
			IF NOT keyword_set(uvalue) THEN uvalue='none'
			CASE ename OF 
				'T_SLIDER' : BEGIN
					IF u.stat.dat THEN BEGIN
						u.stat.time=slider/1.0e3
						widget_control, u.id.t_text,set_value=num2str(u.stat.time,dp=3)
						wspec_plot,u
					ENDIF
				END
				ELSE:
			ENDCASE
		END
		'WIDGET_DRAW' : begin ;process a button event in the draw window
		 case event.id of
		  u.id.draw1 : begin ; in time pane
		   !X=u.plot.x[0] & !Y=u.plot.y[0] & !P=u.plot.p[0] ; restore context
		   if (event.press ne 0) or (event.release ne 0) then begin ; ie not a motion event
	       	    if event.press eq 1 then begin ; LMB click,click - zoom time L,R
		    ;also event.release
		    	click_loc = convert_coord(event.x, event.y, /device, /to_data)
					if u.plot.zoom[0] eq 0 then begin
					 u.plot.zoom[0]=1
					 u.plot.x0[0]=click_loc[0]
					 widget_control, u.id.tau0,set_value=num2str(u.plot.x0[0],dp=3)
					endif else begin
					 u.plot.zoom[0]=0
					 if click_loc[0] lt u.plot.x0[0] then begin
					  u.plot.x1[0]=u.plot.x0[0]
					  u.plot.x0[0]=click_loc[0]
					  widget_control, u.id.tau0,set_value=num2str(u.plot.x0[0],dp=3)
					 endif else u.plot.x1[0]=click_loc[0]
					 widget_control, u.id.tau1,set_value=num2str(u.plot.x1[0],dp=3)
					 wspec_plot_time,u
					endelse
		    endif
		    if event.press eq 2 then begin ; MMB click - set spectrum time
		    	click_loc = convert_coord(event.x, event.y, /device, /to_data)
					u.stat.time=click_loc[0]
					widget_control, u.id.t_slider,set_value=u.stat.time*1000
					widget_control, u.id.t_text,set_value=num2str(u.stat.time,dp=3)
					wspec_plot,u
		    endif
		    if event.press eq 4 then begin ; RMB click - zoom out (ie zoom 0,2)
		    			 u.plot.zoom[0]=0
					 u.plot.x0[0]=0.0
					 u.plot.x1[0]=u.stat.tmax
					 widget_control, u.id.tau0,set_value=num2str(u.plot.x0[0],dp=3)
					 widget_control, u.id.tau1,set_value=num2str(u.plot.x1[0],dp=3)
					 wspec_plot_time,u
		    endif
		   endif else begin ; pointer motion in pane
	  		 ptr_loc = convert_coord(event.x, event.y, /device, /to_data)
			 widget_control, u.id.tim_ind,set_value='t='+num2str(ptr_loc[0],dp=3)
		   endelse
		  end
		  u.id.draw2 : begin ; in wl pane
		   !X=u.plot.x[1] & !Y=u.plot.y[1] & !P=u.plot.p[1] ; restore context
		   if (event.press ne 0) or (event.release ne 0) then begin ; ie not a motion event
	       	    if event.press eq 1 then begin ; LMB click,click - zoom time L,R
		    	;also event.release
		    	click_loc = convert_coord(event.x, event.y, /device, /to_data)
					if u.plot.zoom[1] eq 0 then begin
					 u.plot.zoom[1]=1
					 u.plot.x0[1]=click_loc[0]
					 widget_control, u.id.lam0,set_value=num2str(u.plot.x0[1],dp=3)
					endif else begin
					 u.plot.zoom[1]=0
					 if click_loc[0] lt u.plot.x0[1] then begin
					  u.plot.x1[1]=u.plot.x0[1]
					  u.plot.x0[1]=click_loc[0]
					  widget_control, u.id.lam0,set_value=num2str(u.plot.x0[1],dp=3)
					 endif else u.plot.x1[1]=click_loc[0]
					 widget_control, u.id.lam1,set_value=num2str(u.plot.x1[1],dp=3)
					 wspec_plot_spec,u
					endelse
		    endif
		    if event.press eq 4 then begin ; RMB click - zoom out (ie zoom 1,300)
		    			 u.plot.zoom[1]=0
					 u.plot.x0[1]=u.stat.lmin
					 u.plot.x1[1]=u.stat.lmax
					 widget_control, u.id.lam0,set_value=num2str(u.plot.x0[1],dp=1)
					 widget_control, u.id.lam1,set_value=num2str(u.plot.x1[1],dp=1)
					 wspec_plot_spec,u
		    endif
		   endif else begin ; pointer motion in pane
	  		 ptr_loc = convert_coord(event.x, event.y, /device, /to_data)
			 widget_control, u.id.br_ind,set_value='t='+num2str(ptr_loc[0],dp=3)+' y='+num2str(ptr_loc[1],dp=3)
		   endelse
		  end
		 endcase
		end
   		"WIDGET_TEXT_CH": BEGIN
			widget_control,event.id,get_value=data
			CASE event.id OF
				u.id.shotid : BEGIN
					u.shot=long(data)
					wspec_load_data,u
				END
				u.id.lam0 : BEGIN
					u.plot.x0[1]=float(data)
					wspec_plot_spec,u
				END
				u.id.lam1 : BEGIN
					u.plot.x1[1]=float(data)
					wspec_plot_spec,u
				END			
				u.id.sbr0 : BEGIN
					u.plot.y0[5]=float(data)
					wspec_plot_spec,u
                                END
				u.id.sbr1 : BEGIN
					u.plot.y1[5]=float(data)
					wspec_plot_spec,u
				END
				u.id.tau0 : BEGIN
					u.plot.x0[0]=float(data)
					wspec_plot_time,u
				END
				u.id.tau1 : BEGIN
					u.plot.x1[0]=float(data)
					wspec_plot_time,u
				END	
				u.id.lelem : wspec_plot_spec,u
				u.id.pelem : BEGIN
					u.plot.plabel=int(data)
					wspec_plot_spec,u
				END
				u.id.alab : BEGIN
					u.line.label[0]=data
					wspec_plot_time,u
				END
				u.id.blab : BEGIN
					u.line.label[1]=data
					wspec_plot_time,u
				END
				u.id.clab : BEGIN
					u.line.label[2]=data
					wspec_plot_time,u
				END
				u.id.dlab : BEGIN
					u.line.label[3]=data
					wspec_plot_time,u
                                END
				u.id.t_text : BEGIN
					u.stat.time=float(data)
					widget_control,u.id.t_slider,set_value=u.stat.time*1.0e3
					wspec_plot_time,u
					wspec_plot_spec,u
				END
			ELSE: 
			ENDCASE
		END
		ELSE:
	ENDCASE
	IF button NE 'QUIT' THEN widget_control,event.top,set_uvalue=u		
END

;+
;NAME:
;	W_SPEC
;
;PURPOSE:
;	This procedure launches the W_SPEC widget which is used to
;	analyze the SXR/VUV spectroscopy data from HIREXSR, the
;	McPherson, XEUS and LoWEUS
;
;CALLING SEQUENCE:
;	@run_wspec.pro 		will compile the proper files, set the color
;				tables and launch the widget
;PROCEUDRE:
;	See the C-Mod Wiki page on the W_SPEC widget for instructions
;	(listed under Shot Analysis Tools)
;
;RESTRICTIONS:
;	Due to the limited namespace, many of the GENIE widgets have
;	overlapping function/procedure names so these widgets may not run
;	well in contiguous IDL sessions.
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 1/2012
;
;-

PRO w_spec,shot=shot,time=time
	IF NOT keyword_set(shot) THEN shot=141229
	IF NOT keyword_set(time) THEN time=0.5
	user=logname()
	mdim=get_screen_size()
	IF mdim[0] NE 1600 AND mdim[1] NE 1200 THEN base=widget_base(title='VUV/SXR SPECTROSCOPY',/row,tlb_size_events=1,/scroll,$
		x_scroll_size=mdim[0]*0.95,y_scroll_size=mdim[1]*0.85) ELSE base=widget_base(title='VUV/SXR SPECTROSCOPY',/row,tlb_size_events=1)

	A=widget_base(base,/column)
	B=widget_base(base,/column)
	C=widget_base(base,/column)

	ysz=800

	A0=widget_base(A,/row)
	dum=widget_label(A0,value='TIME EVOLUTION')
	tim_ind=widget_label(A0,value='                        ')
	A1=widget_base(A,frame=2)
	draw1=widget_draw(A1,xsize=400,ysize=ysz,retain=2, /motion_events, /button_events)

	B0=widget_base(B,/row)
	dum=widget_label(B0,value='SPECTRAL BRIGHTNESS')
	br_ind=widget_label(B0,value='                        ')
	B1=widget_base(B,frame=2)
	draw2=widget_draw(B1,xsize=ysz,ysize=ysz,retain=2, /motion_events, /button_events)
	B2=widget_base(B,/row)
	t_text=widget_text(B2,xsize=8,ysize=1,/edit)
	t_slider=widget_slider(B2,xsize=ysz-80,min=0,max=2000,value=time*1000,/drag,/suppress)

	dum=widget_label(C,value='SETUP')
	C1=widget_base(C,/column,xsize=320,ysize=2.3*ysz/3,/frame)
	C2=widget_base(C,/column,xsize=320,ysize=0.8*ysz/3,/frame)
	C1p1=widget_base(C1,/row)
	dum = widget_label(C1p1,value='SHOT: ')
	shotid = widget_text(C1p1,xsize=10,ysize=1,/edit)
	dum = widget_label(C1p1,value='')
	load= widget_button(C1p1,value='LOAD')
	dum = widget_label(C1p1,value='')
	stop= widget_button(C1p1,value='STOP')
	dum = widget_label(C1p1,value='')
	print= widget_button(C1p1,value='PRINT')
	dum = widget_label(C1p1,value='')
	quit= widget_button(C1p1,value='QUIT')
	C1p1x=widget_base(C1,/row)
	setcurr= widget_button(C1p1x,value='GET CURRENT SHOT')
	dum=widget_label(C1p1x,value=' LOAD:')
	C1p1y=widget_base(C1p1x,/row,/nonexcl)
	lspred=widget_button(C1p1y,value=' S')
	llow=widget_button(C1p1y,value=' L')
	lxeus=widget_button(C1p1y,value=' X')

	C1p2=widget_base(C1,/row)
	message=widget_text(C1p2,xsize=45,ysize=6,/scroll)
	C1p3=widget_base(C1,/row,/nonexcl)
	xeus=widget_button(C1p3,value=' XEUS ')
	low=widget_button(C1p3,value=' LOWEUS ')
	spreds=widget_button(C1p3, value=' SPRED (s)')
	spredl=widget_button(C1p3, value=' SPRED (l)')

	space=widget_base(C1,/row,ysize=5)
	C1p5=widget_base(C1,/row)
	dum=widget_label(C1p5,value=' LAM RANGE: ')
 	lam0=widget_text(C1p5,xsize=6,ysize=1,/edit)
	dum=widget_label(C1p5,value=' to ')
 	lam1=widget_text(C1p5,xsize=6,ysize=1,/edit)
	dum=widget_label(C1p5,value='[Ang] ')
	lzoom=widget_button(C1p5,value='ZOOM')
	C1p5a=widget_base(C1,/row,/nonexcl)
	outz=widget_button(C1p5a, value='ZOOM OUT')
	lowz=widget_button(C1p5a, value='LOW-Z He/H-LIKE (B -> Ne)')
	C1p5b=widget_base(C1,/row,/nonexcl)	
	midz=widget_button(C1p5b, value='MID-Z Li-LIKE (Ti -> Cu)')
	highz=widget_button(C1p5b, value='HIGH-Z Mo, W')
	C1p5c=widget_base(C1,/row,/nonexcl)	
	spred=widget_button(C1p5c, value='SPRED')

	C1p6=widget_base(C1,/row)
	dum=widget_label(C1p6,value=' SPECBR RANGE: ')
 	sbr0=widget_text(C1p6,xsize=5,ysize=1,/edit)
	dum=widget_label(C1p6,value=' to ')
 	sbr1=widget_text(C1p6,xsize=5,ysize=1,/edit)
	dum=widget_label(C1p6,value='[ph/s/m^2]')
	C1p8=widget_base(C1,/row)
	dum=widget_label(C1p8,value='')
	C1p8x=widget_base(C1p8,/row,/nonexcl)
	logsbr=widget_button(C1p8x, value='LOG Y')
	autosbr=widget_button(C1p8x, value='AUTOSCALE Y')
	label=widget_button(C1p8x,value='OPLOT LINES')
	C1p8y=widget_base(C1,/row)
	dum=widget_label(C1p8y,value='ELEM: ')
	lelem=widget_text(C1p8y,xsize=15,ysize=1,/edit)
	dum=widget_label(C1p8y,value='P: ')
	pelem=widget_text(C1p8y,xsize=3,ysize=1,/edit)
	C1p8yy=widget_base(C1p8y,/row,/nonexcl)
	aelem=widget_button(C1p8yy, value=' ALL')
	oelem=widget_button(C1p8yy, value='2nd')
	C1p7=widget_base(C1,/row)
	dum=widget_label(C1p7,value=' TIME RANGE: ')
 	tau0=widget_text(C1p7,xsize=5,ysize=1,/edit)
	dum=widget_label(C1p7,value=' to ')
 	tau1=widget_text(C1p7,xsize=5,ysize=1,/edit)
	dum=widget_label(C1p7,value='[sec]')
	tzoom=widget_button(C1p7,value='ZOOM')
	C1p8=widget_base(C1,/row,/nonexcl)
	tline=widget_button(C1p8,value=' OPLOT INDICATOR ON TIME TRACE ')	
	

	dum=widget_label(C2,value='LINE BRIGHTNESS')
	C2p1=widget_base(C2,/row)
	dum=widget_label(C2p1,value='     LINE        LAM1          LAM2           PLOT')
	C2p2=widget_base(C2,/row)
	dum=widget_label(C2p2,value='A')
	alab=widget_text(C2p2,xsize=8,ysize=1,/edit)
	dum=widget_label(C2p2,value=' ')
	alam0=widget_text(C2p2,xsize=5,ysize=1,/edit)
	dum=widget_label(C2p2,value=' to ')
	alam1=widget_text(C2p2,xsize=5,ysize=1,/edit)
	afind=widget_button(C2p2,value=' FIND ')
	C2p2x=widget_base(C2p2,/row,/nonexcl)
	aplot=widget_button(C2p2x,value=' ')
	C2p3=widget_base(C2,/row)
	dum=widget_label(C2p3,value='B')
	blab=widget_text(C2p3,xsize=8,ysize=1,/edit)
	dum=widget_label(C2p3,value=' ')
	blam0=widget_text(C2p3,xsize=5,ysize=1,/edit)
	dum=widget_label(C2p3,value=' to ')
	blam1=widget_text(C2p3,xsize=5,ysize=1,/edit)
	bfind=widget_button(C2p3,value=' FIND ')
	C2p3x=widget_base(C2p3,/row,/nonexcl)
	bplot=widget_button(C2p3x,value=' ')
	C2p4=widget_base(C2,/row)
	dum=widget_label(C2p4,value='C')
	clab=widget_text(C2p4,xsize=8,ysize=1,/edit)
	dum=widget_label(C2p4,value=' ')
	clam0=widget_text(C2p4,xsize=5,ysize=1,/edit)
	dum=widget_label(C2p4,value=' to ')
	clam1=widget_text(C2p4,xsize=5,ysize=1,/edit)
	cfind=widget_button(C2p4,value=' FIND ')
	C2p4x=widget_base(C2p4,/row,/nonexcl)
	cplot=widget_button(C2p4x,value=' ')
	C2p5=widget_base(C2,/row)
	dum=widget_label(C2p5,value='D')
	dlab=widget_text(C2p5,xsize=8,ysize=1,/edit)
	dum=widget_label(C2p5,value=' ')
	dlam0=widget_text(C2p5,xsize=5,ysize=1,/edit)
	dum=widget_label(C2p5,value=' to ')
	dlam1=widget_text(C2p5,xsize=5,ysize=1,/edit)
	dfind=widget_button(C2p5,value=' FIND ')
	C2p5x=widget_base(C2p5,/row,/nonexcl)
	dplot=widget_button(C2p5x,value=' ')

	
	id={base:base,draw1:draw1,draw2:draw2,t_slider:t_slider,t_text:t_text,$
		shotid:shotid,load:load,quit:quit,print:print,stop:stop,message:message,setcurr:setcurr,$
		br_ind:br_ind,tim_ind:tim_ind,$
		lspred:lspred,llow:llow,lxeus:lxeus,$
		xeus:xeus,low:low,spreds:spreds,spredl:spredl,$
		lam0:lam0,lam1:lam1,sbr0:sbr0, sbr1:sbr1,tau0:tau0, tau1:tau1,$
		logsbr:logsbr,autosbr:autosbr,label:label,lelem:lelem,pelem:pelem,aelem:aelem,oelem:oelem,$
		lzoom:lzoom,outz:outz,lowz:lowz,midz:midz,highz:highz,spred:spred,$
		tline:tline,tzoom:tzoom,$
		alab:alab, alam0:alam0,alam1:alam1,afind:afind,aplot:aplot,$
		blab:blab, blam0:blam0,blam1:blam1,bfind:bfind,bplot:bplot,$
		clab:clab, clam0:clam0,clam1:clam1,cfind:cfind,cplot:cplot,$
		dlab:dlab, dlam0:dlam0,dlam1:dlam1,dfind:dfind,dplot:dplot}
	load={s:0,l:1,x:1}
	stat={time:time,tmax:2.0,lmin:0.0,lmax:0.0,dat:0,ps:0,col:!p.color,pscol:0,tline:1,load:load,xeus:0,low:0,spreds:0,spredl:0}
	plot={	asize:[400,ysz],$
		bsize:[ysz,ysz],$
		y0:float([0,0,  0,0,0,0]),$ ;dalpha, sxr, dens, pow, br, specbr
		y1:float([1,100,2,1,1,1]),$
		x0:[0.0,1.0],$ ;time, wavelength
		x1:[2.0,300.0],$
		zoom:[0,0],$
		xeus:1,	low:1,	spreds:0, spredl:0,$
		nt:9,	nm:2,	auto:1,	label:0,plabel:100,order:1,ntau:1000,log:0,$
		col:[!p.color,196,96],pscol:[0,0,0],$
		x:[!X,!X],y:[!Y,!Y],p:[!P,!P]$
		}
	line={	label:['C VI',	'N VII',	'O VIII',	'Fe XXIII'],$
		lam0: [33.49,		24.61,		18.70,		132.4],$
		lam1: [34.07,		25.08,		19.08,		133.5],$
		plot:[1,1,1,1],$
		col:[50,90,130,180],$
		pscol:[30,100,150,200],$
		spec:['xeus','xeus','xeus','loweus']$
		}

	u={id:id,shot:shot,stat:stat,plot:plot,line:line}
	widget_control,base,set_uvalue=u
	widget_control,u.id.shotid,set_value=num2str(u.shot,1)
	widget_control, u.id.lam0,set_value=num2str(u.plot.x0[1],dp=1)
	widget_control, u.id.lam1,set_value=num2str(u.plot.x1[1],dp=1)
	widget_control, u.id.tau0,set_value=num2str(u.plot.x0[0],dp=1)
	widget_control, u.id.tau1,set_value=num2str(u.plot.x1[0],dp=1)
	widget_control, u.id.sbr0,set_value=num2str(u.plot.y0[5],dp=1)
	widget_control, u.id.sbr1,set_value=num2str(u.plot.y1[5],dp=2)
	widget_control, u.id.tline,set_button=u.stat.tline
	widget_control, u.id.t_text,set_value=num2str(u.stat.time,dp=3)
	widget_control, u.id.autosbr,set_button=u.plot.auto
	widget_control, u.id.lspred,set_button=u.stat.load.s
	widget_control, u.id.llow,set_button=u.stat.load.l
	widget_control, u.id.lxeus,set_button=u.stat.load.x
	widget_control, u.id.lelem,set_value='B,C,O,Fe'
	widget_control, u.id.pelem,set_value=num2str(u.plot.plabel,1)
	widget_control, u.id.alab,set_value=u.line.label[0]
	widget_control, u.id.alam0,set_value=num2str(u.line.lam0[0],dp=2)
	widget_control, u.id.alam1,set_value=num2str(u.line.lam1[0],dp=2)
	widget_control, u.id.aplot,set_button=u.line.plot[0]
	widget_control, u.id.blab,set_value=u.line.label[1]
	widget_control, u.id.blam0,set_value=num2str(u.line.lam0[1],dp=2)
	widget_control, u.id.blam1,set_value=num2str(u.line.lam1[1],dp=2)
	widget_control, u.id.bplot,set_button=u.line.plot[1]
	widget_control, u.id.clab,set_value=u.line.label[2]
	widget_control, u.id.clam0,set_value=num2str(u.line.lam0[2],dp=2)
	widget_control, u.id.clam1,set_value=num2str(u.line.lam1[2],dp=2)
	widget_control, u.id.cplot,set_button=u.line.plot[2]
	widget_control, u.id.dlab,set_value=u.line.label[3]
	widget_control, u.id.dlam0,set_value=num2str(u.line.lam0[3],dp=2)
	widget_control, u.id.dlam1,set_value=num2str(u.line.lam1[3],dp=2)
	widget_control, u.id.dplot,set_button=u.line.plot[3]

	!except=0
	widget_control,base,/realize
	;wspec_load_data,u
	xmanager,'w_spec',base, event_handler='w_spec_event'

END


;Main program
w_spec
END

