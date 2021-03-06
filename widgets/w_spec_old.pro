PRO wspec_load_mcp,shot,specbr,lam,time,status=status
	mdsopen,'spectroscopy',shot
	specbr=mdsvalue('_sig=\SPECTROSCOPY::TOP.VUV.ANALYSIS:MCP_BRI_MOD',/quiet,status=status)
	lam=mdsvalue('dim_of(_sig,1)',/quiet)
	time=mdsvalue('dim_of(_sig,0)',/quiet)
	mdsclose,'spectroscopy',shot
	specbr/=10.0
END

PRO wspec_load_vuv,shot,specbr,lam,time,xeus=xeus,loweus=loweus,status=status
	IF keyword_set(xeus) THEN path='\SPECTROSCOPY::TOP.XEUS'
	IF keyword_set(loweus) THEN path='\SPECTROSCOPY::TOP.LOWEUS'	
	mdsopen,'spectroscopy',shot
	specbr=mdsvalue('_sig='+path+':SPEC',/quiet,status=status)
	lam=mdsvalue('dim_of(_sig,0)',/quiet)
	time=mdsvalue('dim_of(_sig,1)',/quiet)
	mdsclose,'spectroscopy',shot
END

PRO wspec_load_xeus,shot,specbr,lam,time,status=status
	wspec_load_vuv,shot,specbr,lam,time,/xeus,status=status
;	IF NOT status THEN BEGIN
;		load_xeus,shot,cnts
;		nframes=n(cnts[0,*])+1
;		xeus_lam_time,lam
;		specbr=cnts/65535.0
;		time=-0.05+indgen(nframes)*1.977/1.0e3
;	ENDIF
;	status=1
END

PRO wspec_load_loweus,shot,specbr,lam,time,status=status
	wspec_load_vuv,shot,specbr,lam,time,/loweus,status=status
;	IF NOT status THEN BEGIN
;		load_loweus,shot,cnts
;		nframes=n(cnts[0,*])+1
;		lowes_lam_time,lam
;		specbr=cnts/65535.0
;		time=-0.05+indgen(nframes)*1.977/1.0e3
;	ENDIF
;	status=1
END

PRO wspec_load_helike,shot,specbr,lam,time,status=status
	mdsopen,'spectroscopy',shot
	specbr=mdsvalue('_sig=\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.SPEC:SPECBR',/quiet,status=status)
	tau=mdsvalue('dim_of(_sig,1)',/quiet)
	lam=mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.SPEC:LAM',/quiet)
	mdsclose,'spectroscopy',shot
	IF status THEN BEGIN
		nch=n(where(specbr[*,0,0] NE -1))+1		;assume CHMAP is not time evolving and take the middle (core) ch
		lam=reform(lam[nch/2,0,*])
		tmp=where(tau NE -1)
		specbr=reform(specbr[nch/2,*,*])
		specbr=reform(specbr[tmp,*])
		time=tau[tmp]
		specbr=rotate(specbr,4)/1.0e3
	ENDIF
END

PRO wspec_load_hlike,shot,specbr,lam,time,status=status
	mdsopen,'spectroscopy',shot
	specbr=mdsvalue('_sig=\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HLIKE.SPEC:SPECBR',/quiet,status=status)
	tau=mdsvalue('dim_of(_sig,1)',/quiet)
	lam=mdsvalue('\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HLIKE.SPEC:LAM',/quiet)
	mdsclose,'spectroscopy',shot

	IF status THEN BEGIN
		nch=n(where(specbr[*,0,0] NE -1))+1		;assume CHMAP is not time evolving and take the middle (core) ch
		lam=reform(lam[nch/2,0,*])
		tmp=where(tau NE -1)
		specbr=reform(specbr[nch/2,*,*])
		specbr=reform(specbr[tmp,*])
		time=tau[tmp]
		specbr=rotate(specbr,4)/1.0e3
	ENDIF
END

PRO wspec_load_2pid,shot,time,prad,status=status
	mdsopen,'spectroscopy',shot
	prad=mdsvalue('_sig=\SPECTROSCOPY::TOP.BOLOMETER:TWOPI_DIODE',/quiet,status=status)*2.75/1.0e3		;define a 2.75 prad-2pi "fudge-factor" [MW]
	time=mdsvalue('dim_of(_sig,0)',/quiet)
	mdsclose,'spectroscopy',shot
END

PRO wspec_load_corextomo,shot,time,br,status=status
	IF NOT keyword_set(ch) THEN ch=19
	IF ch LT 10 THEN ch_str='0' ELSE ch_str=''
	IF NOT keyword_set(array) THEN array=3
	mdsopen,'xtomo',shot
	br=mdsvalue('\XTOMO::TOP.BRIGHTNESSES.ARRAY_'+num2str(array,1)+':CHORD_'+ch_str+num2str(ch,1),/quiet,status=status)
	time=mdsvalue('dim_of(\XTOMO::TOP.BRIGHTNESSES.ARRAY_'+num2str(array,1)+':CHORD_'+ch_str+num2str(ch,1)+')',/quiet)
	mdsclose,'xtomo',shot
	IF status THEN BEGIN
		bl=mean(br[0:ipt(time,-0.01)])
		br-=bl
	ENDIF
END

PRO wspec_load_nl04,shot,time,nl04,status=status
	mdsopen,'electrons',shot
	nl04=mdsvalue('\ELECTRONS::TOP.TCI.RESULTS:NL_04',/quiet,status=status)*1.0e-20
	time=mdsvalue('dim_of(\ELECTRONS::TOP.TCI.RESULTS:NL_04)',/quiet)
	mdsclose,'electrons',shot
END

PRO wspec_load_rfpow,shot,time,rfpow,status=status
	mdsopen, 'rf', shot
	rfpow=mdsvalue('_sig=\rf::RF_POWER_NET',/quiet,status=status)
	time=mdsvalue('dim_of(_sig,0)',/quiet)
	mdsclose,'rf',shot
END

PRO wspec_load_lhpow,shot,time,lhpow,status=status
	mdsopen, 'lh',shot
	lhpow=mdsvalue('_sig=\LH::TOP.RESULTS:NETPOW',/quiet,status=status)
	time=mdsvalue('dim_of(_sig,0)',/quiet)
	mdsclose, 'lh',shot
	if n(lhpow) ne n(time) then begin
	 nlh=min([ n(lhpow), n(time) ])
	 time=time[0:nlh]
	 lhpow=lhpow[0:nlh]
	 print,'fixed unequal data and time vector lengths for lh power'
	endif
END

;6/15/12 - moved line list into GENIE
PRO wspec_load_linestr,line,convert=convert
	;csvpath='/usr/local/cmod/codes/spectroscopy/vuv/vuv_line_list.csv'
	;datpath='/usr/local/cmod/codes/spectroscopy/vuv/vuv_line_list.dat'
	
	;csvpath='/usr/local/cmod/idl/GENIE/IMPSPEC/vuv_line_list.csv'
	;datpath='/usr/local/cmod/idl/GENIE/IMPSPEC/vuv_line_list.dat'
	;IF keyword_set(convert) THEN csv_convert,csvpath,datpath
	;restore,datpath

	GENIE_PATH=getenv('GENIE_PATH') & if GENIE_PATH eq '' then GENIE_PATH='/usr/local/cmod/idl/GENIE/'
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
	IF size(u.dat.prad,/type) EQ 8 THEN BEGIN
		tmp=where(u.dat.prad.time GE xr[0] AND u.dat.prad.time LE xr[1])
		IF tmp[0] NE -1 THEN u.plot.y1[0]=max(u.dat.prad.prad[tmp])*1.05
        ENDIF
	IF size(u.dat.xtomo,/type) EQ 8 THEN BEGIN
		tmp=where(u.dat.xtomo.time GE xr[0] AND u.dat.xtomo.time LE xr[1])
		IF tmp[0] NE -1 THEN u.plot.y1[1]=max(u.dat.xtomo.br[tmp])*1.05
        ENDIF
	IF size(u.dat.nl04,/type) EQ 8 THEN BEGIN
		tmp=where(u.dat.nl04.time GE xr[0] AND u.dat.nl04.time LE xr[1])
		IF tmp[0] NE -1 THEN u.plot.y1[2]=max(u.dat.nl04.nl04[tmp])*1.05
        ENDIF	
	IF size(u.dat.rf,/type) EQ 8 THEN BEGIN
		tmp=where(u.dat.rf.time GE xr[0] AND u.dat.rf.time LE xr[1])
		IF tmp[0] NE -1 THEN u.plot.y1[3]=max(u.dat.rf.pow[tmp])*1.05
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
		IF NOT u.stat.low AND u.stat.mcp THEN spec='mcp'
	ENDIF
	IF lamc GT 3.7 AND lamc LT 3.8 THEN spec='hlike'
	IF lamc GT 3.93 AND lamc LT 4.0 THEN spec='helike'

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
		'hlike' : BEGIN
			ntime=n(u.dat.hrxh.time)+1
			*u.dat.br[4+line]=u.dat.hrxh.time
			tmp=where(u.dat.hrxh.lam GE u.line.lam0[line] AND u.dat.hrxh.lam LE u.line.lam1[line])
			br=sum_array(u.dat.hrxh.specbr[tmp,*],/j)-0.5*(u.dat.hrxh.specbr[tmp[0],*]+u.dat.hrxh.specbr[last(tmp),*])*(n(tmp)+1)
			*u.dat.br[line]=br
                END
		'helike' : BEGIN
			ntime=n(u.dat.hrxhe.time)+1
			*u.dat.br[4+line]=u.dat.hrxhe.time
			tmp=where(u.dat.hrxhe.lam GE u.line.lam0[line] AND u.dat.hrxhe.lam LE u.line.lam1[line])
			br=sum_array(u.dat.hrxhe.specbr[tmp,*],/j)-0.5*(u.dat.hrxhe.specbr[tmp[0],*]+u.dat.hrxhe.specbr[last(tmp),*])*(n(tmp)+1)
			*u.dat.br[line]=br
                 END
		'mcp' : BEGIN
			ntime=n(u.dat.mcp.time)+1
			*u.dat.br[4+line]=u.dat.mcp.time
			tmp=where(u.dat.mcp.lam GE u.line.lam0[line] AND u.dat.mcp.lam LE u.line.lam1[line])
			br=sum_array(u.dat.mcp.specbr[tmp,*],/j)-0.5*(u.dat.mcp.specbr[tmp[0],*]+u.dat.mcp.specbr[last(tmp),*])*(n(tmp)+1)
			*u.dat.br[line]=br
                END
		ELSE : print, 'cannot calc_br for case : '+spec
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
	plot,[0],[0],xr=xr,pos=pos,/xsty,chars=1.0*ls,xticks=xt,xminor=nm,xtickname=replicate(' ',nt+2),yr=yr,ysty=1,ytit='[MW]'
	xyouts,1.05*(xr[1]-xr[0])+xr[0],0.37*(yr[1]-yr[0])+yr[0],'2'+n2g('pi')+' Diode',chars=1.0*ls,orient=90
	IF size(u.dat.prad,/type) EQ 8 THEN oplot,tau,interpol(u.dat.prad.prad,u.dat.prad.time,tau)
	IF u.stat.tline THEN oplot,tau0,yr,linestyle=1.0,color=200
	IF u.stat.ps THEN xyouts,0.8*(xr[1]-xr[0])+xr[0],1.05*(yr[1]-yr[0])+yr[0],num2str(u.shot,1),chars=0.8*ls
	i=3
	yr=[u.plot.y0[4-i],u.plot.y1[4-i]]
	pos=[0.13,low+i*del,0.9,low+(i+1)*del]
	plot,[0],[0],xr=xr,pos=pos,/xsty,chars=1.0*ls,xticks=xt,xminor=nm,xtickname=replicate(' ',nt+2),yr=yr,/ysty,/noerase,ytit='[kW/m!u2!n]'
	xyouts,1.05*(xr[1]-xr[0])+xr[0],0.2*(yr[1]-yr[0])+yr[0],'Core XTOMO',chars=1.0*ls,orient=90
	IF size(u.dat.xtomo,/type) EQ 8 THEN oplot,tau,interpol(u.dat.xtomo.br,u.dat.xtomo.time,tau)
	IF u.stat.tline THEN oplot,tau0,yr,linestyle=1.0,color=200

	i=2
	yr=[u.plot.y0[4-i],u.plot.y1[4-i]]
	pos=[0.13,low+i*del,0.9,low+(i+1)*del]
	plot,[0],[0],xr=xr,pos=pos,/xsty,chars=1.0*ls,xticks=xt,xminor=nm,xtickname=replicate(' ',nt+2),yr=yr,/ysty,/noerase,ytit='[10!u20!n m!u-2!n]'
	xyouts,1.05*(xr[1]-xr[0])+xr[0],0.45*(yr[1]-yr[0])+yr[0],'nl04',chars=1.0*ls,orient=90
	IF size(u.dat.nl04,/type) EQ 8 THEN oplot,tau,interpol(u.dat.nl04.nl04,u.dat.nl04.time,tau)
	IF u.stat.tline THEN oplot,tau0,yr,linestyle=1.0,color=200

	i=1
	yr=[u.plot.y0[4-i],u.plot.y1[4-i]]
	pos=[0.13,low+i*del,0.9,low+(i+1)*del]
	plot,[0],[0],xr=xr,pos=pos,/xsty,chars=1.0*ls,xticks=xt,xminor=nm,xtickname=replicate(' ',nt+2),yr=yr,/ysty,/noerase, ytit='[MW]'
	xyouts,1.05*(xr[1]-xr[0])+xr[0],0.25*(yr[1]-yr[0])+yr[0],'LH',chars=1.0*ls,orient=90,color=196
	xyouts,1.05*(xr[1]-xr[0])+xr[0],0.6*(yr[1]-yr[0])+yr[0],'ICRF',chars=1.0*ls,orient=90
	IF size(u.dat.rf,/type) EQ 8 THEN BEGIN
		oplot,tau,interpol(u.dat.rf.pow,u.dat.rf.time,tau)
	ENDIF
	IF size(u.dat.lh,/type) EQ 8 THEN BEGIN
		oplot,tau,interpol(u.dat.lh.pow,u.dat.lh.time,tau),color=196
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

	;dat={prad:prad,xtomo:xtomo,nl04:nl04,rf:rf,lh:lh}
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.asize[0],u.plot.asize[1],0,0,0]
	ENDIF ELSE device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm

	u.plot.x[0]=!X & u.plot.y[0]=!Y & u.plot.p[0]=!P ; save context
END

FUNCTION wspec_maxsbr,u
	xr=[u.plot.x0[1],u.plot.x1[1]]
	ymax=0.0
	IF u.stat.mcp AND u.plot.mcp THEN BEGIN
		index=ipt(u.stat.time,u.dat.mcp.time)	
		tmp=where(u.dat.mcp.lam GE xr[0] AND u.dat.mcp.lam LE xr[1])
		IF tmp[0] NE -1 THEN ymax=ymax > max(u.dat.mcp.specbr[tmp,index])
        ENDIF
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
	IF u.stat.xeus AND u.plot.xeus THEN BEGIN
		index=ipt(u.stat.time,u.dat.xeus.time)	
		tmp=where(u.dat.xeus.lam GE xr[0] AND u.dat.xeus.lam LE xr[1])
		IF tmp[0] NE -1 THEN ymax=ymax > max(u.dat.xeus.specbr[tmp,index])
        ENDIF
	IF u.stat.hrxh AND u.plot.hrxh THEN BEGIN
		index=ipt(u.stat.time,u.dat.hrxh.time)	
		tmp=where(u.dat.hrxh.lam GE xr[0] AND u.dat.hrxh.lam LE xr[1])
		IF tmp[0] NE -1 AND index NE -1 THEN ymax=ymax > max(u.dat.hrxh.specbr[tmp,index])
        ENDIF
	IF u.stat.hrxhe AND u.plot.hrxhe THEN BEGIN
		index=ipt(u.stat.time,u.dat.hrxhe.time)	
		tmp=where(u.dat.hrxhe.lam GE xr[0] AND u.dat.hrxhe.lam LE xr[1])
		IF tmp[0] NE -1 AND index NE -1 THEN ymax=ymax > max(u.dat.hrxhe.specbr[tmp,index])
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

	IF u.stat.mcp AND u.plot.mcp THEN BEGIN
		index=ipt(u.stat.time,u.dat.mcp.time)
		IF index NE -1 THEN oplot,u.dat.mcp.lam,u.dat.mcp.specbr[*,index],col=col[1]
		u.dat.mcp.index=index
        ENDIF 

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
 
	IF u.stat.hrxhe AND u.plot.hrxhe THEN BEGIN
		index=ipt(u.stat.time,u.dat.hrxhe.time)
		IF index NE -1 THEN oplot,u.dat.hrxhe.lam,u.dat.hrxhe.specbr[*,index],col=col[0]
		u.dat.hrxhe.index=index
        ENDIF 

	IF u.stat.hrxh AND u.plot.hrxh THEN BEGIN
		index=ipt(u.stat.time,u.dat.hrxh.time)
		IF index NE -1 THEN oplot,u.dat.hrxh.lam,u.dat.hrxh.specbr[*,index],col=col[0]
		u.dat.hrxh.index=index
        ENDIF 

	ylab=1.02*(yr[1]-yr[0])+yr[0]
	xlab=0.1*(xr[1]-xr[0])+xr[0]	
	FOR i=0,n(u.line.plot) DO BEGIN
		IF u.line.plot[i] THEN BEGIN
			CASE u.line.spec[i] OF
				'xeus' : BEGIN
					tmp=where(u.dat.xeus.lam GE u.line.lam0[i] AND u.dat.xeus.lam LE u.line.lam1[i])
					index=ipt(u.stat.time,u.dat.xeus.time)
					IF index NE -1 THEN oplot,u.dat.xeus.lam[tmp],u.dat.xeus.specbr[tmp,index],col=lcol[i]
                                END

				'loweus' : BEGIN
					tmp=where(u.dat.low.lam GE u.line.lam0[i] AND u.dat.low.lam LE u.line.lam1[i])
					index=ipt(u.stat.time,u.dat.low.time)
					IF index NE -1 THEN oplot,u.dat.low.lam[tmp],u.dat.low.specbr[tmp,index],col=lcol[i]
                                END
				'hlike' : BEGIN
					tmp=where(u.dat.hrxh.lam GE u.line.lam0[i] AND u.dat.hrxh.lam LE u.line.lam1[i])
					index=ipt(u.stat.time,u.dat.hrxh.time)
					IF index NE -1 THEN oplot,u.dat.hrxh.lam[tmp],u.dat.hrxh.specbr[tmp,index],col=lcol[i]
				END
				'helike' : BEGIN
					tmp=where(u.dat.hrxhe.lam GE u.line.lam0[i] AND u.dat.hrxhe.lam LE u.line.lam1[i])
					index=ipt(u.stat.time,u.dat.hrxhe.time)
					IF index NE -1 THEN oplot,u.dat.hrxhe.lam[tmp],u.dat.hrxhe.specbr[tmp,index],col=lcol[i]
				END
  				'mcp' : BEGIN
					tmp=where(u.dat.mcp.lam GE u.line.lam0[i] AND u.dat.mcp.lam LE u.line.lam1[i])
					index=ipt(u.stat.time,u.dat.mcp.time)
					IF index NE -1 THEN oplot,u.dat.mcp.lam[tmp],u.dat.mcp.specbr[tmp,index],col=lcol[i]
				END
				else  : print, u.line.spec[i]
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
	u.stat.tmax=0	

	;load 2pi-diode prad
	wspec_load_2pid,shot,time,prad,status=status
	IF status THEN prad={time:time,prad:prad} ELSE prad=-1

	;load core xtomo brightness
	wspec_load_corextomo,shot,time,br,status=status
	IF status THEN xtomo={time:time,br:br*1.0e-3} ELSE xtomo=-1

	;load nl04
	wspec_load_nl04,shot,time,nl04,status=status
	IF status THEN nl04={time:time,nl04:nl04} ELSE nl04=-1

	;load ICRF power
	wspec_load_rfpow,shot,time,rfpow,status=status
	IF status THEN rf={time:time,pow:rfpow} ELSE rf=-1
	rft=time

	;load LH power
	wspec_load_lhpow,shot,time,lhpow,status=status
	IF status THEN BEGIN
		time=[rft[0],time,last(rft)]
		lhpow=[0,lhpow,0]
		lh={time:time,pow:lhpow/1.0e3} 
	ENDIF ELSE lh=-1
	widget_control,u.id.message,set_value='SHOT '+num2str(shot,1)+': Time Trace Data Loaded',/append

	;load McPherson spectra
	IF u.stat.load.m THEN wspec_load_mcp,shot,specbr,lam,time,status=status ELSE status=0
	IF status THEN BEGIN
		mcp={specbr:specbr,lam:lam,time:time,index:ipt(time,u.stat.time)} 
		u.stat.mcp=1
		if max(time) gt u.stat.tmax then u.stat.tmax=max(time)
		widget_control,u.id.message,set_value='SHOT '+num2str(shot,1)+': McPherson Loaded',/append
		widget_control,u.id.mcp,set_button=1
        ENDIF ELSE BEGIN
		mcp=-1
		u.stat.mcp=0
		IF u.stat.load.m THEN widget_control,u.id.message,set_value='SHOT '+num2str(shot,1)+': No McPherson Data',/append
		widget_control,u.id.mcp,set_button=0
        ENDELSE

	;load XEUS spectra
	IF u.stat.load.x THEN wspec_load_xeus,shot,specbr,lam,time,status=status ELSE status=0
	IF status THEN BEGIN
		xeus={specbr:specbr,lam:lam,time:time,index:ipt(time,u.stat.time)} 
		u.stat.xeus=1
		if max(time) gt u.stat.tmax then u.stat.tmax=max(time)
		widget_control,u.id.message,set_value='SHOT '+num2str(shot,1)+': XEUS Loaded',/append
		widget_control,u.id.xeus,set_button=1
        ENDIF ELSE BEGIN
		xeus=-1
		u.stat.xeus=0
		IF u.stat.load.x THEN widget_control,u.id.message,set_value='SHOT '+num2str(shot,1)+': No XEUS Data',/append
		widget_control,u.id.xeus,set_button=0
	ENDELSE

	;load LoWEUS spectra
	IF u.stat.load.l THEN wspec_load_loweus,shot,specbr,lam,time,status=status ELSE status=0
	IF status THEN BEGIN
		low={specbr:specbr,lam:lam,time:time,index:ipt(time,u.stat.time)} 
		u.stat.low=1
		if max(time) gt u.stat.tmax then u.stat.tmax=max(time)
		widget_control,u.id.message,set_value='SHOT '+num2str(shot,1)+': LoWEUS Loaded',/append
		widget_control,u.id.low,set_button=1
        ENDIF ELSE BEGIN
		low=-1
		u.stat.low=0
		IF u.stat.load.l THEN widget_control,u.id.message,set_value='SHOT '+num2str(shot,1)+': No LoWEUS Data',/append
		widget_control,u.id.low,set_button=0
        ENDELSE

	;load He-like HIREXSR
	IF u.stat.load.h THEN wspec_load_helike,shot,specbr,lam,time,status=status ELSE status=0
	IF status THEN BEGIN
		hrxhe={specbr:specbr,lam:lam,time:time,index:ipt(time,u.stat.time)} 
		u.stat.hrxhe=1
		if max(time) gt u.stat.tmax then u.stat.tmax=max(time)
		widget_control,u.id.message,set_value='SHOT '+num2str(shot,1)+': HIREXSR (He) Loaded',/append
		widget_control,u.id.hrxhe,set_button=1
        ENDIF ELSE BEGIN
		hrxhe=-1
		u.stat.hrxhe=0
		if max(time) gt u.stat.tmax then u.stat.tmax=max(time)
		IF u.stat.load.h THEN widget_control,u.id.message,set_value='SHOT '+num2str(shot,1)+': No HIREXSR (He) Data',/append
		widget_control,u.id.hrxhe,set_button=0
        ENDELSE

	;load H-like HIREXSR
	IF u.stat.load.h THEN wspec_load_hlike,shot,specbr,lam,time,status=status ELSE status=0
	IF status THEN BEGIN
		hrxh={specbr:specbr,lam:lam,time:time,index:ipt(time,u.stat.time)} 
		u.stat.hrxh=1
		if max(time) gt u.stat.tmax then u.stat.tmax=max(time)
		widget_control,u.id.message,set_value='SHOT '+num2str(shot,1)+': HIREXSR (H) Loaded',/append
		widget_control,u.id.hrxh,set_button=1
        ENDIF ELSE BEGIN
		hrxh=-1
		u.stat.hrxh=0
		IF u.stat.load.h THEN widget_control,u.id.message,set_value='SHOT '+num2str(shot,1)+': No HIREXSR (H) Data',/append
		widget_control,u.id.hrxh,set_button=0
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
	dat={mcp:mcp,low:low,xeus:xeus,hrxhe:hrxhe,hrxh:hrxh,prad:prad,xtomo:xtomo,nl04:nl04,rf:rf,lh:lh,line:line,br:br,sbr:sbr}
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
	widget_control,u.id.hrxhz,set_button=0
	widget_control,u.id.hrxhez,set_button=0
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
				"MCP" : BEGIN
					IF u.stat.dat THEN BEGIN
						IF event.select EQ 1 AND u.stat.mcp THEN u.plot.mcp=1 ELSE u.plot.mcp=0
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
					u.plot.x0[1]=1.0
					u.plot.x1[1]=300.0
					widget_control, u.id.lam0,set_value=num2str(u.plot.x0[1],dp=1)
					widget_control, u.id.lam1,set_value=num2str(u.plot.x1[1],dp=1)
					wspec_plot_spec,u
                                END
				'HRXHEZ' : BEGIN
					wspec_reset_lzoom,u
					u.plot.x0[1]=3.94
					u.plot.x1[1]=4.00
					widget_control, u.id.lam0,set_value=num2str(u.plot.x0[1],dp=1)
					widget_control, u.id.lam1,set_value=num2str(u.plot.x1[1],dp=1)
					wspec_plot_spec,u
                                END
				'HRXHZ' : BEGIN
					wspec_reset_lzoom,u
					u.plot.x0[1]=3.72
					u.plot.x1[1]=3.80
					widget_control, u.id.lam0,set_value=num2str(u.plot.x0[1],dp=1)
					widget_control, u.id.lam1,set_value=num2str(u.plot.x1[1],dp=1)
					wspec_plot_spec,u
                                END
				'LMCP' : IF event.select EQ 1 THEN u.stat.load.m=1 ELSE u.stat.load.m=0
				'LLOW' : IF event.select EQ 1 THEN u.stat.load.l=1 ELSE u.stat.load.l=0
				'LXEUS' : IF event.select EQ 1 THEN u.stat.load.x=1 ELSE u.stat.load.x=0
				'LHRX' : IF event.select EQ 1 THEN u.stat.load.h=1 ELSE u.stat.load.h=0
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
					 u.plot.x0[1]=1.0
					 u.plot.x1[1]=300.0
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
	IF NOT keyword_set(shot) THEN shot=mdscur_shot("cmod")
	IF NOT keyword_set(time) THEN time=1.0
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
	t_slider=widget_slider(B2,xsize=ysz-80,min=0,max=2000,value=1000,/drag,/suppress)

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
	lmcp=widget_button(C1p1y,value=' M')
	llow=widget_button(C1p1y,value=' L')
	lxeus=widget_button(C1p1y,value=' X')
	lhrx=widget_button(C1p1y,value=' H')

	C1p2=widget_base(C1,/row)
	message=widget_text(C1p2,xsize=45,ysize=6,/scroll)
	C1p3=widget_base(C1,/row,/nonexcl)
	xeus=widget_button(C1p3,value=' XEUS ')
	low=widget_button(C1p3,value=' LOWEUS ')
	mcp=widget_button(C1p3, value=' McPherson ')
	C1p4=widget_base(C1,/row,/nonexcl)
	hrxh=widget_button(C1p4, value=' HIREXSR (H) ')
	hrxhe=widget_button(C1p4, value=' HIREXSR (He) ')
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
	hrxhez=widget_button(C1p5c, value='He-LIKE Ar')
	hrxhz=widget_button(C1p5c, value='H-LIKE Ar')

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
		lmcp:lmcp,llow:llow,lxeus:lxeus,lhrx:lhrx,$
		xeus:xeus,low:low,mcp:mcp,hrxhe:hrxhe,hrxh:hrxh,$
		lam0:lam0,lam1:lam1,sbr0:sbr0, sbr1:sbr1,tau0:tau0, tau1:tau1,$
		logsbr:logsbr,autosbr:autosbr,label:label,lelem:lelem,pelem:pelem,aelem:aelem,oelem:oelem,$
		lzoom:lzoom,outz:outz,lowz:lowz,midz:midz,highz:highz,hrxhez:hrxhez,hrxhz:hrxhz,$
		tline:tline,tzoom:tzoom,$
		alab:alab, alam0:alam0,alam1:alam1,afind:afind,aplot:aplot,$
		blab:blab, blam0:blam0,blam1:blam1,bfind:bfind,bplot:bplot,$
		clab:clab, clam0:clam0,clam1:clam1,cfind:cfind,cplot:cplot,$
		dlab:dlab, dlam0:dlam0,dlam1:dlam1,dfind:dfind,dplot:dplot}
	load={m:0,l:1,x:1,h:1}
	stat={time:time,tmax:2,dat:0,ps:0,col:!p.color,pscol:0,tline:1,load:load,xeus:0,low:0,mcp:0,hrxh:0,hrxhe:0}
	plot={	asize:[400,ysz],$
		bsize:[ysz,ysz],$
		y0:float([0,0,  0,0,0,0]),$ ;prad, sxr, nl04, rf, br, specbr
		y1:float([1,100,2,1,1,1]),$
		x0:[0.0,1.0],$ ;time, wavelength
		x1:[2.0,300.0],$
		zoom:[0,0],$
		xeus:1,	low:1,	mcp:1,	hrxh:1,	hrxhe:1,$
		nt:9,	nm:2,	auto:1,	label:1,plabel:100,order:1,ntau:1000,log:0,$
		col:[!p.color,196,96],pscol:[0,0,0],$
		x:[!X,!X],y:[!Y,!Y],p:[!P,!P]$
		}
	line={	label:['Mo XXXII',	'Ar XV',	'Ca XVII',	'Fe XXIII'],$
		lam0: [126.1,		220.4,		192.4,		132.4],$
		lam1: [129.2,		221.6,		193.4,		133.5],$
		plot:[1,1,1,1],$
		col:[50,90,130,180],$
		pscol:[30,100,150,200],$
		spec:['loweus','loweus','loweus','loweus']$
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
	widget_control, u.id.lmcp,set_button=u.stat.load.m
	widget_control, u.id.llow,set_button=u.stat.load.l
	widget_control, u.id.lxeus,set_button=u.stat.load.x
	widget_control, u.id.lhrx,set_button=u.stat.load.h
	widget_control, u.id.lelem,set_value='B,F,Ar,Fe,Mo'
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

