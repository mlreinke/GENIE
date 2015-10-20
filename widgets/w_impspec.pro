PRO w_impspec_plotbr,u
	IF u.stat.ps THEN BEGIN
		xsize=7.0
		ysize=4.5
		ls=0.9
		col=u.plot.pscol 
		tit=num2str(u.shot,1)+' t='+num2str(u.stat.time,dp=3)
		d_old=!d
		device,xsize=xsize,ysize=ysize,/in
	ENDIF ELSE BEGIN
		ls=1.0
		col=u.plot.col
		widget_control,u.id.draw2,get_value=draw_win
		window,0,xsize=u.plot.bsize[0],ysize=u.plot.bsize[1],/pixmap
		tit=''
        ENDELSE	
	ibr=*u.dat.br[u.stat.zindex]
	IF u.plot.bxa THEN BEGIN
		xr=[min(ibr.(u.stat.lindex)[*,1]),max(ibr.(u.stat.lindex)[*,1])]
		u.plot.bxr=xr
		widget_control,u.id.bx0,set_value=num2str(u.plot.bxr[0],dp=2)
		widget_control,u.id.bx1,set_value=num2str(u.plot.bxr[1],dp=2)
        ENDIF 
	xr=u.plot.bxr
	tmp=where(ibr.(u.stat.lindex)[*,1] GE xr[0] AND ibr.(u.stat.lindex)[*,1] LE xr[1])
	IF u.plot.bya THEN BEGIN
		yr=[0,max(ibr.(u.stat.lindex)[tmp,0]+ibr.(u.stat.lindex)[tmp,2])]*1.05
		u.plot.byr=yr
		widget_control,u.id.by0,set_value=num2str(u.plot.byr[0],dp=3)
		widget_control,u.id.by1,set_value=num2str(u.plot.byr[1],dp=3)
        ENDIF 
	yr=u.plot.byr
	plot,[0],[0],xr=xr,yr=yr,/xsty,/ysty,xtit='Time [sec]',ytit='Brightness'

	FOR i=0,ibr.nlines-1 DO BEGIN
		IF i EQ u.stat.lindex THEN color=col[3] ELSE color=col[0]
		oploterror,ibr.(i)[tmp,1],ibr.(i)[tmp,0],ibr.(i)[tmp,2],color=color,errcolor=color
        ENDFOR
	index=ipt(ibr.(u.stat.lindex)[*,1],u.stat.time)
	IF index[0] NE -1 THEN oploterror,[ibr.(u.stat.lindex)[index,1]],[ibr.(u.stat.lindex)[index,0]],[ibr.(u.stat.lindex)[index,2]],color=col[4],errcolor=col[4]
	oplot,u.stat.time*[1,1],yr,linestyle=1.0,color=col[4]	
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.asize[0],u.plot.asize[1],0,0,0]
	ENDIF ELSE device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm

END

PRO w_impspec_plotcoefs,u
	IF u.stat.ps THEN BEGIN
		xsize=7.0
		ysize=4.5
		ls=0.9
		col=u.plot.pscol 
		tit=num2str(u.shot,1)+' t='+num2str(u.stat.time,dp=3)
		d_old=!d
		device,xsize=xsize,ysize=ysize,/in
	ENDIF ELSE BEGIN
		ls=1.0
		col=u.plot.col
		widget_control,u.id.draw2,get_value=draw_win
		window,0,xsize=u.plot.bsize[0],ysize=u.plot.bsize[1],/pixmap
		tit=''
        ENDELSE	
	aplot=u.plot.aplot
	icoefs=*u.dat.coefs[u.stat.zindex]
	ibr=*u.dat.br[u.stat.zindex]
	x=ibr.(u.stat.lindex)[*,1]
	icoefs=icoefs.(u.stat.lindex)
	size=size(icoefs)
	nsub=(size[1]-2)/3
	CASE aplot OF 
		0 : BEGIN
			y=icoefs[indgen(nsub)*3,*,0]
			yerr=icoefs[indgen(nsub)*3,*,1]
			ytit='A(0) - HEIGHT'
			nlines=nsub
			index=u.stat.jindex
                END
		1 : BEGIN
			y=icoefs[indgen(nsub)*3+1,*,0]
			yerr=icoefs[indgen(nsub)*3+1,*,1]
			ytit='A(1) - MEAN'
			nlines=nsub
			index=u.stat.jindex
                END
		2 : BEGIN
			y=icoefs[indgen(nsub)*3+2,*,0]
			yerr=icoefs[indgen(nsub)*3+2,*,1]
			ytit='A(2) - WIDTH'
			nlines=1
			index=0		
                END
		3 : BEGIN
			y=icoefs[nsub*3,*,0]
			yerr=icoefs[nsub*3,*,1]
			ytit='A(3) - SLOPE'
			nlines=1
			index=0
                END
		4 : BEGIN
			y=icoefs[nsub*3+1,*,0]
			yerr=icoefs[nsub*3+1,*,1]
			ytit='A(4) - DC'
			nlines=1
			index=0
                    END
		ELSE :		
        ENDCASE
	
	IF u.plot.bxa THEN BEGIN
		xr=[min(x),max(x)]
		u.plot.bxr=xr
		widget_control,u.id.bx0,set_value=num2str(u.plot.bxr[0],dp=2)
		widget_control,u.id.bx1,set_value=num2str(u.plot.bxr[1],dp=2)
        ENDIF 
	xr=u.plot.bxr
	tmp=where(x GE xr[0] AND x LE xr[1])
	IF u.plot.bya THEN BEGIN
		yr=[min(y[index,tmp]-yerr[index,tmp]),max(y[index,tmp]+yerr[index,tmp])]
		IF aplot EQ 0 OR aplot EQ 4 THEN yr[0]=0	
		u.plot.byr=yr
		widget_control,u.id.by0,set_value=num2str(u.plot.byr[0],dp=3)
		widget_control,u.id.by1,set_value=num2str(u.plot.byr[1],dp=3)
        ENDIF
	yr=u.plot.byr

	plot,[0],[0],xr=xr,yr=yr,/xsty,/ysty,xtit='Time [sec]',ytit=ytit
	FOR i=0,nlines-1 DO BEGIN
		IF i EQ index THEN color=col[3] ELSE color=col[0]
		oploterror,x[tmp],y[i,tmp],yerr[i,tmp],color=color,errcolor=color
        ENDFOR
	tindex=ipt(x,u.stat.time)
	IF tindex[0] NE -1 THEN oploterror,[x[index]],[y[index,tindex]],[yerr[index,tindex]],color=col[4],errcolor=col[4]
	oplot,u.stat.time*[1,1],yr,linestyle=1.0,color=col[4]
	oplot,xr,[0,0],linestyle=1.0,color=col[0]	

	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.asize[0],u.plot.asize[1],0,0,0]
	ENDIF ELSE device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm

END

PRO w_impspec_plotspec,u
	IF u.stat.ps THEN BEGIN
		xsize=7.0
		ysize=4.5
		ls=0.9
		col=u.plot.pscol 
		tit=num2str(u.shot,1)+' t='+num2str(u.stat.time,dp=3)
		d_old=!d
		device,xsize=xsize,ysize=ysize,/in
	ENDIF ELSE BEGIN
		ls=1.0
		col=u.plot.col
		widget_control,u.id.draw1,get_value=draw_win
		window,0,xsize=u.plot.asize[0],ysize=u.plot.asize[1],/pixmap
		tit=''
        ENDELSE

	i=u.stat.zindex
	j=u.stat.lindex
	ispec=*u.dat.ispec[i]
	icoefs=*u.dat.coefs[i]
	ibr=*u.dat.br[i]
	ispec=*u.dat.ispec[u.stat.zindex]
	IF size(ispec,/type) EQ 8 THEN BEGIN
		spec=ispec.(j).spec
	 	coefs=icoefs.(j)[*,ipt(ibr.(j)[*,1],u.stat.time),0]
		index=where(u.dat.spec EQ spec)
		data=*u.dat.data[index[0]]
		k=ipt(u.stat.time,data.time)
		tmp=where(data.lam GE ispec.(j).dlam[0] AND data.lam LE ispec.(j).dlam[1])
		x=data.lam[tmp]
		y=data.specbr[tmp,k]
		yerr=data.sig[tmp,k]
      
		yfit=gaussian_fits(x,coefs)
	 	xr=[min(x), max(x)]
		xplot=make(x[0],last(x),u.plot.nx)
		gfit=gaussian_fits(xplot,coefs)
		yr=[min(y-yerr) < 0,(max(y+yerr) > max(gfit))*1.05]
		plot,[0],[0],xr=xr,yr=yr,/xsty,/ysty,ytit='Spec. Bright.',xtit='Wavelength [Ang]',chars=1.0*ls,pos=[0.1,0.4,0.98,0.98]
		oploterror,x,y,yerr,psym=8,color=col[0]
		oplot,xr,[0,0],linestyle=1
		oplot,xplot,gfit,color=col[2]
		ilines=(n(coefs)-1)/3
		FOR i=0,ilines-1 DO BEGIN
			IF i EQ u.stat.jindex THEN BEGIN
				color=col[3] 
				thick=2.0
                        ENDIF ELSE BEGIN
				color=col[1]
				thick=1.0
			ENDELSE
			oplot,xplot,gaussian_fits(xplot,coefs[3*i:3*(i+1)-1]),color=color,linestyle=2.0,thick=thick
		ENDFOR
		oplot,xplot,coefs[3*ilines]*(xplot-xplot[0])+last(coefs),color=col[1],linestyle=3.0
	
		yr=max(abs(yfit-y)+yerr)
		plot,[0],[0],xr=xr,yr=[-yr,yr],/xsty,/ysty,ytit='Spec. Bright.',xtit='Wavelength [Ang]',chars=1.0*ls,pos=[0.1,0.08,0.98,0.32],/noerase
		oploterror,x,y-yfit,yerr,psym=8
		oplot,xr,[0,0],linestyle=1
	
	 	IF NOT u.stat.ps THEN BEGIN
			wset,draw_win
			device,copy=[0,0,u.plot.asize[0],u.plot.asize[1],0,0,0]
		ENDIF ELSE device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm
	ENDIF
END

PRO w_impspec_load,u
	tree=getenv('IMPSPEC_MDS_TREE')
	node=getenv('IMPSPEC_MDS_PATH')

	WIDGET_CONTROL,/hourglass
	IF u.stat.dat THEN BEGIN
		heap_free,u.dat.ispec
		heap_free,u.dat.br
		heap_free,u.dat.coefs
	ENDIF
	shot=u.shot
	
	zlist=u.z
	nz=u.stat.nz
	zchk=intarr(nz)
	ispec=ptrarr(nz,/allocate_heap)
	br=ptrarr(nz,/allocate_heap)
	coefs=ptrarr(nz,/allocate_heap)

	mdsopen,tree,shot
	widget_control,u.id.message,set_value='LOADING SHOT '+num2str(shot,1),/append
	FOR i=0,nz-1 DO BEGIN
		zstr=read_atomic_name(zlist[i])
		nlines=mdsvalue(node+'.'+strupcase(zstr[0])+':NLINES',/quiet,status=status)
		IF status THEN BEGIN
			zchk[i]=nlines
			widget_control,u.id.message,set_value='LOADING '+num2str(nlines,1)+' LINES FOR Z='+num2str(zlist[i],1) ,/append
			xbr={elem:zstr,z:zlist[i],nlines:nlines}
			xcoefs={elem:zstr,z:zlist[i],nlines:nlines}
			FOR k=0,nlines-1 DO BEGIN
				j=nlines-1-k
				path=node+'.'+strupcase(zstr[0])+'.LINE'+num2str(j,1)
				jbr=mdsvalue('_sig='+path+':BR',/quiet)
				jtime=mdsvalue('dim_of(_sig,0)',/quiet)
				jbsig=mdsvalue('dim_of(_sig,1)',/quiet)
				jcoefs=mdsvalue('_sig='+path+':COEFS',/quiet)
				jcsig=mdsvalue('dim_of(_sig,1)',/quiet)
				addstr="'line"+num2str(j,1)+"',[[[jcoefs]],[[jcsig]]]"
				result=execute("xcoefs=create_struct("+addstr+",xcoefs)")	
				addstr="'line"+num2str(j,1)+"',[[jbr],[jtime],[jbsig]]"
				result=execute("xbr=create_struct("+addstr+",xbr)")			
                        ENDFOR
			*br[i]=xbr
			*coefs[i]=xcoefs
                ENDIF ELSE BEGIN
			widget_control,u.id.message,set_value='No Z='+num2str(zlist[i],1)+' DATA' ,/append
			*ispec[i]=-1
			*br[i]=-1
			*coefs[i]=-1
                ENDELSE

        ENDFOR
	mdsclose,tree,shot
	widget_control,u.id.message,set_value='SHOT '+num2str(shot,1)+' LOADED',/append
	FOR i=0,nz-1 DO IF zchk[i] NE 0 THEN *ispec[i]=read_ispec_tree(shot,zlist[i])
	instr=''	;create list spectrometers for all lines
	FOR i=0,n(ispec) DO BEGIN
		xspec=*ispec[i]
		IF size(xspec,/type) EQ 8 THEN FOR j=0,xspec.nlines-1 DO instr=[instr,xspec.(j).spec]
	ENDFOR
	instr=instr[1:*]
	spec=instr[0]	;create list of unique spectrometers
	FOR i=0,n(instr) DO BEGIN
		tmp=where(spec EQ instr[i])
		IF tmp[0] EQ -1 THEN spec=[spec,instr[i]]
        ENDFOR
	data=ptrarr(n(spec)+1,/allocate_heap)	;load each spectrometer data structure into a pointer
	FOR i=0,n(spec) DO BEGIN
		vuv_load_spec,shot,spec[i],specbr,lam,time,sigbr=sigbr
		idat={specbr:specbr,sig:sigbr,lam:lam,time:time,nt:n(time)+1}
		*data[i]=idat
        ENDFOR

	u.stat.dat=1 
	dat={data:data,spec:spec,zchk:zchk,ispec:ispec,br:br,coefs:coefs}
	u={id:u.id,shot:u.shot,stat:u.stat,plot:u.plot,z:u.z,dat:dat}
	widget_control,u.id.base, set_uvalue=u	

	ispec=*u.dat.ispec[u.stat.zindex]
	if size(ispec,/type) eq 8 then begin
	 w_impspec_write_ltable,u
	 w_impspec_write_itable,u
	 w_impspec_plot,u
	endif
END

PRO w_impspec_write_ztable,u
	widget_control,u.id.z_table,set_value=transpose(u.z)
	widget_control,u.id.z_table,background_color=[255,255,255]
	widget_control,u.id.z_table,background_color=[0,200,200],use_table_select=[0,u.stat.zindex,0,u.stat.zindex]
END

PRO w_impspec_write_ltable,u
	widget_control,u.id.l_table,get_value=table
	ispec=*u.dat.ispec[u.stat.zindex]
	;if ispec[0] ne -1 then begin
	 FOR i=0,ispec.nlines-1 DO BEGIN
		table[0,i]=num2str(i,1)
		table[1,i]=num2str(ispec.(i).lam[0],dp=3)
		table[2,i]=ispec.(i).label[0]
		table[3,i]=ispec.(i).iso[0]
         ENDFOR
	 FOR i=ispec.nlines,n(table[0,*]) DO table[*,i]=['','','','']
	 widget_control,u.id.l_table,set_value=table
	widget_control,u.id.l_table,background_color=[255,255,255]
	widget_control,u.id.l_table,background_color=[0,200,200],use_table_select=[0,u.stat.lindex,3,u.stat.lindex]
	 u.stat.nl=ispec.nlines
	;endif
END

PRO w_impspec_write_itable,u
	widget_control,u.id.i_table,get_value=table
	ispec=*u.dat.ispec[u.stat.zindex]
	;if ispec[0] ne -1 then begin
	 iline=ispec.(u.stat.lindex)
	 FOR i=0,iline.ilines-1 DO BEGIN
		table[0,i]=num2str(i,1)
		table[1,i]=num2str(iline.z[i],1)
		table[2,i]=num2str(iline.lam[i],dp=3)
		table[3,i]=iline.label[i]
		table[4,i]=iline.iso[i]
         ENDFOR
	 FOR i=iline.ilines,n(table[0,*]) DO table[*,i]=['','','','','']
	 widget_control,u.id.i_table,set_value=table
	widget_control,u.id.i_table,background_color=[255,255,255]
	widget_control,u.id.i_table,background_color=[0,200,200],use_table_select=[0,u.stat.jindex,4,u.stat.jindex]
	 u.stat.ni=iline.ilines
	;endif
END

PRO w_impspec_plot,u
	ispec=*u.dat.ispec[u.stat.zindex]
	if size(ispec,/type) eq 8 then begin
	 w_impspec_plotspec,u
	 w_impspec_plotb,u
	endif
END

PRO w_impspec_plotb,u
	IF u.plot.br THEN w_impspec_plotbr,u ELSE w_impspec_plotcoefs,u
END

PRO w_impspec_event,event
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
				"LOAD" : w_impspec_load,u
				"QUIT" : BEGIN
					IF u.stat.dat THEN BEGIN
						heap_free,u.dat.ispec
						heap_free,u.dat.br
						heap_free,u.dat.coefs
						heap_free,u.dat.spec
					ENDIF
					widget_control,event.top,/destroy
					!except=1
					heap_gc
				END
				"PRINT" : BEGIN
					psplot
					u.stat.ps=1

					psc
					u.stat.ps=0
					xwplot
					set_plot,'x'	
				END
				"STOP" : BEGIN
					stop
                                END
        			"SETCURR": BEGIN
            				u.shot = mdsvalue('current_shot("cmod")')
            				widget_control,u.id.shotid,set_value=num2str(u.shot,1)
					w_impspec_load,u
                                END
				"LM" : BEGIN
					IF u.stat.lindex NE 0 THEN u.stat.lindex-=1
					widget_control,u.id.litext,set_value=num2str(u.stat.lindex,1)
					IF u.stat.dat THEN BEGIN
						w_impspec_write_itable,u
						w_impspec_write_ltable,u
						w_impspec_plot,u
					ENDIF
                                END
				"LP" : BEGIN
					IF u.stat.lindex NE u.stat.nl-1 THEN u.stat.lindex+=1
					widget_control,u.id.litext,set_value=num2str(u.stat.lindex,1)
					IF u.stat.dat THEN BEGIN
						w_impspec_write_itable,u
						w_impspec_write_ltable,u
						w_impspec_plot,u
					ENDIF
                                END
				"IM" : BEGIN
					IF u.stat.dat AND NOT u.plot.br THEN BEGIN	
						IF u.stat.jindex NE 0 THEN u.stat.jindex-=1
						w_impspec_plot,u
						w_impspec_write_itable,u
						widget_control,u.id.itext,set_value=num2str(u.stat.jindex,1)
					ENDIF
                                END
				"IP" : BEGIN
					IF u.stat.dat AND NOT u.plot.br THEN BEGIN	
						IF u.stat.jindex NE u.stat.ni-1 THEN u.stat.jindex+=1
						w_impspec_plot,u
						w_impspec_write_itable,u
						widget_control,u.id.itext,set_value=num2str(u.stat.jindex,1)
					ENDIF
                                END
				"ZM" : BEGIN
					IF u.stat.zindex NE 0 THEN u.stat.zindex-=1
					widget_control,u.id.zitext,set_value=num2str(u.z[u.stat.zindex],1)
					IF u.stat.dat THEN BEGIN
						w_impspec_write_ltable,u
						w_impspec_write_itable,u
						w_impspec_write_ztable,u
						w_impspec_plot,u
					ENDIF
                                END
				"ZP" : BEGIN
					IF u.stat.zindex NE u.stat.nz-1 THEN u.stat.zindex+=1
					widget_control,u.id.zitext,set_value=num2str(u.z[u.stat.zindex],1)
					IF u.stat.dat THEN BEGIN
						w_impspec_write_ltable,u
						w_impspec_write_itable,u
						w_impspec_write_ztable,u
						w_impspec_plot,u
					ENDIF
                                END
				"BPLOT" : BEGIN
					u.plot.br=1
					u.plot.aplot=0
					w_impspec_plotbr,u
                                END
				"A0PLOT" : BEGIN
					u.plot.br=0
					u.plot.aplot=0
					w_impspec_plotcoefs,u
                                END
				"A1PLOT" : BEGIN
					u.plot.br=0
					u.plot.aplot=1
					w_impspec_plotcoefs,u
                                END
				"A2PLOT" : BEGIN
					u.plot.br=0
					u.plot.aplot=2
					w_impspec_plotcoefs,u
                                END
				"A3PLOT" : BEGIN
					u.plot.br=0
					u.plot.aplot=3
					w_impspec_plotcoefs,u
                                END
				"A4PLOT" : BEGIN
					u.plot.br=0
					u.plot.aplot=4
					w_impspec_plotcoefs,u
                                END
				"BXAUTO" : BEGIN
					IF event.select EQ 1 THEN u.plot.bxa=1 ELSE u.plot.bxa=0
					w_impspec_plotb,u
				END
				"BYAUTO" : BEGIN
					IF event.select EQ 1 THEN u.plot.bya=1 ELSE u.plot.bya=0
					w_impspec_plotb,u
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
						w_impspec_plot,u
					ENDIF
				END
				ELSE:
			ENDCASE
		END
		'WIDGET_DRAW' : BEGIN ;process a button event in the draw window
		 	CASE event.id OF
		  	u.id.draw2 : BEGIN ; in time pane
		   		IF (event.press NE 0) OR (event.release NE 0) THEN BEGIN ; ie not a motion event
	       	   			IF event.press EQ 1 THEN BEGIN ; LMB click,click - zoom time L,R
		    				click_loc = convert_coord(event.x, event.y, /device, /to_data)
						IF u.plot.zoom[0] EQ 0 THEN BEGIN
						 	u.plot.zoom[0]=1
							u.plot.bxr[0]=click_loc[0]
							widget_control, u.id.bx0,set_value=num2str(u.plot.bxr[0],dp=3)
                                                ENDIF ELSE BEGIN
							u.plot.zoom[0]=0
						 	IF click_loc[0] LT u.plot.bxr[0] THEN BEGIN
								u.plot.bxr[1]=u.plot.x0[0]
							  	u.plot.bxr[0]=click_loc[0]
							  	widget_control, u.id.bx0,set_value=num2str(u.plot.bxr[0],dp=3)
                                               	 	ENDIF ELSE  u.plot.bxr[1]=click_loc[0]
							widget_control, u.id.bx1,set_value=num2str(u.plot.bxr[1],dp=3)
						 	w_impspec_plotb,u
						ENDELSE	
		    			ENDIF
		    			if event.press eq 2 then begin ; MMB click - set spectrum time
		    				click_loc = convert_coord(event.x, event.y, /device, /to_data)
						u.stat.time=click_loc[0]
						widget_control, u.id.t_slider,set_value=u.stat.time*1000
						widget_control, u.id.t_text,set_value=num2str(u.stat.time,dp=3)
						w_impspec_plotb,u
		    			ENDIF
		    			if event.press eq 4 then begin ; RMB click - zoom out (ie zoom 0,2)
		    				 u.plot.zoom[0]=0
						 u.plot.bxr[0]=0.0
						 u.plot.bxr[1]=u.stat.tmax
						 widget_control, u.id.bx0,set_value=num2str(u.plot.bxr[0],dp=3)
						 widget_control, u.id.bx1,set_value=num2str(u.plot.bxr[1],dp=3)
						 w_impspec_plotb,u
		    			endif
		   		endif else begin ; pointer motion in pane
	  				ptr_loc = convert_coord(event.x, event.y, /device, /to_data)
					widget_control, u.id.tim_ind,set_value='t='+num2str(ptr_loc[0],dp=3)
		  		endelse
                        END
			ELSE :
		 	ENDCASE
		END
   	  
 		"WIDGET_TEXT_CH": BEGIN
			widget_control,event.id,get_value=data
			CASE event.id OF
				u.id.shotid : BEGIN
					u.shot=long(data)
					w_impspec_load,u
				END
				u.id.litext : BEGIN
					IF int(data) LE u.stat.nl-1 THEN BEGIN
						u.stat.lindex=int(data)
						IF u.stat.dat THEN w_impspec_plot,u
                                        ENDIF ELSE widget_control,u.id.litext,set_value=num2str(u.stat.lindex,1)
                                END
				u.id.zitext : BEGIN
					tmp=where(u.z EQ int(data[0]))
					IF tmp[0] NE -1  THEN BEGIN
						u.stat.zindex=tmp[0]
						IF u.stat.dat THEN BEGIN
							w_impspec_write_ltable,u
							w_impspec_plot,u
						ENDIF
                                        ENDIF ELSE widget_control,u.id.zitext,set_value=num2str(u.z[u.stat.zindex],1)
                                END
				u.id.bx0 : BEGIN
					u.plot.bxr[0]=float(data)
					IF u.stat.dat THEN BEGIN
						w_impspec_plotb,u
					ENDIF
       	                        END
				u.id.bx1 : BEGIN
					u.plot.bxr[1]=float(data)
					IF u.stat.dat THEN BEGIN
						w_impspec_plotb,u
					ENDIF
                                END
				u.id.by0 : BEGIN
					u.plot.byr[0]=float(data)
					IF u.stat.dat THEN BEGIN
						w_impspec_plotb,u
					ENDIF
       	                        END
				u.id.by1 : BEGIN
					u.plot.byr[1]=float(data)
					IF u.stat.dat THEN BEGIN
						w_impspec_plotb,u
					ENDIF
       	                        END		
			ELSE: 
			ENDCASE
                    END
		"WIDGET_TABLE_CELL_SEL" : BEGIN
			CASE ename OF 
				"Z_TABLE" : BEGIN
					IF event.sel_top NE -1 THEN u.stat.zindex=event.sel_top
					IF u.stat.dat THEN BEGIN
						ispec=*u.dat.ispec[u.stat.zindex]
						if size(ispec,/type) eq 8 then begin
						 w_impspec_write_ztable,u
						 w_impspec_write_ltable,u
						 w_impspec_write_itable,u
						 w_impspec_plot,u
						endif
						widget_control,u.id.zitext,set_value=num2str(u.z[u.stat.zindex],1)
					ENDIF
                                END
				"L_TABLE" : BEGIN
					IF event.sel_top NE -1 THEN u.stat.lindex=event.sel_top
					IF u.stat.dat THEN BEGIN
						ispec=*u.dat.ispec[u.stat.zindex]
						if size(ispec,/type) eq 8 then begin
						 w_impspec_write_ltable,u
						 w_impspec_write_itable,u
						 w_impspec_plot,u
						endif
						widget_control,u.id.litext,set_value=num2str(u.stat.lindex,1)
					ENDIF
                                     END
				"I_TABLE" : BEGIN
					IF event.sel_top NE -1 THEN u.stat.jindex=event.sel_top
					IF u.stat.dat THEN BEGIN
						ispec=*u.dat.ispec[u.stat.zindex]
						if size(ispec,/type) eq 8 then begin
						 w_impspec_write_itable,u
						 w_impspec_plot,u
						endif
						widget_control,u.id.itext,set_value=num2str(u.stat.jindex,1)
					ENDIF
                                     END

				ELSE : 
                        ENDCASE

		END
		ELSE:
	ENDCASE
	IF button NE 'QUIT' THEN widget_control,event.top,set_uvalue=u		
END

;+
;NAME:
;	W_IMPSPEC
;
;PURPOSE:
;	This procedure launches the W_IMPSPEC widget which is used to
;	check/modify the SXR/VUV/visible spectroscopy fits that
;	populate time-evolving line-brightness tree nodes
;
;CALLING SEQUENCE:
;	@run_wimpspec.pro 	will compile the proper files, set the color
;				tables and launch the widget
;PROCEUDRE:
;	See the C-Mod Wiki page on the W_IMSPEC widget for instructions
;	(listed under Shot Analysis Tools)
;
;RESTRICTIONS:
;	Due to the limited namespace, many of the GENIE widgets have
;	overlapping function/procedure names so these widgets may not run
;	well in contiguous IDL sessions.
;
;	***W_IMPSPEC and the IMPSPEC infrastructure is still
;	evolving.  Please see M.L. Reinke if a change in functionality occurs
;	that breaks established analysis routines.***
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 3/2012
;
;-

PRO w_impspec,shot=shot,time=time,zlist=zlist
	get_genie_env,MACHINE=mach
	CASE mach OF 
		'cmod' : BEGIN
			IF NOT keyword_set(shot) THEN shot=1120210003
			IF NOT keyword_set(time) THEN time=1.0
			IF NOT keyword_set(zlist) THEN zlist=[5,7,8,9,10,18,42,74]
			machtr=[0.01,1.65]
		END
		'nstxu' : BEGIN
			IF NOT keyword_set(shot) THEN shot=141229
			IF NOT keyword_set(time) THEN time=0.3
			IF NOT keyword_set(zlist) THEN zlist=[7,8]
			machtr=[0.01,1.0]
		END
	ENDCASE


	loadct,12,/silent
	base=widget_base(title='IMPSPEC ANALYSIS',/row,tlb_size_events=1)
	
	B=widget_base(base,/column)
	A=widget_base(base,/column)

	ysz=900
	axsz=500
	A1=widget_base(A,/column,xsize=axsz,ysize=ysz*0.4,/frame)
	A2=widget_base(A,/column,xsize=axsz,ysize=ysz*0.55,/frame)

	bxsz=750
	draw1=widget_draw(B,xsize=bxsz,ysize=0.65*ysz,/frame)
	B1=widget_base(B,/row)
	B1p1=widget_base(B1,/row)
	i_table=widget_table(B1p1,xsize=5,ysize=15,column_lab=['#','Z','LAM','LABEL','ISO-E'],value=strarr(5,15),/no_row_headers,y_scroll_size=10,$
		column_width=[35,35,75,150,60],/all_events,alignment=1)

	B2=widget_base(B1,/column,xsize=bxsz*0.25,ysize=0.295*ysz,/frame)
	B2a=widget_base(B2,/row)
	bplot=widget_button(B2,value='BRIGHT')
	B2b=widget_base(B2,/row)
	dum = widget_label(B2b,value='          FIT COEFS')
	B2bx=widget_base(B2,/row)
	space=widget_base(B2bx,/row,xsize=4)
	im=widget_button(B2bx,value=' < ')
	dum = widget_label(B2bx,value='ILINE# ')
	itext = widget_text(B2bx,xsize=3,ysize=1,/edit)
	ip=widget_button(B2bx,value=' > ')
	B2c=widget_base(B2,/row)
	dum = widget_label(B2c,value='')
	a0plot=widget_button(B2c,value='A0')
	dum = widget_label(B2c,value='')
	a1plot=widget_button(B2c,value='A1')
	dum = widget_label(B2c,value='')
	a2plot=widget_button(B2c,value='A2')
	dum = widget_label(B2c,value='')
	a3plot=widget_button(B2c,value='A3')
	dum = widget_label(B2c,value='')
	a4plot=widget_button(B2c,value='A4')
	B2d=widget_base(B2,/row)
	by0=widget_text(B2d,xsize=4,ysize=1,/edit)
	dum = widget_label(B2d,value='< Y <')
	by1=widget_text(B2d,xsize=4,ysize=1,/edit)
	B2dx=widget_base(B2d,/row,/nonexcl)
	byauto=widget_button(B2dx,value='AUTO')
	B2e=widget_base(B2,/row)
	bx0=widget_text(B2e,xsize=4,ysize=1,/edit)
	dum = widget_label(B2e,value='< X <')
	bx1=widget_text(B2e,xsize=4,ysize=1,/edit)
	B2ex=widget_base(B2e,/row,/nonexcl)
	bxauto=widget_button(B2ex,value='AUTO')


	B3=widget_base(B,/row)
	t_text=widget_text(B3,xsize=8,ysize=1,/edit)
	t_slider=widget_slider(B3,xsize=bxsz-120,min=0,max=2000,value=time*1.0e3,/drag,/suppress)


	A1p1=widget_base(A1,/row)
	dum = widget_label(A1p1,value='SHOT: ')
	shotid = widget_text(A1p1,xsize=10,ysize=1,/edit)
	dum = widget_label(A1p1,value='')
	load= widget_button(A1p1,value='LOAD')
	dum = widget_label(A1p1,value='')
	save= widget_button(A1p1,value='SAVE')
	dum = widget_label(A1p1,value='')
	quit= widget_button(A1p1,value='QUIT')
	dum = widget_label(A1p1,value='')
	print= widget_button(A1p1,value='PRINT')
	dum = widget_label(A1p1,value='')
	stop= widget_button(A1p1,value='STOP')
	dum = widget_label(A1p1,value='')
	setcurr= widget_button(A1p1,value='SET CURRENT SHOT')

	A1p2=widget_base(A1,/row)
	message=widget_text(A1p2,xsize=45,ysize=6,/scroll)

	A1p3=widget_base(A1,/row)
	nz=n(zlist)+1
	space=widget_base(A1p3,/row,xsize=25)
	z_table=widget_table(A1p3,xsize=1,ysize=nz,column_lab=['Z'],value=intarr(nz),/no_row_headers,y_scroll_size=4,column_width=35,/all_events,alignment=1)
	space=widget_base(A1p3,/row,xsize=55)
	l_table=widget_table(A1p3,xsize=4,ysize=15,column_lab=['#','LAM','LABEL','ISO-E'],value=strarr(4,15),/no_row_headers,y_scroll_size=4,$
		column_width=[35,75,100,60],/all_events,alignment=1)
	A1p4=widget_base(A1,/row)
	zm=widget_button(A1p4,value=' < ')
	dum = widget_label(A1p4,value='Z=')
	zitext = widget_text(A1p4,xsize=3,ysize=1,/edit)
	zp=widget_button(A1p4,value=' > ')
	space=widget_base(A1p4,/row,xsize=85)
	lm=widget_button(A1p4,value=' < ')
	dum = widget_label(A1p4,value='LINE# ')
	litext = widget_text(A1p4,xsize=3,ysize=1,/edit)
	lp=widget_button(A1p4,value=' > ')

	draw2=widget_draw(A2,xsize=axsz,ysize=0.55*ysz,/frame,/button_events)

	id={base:base,draw1:draw1,draw2:draw2,t_slider:t_slider,t_text:t_text,$
		bplot:bplot,im:im,itext:itext,ip:ip,a0plot:a0plot,a1plot:a1plot,a2plot:a2plot,a3plot:a3plot,a4plot:a4plot,$
		by0:by0,by1:by1,byauto:byauto,bx0:bx0,bx1:bx1,bxauto:bxauto,$
		shotid:shotid,load:load,save:save,quit:quit,print:print,stop:stop,message:message,setcurr:setcurr,$
		z_table:z_table,l_table:l_table,zm:zm,zitext:zitext,zp:zp,lm:lm,litext:litext,lp:lp,$
		i_table:i_table}

	stat={time:time,dat:0,ps:0,lindex:0,zindex:0,jindex:0,ni:0,nl:0,nz:n(zlist)+1,tmax:machtr[1]}
	plot={asize:[bxsz,0.65*ysz],bsize:[axsz,0.55*ysz],col:[255,50,180,90,150],pscol:[0,30,200,100,150],$
		nx:100,br:1,bxr:machtr,byr:[0.0,0.0],bya:1,bxa:0,aplot:0,zoom:0}
	u={id:id,shot:shot,stat:stat,plot:plot,z:zlist}
	widget_control,base,set_uvalue=u
	widget_control,u.id.shotid,set_value=num2str(u.shot,1)
	widget_control,u.id.z_table,set_value=transpose(u.z)
	widget_control,u.id.z_table,background_color=[255,255,255]
	widget_control,u.id.z_table,background_color=[0,200,200],use_table_select=[0,u.stat.zindex,0,u.stat.zindex]
	widget_control,u.id.zitext,set_value=num2str(u.z[u.stat.zindex],1)
	widget_control,u.id.litext,set_value=num2str(u.stat.lindex,1)
	widget_control,u.id.itext,set_value=num2str(u.stat.jindex,1)
	widget_control,u.id.bx0,set_value=num2str(u.plot.bxr[0],dp=2)
	widget_control,u.id.bx1,set_value=num2str(u.plot.bxr[1],dp=2)
	widget_control,u.id.by0,set_value=num2str(u.plot.byr[0],dp=3)
	widget_control,u.id.by1,set_value=num2str(u.plot.byr[1],dp=3)
	widget_control,u.id.byauto,set_button=u.plot.bya
	widget_control,u.id.bxauto,set_button=u.plot.bxa

	!except=0
	widget_control,base,/realize
	w_impspec_load,u
	xmanager,'w_impspec',base
END
