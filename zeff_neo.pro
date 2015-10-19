FUNCTION make,a,b,n
	output=findgen(n)/(n-1)*(b[0]-a[0])+a[0]
	RETURN,output
END

PRO get_ts_local,shot,efit,ts,verb=verb

	mdsopen,'electrons',shot
	te_tsc_n = mdsvalue('\yag_new.results.profiles:te_rz',/quiet,status=status)
	IF status THEN BEGIN
  		te_tsc_n = 1000.*te_tsc_n
  		ne_tsc_n = mdsvalue('\yag_new.results.profiles:ne_rz',/quiet,status=status)
  		r_tsc_n = mdsvalue('yag_new.results.profiles:r_mid_t')
  		t = mdsvalue('dim_of(\yag_new.results.profiles:te_rz)')
  		nt = n_elements(t)
  		te_tsc_n = median(te_tsc_n,3,dim=2)  ;median filter to remove channel drop-outs
  		ne_tsc_n = median(ne_tsc_n,3,dim=2)

  		IF shot LT 1070807001 THEN BEGIN 		;skip after 10708 where "old" data doesn't exist
    			te_tsc_o = mdsvalue('1000.*\electrons::top.yag.results.global.profile:te_rz_t',/quiet)
    			ne_tsc_o = mdsvalue('\electrons::top.yag.results.global.profile:ne_rz_t',/quiet)
    			r_tsc_o = mdsvalue('\electrons::top.yag.results.global.profile:r_mid_t',/quiet)
  		ENDIF

  		;get edge data
  		te_tse = mdsvalue('\ts_te',/quiet,status=status)
  		IF status THEN BEGIN
   			ne_tse = mdsvalue('\ts_ne',/quiet,status=status)
    			r_tse = mdsvalue('\ts_rmid',/quiet,status=status)
    			t_tse = mdsvalue('dim_of(\ts_te)')
    			tse_stat = 1
    			IF keyword_set(verb) THEN print, 'Got Edge TS'
  		ENDIF ELSE BEGIN
			IF keyword_set(verb) THEN print,'Error Getting Edge TS'
     			tse_stat = 0
     			r_tse = fltarr(nt,1) + .91
     			te_tse = fltarr(nt,1)
     			ne_tse = fltarr(nt,1)
		ENDELSE
	ENDIF ELSE BEGIN
	  	IF keyword_set(verb) THEN print, 'Error Getting  TS'
  		ts = {te:fltarr(efit.nt,efit.nr),ner:fltarr(efit.nt,efit.nr),te0:fltarr(efit.nt),ne0:fltarr(efit.nt)}
 	 	RETURN
	ENDELSE
	mdsclose,'electrons',shot


  	;now join and sort data
 	IF shot LT 1070807001 THEN BEGIN
    		r_temp = transpose([transpose(r_tsc_o),transpose(r_tsc_n),transpose(r_tse)]) > 0
		te_temp = transpose([transpose(te_tsc_o),transpose(te_tsc_n),transpose(te_tse)])
    		ne_temp = transpose([transpose(ne_tsc_o),transpose(ne_tsc_n),transpose(ne_tse)])
  	ENDIF ELSE BEGIN
    		r_temp = transpose([transpose(r_tsc_n),transpose(r_tse)]) > 0
    		te_temp = transpose([transpose(te_tsc_n),transpose(te_tse)])
    		ne_temp = transpose([transpose(ne_tsc_n),transpose(ne_tse)])
	ENDELSE

  	r_raw = r_temp*0
  	te_raw = te_temp*0
  	ne_raw = ne_temp*0

  	s = size(r_raw)
  	r_pad = fltarr(s(1)) + .5
  	te_pad =  fltarr(s(1))
  	ne_pad =  fltarr(s(1))

  	FOR j = 0,n_elements(t)-1 DO BEGIN  ;loop over times
    		i = sort(r_temp(j,*))
    		r_raw(j,*) = r_temp(j,i) 
    		te_raw(j,*) =  te_temp(j,i)
    		ne_raw(j,*) =  ne_temp(j,i)
    		te_pad(j) = te_raw(j,0)
    		ne_pad(j) = ne_raw(j,0)
  	ENDFOR

  	;pad array by duplicating value at smallest R - this helps with 2d interpolation
  	r_raw = transpose([transpose(r_pad),transpose(r_raw)])
  	te_raw = transpose([transpose(te_pad),transpose(te_raw)])
  	ne_raw = transpose([transpose(ne_pad),transpose(ne_raw)])

  	;now truncate to include only actual data
  	rtest = total(r_raw,2)
  	j = where(rtest gt 0)
 	r_raw = r_raw(j,*)
  	te_raw = te_raw(j,*)
  	ne_raw = ne_raw(j,*)
  	t = t[j]

  	;try some smoothing in addition to the median filtering done above
  	width = [1,5] 
  	nex = smooth(ne_raw,width)
  	tex = smooth(te_raw,width)

  	;map data onto efit time and space grid
  	s = size(r_raw)
  	t_2d = t#(fltarr(s(2))+1.)
  	triangulate,r_raw,t_2d,tr
  	te_fine = trigrid(r_raw,t_2d,tex,tr,xout=efit.rgrid,yout=efit.tgrid)
  	ne_fine = trigrid(r_raw,t_2d,nex,tr,xout=efit.rgrid,yout=efit.tgrid)

  	ts = {te:transpose(te_fine),ner:transpose(ne_fine),r_raw:r_raw,t_raw:t,te_raw:te_raw,ne_raw:ne_raw}
	IF keyword_set(verb) THEN print, 'Got Core TS'
END

;+
;NAME:
;	ZEFF_NEO
;
;PURPOSE:
;	This procedure is used to estimate the time evolution of zeff by comparing the measured current
;	with the current calculated using neoclassical conductivity.
;
;CALLING SEQUENCE:
;	ZEFF_NEO,shot,zeff,times
;
;INPUTS:
;	shot:		LONG	shot number
;
;OPTIONAL INPUTS:
;	zion:		INT	charge of the plasma main ion DEFAULT: 1.0
;	zeff_max	FLOAT	of the maximum zeff possible DEFAULT: 4.0
;	n_zeff:		INT	number of points to grid zeff from zion -> zeff_max DEFAULT: 8
;	dt:		FLOAT	time interval [sec] to averge over when calculating zeff DEFAULT: 0.1
;	trange:		FLTARR	[t1,t2] of the time window [sec] to do calculation DEFAULT: [0.5,1.5]
;
;KEYWORD PARAMETERS:
;	plot:		/plot will plot output to the active window for debugging purposes
;	verb:		/verb will print output to the terminal
;	qfit:		/qfit will use the quick fits to the Thomson data instedy of GET_TS_LOCAL
;
;OUTPUTS:
;	zeff:		FLTARR	[(t2-t1)/dt] of the zeff
;	times:		FLTARR	[(t2-t1)/dt] of the time points for the zeff [sec]
;
;OPTIONAL OUTPUTS:
;	ts:		STRUC of the Thomson scattering (raw or QFIT) data used in the calculations
;	efit:		STRUC of the efit equilibrium data used
;
;PROCEDURE:
;	This procedure was adapted from Martin Greenwald's confinement widget in
;	/home/g/transport/wconfinement/confinement_get.pro.  It uses the formulas from 
;	O. Sauter et al. PoP v6 pg2834.  Bootstrap current is currently not included.
;
;	Plasma current is calculated using E|| from measurements and conductivity from measurements
;	and a range of assumed zeff.  Ip_neo vs. Zeff is compared to Ip_measured to determine Zeff.
;
;MODIFICATION HISTORY:
;	Adapted from: CONFINEMENT_GET.PRO (g@psfc.mit.edu)
;	Written by:	ML Reinke 1/6/09
;	12/10/2010	M.L. Reinke - added the ability to use QFIT for the density and temperature profiles
;	6/8/2011	M.L. Reinke - updated the documentation
;
;-

PRO zeff_neo,shot,zeff,times,zion=zion,zeff_max=zeff_max,n_zeff=n_zeff,dt=dt,trange=trange,plot=plot,verb=verb,ts=ts,efit=efit,qfit=qfit

	IF NOT keyword_set(zion) THEN zion=1.0
	IF NOT keyword_set(zeff_max) THEN zeff_max=4.0
	IF NOT keyword_set(n_zeff) THEN n_zeff=8
	IF NOT keyword_set(dt) THEN dt=0.1
	IF NOT keyword_set(trange) THEN trange=[0.5,1.5]	
	IF NOT keyword_set(bs) THEN bs=0.0
	zeff_grid=make(zion,zeff_max,n_zeff)

	;load EFIT data from tree
	mdsopen,'analysis',shot
  	t = mdsvalue('dim_of(\efit_geqdsk:pcurrt,2)')	;sec
	rout = mdsvalue('\efit_aeqdsk:rout')/100.0	;geometric center meters
   	rmid = mdsvalue('\EFIT_RMID') 			;outer midplane intercept of normalized rmid(time,psi)
  	ip = abs(mdsvalue('\efit_aeqdsk:pasmat'))	;plasma current - amps
  	ip2 = mdsvalue('\efit_aeqdsk:cpasma')		;alt plasma current - amps?
  	li = mdsvalue('\efit_aeqdsk:ali')
  	wmhd = mdsvalue('\efit_aeqdsk:wplasm')		;plasma energy - joules       	
	psurf = mdsvalue('\efit_aeqdsk:sibdry*$2pi')
  	nbbbs = mdsvalue('\EFIT_GEQDSK:NBBBS')
  	rbbbs = mdsvalue('\efit_geqdsk:rbbbs')		;r boundary cm
  	zbbbs = mdsvalue('\efit_geqdsk:zbbbs')		;z boundary cm 
  	volume = mdsvalue('\efit_aeqdsk:vout')/1.0e6	;volume inside lcfs (cm3)
  	volp = mdsvalue('\efit_fitout:volp')   		;volume of each flux surface (cm3)
  	qpsi = mdsvalue('\efit_fitout:qpsi')		;qpsi(t,psi) 
	mdsclose,'analysis',shot
  	IF keyword_set(verb) THEN print, 'EFIT data loaded'
 	s = size(rmid)
 	ntime = s(1)
  	nrad = s(2)

	;calculate p_ohmic 
	vsurf = deriv(t,psurf)
	len2 = fltarr(ntime)
        FOR i=0,ntime-1 DO BEGIN
                nb = nbbbs(i)
                IF nb NE 0 THEN BEGIN
                        rb = reform(rbbbs(i,0:nb-1),/over)
                        zb = reform(zbbbs(i,0:nb-1),/over)
                        len2(i) = total(sqrt( (rb-shift(rb,-1))^2 + (zb-shift(zb,-1))^2) )^2
                ENDIF ELSE len2(i) = 1.e10
        ENDFOR
	poynt = vsurf * ip
	wmag = 2.e-7*!pi * li * ip2^2 * volume / len2
	dwmag = deriv(t,wmag)
	dwdt = deriv(t,wmhd)
	poh = abs(poynt - dwmag)

	;create regular array for interpolation  
	nt = fix(1000.*(max(t)-min(t)))  	;standard mapping to 1 msec resolution
	tgrid = min(t) + .001*findgen(nt)   
	nr = fix(500.*(.91 - min(rmid)))
	rgrid = min(rmid) + .002*findgen(nr)  	;5 mm resolution

	;interpolate 1D data
	rout = interpol(rout,t,tgrid)	
	ip = interpol(ip,t,tgrid)
	poh = interpol(poh,t,tgrid)

	;interpolate 2D data
	t_2d = t#(fltarr(s(2))+1.)
	triangulate,rmid,t_2d,tr
	qpsi = transpose(trigrid(rmid,t_2d,qpsi,tr,xout=rgrid,yout=tgrid))
	volp = transpose(trigrid(rmid,t_2d,volp,tr,xout=rgrid,yout=tgrid)) 

	;calculate data
	dv = (volp - shift(volp,0,1)) > 0	;get differential volumes
	vres=poh/ip
	efit={nt:nt,nr:nr,rout:rout,rgrid:rgrid,tgrid:tgrid,qpsi:qpsi,vres:vres,dv:dv}	;efit structure
	IF keyword_set(verb) THEN print, 'EFIT data interpolated'
 
	;hardcode use of TS for density and temperature
	IF NOT keyword_set(qfit) THEN BEGIN
	 	get_ts_local,shot,efit,ts,verb=verb
		ner = ts.ner		;dens on the [tgrid,rgrid]
		te = ts.te
	ENDIF ELSE BEGIN
		qfit,shot,dens,temp,rmaj,tau
		ner=fltarr(efit.nt,efit.nr)
		te=fltarr(efit.nt,efit.nr)
		FOR i=0,efit.nt-1 DO BEGIN
			index=ipt(tau,efit.tgrid[i])
			IF index[0] NE -1 THEN BEGIN
				ner[i,*]=interpol(dens[*,index],rmaj[*,index],efit.rgrid)*1.0e20
				te[i,*]=interpol(temp[*,index],rmaj[*,index],efit.rgrid)*1.0e3
			ENDIF
		ENDFOR
		ts={ner:ner,te:te}
	ENDELSE

	;calculate neoclassical plasma current for various assumed zeff
	ip_neo=fltarr(efit.nt,n_zeff)
	rvec = fltarr(efit.nr)+1
 	tvec = fltarr(efit.nt)+1   
   	majr = efit.rout#rvec
    	minr = tvec#efit.rgrid - majr
    	eps = minr/majr
    	ft = sqrt(eps)   ;approximate trapping fraction
	;ftu=1.0-1.0/(1.0-eps^2)^0.5*(1-1.5*sqrt(eps)+0.5*eps^1.5)
	;ftl=3.0*sqrt(2.0)/!pi*sqrt(eps)
	;ft=0.75*ftu+0.25*ftl
    	lambdae = 31.3 - alog(sqrt(ner)/te)  ;2d variables
    	lambdai = 30. - alog(zion^3.0*sqrt(ner)/te^1.5)
    	ephi = (efit.vres/(2.*!pi*efit.rout))#rvec  ;toroidal electric field
	IF keyword_set(verb) THEN print, 'Calculating Neoclassical Ip for Zeff'
	FOR i=0,n_zeff-1 DO BEGIN
		IF keyword_set(verb) THEN print, 'Zeff ='+string(zeff_grid[i])
		zneo = zeff_grid[i] + fltarr(efit.nt,efit.nr)
    		nz = 0.58 + 0.74/(0.76+zneo) ;2d
    		sig_spitz = 1.9e4*te^1.5/(zneo*nz*lambdae)
    		nustare  = 6.92e-18 *efit.qpsi*majr*ner*zneo*lambdae/(te^2*eps^1.5)
    		ft_eff_33 = ft/(1. + (0.55-0.1*ft)*nustare^.5 + 0.45*(1.-ft)*nustare/zneo^1.5)
    		sig_ratio = 1. - (1. + 0.36/zneo)*ft_eff_33 + (0.59/zneo)*ft_eff_33^2 - (0.23/zneo)*ft_eff_33^3
    		sig_neo = sig_spitz*sig_ratio
    		j = abs(sig_neo*ephi)

    		;now deal with nan's that show up in this calculation
    		mask = fltarr(efit.nt,efit.nr) 
    		ij = where(finite(j),count)
    		if count gt 0 then mask(ij) = 1.
    		ij = where(mask eq 0,count)
    		if count gt 0 then j(ij) = 0
    		ip_neo[*,i] = total(j*efit.dv/(6.28*majr),2)
	ENDFOR
	ip*=1.0e-6
	ip_neo*=1.0e-6

	times=fltarr(floor((trange[1]-trange[0])/dt))
	zeff=fltarr(n_elements(times))
	FOR i=0,n_elements(times)-1 DO BEGIN
		times[i]=trange[0]+(i+0.5)*dt
		tmp=where(tgrid GE times[i]-dt/2.0 AND tgrid LT times[i]+dt/2.0)
		IF tmp[0] NE -1 THEN zeff[i]=interpol(zeff_grid,total(ip_neo[tmp,*],1)/n_elements(tmp),mean(ip[tmp]))
	ENDFOR
			
	IF keyword_set(plot) THEN BEGIN
		!p.multi=[0,0,2]
		colors=make(0,200,n_zeff)
		plot, tgrid,ip,xr=[0.0,2.0],/xsty,yr=[0,max(ip)*2.0],/ysty,xtit='Time [sec]',ytit='Plasma Current [MA]',chars=1.2,tit=string(shot)+' BS = '+num2str(bs,dp=2)
		FOR i=0,n_zeff-1 DO BEGIN
			oplot, tgrid,ip_neo[*,i],color=colors[i]
			xyouts,1.6,max(ip)*2.0-(i+1)*0.15*max(ip),string(zeff_grid[i]),color=colors[i]
		ENDFOR

		plot, times,zeff,psym=-6,xtit='Time[sec]',xr=[0.0,2.0],/xsty,yr=[0,max(zeff)*1.1],/ysty,ytit='Zeff',chars=1.2
		!p.multi=[0]
	ENDIF

 END
