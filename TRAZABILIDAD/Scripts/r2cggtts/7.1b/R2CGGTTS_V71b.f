      program cctf
      
c !!! fr integer pred freqg real
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  VERSION 7.1a : April 2016
c
c This code generates the CGGTTS files in a MULTI-CHANNEL approach
c     from RINEX files at a sampling rate of 30 sec.
c     using broadcast orbits and dual-frequency code
c
c the files needed are:
c    Rinex 3.02 observation file of the day --> 'rinex_obs'
c    Rinex 3.02 observation file of the day after --> 'rinex_obs_p'
c    GPS Rinex 3.02 Navigation file of the day --> 'rinex_nav_gps'
c    GPS Rinex 3.02 Navigation file of the day after --> 'rinex_nav_p_gps'
c    GAL Rinex 3.02 Navigation file of the day --> 'rinex_nav_gal'
c    GAL Rinex 3.02 Navigation file of the day after --> 'rinex_nav_p_gal'
c    GLONASS Rinex 3.02 Navigation file of the day --> 'rinex_nav_glo'
c    GLONASS Rinex 3.02 Navigation file of the day after --> 'rinex_nav_p_glo'
c    coordinates, calibration delays --> 'paramCGGTTS.dat'
c    If the receiver measrues GPS C1 rather than P1, use IGS Biases of the month --> 'biasC1P1.dat'
c
c   an optional inputFile.dat contains the name of the input files, the last line must be the mjd.
c
c the soft asks for the mjd of the day to be treated
c
c The output is in the files 'CGGTTS.GPS CGGTTS.GAL, CGGTTS.GLO'
c
c 2016:
c April 13: change n1(3) into n1(5) in variable definitions of the main
c           add sec in real*8 in subroutine read_head_obs
c April 14: add "BACKSPACE(LUOUT)" in the subroutine Header
c           put Trel=0 for GLONASS satellites
c April 24: allows to work without the observation file for D+1
c           Elevation computed for an ellipsoidal Earth
c           Allows using GLONASS with C/P codes
c           Allows using C1X and C5X for Galileo
c  May 7:   change format 625
c           add "GAL" in the reading of RECDLY for E1 and E5a 
c
c Changes to V7.1b (jan 2017):
c  - increase the size of brorb from 5000 to 15000
c  - modify the format 175 in read_head_obs, i2 was missing
c  - change "if(systime ne 'GPS') call leapscorr" into "if(systime eq 'GLO').."
c  - add the missing declaration of colGLOC1 and colGLOC2
c  - correct col1(iii) into col2(iii) in secion 'reading of observations' 
c  - add 'if(.not.status)goto 150' after 'call satpos'
c  - change "mod" in "dmod"  for secg in subroutine leapscorr
c  - remove the test on satellite health
c  - choose the ephem on the TOC, with a validity of 4h.
c
c contact: Pascale DEFRAIGNE (p.defraigne@oma.be)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      integer*4
     1       LUFIL, LUSCRN, LUPARA,  LUINP,  LULOG,  LUOUT, LUBIAS, 
     2      LUFOBS, LUFOBP,  LUNAV, LUNAVP,   IERR, istop ,IER(20:25)

      integer TRKL(45),yy,mm,dd,hh,min,he,prn,flag(45),mi,satnum(45),
     1  satn(45),secsch,leapsec,col1(45),col2(45),yn,mn,dn,hn,minn,
     1  schsat(300),schsatm(300),schh(300),schmin(300),ic(45),WN,
     2  ichoice_sys,mjdtrk,IOE,ck,elv,azth,ch,fr,hc,pprn_gal(15000),
     2  pprn_gps(15000),colGPSP1,colGPSP2,colGPSC1,colGLOP1,colGLOP2,
     2  colGALE1,colGALE5,fact,snr(45),nb_obs(5),secn,heure,binr(3,2),
     1  iGP1,iGP2, iGC1,iRP1,iRP2,iRC1,iE1,iE5,n1(5),colGLOC2,colGLOC1

     
      real*8 pi,vitlum,schtime(300),sec,mjd,mjddep, a0,a1,corrgeom,     
     1  brorb_gps(15000,7,4),corprev,alpha,we,ds,djul,last,trel,tau,     
     1  deltatk,tgd,brorb_gal(15000,7,4),antpos(3),x(3),xs(3),lla(3),
     1  Az,El,moment,timefirst,timestop,timenav,toc,trc,Ttr,tropo(26),           
     1  mdtr,smdt,stdv, brorb(15000,7,4),svclk_gps(15000,3),iono(26),    
     1  mdio,smdi,measiono(26,45),isg,svclk_gal(15000,3), refsv(26),    
     1  refsvmil,srsv,refgps(26),refgpsmil,srgps,dsg,resquad(26,45),     
     1  pcor,cr,clocksat,svclk(3),tclock,recdelP1,recdelP2,refdel,      
     2  cabdel,deltaT(26),p2valobs(45),epoch(26),valobs(40,45), 
     4  p1valobs(45),datenav_gps(15000),bias(40),NS,DeltaN,mf,zpd,
     4  NSlog,datenav_gal(15000),factor(3,2),freq(3,2),ffactor(3,2),
     4  lambda(3,2),k1(3),k2(3),timenav1,datenav1(15000),diftime,
     4  recdelE5,recdelE1,difsave,tobs_first
      
C    GLONASS ADDITION 	 
      real*8 recdelP1glo,recdelP2glo,tauCa,tauCb,
     1    r,datenavglo(1000),clkoffset,temp
      real*8 coef1,coef2,trcprev,dte,tauNtmp,gammaNtmp,tktmp,
     1    timenavprev,xsp(3),xsm(3),vit(3),delta,x1,x2,x3
      real*8 sp3(1000,30,3),tocglo(1000),clock(1000,30,3),z,
     1    v(1000,30,3),g(1000,30,3),r2,r3,xsn,ysn,zsn,vxsn,vysn,vzsn,
     1    gxsn,gysn,gzsn,health(1000,30),freqg(1000,30),secg,
     1    recdelP1g,recdelP2g
      integer i,prng,yg,mg,dg,hg,ming,nbglo,xx,nl,a
      
      

      logical status, stat , LOGOP,ofset,exists,oofset
      character*2 class,schcl(15000),schclm(15000), brol, year,month,day
      character*3 frc,systime,obs(40),obsgal,ssystime
      CHARACTER*7 CHSTAT, CHOLD, CHUNKO,lengC1P1,sversion
      character*30 param,REF,LAB,RCVR,REVDATE,COMMENTS 
      character*20 comment,dateC1P1,calID
      character*500 string
      character*125 name
      character*512 FNAVGPS, FNAVGAL, FNAVGLO 
      character*512 FNAVPGPS, FNAVPGAL, FNAVPGLO
      character*512 FOBS, FOBSP,  FLLOG, NAMFIL
      character*512 FLOUT_GPS, FLOUT_GAL,FLOUT_GLO
      character*1 param1(30),syst
      CHARACTER*1200 strerr

      equivalence(param1(1),param)

c  HEADER common:
      common/headint/ ch,ICHOICE_SYS
      common/headrel8/ antpos, recdelP1,recdelP2,recdelE1,
     1                 recdelE5,recdelP1g,recdelP2g,cabdel,refdel
      common/headchar/ LAB,REF,RCVR,REVDATE,COMMENTS,calID
      common/leap/leapsec

      freq(1,1)=1575.42d0   !GPS L1
      freq(1,2)=1227.6d0    !GPS L2
      freq(2,1)=1602.d0        !GLONASS L1
      freq(2,2)=1246.d0        !GLONASS L2
      freq(3,1)=1575.42d0   !GALILEO L1
      freq(3,2)=1176.45d0   !GALILEO L5a
      
      
      WE=7.292115d-5
      vitlum=299792458.0d0
      pi=3.1415926535898d0     
       c2=vitlum*vitlum
      do i=1,3
       do ii=1,2
       lambda(i,ii)=vitlum/freq(i,ii)
       enddo
       k1(i)=(freq(i,1)**2/(freq(i,1)**2-freq(i,2)**2))
       k2(i)=(freq(i,2)**2/(freq(i,1)**2-freq(i,2)**2))
      enddo

c
      DATA     LUKEYB,      LUSCRN,      LUPARA
     +/          5,            6,          9/
      DATA      LUPRT,      LULOG
     +/           92,         94/
      DATA      LUINP,      LUFOBS,      LUFOBP
     +/         93,            7,         71/
      DATA      LUNAVGPS,      LUNAVPGPS,      LUBIAS
     +/         20,           21,         26/
      DATA      LUNAVGAL,      LUNAVPGAL  
     +/         22,           23/
      DATA      LUNAVGLO,      LUNAVPGLO  
     +/         24,           25/
      DATA      LUOUT_GAL,      LUOUT_GPS,  LUOUT_GLO
     +/         14,           15,              16/

      DATA CHOLD, CHUNKO /'OLD    ','UNKNOWN'/


c     **** start of executable code ****

      hc=0
      exists=.false.

c     Define the version of the Program (variable so it is modified only once)
c     This should be 3 characters maximum otherwise modify the defined size
      sversion='7.1b'

c  default file names
      FNAVGPS = 'rinex_nav_gps'
      FNAVGAL = 'rinex_nav_gal'
      FNAVGLO = 'rinex_nav_glo'
      FNAVPGPS = 'rinex_nav_p_gps'
      FNAVPGAL = 'rinex_nav_p_gal'
      FNAVPGLO = 'rinex_nav_p_glo'
      FOBS = 'rinex_obs'
      FOBSP = 'rinex_obs_p'
      FLOUT_gal = 'CGGTTS.GAL'
      FLOUT_gps = 'CGGTTS.GPS'
      FLOUT_glo = 'CGGTTS.GLO'
      FLLOG = 'CGGTTS.log'
      IER=0
      
c  default set flag modified julian day is not defined
      mjd = -1.0
                   
c    open input file if it exists
c    subroutine opfil cannot be used as the inputFile may not exist then ierr is not 0
      LUFIL = LUINP
      NAMFIL = 'inputFile.dat'
      CHSTAT = CHOLD
      OPEN(UNIT=LUFIL, FILE= NAMFIL,
     1     STATUS= CHSTAT, IOSTAT= IERR)

       IF(IERR .EQ. 0)THEN         
      REWIND LUINP

 5001 CONTINUE
        READ(LUINP,'(a30)', END=5003) PARAM
        IF(PARAM .EQ. 'FILE_RINEX_NAV_GPS') READ(LUINP,5000) FNAVGPS
        IF(PARAM .EQ. 'FILE_RINEX_NAV_P_GPS') READ(LUINP,5000) FNAVPGPS         
        IF(PARAM .EQ. 'FILE_RINEX_NAV_GLO') READ(LUINP,5000) FNAVGLO
        IF(PARAM .EQ. 'FILE_RINEX_NAV_P_GLO') READ(LUINP,5000) FNAVPGLO
        IF(PARAM .EQ. 'FILE_RINEX_NAV_GAL') READ(LUINP,5000) FNAVGAL
        IF(PARAM .EQ. 'FILE_RINEX_NAV_P_GAL') READ(LUINP,5000) FNAVPGAL
        IF(PARAM .EQ. 'FILE_RINEX_OBS') READ(LUINP,5000) FOBS  
        IF(PARAM .EQ. 'FILE_RINEX_OBS_P') READ(LUINP,5000) FOBSP
        IF(PARAM .EQ. 'FILE_CGGTTS_LOG') READ(LUINP,5000) FLLOG
        IF(PARAM .EQ. 'FILE_CGGTTS_GPS') READ(LUINP,5000) FLOUT_GPS
	IF(PARAM .EQ. 'FILE_CGGTTS_GLO') READ(LUINP,5000) FLOUT_GLO
 	IF(PARAM .EQ. 'FILE_CGGTTS_GAL') READ(LUINP,5000) FLOUT_GAL                       
        IF(PARAM .EQ. 'MODIFIED_JULIAN_DAY') READ(LUINP,*) MJD
        IF(PARAM .EQ. 'MODIFIED_JULIAN_DAY') goto 5003

 5000   FORMAT(A512)

      GOTO 5001
                         
      ENDIF

 5003 CONTINUE
 
c  open LOG file 
          LUFIL = LULOG
          NAMFIL = FLLOG
          CHSTAT = CHUNKO
          CALL OPFIL
     1        (LUFIL, LUSCRN, LULOG,IERR,CHSTAT,LOGOP,NAMFIL)
          IF(IERR .EQ. 0) LOGOP = .TRUE.   
 
c Write the version of the program at start-up
      write(LUSCRN,*)' Program rin2cgg version ',sversion
      IF(LOGOP) WRITE(LULOG,*) ' Program rin2cgg version ',sversion

      recdelP1 = 0
      recdelP2 = 0
      recdelP1g =0
      recdelP2g =0
      recdelE5 = 0
      recdelE1 = 0
      cabdel=0
      refdel=0
      calID='Unknown'
      COMMENTS='NO COMMENT'
      
c  open parameter file must exist
      LUFIL = LUPARA
      NAMFIL = 'paramCGGTTS.dat'
      CHSTAT = CHOLD
      CALL OPFIL
     1    (LUFIL, LUSCRN, LULOG,  IERR, CHSTAT,  LOGOP, NAMFIL)

c     READING OF PARAMETERS:
      recdelE5=0.d0
      do  ipar=1,25
         read(lupara,'(a30)',end=21)param
         if(param.eq.'REV DATE')     read(LUPARA,'(a30)')revdate
         if(param.eq.'RCVR')         read(LUPARA,'(a30)')rcvr
         if(param.eq.'CH')           read(LUPARA,*)CH
         if(param.eq.'LAB NAME')     read(LUPARA,'(a30)')LAB
         if(param.eq.'X COORDINATE') read(LUPARA,*)antpos(1)
         if(param.eq.'Y COORDINATE') read(LUPARA,*)antpos(2)
         if(param.eq.'Z COORDINATE') read(LUPARA,*)antpos(3)
         if(param.eq.'COMMENTS')     read(LUPARA,'(a30)')COMMENTS
         if(param.eq.'REF')          read(LUPARA,'(a30)')REF
         if(param.eq.'INT DELAY P1 GPS (in ns)')
     1                               read(LUPARA,*)recdelP1
         if(param.eq.'INT DELAY C1 GPS (in ns)')
     1                               read(LUPARA,*)recdelP1   ! the same HW delay is used for GPS C1, L1 and Galileo E1
         if(param.eq.'INT DELAY P2 GPS (in ns)')
     1                               read(LUPARA,*)recdelP2
         if(param.eq.'INT DELAY P1 GLO (in ns)')
     1                               read(LUPARA,*)recdelP1g
         if(param.eq.'INT DELAY C1 GLO (in ns)')
     1                               read(LUPARA,*)recdelP1g   ! the same HW delay is used for GPS C1, L1 and Galileo E1
         if(param.eq.'INT DELAY P2 GLO (in ns)')
     1                               read(LUPARA,*)recdelP2g
         if(param.eq.'INT DELAY E1 GAL (in ns)')
     1                               read(LUPARA,*)recdelE1   ! the same HW delay is used for GPS C1, L1 and Galileo E1
         if(param.eq.'INT DELAY E5a GAL (in ns)')
     1                               read(LUPARA,*)recdelE5
         if(param.eq.'ANT CAB DELAY (in ns)')
     1                               read(LUPARA,*)cabdel
         if(param.eq.'CLOCK CAB DELAY XP+XO (in ns)')
     1                               read(LUPARA,*)refdel
         if(param.eq.'CALIBRATION REFERENCE')
     1                               read(LUPARA,'(A20)')calID
         if(param.eq.'LEAP SECOND')  read(LUPARA,*)leapsec
       enddo
             
c     open outputfile new or existing file will be overwritten      
 21   LUFIL = LUOUT_gal
      NAMFIL = FLOUT_gal
      CHSTAT = CHUNKO
      CALL OPFIL
     1    (LUFIL, LUSCRN, LULOG,  IERR, CHSTAT,  LOGOP, NAMFIL)
      LUFIL = LUOUT_GPS
      NAMFIL = FLOUT_GPS
      CHSTAT = CHUNKO
      CALL OPFIL
     1    (LUFIL, LUSCRN, LULOG,  IERR, CHSTAT,  LOGOP, NAMFIL)
      LUFIL = LUOUT_GLO
      NAMFIL = FLOUT_GLO
      CHSTAT = CHUNKO
      CALL OPFIL
     1    (LUFIL, LUSCRN, LULOG,  IERR, CHSTAT,  LOGOP, NAMFIL)
                                   
c tropo parameters       
       NS=324.8d0
       DeltaN=-7.32d0*dexp(0.005577d0*Ns)  
       NSLOG= dlog((NS+DeltaN)/105.d0 ) 
                
c  internal delays in m:
      recdelP1 = recdelP1*vitlum*1.0d-9
      recdelP2 = recdelP2*vitlum*1.0d-9
      recdelP1g = recdelP1g*vitlum*1.0d-9
      recdelP2g = recdelP2g*vitlum*1.0d-9
      recdelE5 = recdelE5*vitlum*1.0d-9
      recdelE1 = recdelE1*vitlum*1.0d-9
      
c  cable delays in 0.1 ns:
      cabdel = cabdel*10.
      refdel = refdel*10.      
      
c     DETERMINATION OF THE TRACKS OF THE DAY
c     MJDdep is the reference date for the schedules  (1 October 1997)

      nsch=89
      ileaps=leapsec/30
      ileaps=(ileaps+1)*30

c  check if modified julian date has been read from input file
      if(mjd .eq. -1.0)then
         write(LUSCRN,*) ' enter mjd:'
         read(LUKEYB,*)mjd
      endif
      
      mjddep=50722.0
      do 22 i=1,nsch
       its=2+16*(i-1)
       mi=mod(its,60)
       he=its/60
       call decal(mjd,mjddep,he,mi,schh(i),schmin(i))
       schtime(i)=mjd+schh(i)/24.d0+schmin(i)/24.d0/60.d0
     1    + dble(ileaps)/86400.d0    ! in GPS TIME
   22 continue     
    
      call classement(schh,schmin,schtime,89)
       
      if(schmin(89).lt.43)then
      nsch=90
      schh(90)=schh(89)
      schmin(90)=schmin(89)+16
      schtime(90)=schtime(89)+dble(16.d0/1440.d0)
      endif
      
c  open observation file must exist      
      LUFIL = LUFOBS
      NAMFIL = FOBS
      CHSTAT = CHOLD
      CALL OPFIL
     1    (LUFIL, LUSCRN, LULOG,  IERR, CHSTAT,  LOGOP, NAMFIL)
      
c  open observation file next day must exist      
      LUFIL = LUFOBP
      NAMFIL = FOBSP
      CHSTAT = CHOLD
      OPEN(UNIT=LUFIL, FILE= NAMFIL,
     1     STATUS= CHSTAT, IOSTAT= IERR)     
      if(IERR.ne.0)write(6,*)'NO OBS FILE FOR DAY+1'
      if(IERR.ne.0)IER(LUFIL)=1
      
c  open GPS navigation file must exist      
      LUFIL = LUNAVGPS
      NAMFIL = FNAVGPS
      CHSTAT = CHOLD
      OPEN(UNIT=LUFIL, FILE= NAMFIL,
     1     STATUS= CHSTAT, IOSTAT= IERR)     
      if(IERR.ne.0)write(6,*)'NO NAV FILE FOR GPS'
      if(IERR.ne.0)IER(LUFIL)=1

c  open navigation fil next day must exist      
      LUFIL = LUNAVPGPS
      NAMFIL = FNAVPGPS
      CHSTAT = CHOLD
      OPEN(UNIT=LUFIL, FILE= NAMFIL,
     1     STATUS= CHSTAT, IOSTAT= IERR)
      if(IERR.ne.0)write(6,*)'NO NAV FILE FOR GPS-P'
      if(IERR.ne.0)IER(LUFIL)=1

c  open GLO navigation file must exist      
      LUFIL = LUNAVGLO
      NAMFIL = FNAVGLO
      CHSTAT = CHOLD
      OPEN(UNIT=LUFIL, FILE= NAMFIL,
     1     STATUS= CHSTAT, IOSTAT= IERR)     
      if(IERR.ne.0)write(6,*)'NO NAV FILE FOR GLONASS'
      if(IERR.ne.0)IER(LUFIL)=1

c  open navigation fil next day must exist      
      LUFIL = LUNAVPGLO
      NAMFIL = FNAVPGLO
      CHSTAT = CHOLD
      OPEN(UNIT=LUFIL, FILE= NAMFIL,
     1     STATUS= CHSTAT, IOSTAT= IERR)
      if(IERR.ne.0)write(6,*)'NO NAV FILE FOR GLONASS-P'
      if(IERR.ne.0)IER(LUFIL)=1
            
c  open GAL navigation file must exist      
      LUFIL = LUNAVGAL
      NAMFIL = FNAVGAL
      CHSTAT = CHOLD
      OPEN(UNIT=LUFIL, FILE= NAMFIL,
     1     STATUS= CHSTAT, IOSTAT= IERR)
      if(IERR.ne.0)write(6,*)'NO NAV FILE FOR GAL'
      if(IERR.ne.0)IER(LUFIL)=1

c  open GAL navigation fil next day must exist      
      LUFIL = LUNAVPGAL
      NAMFIL = FNAVPGAL
      CHSTAT = CHOLD
      OPEN(UNIT=LUFIL, FILE= NAMFIL,
     1     STATUS= CHSTAT, IOSTAT= IERR)
      if(IERR.ne.0) write(6,*)'NO NAV FILE FOR GAL-P'
      if(IERR.ne.0)IER(LUFIL)=1

c   ----- -----  reading of ephemerides for GPS and GALILEO -- --------
c HEADER 
       heure=0
       jkl_gps=0
       jkl_gal=0
       
       DO 75 IFIL=20,23       ! loop on nav files for GPS and GALILEO
  779  format(40x,i4,2i2,1x,3i2)    
  6    if(IER(IFIL).eq.0)then
        read(IFIL,'(60x,a20)')comment
        if(comment.eq.'RINEX VERSION / TYPE')then
	backspace(IFIL)
         read(IFIL, '(i6)')iv
         if(iv.ne.3)then
          write(6,*)'navigation file not in version 3'
          stop
         else
          goto 6
         endif
        elseif(comment.eq.'END OF HEADER')then	
        yn=yn-2000
        call datemjd(yn,mn,dn,heure,0,0.d0,Timenav1)
         goto 7
        elseif(comment.eq.'PGM / RUN BY / DATE')then
	backspace(IFIL)
         read(IFIL, 779)yn,mn,dn,hn,minn,secn
	goto 6	
	else
         goto 6
        endif
      
  7    continue

c ephemerides 
       do 9 jk=1,15000
        read(IFIL,90,end=12)syst,prn,yn,mn,dn,hn,minn,secn,
     1                            (svclk(i),i=1,3)

        yn=yn-2000
        call datemjd(yn,mn,dn,hn,minn,dble(secn),Timenav)
	if(syst.eq.'E')then	!  GALILEO 
	jkl_gal=jkl_gal+1
	datenav_gal(jkl_gal)=timenav   
        pprn_gal(jkl_gal)=prn+200
	svclk_gal(jkl_gal,1)=svclk(1)
	svclk_gal(jkl_gal,2)=svclk(2)
	svclk_gal(jkl_gal,3)=svclk(3)
         do  j=1,7
         read(IFIL,91)(brorb_gal(jkl_gal,j,i),i=1,4)   
         enddo
	elseif(syst.eq.'G')then  
	jkl_gps =  jkl_gps +1           !  GPS
	datenav_gps(jkl_gps)=timenav  
        pprn_gps(jkl_gps)=prn
	svclk_gps(jkl_gps,1)=svclk(1)
	svclk_gps(jkl_gps,2)=svclk(2)
	svclk_gps(jkl_gps,3)=svclk(3)
         do  j=1,7
         read(IFIL,91)(brorb_gps(jkl_gps,j,i),i=1,4)   
         enddo
        endif
   9   continue
 
  12   continue
       endif
 75    CONTINUE        !end of navigation files loop
   

   90 format(a1,i2,1x,i4,5i3,3d19.12)
   91 format(4x,4d19.12)
   
c -----Reading of GLONASS ephemerides -- ----- --------

      health=0
      freqg=0
      v=0
      g=0
      sp3=0
      clock=0
      nbephem=1
      Timenavprev=0.
      jkl=0
      k=0

C Header

       do 16 ifil=24,25
       comment='NC     '
       IF(IER(ifil).eq.0)THEN
       do while (comment.ne.'END OF HEADER')
	   read(ifil,'(60x,a20)')comment

         if(comment.eq.'CORR TO SYSTEM TIME') then
	   backspace ifil
	   read(ifil,'(21X,d19.12)')tauCa
	   tauCa=-tauCa
         endif	   
       enddo
	   
      
C Read content of current file  -  Extract only 30min interval ephemerides
       do i=1,1500   !**
       read(ifil,95,end=16)prng,yg,mg,dg,hg,ming,is,
     1    TauNtmp,gammaNtmp,tktmp
 95     format(1x,i2,1x,i4,5(1x,i2),3(E19.12))
        secg=dble(is)
	yg=yg-2000

       if(ming.ne.15.and.ming.ne.45)then
	 read(IFIL,*)
	 read(IFIL,*)
	 read(IFIL,*)
	
       else

        call leapscorr(yg,mg,dg,hg,ming,secg,leapsec)
	call datemjd(yg,mg,dg,hg,ming,secg,Timenav)

	if(ifil.eq.24)jkl=hg*2+(ming-15)/30+1
	if(ifil.eq.25)jkl=(hg+24)*2+(ming-15)/30+1
	datenavglo(jkl)=Timenav
	clock(jkl,prng,1)=TauNtmp
	clock(jkl,prng,2)=gammaNtmp			
	clock(jkl,prng,3)=tktmp
        read(IFIL,'(4X,4d19.12)')sp3(jkl,prng,1),
     1        v(jkl,prng,1),g(jkl,prng,1),health(jkl,prng)
        read(IFIL,'(4X,4d19.12)')sp3(jkl,prng,2),
     1        v(jkl,prng,2),g(jkl,prng,2),freqg(jkl,prng)
        read(IFIL,'(4X,3d19.12)')sp3(jkl,prng,3),
     1        v(jkl,prng,3),g(jkl,prng,3)
       endif
       
       enddo   !**
       ENDIF
  16   continue

      nbglo=jkl
 		
      do 660 jkl=1,nbglo
         djul=datenavglo(jkl)-44244.D0
         djul=dmod(djul,7.D0)  
         tocglo(jkl)=djul*86400.D0
  660  continue

c   ----- -----  reading of the header of the observation file -- --------
      ofset=.false.
      
      call read_head_obs(nb_obs,LUFOBS,colGPSP1,colGPSP2,
     1     colGPSC1,colGLOP1,colGLOP2,colGLOC1,colGLOC2,colGALE1,
     1     colGALE5,factor,systime,ofset,tobs_first)

c   reading of the header of the file of the day after
      if(IER(LUFOBP).eq.0)
     1 call read_head_obs(n1,LUFOBP,iGP1,iGP2,
     1     iGC1,iRP1,iRP2,iRC1,irC2,iE1,
     1     iE5,ffactor,ssystime,oofset,moment)
      iobs_P=0

      if(iGP2.ne.colGPSP2.or.iGP1.ne.colGPSP1.or.
     1   iGC1.ne.colGPSC1.or.iRP1.ne.colGLOP1)iobs_P=1
      if(iRP2.ne.colGLOP2.or.iRC1.ne.colGLOC1)iobs_P=1
      if(iE1.ne.colGALE1.or.iE5.ne.colGALE5)iobs_P=1
      if(ssystime.ne.systime)iobs_P=1
      if(oofset.neqv.ofset)iobs_P=1
      do ii=1,3
      if(n1(II).ne.nb_obs(II))iobs_P=1
      do jj=1,2
       if(ffactor(ii,jj).ne.factor(ii,jj))iobs_P=1
      enddo
      enddo
             
58    ipass=0
      ipassv=0
      itestC1=0
      if(colGPSP1.eq.0)itestC1=1
      if(colGLOP1.eq.0)itestC1glo=1
      
c    last reading in the rinex file of the day     
      last=mjd+1.d0-30.d0/86400.d0
      
c    ds= security delay for time comparisons     (moment and timefirst or timestop)                
      ds=0.1d0/86400.d0           

      IF(itestC1.eq.1)THEN
      bias=99.
      
c open BIAS file if necessary

      LUFIL = LUBIAS
      NAMFIL = 'biasC1P1.dat'
      CHSTAT = CHOLD
      OPEN(UNIT=LUFIL, FILE= NAMFIL,
     1     STATUS= CHSTAT, IOSTAT= IERR)
      if(IERR.ne.0)write(6,*)'NO BIAS C1P1 FILE FOR GPS'
      if(IERR.ne.0)IER(LUFIL)=1
      
c READING OF BIASES C1-P1:

      if(IER(LUBIAS).eq.0)then
       read(lubias,801)lengC1P1,dateC1P1
 802   read(lubias,'(a30)')param
       if(param(1:1).ne.'G')goto 802
       backspace(lubias)     
       do iprn=1,40
        read(lubias,803,end=627)inp,bias(inp)
       enddo
 627   continue
      endif

 803  format(1x,i2,20x,f12.1)
 801  format(7x,a7,24x,a20)
      ENDIF
      
      
c  CREATE THE HEADER OF THE OUTPUT FILE 
628   call header(sversion,itestc1,lengC1P1,dateC1P1,luout_gps)
      write(LUOUT_gps,105)
      write(LUOUT_gps,102)

      call header(sversion,itestc1,lengC1P1,dateC1P1,luout_gal)
      write(LUOUT_gal,105)
      write(LUOUT_gal,102)

      call header(sversion,itestc1glo,lengC1P1,dateC1P1,luout_glo)
      write(LUOUT_glo,105)
      write(LUOUT_glo,102)
   
c같같같같같같같같같같같같같같같같같같같      
c    begin the run for each track
c같같같같같같같같같같같같같같같같같같같        
      iv=0
      moment=0
      
      do 100 isch=1,nsch
       if(schtime(isch).lt.tobs_first.and.tobs_first.ne.0)goto 100
       if(schtime(isch).lt.moment)goto 100
       if(iv.eq.1)goto 100

      ipk=0
      ksec=0
      timefirst=schtime(isch)
      do isat=1,45
       TRKL(isat)=0  
        do  kkk=1,26
          measiono(kkk,isat)=0.d0
          resquad(kkk,isat)=0.d0
        enddo
      enddo
       
      moment=0.d0
      timestop=schtime(isch)+13.d0/60.d0/24.d0-1.d0/86400.d0
      
c reading of observations     
 71    format(2x,i4,4(1x,i2.2),f11.7,2x,i1,i3,12x,f15.12)

       LU=LUFOBS
  89   if(iv.eq.1)LU=LUFOBP
       if(iv.eq.1.and.IER(LUFOBP).ne.0)goto 100
       read(LU,'(a500)',end=1900)string
       clkoffset=0.d0
        if(string(1:1).ne.'>')goto 89       
       read(string,'(29x,i3)')lflag
       if(lflag.eq.0)then
        read(string,71)yy,mm,dd,hh,min,sec,lflag,nnsat,clkoffset
        elseif(lflag.eq.4)then
        read(string,'(29x,2i3)')lflag,nline
        do  ii=1,nline
        read(LU,*)
        enddo
        goto 89
       elseif(lflag.eq.2)then        
        write(strerr,*)'Moving antenna--stop'
        istop = 103
        call erstop (istop, strerr, logop, luscrn, lulog)
       endif 
       goto 897

1900   iv=1             ! if end of file in the first file, needed when used with hourly files. 
       ipass=0
       goto 89   
   
 897   deltaTK=0.
       
       if(dabs(sec-60.).le.0.003)deltaTK=60.-sec
       if(dabs(sec-30.).le.0.003)deltaTK=30.-sec
       if(dabs(sec).le.0.003)deltaTK=-sec
       sec=sec+deltaTK
       if(ofset)deltaTK= clkoffset
        
       if(systime.eq.'GLO')call leapscorr(yy,mm,dd,hh,min,sec,leapsec)
       yy=yy-2000
       call datemjd(yy,mm,dd,hh,min,sec,moment)

c  Verify the DATE at first reading
       if(ipass.eq.0.)then
        ipass=1
        idate=int(moment-mjd)
	if(iv.eq.0)then
         if(idate.ne.0)
     1   then
          write(LUSCRN,*)
     1     'MJD different from dates in the file "rinex_obs" '      
          IF(LOGOP) WRITE(LULOG,*)
     1     'MJD different from dates in the file "rinex_obs" '      
          call exit(5)
         endif
	else
         if(idate.ne.1)
     1   then
          write(LUSCRN,*)
     1     'please verify the dates in the file "rinex_obs_p" '      
          IF(LOGOP) WRITE(LULOG,*)
     1     'please verify the dates in the file "rinex_obs_p" '      
          call exit(5)
         endif
	endif
       endif
       
c  reading of  observations
      do 4 iii=1,nnsat
 887    read(LU,'(a500)')string
        if(string.eq.' ')goto 887
	
	read(string,'(a1)')syst
	if(syst.eq.'G')then
	isys=1
	elseif(syst.eq.'R')then
	isys=2
	elseif(syst.eq.'E')then
	isys=3	
	elseif(syst.eq.'C')then
	isys=4	
	elseif(syst.eq.'S')then
	isys=5	
	endif
	
       read(string,73)syst,satnum(iii),
     1              (valobs(k,iii),snr(k),k=1,nb_obs(isys))
	if(syst.eq.'G')then	
        col1(iii)=colGPSP1
        if(colGPSP1.eq.0)col1(iii)=colGPSC1
        col2(iii)=colGPSP2
       elseif(syst.eq.'R')then	
        col1(iii)=colGLOP1
        if(colGLOP1.eq.0)col1(iii)=colGLOC1
        col2(iii)=colGLOP2
        if(colGLOP2.eq.0)col2(iii)=colGLOC2   !! was col1(iii) corrected into col2 on Jan 10, 2017
	satnum(iii)=satnum(iii)+100
       elseif(syst.eq.'E')then	
        col1(iii)=colGALE1
        col2(iii)=colGALE5
	satnum(iii)=satnum(iii)+200  
       elseif(syst.ne.'G'.and.syst.ne.'R'.and.syst.ne.'E')then
       	col1(iii)=0
	col2(iii)=0
	satnum(iii)=satnum(iii)+500       
       endif

   73  format(a1,i2.2,30(f14.3,1x,i1))
    4 continue
 

c construction of the 13 min vectors 
      if(ksec.eq.0.and.moment.gt.timestop+ds)goto 100
      IF(moment.ge.timefirst-ds.and.moment.le.timestop+ds)THEN
       ksec=ksec+1
       deltat(ksec)=deltaTK
       
c in order to have the same satellite list at each step
       if(ipk.eq.1)then
       p1valobs=0.d0
       p2valobs=0.d0
       do  i=1,nnsat
        do  j=1,nbsat
         if(satnum(i).eq.satn(j))then
           if(col1(i).ne.0)p1valobs(j)=valobs(col1(i),i)
           if(col2(i).ne.0)p2valobs(j)=valobs(col2(i),i)
         endif
        enddo
        enddo
      endif
	
c starting a new track
      if(ipk.eq.0)then
       nbsat=nnsat
       do  i=1,nbsat
        satn(i)=satnum(i)
	p1valobs(i)=valobs(col1(i),i)
	p2valobs(i)=valobs(col2(i),i)
       enddo
       ipk=1
      endif
            
c construction of the vectors 
       epoch(ksec)=moment
       do 5 i=1,nbsat
        if(p1valobs(i).ne.0.d0.and.p2valobs(i).ne.0.d0)then
	 if(satn(i).lt.100)then						!GPS
          if(itestC1.eq.1) p1valobs(i)=
     1             p1valobs(i)+bias(satn(i))*vitlum*1.0d-9
           resquad(ksec,i)=k1(1)*(p1valobs(i)/factor(1,1)-recdelP1)
     1                 -k2(1)*(p2valobs(i)/factor(1,2)-recdelP2)
           measiono(ksec,i)=(p1valobs(i)/factor(1,1)-recdelP1)
     1                      -resquad(ksec,i)
	 elseif(satn(i).lt.200)then					! GLONASS
           resquad(ksec,i)=k1(2)*(p1valobs(i)/factor(2,1)-recdelP1g)
     1                 -k2(2)*(p2valobs(i)/factor(2,2)-recdelP2g)            
           measiono(ksec,i)=(p1valobs(i)/factor(2,1)-recdelP1g)
     1                     -resquad(ksec,i)
         elseif(satn(i).lt.300)then					! Galileo
	  resquad(ksec,i)=k1(3)*(p1valobs(i)/factor(3,1)-recdelE1)
     1                 -k2(3)*(p2valobs(i)/factor(3,2)-recdelE5)           
          measiono(ksec,i)=(p1valobs(i)/factor(3,1)-recdelE1)
     1                     -resquad(ksec,i)
	 endif  
         TRKL(i)=TRKL(i)+1
         ic(i)=ic(i)+1

        endif
   5   continue
      if(ksec.eq.26)goto 20
      ENDIF 
        
      if(moment.lt.timestop-ds)goto 89      

c considering only GPS and/or GALILEO
20    iprnmax=300
      iprnmin=0
    
      do 150 isat=1,nbsat
        
c     keep only desired satellites 
       IF(TRKL(isat).eq.26.and.satn(isat).lt.iprnmax.
     1                    and.satn(isat).gt.iprnmin)THEN
            
c  choice of the good ephemeris
      jklsave=0
      tgd=0.
      fr=0

      IF(satn(isat).lt.100)THEN   			  ! GPS
         difsave=100000.d0
	 WN=int((epoch(13)-44244.d0)/7.d0)
         djul=epoch(13)-44244.D0
         djul=dmod(djul,7.D0)  
         trc=djul*86400.D0 
       do jkl=1,jkl_gps
         if(pprn_gps(jkl).eq.satn(isat).and.brorb_gps(jkl,6,2).eq.0)then
          diftime=schtime(isch)-datenav_gps(jkl)
	  if(dabs(diftime).lt.dabs(difsave).and.diftime.le.0.and.
     1                   dabs(diftime).le.0.1)then
          jklsave=jkl	
	  difsave=diftime
          djul=datenav_gps(jkl)-44244.D0
          endif	
         endif	
      enddo
         if(jklsave.eq.0)goto 150
        svclk(1)=svclk_gps(jklsave,1)
        svclk(2)=svclk_gps(jklsave,2)
        svclk(3)=svclk_gps(jklsave,3)	
       brorb=brorb_gps
       tgd=brorb(jklsave,6,3)
       IOE=int(brorb(jklsave,1,1))	 
	
       ELSEIF(satn(isat).lt.200)THEN  			  !GLONASS      
       
       	 do jkl=1,1000
           if(jkl==1) then
	        if(datenavglo(1).gt.schtime(isch)+6.5/24./60.)goto 6600
	    endif
           if(datenavglo(jkl).gt.schtime(isch)+6.5/24./60.)cycle
	   if((schtime(isch)+6.5/24./60.)-datenavglo(jkl).
     1         le.(30.d0/(24*60.d0))) goto 6600
	 enddo
	 if(jkl.eq.1001)goto 150
 6600    jklsave=jkl
         IOE=jklsave*2
	 fr=freqg(jklsave,satn(isat)-100)
 
       ELSEIF(satn(isat).gt.200)THEN    		 !Galileo
       
         ict=258	! E5a
         difsave=100000.d0
         WN=int((epoch(13)-44244.d0)/7.d0)
         djul=epoch(13)-44244.D0
         djul=dmod(djul,7.D0)  
         trc=djul*86400.D0 
       do 110 jkl=1,jkl_gal
         if(pprn_gal(jkl).eq.satn(isat)
     1       .and.brorb_gal(jkl,5,2).eq.ict)then
     
c   Discard the block either DVS=1 or HS=1       
          call int2bin(int(brorb_gal(jkl,6,2)),binr)   	
c	  if(binr(2,2).eq.1)goto 110
     
          diftime=schtime(isch)-datenav_gal(jkl)
          if(dabs(diftime).lt.dabs(difsave).and.diftime.gt.0.and.
     1                   dabs(diftime).le.0.167)then			! validity of nav message: 4 h

          difsave=diftime
	  jklsave=jkl
          djul=datenav_gal(jkl)-44244.D0
	  endif
         endif	
 110   continue
         if(jklsave.eq.0)goto 150
       
        svclk(1)=svclk_gal(jklsave,1)
        svclk(2)=svclk_gal(jklsave,2)
        svclk(3)=svclk_gal(jklsave,3)
	brorb=brorb_gal
        tgd=brorb(jklsave,6,3)
        IOE=int(brorb(jklsave,1,1))
		 
       ENDIF

 	 jkl=jklsave 
         djul=dmod(djul,7.D0)  
         toc=djul*86400.D0

c computation of the corrected 26 pseudoranges for the 13 min track       

  659    do 120 ksec=1,TRKL(isat)       
         if(resquad(ksec,isat).ne.0.)then
         iii=0
  	 jjj=0
         refgps(ksec)=0.d0
	 refsv(ksec)=0.d0
         corprev=0.d0
         djul=epoch(ksec)-44244.D0
         djul=dmod(djul,7.D0)  
         trc=djul*86400.D0           			 ! reception time	 
	 
c geometric distance	 
  54     if(iii.eq.0)
     1     tau=resquad(ksec,isat)/vitlum-refgps(ksec)*1.0d-10
         if(iii.eq.1)tau=corrgeom/vitlum-refgps(ksec)*1.0d-10
         ttr=trc-tau      				   ! emission time
         
	 if (satn(isat).lt.100.or.satn(isat).gt.200)then	! GPS/Galileo
	  call satpos(brorb,Ttr,Trel,X,status,jkl,
     1         LUSCRN, LULOG, LOGOP,satn(isat))
          if(.not.status)goto 150
          alpha=tau*WE				! Sagnac (WE=Earth rotation)
          xs(1)=x(1)*dcos(alpha)+x(2)*dsin(alpha)
          xs(2)=-x(1)*dsin(alpha)+x(2)*dcos(alpha)
          xs(3)=x(3)
          corrgeom=dsqrt((antpos(1)-Xs(1))**2
     1              +(antpos(2)-Xs(2))**2+(antpos(3)-Xs(3))**2)
          if(dabs(corrgeom-corprev).lt.1.0d-5)goto 64
          corprev=corrgeom
          iii=1
          goto 54
 64       continue
         else							!  GLONASS
c  Xs and vit by R.Kutta integration
            call integ(ttr,tocglo(JKL),sp3,v,g,tau,Xs,
     1      satn(isat)-100,vit,exists,LULOG,jkl) 
            if(health(jkl,satn(isat)-100).ne.0.)goto 150
            if(.not.exists)goto 150     

            corrgeom=dsqrt((antpos(1)-Xs(1))**2
     1       +(antpos(2)-Xs(2))**2+(antpos(3)-Xs(3))**2)
            if(dabs(corrgeom-corprev).lt.1.0d-5)goto 164
            corprev=corrgeom
            iii=1
            goto 54
 164        continue 
C  no Trel for GLONASS as already included in the satellite clock
	   trel=0

	 endif
	            
c  measured iono (unit= m)
         iono(ksec)=measiono(ksec,isat) 	            
	  
c  tropo (unit= m)
         call XYZLLA(antpos,lla)
         call AZELd(Xs,antpos,lla,Az,El)
	 if(El.lt.0.174533)goto 150
         mf = 1.d0/(dsin(El)+0.00143d0/(dtan(El)+0.0455d0))
         if(lla(3).le.1000.)then
            zpd = 2162.d0 + Ns*(1.-lla(3)/1000.d0) 
     1                 + 0.5d0*DeltaN*(1.-(lla(3)/1000.d0)**2)
         else	
  	    zpd = 732.d0 - 8.d0 * (NS+DeltaN)/NSlog *
     1         (dexp(-NSlog)-dexp(0.125d0*(1.d0-lla(3)/1000.d0)*NSlog))
         endif
         tropo(ksec)=mf*zpd/1000.d0      
            
	    
c  corrected pseudorange            
         pcor = resquad(ksec,isat) + trel*vitlum - tropo(ksec)  
         CR = Pcor-corrgeom
            
c  computation of REFSV (unit= 0.1 ns)
         if(.not.ofset)refsv(ksec)= (CR/vitlum)*1.0d+10- cabdel+refdel
         if(ofset)refsv(ksec)= (CR/vitlum)*1.0d+10 
     1           - cabdel + refdel  +deltaT(ksec)*1.0d+10
           
c  satellite clock correction (unit= 0.1 ns)
	 if (satn(isat).lt.100.or.satn(isat).gt.200)then	! GPS/Galileo
         tclock=ttr-toc
         if(tclock.gt.302400.d0)tclock=tclock-604800.d0
         if(tclock.lt.-302400.d0)tclock=tclock+604800.d0
         clocksat=(svclk(1)+svclk(2)*tclock
     1        +svclk(3)*tclock*tclock)*1.0d+10
         else                                           ! GLONASS
            djul=datenavglo(jkl)-44244.D0
            djul=dmod(djul,7.D0)
            tocg=djul*86400.D0
           tclock=ttr-tocg
           if(tclock.gt.302400.d0)tclock=tclock-604800.d0
           if(tclock.lt.-302400.d0)tclock=tclock+604800.d0
           clocksat=clock(jkl,satn(isat)-100,1)*1.0d+10+
     1          clock(jkl,satn(isat)-100,2)*tclock*1.0d+10

         endif

c  computation of  REFSYS (called refgps for historical reason) (unit= 0.1 ns)
         refgps(ksec)=refsv(ksec)+clocksat   
	 if(dabs(refgps(ksec)).gt.1d+4.and.jjj.eq.0)then
	 jjj=1
	 goto 54
	 endif
      
       else
        refgps(ksec)=0.
        refsv(ksec)=0.
        tropo(ksec)=0.
        iono(ksec)=0.
       endif
c      end of do loop   
 120   continue
      

c  AZth et ELv at the midpoint of the track
       corprev=0.
       iii=0
       djul=schtime(isch)+6.5/24./60.-44244.0
     1       -dble(ileaps)/86400.+dble(leapsec)/86400.
       djul=dmod(djul,7.D0)  
       trc=djul*86400.D0
  53   if(iii.eq.0)tau=resquad(13,isat)/vitlum
       if(iii.eq.1)tau=corrgeom/vitlum
       ttr=trc-tau
	 if (satn(isat).lt.100.or.satn(isat).gt.200)then	! GPS/Galileo
       call satpos(brorb,Ttr,Trel,X,status,jkl,
     1         LUSCRN, LULOG, LOGOP,satn(isat))        
       alpha=tau*WE
       xs(1)=x(1)*dcos(alpha)+x(2)*dsin(alpha)
       xs(2)=-x(1)*dsin(alpha)+x(2)*dcos(alpha)
       xs(3)=x(3)
       corrgeom=dsqrt((antpos(1)-Xs(1))**2
     1          +(antpos(2)-Xs(2))**2+(antpos(3)-Xs(3))**2)
       if(dabs(corrgeom-corprev).lt.1.0d-5)goto 63
       corprev=corrgeom
       iii=1
       goto 53
 63    continue
        else							! GLONASS
	call integ(ttr,tocglo(JKL),sp3,v,g,tau,Xs,
     1 satn(isat)-100,vit,exists,LULOG,jkl)    
        if(.not.exists)cycle    
       corrgeom=dsqrt((antpos(1)-Xs(1))**2
     1          +(antpos(2)-Xs(2))**2+(antpos(3)-Xs(3))**2)
       if(dabs(corrgeom-corprev).lt.1.0d-5)goto 763
       corprev=corrgeom
       iii=1
       goto 53
 763    continue
       endif
        
         call AZELd(Xs,antpos,lla,Az,El)
       el=el/pi*1800.d0
       az=az/pi*1800.d0
       elv=nint(el)
       azth=nint(az)
       if(el.lt.100)goto 150
           
c  linear fits
       call fitlin(epoch,trkl(isat),refsv,refsvmil,srsv,stdv)
       call fitlin(epoch,trkl(isat),refgps,refgpsmil,srgps,dsg)
       call fitlin(epoch,trkl(isat),tropo,mdtr,smdt,stdv)
       call fitlin(epoch,trkl(isat),iono,mdio,smdi,isg)
            
c  iono and tropo in 0.1 ns
       mdtr=mdtr/vitlum*1.0d+10
       mdio=(mdio/vitlum- tgd)*1.0d+10
       isg=isg/vitlum*1.0d+10
      
c  drifts in  0.1ps/s
       srsv=srsv*1000.d0/86400.d0
       srgps=srgps*1000.d0/86400.d0
       smdt=smdt/vitlum*1.0d+10 *1000.d0 /86400.d0
       smdi=smdi/vitlum*1.0d+10 *1000.d0 /86400.d0     

       class='FF'
               
c  writing of results241
  154  trkl(isat)=trkl(isat)*30
       secsch=0
       mjdtrk=int(mjd)
       ick=0

        if(satn(isat).lt.100)then
	syst='G'
	frc='L3P'
	ifile=LUOUT_GPS
	elseif(satn(isat).lt.200)then
	syst='R'
	frc='L3P'
	satn(isat)=satn(isat)-100
	ifile=LUOUT_GLO
	elseif(satn(isat).gt.200)then
	syst='E'
	frc='L3E'
	satn(isat)=satn(isat)-200
	ifile=LUOUT_GAL
	endif
	
       write(name,104)syst,satn(isat),class,mjdtrk,schh(isch),
     1  schmin(isch),secsch,trkl(isat),elv,azth,nint(refsvmil),
     2  nint(srsv),nint(refgpsmil),nint(srgps),nint(dsg),
     3  ioe,nint(mdtr),nint(smdt),nint(mdio),nint(smdi),
     4  nint(mdio),nint(smdi),nint(isg),fr,hc,frc
        do  kk=1,125
        ick=ick+ichar(name(kk:kk))
        enddo
        ck=mod(ick,256)
        
        write(ifile,103)syst,satn(isat),class,mjdtrk,schh(isch),
     1  schmin(isch),secsch,trkl(isat),elv,azth,nint(refsvmil),
     1  nint(srsv),nint(refgpsmil),nint(srgps),nint(dsg),
     1  ioe,nint(mdtr),nint(smdt),nint(mdio),nint(smdi),
     1  nint(mdio),nint(smdi),nint(isg),fr,hc,frc,ck

 105  format('SAT',1x,'CL',2x,'MJD',2x,'STTIME',1x,'TRKL',1x,
     1 'ELV',1x,'AZTH',3x,'REFSV',6x,'SRSV',5x,'REFSYS',4x,'SRSYS',
     2  2x,'DSG ','IOE ','MDTR ','SMDT ','MDIO ',
     3 'SMDI MSIO SMSI ISG FR HC FRC CK')
 102  format(13x,'hhmmss',2x,'s',2x,'.1dg ','.1dg',4x,'.1ns',5x,
     1  '.1ps/s',5x,'.1ns',4x,'.1ps/s ','.1ns',
     2  5x,'.1ns.1ps/s.1ns.1ps/s.1ns.1ps/s.1ns')
 104  format(a1,i2,1x,a2,1x,i5,1x,3i2.2,1x,i4,1x,i3,1x,i4,sp,1x,i11,
     1   1x,i6,ss,1x,i11,sp,1x,i6,1x,ss,i4,1x,i3,1x,i4,1x,
     2   sp,i4,ss,1x,i4,1x,sp,i4,ss,1x,i4,1x,i4,1x,i3,1x,i2,
     3  1x,i2,1x,a3,1x)    
 103  format(a1,i2.2,1x,a2,1x,i5,1x,3i2.2,1x,i4,1x,i3.3,1x,i4,sp,1x,
     1   i11,1x,i6,ss,1x,i11,sp,1x,i6,1x,ss,i4,1x,i3,1x,i4,1x,
     2   sp,i4,ss,1x,i4,1x,sp,i4,ss,1x,i4,1x,i4,1x,i3,1x,i2,
     3  1x,i2,1x,a3,1x,z2.2)    

      ENDIF   

  150 continue
  
  100 continue
      
      WRITE(LUSCRN,*) ' Program executed successfully.'
      IF(LOGOP) WRITE(LULOG,*) ' Program executed successfully.'

c     close all files
      LUFIL = LUPARA
            CALL CLFIL
     1    (LUFIL, LUSCRN, LULOG, IERR, LOGOP)

      LUFIL = LUOUT 
            CALL CLFIL
     1    (LUFIL, LUSCRN, LULOG, IERR, LOGOP)

      LUFIL = LULOG 
            CALL CLFIL
     1    (LUFIL, LUSCRN, LULOG, IERR, LOGOP)

      LUFIL = LUINP 
            CALL CLFIL
     1    (LUFIL, LUSCRN, LULOG, IERR, LOGOP)

      LUFIL = LUFOBS
            CALL CLFIL
     1    (LUFIL, LUSCRN, LULOG, IERR, LOGOP)

      LUFIL = LUFOBP 
            CALL CLFIL
     1    (LUFIL, LUSCRN, LULOG, IERR, LOGOP)

      LUFIL = LUNAV 
            CALL CLFIL
     1    (LUFIL, LUSCRN, LULOG, IERR, LOGOP)

      LUFIL = LUNAVP
            CALL CLFIL
     1    (LUFIL, LUSCRN, LULOG, IERR, LOGOP)


c     stop the program
c     normal stop when no errors were detected
999   stop 

      end

c     ********************************

      subroutine satpos
     1 (brorb, Ttr, Trel, X, status, jkl, LUSCRN, LULOG, LOGOP,iprn)

C       brorb : broadcast ephemeris
C       Ttr : satellite GPS time
C       Trel: relativistic correction term
C       X   : satellite position

      real*8 brorb(15000,7,4),X(3)
      real*8 Ttr,Trel,M0, dn, ec, A, W0, i0, w, Wdot,
     |       Cuc, Cus, Crc, Crs, Cic, Cis, Toe, Idot,
     |       T, n0, n, M, E, Eold, snu, cnu, nu, phi, 
     |       du, dr, di, u, r, i, Xdash, Ydash, Wc,pi,mu,Wedot,F,
     |       mugal,fgal,mugps,fgps
      
      logical status, LOGOP      
      integer it,iprn
      integer*4
     1  LUSCRN, LULOG

      data pi/3.1415926535898d0/,mugps/3.986005d+14/,
     |       Wedot/7.2921151467d-5/,Fgps/-4.442807633d-10/,
     |       mugal/3.986004418d+14/,Fgal/-4.442807309d-10/
     
      if(iprn.lt.100)then
      F=Fgps
      mu=mugps
      elseif(iprn.gt.200)then
      F=Fgal
      mu=mugal
      endif
     
      
      
     
C     for RINEX files

      M0  = brorb(jkl,1,4)
      dn  = brorb(jkl,1,3)
      W0  = brorb(jkl,3,3)
      i0  = brorb(jkl,4,1)
      w   = brorb(jkl,4,3)
      Wdot= brorb(jkl,4,4)
      idot= brorb(jkl,5,1) 

C     **************
       
      Crs = brorb(jkl,1,2)
      Cuc = brorb(jkl,2,1)
      ec  = brorb(jkl,2,2)
      Cus = brorb(jkl,2,3)
      A   = brorb(jkl,2,4)*brorb(jkl,2,4)    
      Toe = brorb(jkl,3,1)
      Cic = brorb(jkl,3,2)
      Cis = brorb(jkl,3,4)
      Crc = brorb(jkl,4,2) 
      
      status=.true.
      it=0
       
      T= Ttr - Toe
      
      if (T.gt.302400.d0)  T = T - 604800.d0
      if (T.lt.-302400.d0)  T = T + 604800.d0
      n0 = dsqrt(mu / (A*A*A))
      n = n0 + dn
      M = M0 + n*T
      E = M 

   10 it=it+1
      Eold = E
      E = M + ec * dsin(E)
c      write(6,66)it,Eold,E,M,ec
66    format(i3,4f18.5      )
      if ((it.eq.10).or.(dabs(E - Eold).le.(1.0d-8))) goto 12
      goto 10
   12 if (it.eq.10)then
        status = .false. 
        write(LUSCRN,*)'no convergence for E'
        if(LOGOP) write(LULOG,*)'no convergence for E'
        goto 15
      endif     
      snu = dsqrt(1.d0 - ec*ec) * dsin(E) / (1.d0 - ec*dcos(E))
      cnu = (dcos(E) - ec) / (1.d0 - ec*dcos(E))
      nu=datan2(snu,cnu)
      
      phi = nu + w

      du = Cuc*dcos(2.d0*phi) + Cus*dsin(2.d0*phi)
      dr = Crc*dcos(2.d0*phi) + Crs*dsin(2.d0*phi)
      di = Cic*dcos(2.d0*phi) + Cis*dsin(2.d0*phi)
      
      u = phi + du
      r = A*(1.d0 - ec*dcos(E)) + dr
      i = i0 + idot*T +di

      Xdash = r*dcos(u)
      Ydash = r*dsin(u)
     

      Wc= W0 + (Wdot - Wedot)*T - Wedot*Toe
                        
      X(1) = Xdash*dcos(Wc) - Ydash*dcos(i)*dsin(Wc)
      X(2) = Xdash*dsin(Wc) + Ydash*dcos(i)*dcos(Wc)
      X(3) = Ydash*dsin(i)

C     relativistic correction term
      Trel = F * ec * brorb(jkl,2,4) * dsin(E) 

   15 return 
      end

c     ********************************

      subroutine AZELd(Xsat,Xr,lla, az,el)


         REAL*8  X1,Y1,Z1,X2,Y2,Z2, PHI,PLAM, AZ,EL,DIST,  
     1   RLAT , RLON ,PI, SRLON , CRLON,  SRLAT , CRLAT
     1   DX,DY,DZ, DU,DV,DW,xsat(3),xr(3),lla(3)

!@ Initial values

      x2=xsat(1)
      y2=xsat(2)
      z2=xsat(3)
      
      x1=xr(1)
      y1=xr(2)
      z1=xr(3)

      phi=lla(1) 
      plam =lla(2)

!@@@@@@@@@@@@@@@@@@@@@
!@ Compute elevation @
!@@@@@@@@@@@@@@@@@@@@@
      
!
      PI = (DATAN(1.d0)) * 4.D0
      RLAT = PHI
      RLON = PLAM
      SRLAT= DSIN(RLAT)
      CRLAT= DCOS(RLAT)
      SRLON= DSIN(RLON)
      CRLON= DCOS(RLON)
!
       DX = X2-X1
       DY = Y2-Y1
       DZ = Z2-Z1

       DU = -SRLAT*CRLON*DX - SRLAT*SRLON*DY + CRLAT*DZ
       DV = -SRLON*DX + CRLON*DY
       DW = CRLAT*CRLON*DX + CRLAT*SRLON*DY + SRLAT*DZ

        DIST = DSQRT(DX*DX+DY*DY+DZ*DZ)
        AZ = DATAN2(DV,DU)
        EL = DASIN(DW/DIST)
        if (el.lt.0.D0) el=el+pi/2.D0
        if (az.lt.0.D0) az=az+2*pi  

       end subroutine AZELD



c     ********************************

      subroutine XYZLLA(xi,lla)
      
      implicit real*8(d)
      
      real*8 xi(3),lla(3),p,t,st,ct,rn,a,b,pi
      
      pi=3.1415926535898d0 
      a=6378137.0d0
      da=a
      df=0.00335281066d0
      dx=xi(1)
      dy=xi(2)
      dz=xi(3)
      DA2=a*a                                                           
      DB=DA*(1.D0-DF)                                                     
      DE=DSQRT((DA2-DB*DB)/DA2)                                           
      DE2=DE*DE                                                           
      DL=(DATAN(DY/DX))                                              
      IF(DL)6,3,3                                                         
    3 IF(DX)5,8,8                                                         
    5 DL=DL+pi                                                        
      GO TO 8                                                             
    6 IF(DX)5,8,7                                                         
    7 DL=DL+2.*pi                                                        
    8 DXY=DSQRT(DX**2+DY**2)                                              
      DP=DATAN(DZ/DXY)                                                    
      L=1                                                                 
   10 DN=DA/DSQRT(1.D0-(DE*DSIN(DP))**2)                                  
      GO TO (11,20),L                                                     
   11 DZP=DZ+DE2*DN*DSIN(DP)                                              
      DPH=DATAN(DZP/DXY)                                                  
      DRES=DABS(DP-DPH)                                    
      IF(DRES.LT.0.0001D0)GO TO 15                                        
      DP=DPH                                                              
      GO TO 10                                                            
   15 L=2                                                                 
      DCP=DCOS(DPH)                                                       
      DCL=DCOS(DL)                                                   
      GO TO 10                                                            
   20 DH=(DX/(DCP*DCL))-DN  
      lla(1)=dph
      lla(2)=dl
      lla(3)=dh                                              
      RETURN                                                              
      END                    


c     ********************************

      subroutine fitlin (time,len,vect,valmil,deriv,stdv)
      
      real*8 vect(26),valmil,deriv,stdv,abscisse,time(26),
     1    s,sx,sxx,sy,sxy,delta,a,b,ord,cal,estim,tmil,t1
      common/leap/leapsec

      stdv=0.d0
      s=0.d0
      sx=0.d0
      sxy=0.d0
      sxx=0.d0
      sy=0.d0
      cal=1.0d3
      tmil=time(1)+6.5d0/60.d0/24.d0
     1       -30.d0/86400.d0+dble(leapsec)/86400.d0
      t1=time(1)
      do 15 j=1,len
        abscisse = time(j)-t1
        if(vect(j).ne.0.d0) then
         ord=(vect(j)-vect(1))/cal
         s = s+ 1.0d0
         sx = sx + abscisse
         sxx = sxx +  abscisse * abscisse
         sy = sy + ord
         sxy = sxy + ord * abscisse
        endif
  15  continue
      if( s .gt. 1.d0 )then
        delta = (sxx * s) - (sx * sx)
        A = ((sy * sxx) - (sx * sxy)) / delta
        B = ((s * sxy) - (sx * sy)) / delta
      else
        A = 98.0d0
        B = 98.0d0
      endif
      
      if(A.ne.98.0d0)then      
        valmil=(a+b*(tmil-t1))*cal + vect(1)
        deriv=B*cal
        stdv=0.
        do j=1,len
          estim=(a+b*(time(j)-t1))*cal+ vect(1)
          stdv=stdv+(vect(j)-estim)**2
        enddo
        stdv=dsqrt(stdv/(s-2.))
      else
        valmil=9999999
        deriv=99999999
        stdv=9999999
      endif      
      
      return
      end

c     ********************************

      subroutine datemjd ( IAN, MO, JO, hh,min,sec, dpeju )
c     output : dpeju = mean julian day

      IMPLICIT REAL*8( d )
      REAL*8 sec
      integer hh, NJAO( 12 ), NJAB( 12 ), NJAR( 12 )
      
      DATA NJAO /0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334/
      DATA NJAB /0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335/
      DATA NJAR /0, 31, 59, 90, 120, 151, 181, 212, 243, 263, 294, 324/
      
      NOPT=1     
      dfj= (hh + dble(min)/60. + sec/3600.) / 24.
      ian=ian+2000

    1 I1  =  IAN
      IF ( MO .LE. 2 ) I1  =  IAN - 1

      I2  =  MO  +  1
      IF ( MO .LE. 2 ) I2  =  MO  +  13

      D1 = JO + DFJ

      IF ( I1 .GE. 0 ) IF ( IAN - 1582 ) 2, 3, 5

c     ante christum natum date
      DPEJU =   dINT( 365.25*dble(I1) - 0.75 ) 
     *        + dINT( 30.6001*dble(I2) )+ D1 + 1720994.5
     *         -2400000.5d0
      RETURN

c     date in the julian calendar ( -> 4 oct. 1582 )
    2 DPEJU = dINT( 365.25*dble(I1) ) 
     *  + dINT( 30.6001*dble(I2) ) + D1 + 1720994.5
     *        -2400000.5d0
      RETURN

c     date in 1582 
    3 IF ( MO - 10 ) 2, 4, 5
     
c     date in oct 1582
    4 IF ( NOPT .EQ. 2 ) GO TO 6
c     date in oct 1582 ( yyyy, mm, dd )
      IF ( JO .LT. 15 ) GO TO 2

c     date in the gregorian calendar ( 15 oct 1582 -> )
    5 DPEJU = dINT( 365.25*dble(I1) ) + dINT( 30.6001*dble(I2) ) 
     *         + D1 + 1720994.5 + dble(2-I1/100 + I1/400)
     *        -2400000.5d0
      RETURN

c     date in oct 1582 ( yyyy, ddd )
    6 continue

      END
      

c     *****************************************************************
c     Subroutine for calculatin MJD values
c     *****************************************************************
      subroutine decal(mjd,mjddep,he,mi,schh,schmin)
      
      integer hdecal,schmin,schh,he,mi,t0,tt
      real*8 mjd,mjddep
      
      t0=he*60+mi
      tt=t0-4*int(mjd-mjddep)
      do 15 while (tt.lt.0)
        tt=tt+1436
  15  continue
      schh= int(tt/60)
      schmin=mod(tt,60)
     
      return
      end

c     *****************************************************************
c     Subroutine for handling the time
c     *****************************************************************
      subroutine classement(schh,schm,schtime,nn)
      
      real*8 schtime(nn),provt(nn)
      integer schh(nn),schm(nn),provh(nn),provm(nn)
      
      ndep=1
      do 1000 i=2,nn
       if(schtime(i).lt.schtime(ndep)) ndep=i
 1000 continue
  
      do 1100 i=1,nn-ndep+1
        provt(i)=schtime(ndep+i-1)
        provh(i)=schh(ndep+i-1)
        provm(i)=schm(ndep+i-1)
 1100 continue

      do 1200 i=nn-ndep+2,nn
        provt(i)=schtime(i-(nn-ndep+1))
        provh(i)=schh(i-(nn-ndep+1))
        provm(i)=schm(i-(nn-ndep+1))
 1200 continue
    
      do 1300 i=1,nn
        schtime(i)=provt(i)
        schh(i)=provh(i)
        schm(i)=provm(i)
 1300 continue

      return
      end

c     *****************************************************************
c       Subroutine to create the header of the CGGTTS file
c     *****************************************************************
      subroutine HEADER(sversion,itestc1,lengC1P1,dateC1P1,LUOUT)

      INTEGER*4
     1    LUOUT, ICK, K, I

      real*8 antpos(3),recdelP1,recdelP2,cabdel,refdel,m_ns,recdelE5,
     1       recdelE1,recdelP1g,recdelP2g
      integer ck,ch
      character*30 LAB,REF,RCVR,REVDATE,COMMENTS
      character*80 abstxt
      character*7 lengC1P1,sversion,ADEL
      character*20 dateC1P1,calID
      character*1 C1p1,E5
     
      common/headint/ ch,ICHOICE_SYS
      common/headrel8/ antpos, recdelP1,recdelP2,recdelE1,
     1                 recdelE5,recdelP1g,recdelP2g,cabdel,refdel
      common/headchar/ LAB,REF,RCVR,REVDATE,COMMENTS,calID
      common/leap/leapsec
             
      m_ns=299792458.0d0*1.0d-9
         
      c1p1='P'
      if(itestC1.eq.1)c1p1='C'
      E5='a'
      
      write(LUOUT,612)
      write(LUOUT,613)revdate
      write(LUOUT,614)rcvr,sversion
      if(itestC1.eq.0)write(LUOUT,615)ch
      if(itestC1.eq.1.and.luout.eq.15)
     1   write(LUOUT,6155)ch,lengC1P1,dateC1P1
      if(itestC1.eq.1.and.luout.eq.16)
     1   write(LUOUT,615)ch
      write(LUOUT,616)rcvr
      write(LUOUT,601)LAB
      write(LUOUT,602)antpos(1)
      write(LUOUT,603)antpos(2)
      write(LUOUT,604)antpos(3)
      if(luout.ne.16)write(LUOUT,610)
      if(luout.eq.16)write(LUOUT,910)
      write(LUOUT,611)comments
      
      ADEL='INT DLY'
      if(cabdel.eq.0.and.refdel.ne.0)ADEL='SYS DLY'
      if(cabdel.eq.0.and.refdel.eq.0)ADEL='TOT DLY'
      if(luout.eq.15)
     1  write(LUOUT,605)ADEL,recdelP1/m_ns,c1p1,recdelP2/m_ns,calID
      if(luout.eq.14)
     1  write(LUOUT,625)ADEL,recdelE1/m_ns,recdelE5/m_ns,E5,calID
      if(luout.eq.16)
     1  write(LUOUT,705)ADEL,recdelP1g/m_ns,c1p1,recdelP2g/m_ns,calID
      
      if(cabdel.ne.0)write(LUOUT,606)cabdel/10.0d0
      if(refdel.ne.0)write(LUOUT,607)refdel/10.0d0
      write(LUOUT,608)REF
      
      
 612  format('CGGTTS     GENERIC DATA FORMAT VERSION = 2E')
 613  format('REV DATE = ',a30)
 614  format('RCVR = ',a30,'R2CGGTTS v',a)
 615  format('CH = ',i2,' ')
 6155 format('CH = ',i2,'           ','C1P1bias: ',a7,2x,a20)
 616  format('IMS = ',A30)
 601  format('LAB = ',A30)
 602  format('X =',sp,f12.2,' m ')
 603  format('Y =',sp,f12.2,' m ')
 604  format('Z =',sp,f12.2,' m ')
 605  format(A7,' = ',f6.1,' ns (GPS ',a1,'1), ',
     1        f6.1,' ns (GPS P2)',5x,'CAL_ID = ',A20)
 705  format(A7,' = ',f6.1,' ns (GLO ',a1,'1), ',
     1        f6.1,' ns (GLO P2)',5x,'CAL_ID = ',A20)
625   format(A7,' = ',f6.1,' ns (GAL E1), ',
     1        f6.1,' ns (GAL E5',a1,')',5x,'CAL_ID = ',A20)
 606  format('CAB DLY = ',f6.1,' ns ')
 607  format('REF DLY = ',f6.1,' ns')
 608  format('REF = ',a30)
 609  format('CKSUM = ',z2.2)
 610  format('FRAME = ITRF')
 910  format('FRAME = ITRF, PZ-90->ITRF Dx = 0.0 m, Dy = 0.0 m,',
     1 ' Dz = 0.0 m, ds = 0.0, Rx = 0.0, Ry = 0.0, Rz = 0.000000')
 611  format('COMMENTS = ',a30)
 
      rewind(LUOUT)
      ick=0
      
      do  k=1,25
        read(LUOUT,'(a50)',end=1)abstxt
        do 15 i=1,50
          if(abstxt(i:i+30).ne.'   ')
     1        ick=ick+ichar(abstxt(i:i))
 15     continue
      enddo
 1    continue
          
      abstxt='CKSUM = '
      do 16 i=1,50
        if(abstxt(i:i+2).eq.'   ')
     1  then
          goto 9
        else

          if(i .NE. 1) then            
            ick=ick+ichar(abstxt(i:i))
          endif          
        endif
 16   continue
 
 9    ick=ick+32
      ck=mod(ick,256)
      
      BACKSPACE(LUOUT)
      write(LUOUT,609)ck           
      write(LUOUT,*)
      
      return
      end

      
c     *****************************************************************
c       Subroutine to open a file and perform proper validation for file handling
c     *****************************************************************
      subroutine OPFIL
     1    (LUFIL, LUSCRN, LULOG,  IERR, CHSTAT, LOGOP, NAMFIL)

      IMPLICIT NONE

      INTEGER*4
     1    LUFIL,  LUSCRN, LULOG,  IERR, ISTOP

      CHARACTER*7
     1   CHSTAT

      CHARACTER*1200
     1      STRERR
           
      LOGICAL
     1    LOGOP
             
      CHARACTER*512
     1   NAMFIL

c     open file namfil & error trapping
      OPEN(UNIT=LUFIL, FILE= NAMFIL,
     1     STATUS= CHSTAT, IOSTAT= IERR)

c     open file successful?     
      IF(IERR .NE. 0)THEN         
         WRITE(STRERR,1000)  LUFIL, IERR, NAMFIL         
             
 1000    FORMAT('Emergency stop\n',
     1          ' fail to open file unit',I4,' ierr',i4,'\n',
     2          1x,A512)
         ISTOP = 106         
         CALL erstop (istop, strerr, logop, luscrn, lulog)
      ENDIF

      RETURN
      END

c     *****************************************************************
c       Subroutine to CLOSE a file and perform proper validation for file handling
c     *****************************************************************
      subroutine CLFIL
     1     (LUFIL, LUSCRN, LULOG, IERR, LOGOP)

      IMPLICIT NONE

      INTEGER*4
     1      LUFIL, LUSCRN, LULOG,  IERR, ISTOP

      LOGICAL
     1      LOGOP

      CHARACTER*1200
     1     STRERR      
c     close file lufil
      
      CLOSE (UNIT=LUFIL, IOSTAT=IERR)

c     file successfully closed?      
      IF(IERR .NE. 0)
     1 THEN
c        close failed
         WRITE(STRERR, 1000) LUFIL, IERR
 1000    FORMAT('Fail to close unit ',I3,' ierr= ', I3)
 
         ISTOP = 107
         CALL erstop (istop, strerr, logop, luscrn, lulog)
      ENDIF

      RETURN

      END

c     *****************************************************************
c       Subroutine to exit on an error stop passed by the user.
c     *****************************************************************
      subroutine erstop (istop, strerr, logop, luscrn, lulog)

      implicit none
      integer*4 istop, lulog, luscrn
      CHARACTER*1200 strerr
      logical logop      

      write(luscrn, 1000) istop, strerr
      IF(LOGOP) write(LULOG, 1000) istop, strerr
 1000 format(' istop',i4,a1200)
c      exit is not a Standard Fortran function. However with g77 and f77 the function works ok.
c      If exit is not allowed by your compiler please use stop istop instead 
      call exit(istop)
c      stop istop
      return
      end
c     last line of file

c******************************************
      subroutine read_head_obs(nb_obs,LUFOBS,colGPSP1,colGPSP2,
     1           colGPSC1,colGLOP1,colGLOP2,colGLOC1,colGLOC2,
     1           colGALE1,colGALE5,factor,systime,ofset,tobs_first)
     
       integer IRCV_CLOCK,nb_obs(5),yy,mm,dd,min,hh,colGPSP1,
     2  colGPSP2,colGPSC1,colGLOP1,colGLOP2,colGALE1,colGALE5,
     2  colGLOC1,colGLOC2
      real*8 form,fact,factor(3,2),tobs_first,sec
      
      character*20  comment
      character*3 systime,obs(26),obsgal
      character*1 syst
      
      logical ofset
     
      SYSTIME='GPS'
      IRCV_CLOCK=0
      colGPSP1=0
      colGPSP2=0
      colGPSC1=0
      colGLOP1=0
      colGLOC1=0
      colGLOC2=0
      colGLOP2=0
      colGALE1=0
      colGALE5=0
      tobs_first=0
      
      factor=1.d0
           
  57  read(LUFOBS,'(60x,a20)') comment
  
      if(comment.eq.'RINEX VERSION / TYPE')then
        backspace LUFOBS
        read(LUFOBS,'(f9.2)')form
        if(dint(form).ne.3)then
	write(6,*)'observations NOT in RINEX 3, use another SW'
	stop
	endif
        goto 57
      elseif(comment.eq.'SYS / # / OBS TYPES ')then
        backspace LUFOBS
        read(LUFOBS,70)syst,nbo,(obs(i),i=1,min0(nbo,13))
	if(nbo.gt.13)read(LUFOBS,170)(obs(i),i=14,min0(nbo,26))
	if(nbo.gt.26)read(LUFOBS,170)(obs(i),i=27,nbo)
	if(SYST.eq.'G')then
	nb_obs(1)=nbo
        do  i=1,nbo
         if(obs(i).eq.'C1P'.or.obs(i).eq.'C1D'.or.obs(i).eq.'C1W')
     1              colGPSP1=i
         if(obs(i).eq.'C2P'.or.obs(i).eq.'C2D'.or.obs(i).eq.'C2W')
     1              colGPSP2=i
         if(obs(i).eq.'C1C'.or.obs(i).eq.'C1L'.or.obs(i).eq.'C1X')
     1              colGPSC1=i
        enddo
	elseif(SYST.eq.'R')then
	nb_obs(2)=nbo
        do  i=1,nbo
         if(obs(i).eq.'C1P')colGLOP1=i
         if(obs(i).eq.'C1C')colGLOC1=i
         if(obs(i).eq.'C2P')colGLOP2=i
         if(obs(i).eq.'C2C')colGLOC2=i
        enddo	
	elseif(SYST.eq.'E')then
	nb_obs(3)=nbo
        do  i=1,nbo
         if(obs(i).eq.'C1X'.or.obs(i).eq.'C1C'.or.obs(i).eq.'C1I')
     1                   colGALE1=i            
         if(obs(i).eq.'C5X'.or.obs(i).eq.'C5Q'.or.obs(i).eq.'C5I')
     1                   colGALE5=i
        enddo	
	elseif(SYST.eq.'C')then   !! COMPASS
	nb_obs(4)=nbo
	elseif(SYST.eq.'S')then    !! SBAS
	nb_obs(5)=nbo
	endif	
        goto 57
       elseif(comment.eq.'TIME OF FIRST OBS')then
        backspace 7  
        read(LUFOBS,'(5i6,18x,a3)')yy,mm,dd,hh,min,systime
        call datemjd(yy-2000,mm,dd,hh,min,sec,tobs_first)
        goto 57
      elseif(comment.eq.'RCV CLOCK OFFS APPL ')then
        backspace 7	
        read(LUFOBS,'(i6)')IRCV_CLOCK
       if(IRCV_CLOCK.eq.1)ofset=.true.
c if IRCV_CLOCK=1 : read the rcv clock applied each time, and correct for it.
        goto 57
      elseif(comment.eq.'SYS / SCALE FACTOR  ')then
        backspace 7	
        read(LUFOBS,175)syst,fact,n,(obs(i),i=1,n)
	if(syst.eq.'G')then
	 if(n.eq.0)then
	 factor(1,1)=dble(fact)
	 factor(1,2)=dble(fact)
	 else 
	 do i=1,n
	 obsgal=obs(i)
	 if(obsgal(1:2).eq.'C1')factor(1,1)=dble(fact)
	 if(obsgal(1:2).eq.'C2')factor(1,2)=dble(fact)
	 enddo
	 endif
	elseif(syst.eq.'R')then
	 if(n.eq.0)then
	 factor(2,1)=dble(fact)
	 factor(2,2)=dble(fact)
	 else 
	 do i=1,n
	 obsgal=obs(i)
	 if(obsgal(1:2).eq.'C1')factor(2,1)=dble(fact)
	 if(obsgal(1:2).eq.'C2')factor(2,2)=dble(fact)
	 enddo
	 endif
	elseif(syst.eq.'E')then
	 if(n.eq.0)then
	 factor(3,1)=dble(fact)
	 factor(3,2)=dble(fact)
	 else 
	 do i=1,n
	 obsgal=obs(i)
	 if(obsgal(1:2).eq.'C1')factor(3,1)=dble(fact)
	 if(obsgal(1:2).eq.'C5')factor(3,2)=dble(fact)
	 enddo
	 endif
	endif
		
        goto 57
      elseif(comment.eq.'END OF HEADER       ')then
      return
      else
          goto 57
      endif
      
  170 format(6x,13(1x,a3))     
   70 format(a1,2x,i3,13(1x,a3))
  175 format(a1,1x,i4,2x,i2,12(1x,a3))
      
  
      
      end
      
c************************************************************************      
      subroutine int2bin(intr,binr)
	integer intr,power,toBinary,bineq,tempbin(9),binr(3,2),temp,i
	character*9 charr
        
	if(intr.gt.511)write(6,*)'Error: only 9 bits binary'
	power=0
	temp=intr
	bineq=0
     
        DO WHILE(temp > 0)
          bineq = bineq + (MOD(temp, 2)*(MOD(temp,2)*(10**power)))
          power = power + 1
          temp = temp / 2
        END DO
	toBinary = bineq

        write(charr,'(i9)')toBinary
	read(charr,'(9(i1))')(tempbin(i),i=1,9)

c     E1B Data Validity Status and signal Health Status	
	binr(1,1)=tempbin(9)
	binr(1,2)=tempbin(8)+tempbin(7)*2

c     E5a Data Validity Status and signal Health Status
	binr(2,1)=tempbin(6)
	binr(2,2)=tempbin(5)+tempbin(4)*2

c     E5b Data Validity Status and signal Health Status
        binr(3,1)=tempbin(3)
	binr(3,2)=tempbin(2)+tempbin(1)*2	
      end
      
      

c******************************************************************

c     ******************************************************************
c     Subroutine Runge Kutta 4th Order for Satelitte Orbit integration 
c     ******************************************************************	
	

C       ttr = tps auquel on doit trouver la position
C  		delta = tps entre 2 epehmerides
C 		nbeph = nb epoch avec ephemerides (removed)
C		toc = date de l'ephemeride de l'IOE recherch?(ce n'est plus un vecteur)
C		sp3 = vecteur des eph.
c       xs : position (sortie)
c       vit: vitesse (sortie)
c       exists = test logique qui v?ifie que les positions ont bien ??calcul?s
c       LULOG = code fichier log
c       jkl = indice de l'ephemeride ?utiliser, d?ermin?avec le code suivant:
c             ( datenavglo (jkl) = date en MJD de l'?h??ide, (= toc qui est exprim?en GPSweek))
ccc	             do jkl=1,1000
ccc                if(jkl==1) then
ccc	                 if(datenavglo(jkl).gt.schtime(isch)+6.5/24./60.)goto 6600
ccc	               endif
ccc                if(datenavglo(jkl).gt.schtime(isch)+6.5/24./60.)cycle
ccc	               if((schtime(isch)+6.5/24./60.)-datenavglo(jkl)
ccc           1   .le.(30.d0/(24*60.d0)))then
ccc	                    goto 6600
ccc                endif
ccc	             enddo
ccc        6600     continue
c******************************************************************
        subroutine integ(ttr,toc,sp3,v,g,tau,Xs,
     1   prn,vit,exists,LULOG,jkl)
     
        implicit none
          real*8 omegae,graveart,c20bis
	  parameter (omegae=7.2921151467d-5,graveart=398600.44d0)
	  parameter (c20bis=-1.5d0*1.08263d-3*6378.136d0**2)
	 
	  real*8 ttx,tty,ttz,tsc,trx,try,trz,theta
	  integer prn,LULOG,i
	  
	  real*8 epoch,ttr,toc,tau,r2,r3,dt,xs(3),vit(3),temp(3),
     1   gi(3) 

	  real*8 sp3(1000,30,3),v(1000,30,3),g(1000,30,3),alpha,gs(3)
	  real*8 reste,pas 
	  integer nbcycle,jkl,j
	  logical exists
	  	        
	   xs=0.d0
	   vit=0.d0
	   gi=0.d0
	
	exists=.false.
	alpha=tau*omegae
	 
          dt=ttr-toc
	  if (dt.gt.302400.d0)  dt = dt - 604800.d0
          if (dt.lt.-302400.d0) dt = dt + 604800.d0
     	  pas=10.d0

c    Special case: dt<0, modifications to correct problem with negative number modulo
          if(dt.lt.0)then
              reste=-dmod(-dt,pas)
              pas=-pas
              nbcycle=ABS((dt-reste)/pas)
           else
              reste=dmod(dt,pas)
              nbcycle=(dt-reste)/pas
           endif

c    Extract initial parameters (for toc(jkl))
          if(sp3(jkl,prn,1)==0.d0.and.sp3(jkl,prn,2)==0.d0.and.
     1             sp3(jkl,prn,3)==0.d0) goto 703	  
          xs(1)=sp3(jkl,prn,1)
	  xs(2)=sp3(jkl,prn,2)
	  xs(3)=sp3(jkl,prn,3)
	  vit(1)=v(jkl,prn,1)
	  vit(2)=v(jkl,prn,2)
	  vit(3)=v(jkl,prn,3)
	  gi(1)=g(jkl,prn,1)
	  gi(2)=g(jkl,prn,2)
	  gi(3)=g(jkl,prn,3)

 
c    Determine initial acceleration parameters
          call f(xs,vit,gi,gs)	 
	 
c    Start integrating data with the step "pas"	 ==> xs,vit
            do j=1,nbcycle
	     call kutta(ttr,xs,vit,gs,gi,pas)  
	
c    Calcul corresponding acceleration gs from xs, vit and gi	
	     call f(xs,vit,gi,gs)
            enddo
c    Calcul for the last step lower than "pas" value
		call kutta(ttr,xs,vit,gs,gi,reste)
C     Rotation from PZ-90 to WGS-84 (not needed)

c    Correct positions from the earth rotation		 
	  temp=xs
	  xs(1)=(temp(1)*dcos(alpha)+temp(2)*dsin(alpha))*1.D3
	  xs(2)=(-temp(1)*dsin(alpha)+temp(2)*dcos(alpha))*1.D3
	  xs(3)= temp(3)*1.D3    
	  temp=vit
	  vit=vit*1.D3

  10   if(xs(1)==0..and.xs(2)==0..and.xs(3)==0.)then
	exists=.false.
	write(LULOG,702)prn,ttr
	else
	exists=.true.
	endif
 702   format('In BRDC, missing data for glonass satellite ',i3,' on ',
     1   f10.3,' (in GPS week) -> line skipped')

 703   continue
       end subroutine integ
	 
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	 
      subroutine kutta(ttr,x,v,g,gi,h)   
      implicit none
      real*8 ttr,h,omegae,graveart,c20bis,r2,r3   
      real*8 x(3),v(3),g(3),K(6,4),valu,xc(3),vc(3)
      real*8 gi(3),temp(3)
		 
      parameter (omegae=7.2921151467d-5,graveart=398600.44d0)
      parameter (c20bis=-1.5d0*1.08263d-3*6378.136d0**2)
		 
C     K [x x x x]
C       [x x x x]
C       [x x x x]
C       [x x x x]
C       [x x x x]
C       [x x x x]
		
C        Definition des constantes d'integration
        K=0

        K(1:3,1)=h*g
        K(4:6,1)=h*v

        vc(1)=v(1)+K(1,1)/2.D0
        vc(2)=v(2)+K(2,1)/2.D0
        vc(3)=v(3)+K(3,1)/2.D0
        xc(1)=x(1)+K(4,1)/2.D0
        xc(2)=x(2)+K(5,1)/2.D0
        xc(3)=x(3)+K(6,1)/2.D0
	 
        call f(xc,vc,gi,g)
        K(1:3,2)=h*g
        K(4:6,2)=h*vc
		
        vc(1)=v(1)+K(1,2)/2.D0
        vc(2)=v(2)+K(2,2)/2.D0
        vc(3)=v(3)+K(3,2)/2.D0
        xc(1)=x(1)+K(4,2)/2.D0
        xc(2)=x(2)+K(5,2)/2.D0
        xc(3)=x(3)+K(6,2)/2.D0

        call f(xc,vc,gi,g)
        K(1:3,3) =h*g
        K(4:6,3)=h*vc
		
          vc(1)=v(1)+K(1,3)
          vc(2)=v(2)+K(2,3)
          vc(3)=v(3)+K(3,3)
          xc(1)=x(1)+K(4,3)
          xc(2)=x(2)+K(5,3)
          xc(3)=x(3)+K(6,3)

         call f(xc,vc,gi,g)
	     K(1:3,4) =h*g        
		 K(4:6,4)=h*vc
		 
C 		Positions and velocities determined by RK4:

         vc(1)=v(1)+1./6.D0*(K(1,1)+2*K(1,2)+2*K(1,3)+K(1,4))
         vc(2)=v(2)+1./6.D0*(K(2,1)+2*K(2,2)+2*K(2,3)+K(2,4))
         vc(3)=v(3)+1./6.D0*(K(3,1)+2*K(3,2)+2*K(3,3)+K(3,4))
         xc(1)=x(1)+1./6.D0*(K(4,1)+2*K(4,2)+2*K(4,3)+K(4,4))
         xc(2)=x(2)+1./6.D0*(K(5,1)+2*K(5,2)+2*K(5,3)+K(5,4))
         xc(3)=x(3)+1./6.D0*(K(6,1)+2*K(6,2)+2*K(6,3)+K(6,4))				  
         x=xc
         v=vc 
         end subroutine kutta
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      subroutine f(x,v,gi,g)
      real *8 x(3),v(3),c20bis,graveart,omegae,r2,r3,cinqz2r2,
     1   gi(3),g(3)
      parameter (omegae=7.2921151467d-5,graveart=398600.44d0)
      parameter (c20bis=-1.5d0*1.08263d-3*6378.136d0**2) 
		 
      r2 = x(1)**2+x(2)**2+x(3)**2
      r3 = r2*dsqrt(r2)
      cinqz2r2 = 5.d0*x(3)**2/r2

c       determine acceleration from x,v and g		
      g(1) =graveart*x(1)*
     1    (-1.d0+c20bis*(1.d0-cinqz2r2)/r2)/r3
     1    + omegae*(omegae*x(1)+2.d0*v(2))+gi(1)
      g(2) = graveart*x(2)*
     1    (-1.d0+c20bis*(1.d0-cinqz2r2)/r2)/r3
     1   + omegae*(omegae*x(2)-2.d0*v(1))+gi(2)
      g(3) = graveart*x(3)*
     1   (-1.d0+c20bis*(3.d0-cinqz2r2)/r2)/r3+gi(3)
       end subroutine f
		 
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@		 
        subroutine leapscorr(yg,mg,dg,hg,ming,secg,gleap)
        implicit none
        integer yg,mg,dg,hg,ming,gleap
        real*8 secg
        if(gleap.eq.0)goto 10
        secg=secg+gleap
        if(secg.ge.60)then
           ming=ming+1
           secg=dmod(60.d0,secg)
        endif
        if(ming.ge.60)then
           hg=hg+1
           ming=mod(60,ming)
        endif
        if(hg.ge.24)then
           dg=dg+1
           hg=mod(24,hg)
        endif

 10   continue	
	return
	     end subroutine leapscorr

