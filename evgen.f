c-------------------------------------------------------------------------
c
c   evgen.f,   FORTRAN version 1.0
c
c
c Adrian Dumitru, July 1998
c perform Cooper-Frye FreezeOut with quantum statistics
c and write particle mom./coord. in OSCAR standard
c-------------------------------------------------------------------------


      IMPLICIT REAL*8 (A-H,O-Z)                                         
      INTEGER THET(104),GCLUST(104),NQI(104),NSI(104)
      REAL*8 MCLUST(104),NQU(104),NAQ(104),NST(104),NAS(104)
      REAL*8 MSTR,dpt
      REAL*8 xmax(104)
      CHARACTER*8 reffram
      real ran1

      parameter (i2maxbin=38)	   ! # of pT-bins
      parameter (imaxbin=20)       ! # of rT & tau - bins
      parameter (Nid=10)	   ! max. id to be considered, <=104
      real*8 dNdpt(Nid,i2maxbin,imaxbin,imaxbin,imaxbin)
      real*8 xnorm(Nid)
c
C
      PI=ACOS(-1D0)
      PI2=PI*PI
      HC3=197.32858D0**3
C
      OPEN(13,FILE='evgen_par.dat',STATUS='OLD')
      read(13,*) id_enh1           ! id-nr. to be enhanced
      read(13,*) id_enh2
      read(13,*) id_enh3
      read(13,*) id_enh4
      read(13,*) id_enh5
      read(13,*) rapgap            ! total rapidity gap
      read(13,*) bimp              ! impact parameter (irrelevant; appears only in output)
      read(13,*) RAu               ! initial Radius
      read(13,*) Nevents	   ! # of events
      read(13,*) dpt		   ! pt bin width, in MeV
      read(13,*) drtbin		   ! rt bin width, in units of RAu
      read(13,*) dtbin		   ! t bin width, in units of RAu
      read(13,*) idum		   ! random number seed (<0)
      CLOSE(13)


      dphi=2d0*PI/imaxbin
      xx=ran1(idum)

      READ(5,*) NCLUST,B,MSTR                                           
      READ(5,*) (MCLUST(I),I=1,NCLUST)                                  
      READ(5,*) (GCLUST(I),I=1,NCLUST)                                  
      READ(5,*) (NQU(I),I=1,NCLUST)                                     
      READ(5,*) (NAQ(I),I=1,NCLUST)                                     
      READ(5,*) (NST(I),I=1,NCLUST)                                     
      READ(5,*) (NAS(I),I=1,NCLUST)                                     
      READ(5,*) (THET(I),I=1,NCLUST)                                    
      if (NCLUST.gt.Nid) NCLUST=Nid

C                                                                       
      do 3 i=1,Nid
         do 3 j=1,i2maxbin
            do 3 j2=1,imaxbin
               do 3 j3=1,imaxbin
                  do 3 j4=1,imaxbin
                  dNdpt(i,j,j2,j3,j4)=0d0
 3    continue
C                                                                       
      DO 500 I=1,NCLUST                                                 
        NQI(I)=NQU(I)-NAQ(I)                                            
        NSI(I)=NST(I)-NAS(I)                                            
	if (ABS(NQI(i)-NQU(i)+NAQ(i)).gt.1d-5) then
	  write(*,*) ' NQI, NSI wrong, STOP'
	  STOP
	end if
	if (ABS(NSI(i)-NST(i)+NAS(i)).gt.1d-5) then
	  write(*,*) ' NST, NAS wrong, STOP'
	  STOP
	end if
  500 CONTINUE                                                          
C                                                                       
C
C                                                                       
      line=0
      r0=1d0
      time0=0d0
      open(40,FILE='hyper.dat',STATUS='OLD')
 2    read(40,*,END=6,ERR=99) r,time,eps,vel,tt,xmuq,xmus
c      write(*,*) r,time,eps,vel,tt,xmuq,xmus
        if(r.lt.0d0) goto 2
	gamma  = 1d0/sqrt(1d0-vel**2)
	dr = r-r0
	dt = time-time0
	r0 = r
	time0 = time
	if (line.eq.0) then
	  line=1
	  goto 2
	end if
        if (tt.lt.10d0) goto 2	! Temp>10 MeV cut for security

	do 10 i=1,Nid
c           write(*,*) i
	  fac =(xmuq*NQI(I)+xmus*NSI(I))/tt
	  do 20 j=1,i2maxbin
	  do 20 j4=1,imaxbin
            phi=dble(j4-1)*dphi
	    pt=j*dpt
            xmt=sqrt(pt**2+MCLUST(I)**2)
            s2=cos(phi)
	    a=gamma*vel*pt*s2/tt
	    a2=gamma*xmt/tt
c integral over y-eta from 0 to rapgap
            sum1=0d0
            sum2=0d0
	    h=rapgap/64d0          ! 64 is number of integration points
	    h1=0.77459667d0*h
	    do k=1,64,2            ! 64 is number of integration points
               cheta1=cosh(dble(k)*h-h1)
               cheta2=cosh(dble(k)*h)
               cheta3=cosh(dble(k)*h+h1)
               e1=a2*cheta1-a-fac
               e2=a2*cheta2-a-fac
               e3=a2*cheta3-a-fac
               if (e1.lt.1d2) then
                  ffac1=1d0/(exp(e1)+THET(i))
               else
                  ffac1=0d0
               end if
               if (e2.lt.1d2) then
                  ffac2=1d0/(exp(e2)+THET(i))
               else
                  ffac2=0d0
               end if
               if (e3.lt.1d2) then
                  ffac3=1d0/(exp(e3)+THET(i))
               else
                  ffac3=0d0
               end if
c               write(*,*) ffac1,ffac2,ffac3
               sum1=sum1+5d0*ffac1+8d0*ffac2+5d0*ffac3
               sum2=sum2+
     &            5d0*ffac1*cheta1+8d0*ffac2*cheta2+5d0*ffac3*cheta3
            end do
c h/9 from Gauss Integration, factor 2  because y-eta integration was done only for y-eta>0
            sum1=sum1*h/9d0*2d0
            sum2=sum2*h/9d0*2d0
	    a3 = pt*dt*s2*sum1-xmt*dr*sum2

c negative contribution from timelike part of hypersurface could be cut out
c	    if (a3.lt.0d0) goto 20
            j2=r/drtbin+1
            j3=time/dtbin+1
            dNdpt(i,j,j2,j3,j4)=dNdpt(i,j,j2,j3,j4) + a3*pt*r*time
 20	  continue
 10	continue
	goto 2

  6	continue
        do 8 i=1,Nid
          MCLUST(i)=MCLUST(i)/1d3     ! Conversion MeV -> GeV
          xnorm(i)=0d0	! tot. # of particles of species i
          xmax(i)=0d0	! max. of probab. distrib.
          do 11 j=i2maxbin,1,-1
             sum0=0d0
             do 9 j2=1,imaxbin
                do 9 j3=1,imaxbin
                   do 9 j4=1,imaxbin
            sum0=sum0+dNdpt(i,j,j2,j3,j4)
            if(dNdpt(i,j,j2,j3,j4).gt.xmax(i)) 
     &           xmax(i)=dNdpt(i,j,j2,j3,j4)
 9        continue
          xnorm(i)=xnorm(i)+sum0
 11       continue
c 39.47842 = (2*PI)**2
          xnorm(i)=xnorm(i)*dpt*dphi*RAu**3/HC3*GCLUST(I)/39.47842d0
 8      continue
        dpt=dpt/1d3     ! Conversion MeV -> GeV

c ENHANCE some hadron species
        if(id_enh1.lt.104) xnorm(id_enh1)=xnorm(id_enh1)*1d1
        if(id_enh2.lt.104) xnorm(id_enh2)=xnorm(id_enh2)*1d1
        if(id_enh3.lt.104) xnorm(id_enh3)=xnorm(id_enh3)*1d1
        if(id_enh4.lt.104) xnorm(id_enh4)=xnorm(id_enh4)*1d1
        if(id_enh5.lt.104) xnorm(id_enh5)=xnorm(id_enh5)*1d1

        do 44 i=1,Nid
           do 44 j=1,i2maxbin
              do 44 j2=1,imaxbin
                 do 44 j3=1,imaxbin
                    do 44 j4=1,imaxbin
c normalize probab. distrib. to maximum of 1
                    dNdpt(i,j,j2,j3,j4)=dNdpt(i,j,j2,j3,j4)/xmax(i)
 44     continue

        itotp=0
        do i=1,Nid
           do npart=1,NINT(xnorm(i)*rapgap)
	      id=id_clust(i,idum)
              if(id.ne.0) itotp=itotp+1
           end do
        end do

C %%%%%%%%%%%%%%%% write OSCAR header %%%%%%%%%%%%%%
       OPEN(unit=19,FILE='events.dat')

c define output format
      write (19,901) 'OSC1997A    '
      write (19,901) 'final_id_p_x'
 901  format (a12)

c eqsp = equal speed
c tar  = lab frame
c pro  = proj. frame
      reffram='eqsp'
c     reffram='tar'
c     reffram='pro'

      iapp=197
      iatt=197
      izpp=79
      iztt=79
      ebeam=1E2
      write(19,902) 'BjHydro', '1.0', iapp, izpp, iatt, iztt,
     .     reffram, ebeam, 1
 902  format (2(a8,2x),'(',i3,',',i6,')+(',i3,',',i6,')',2x,a4,2x,
     &     e10.4,2x,i8)


C **************** loop over events **************
      do ievent=1,Nevents
        write(19,903) ievent, itotp, bimp, 0D0
 903  format (i10,2x,i10,2x,f8.3,2x,f8.3)

C ++++++++++++++++ loop over particle species ++++++++++++++++++
        icount=0
        do i=1,Nid
C ---------- loop over all particles of given species in given event --------
           do npart=1,NINT(xnorm(i)*rapgap)
c determine OSCAR particle id
	      id=id_clust(i,idum)
c this particle species might not be interesting ...
	      if (id.eq.0) goto 77
              icount=icount+1
 222          ipTbin=1+INT(real(i2maxbin)*ran1(idum))
              irTbin=1+INT(real(imaxbin)*ran1(idum))
              itbin =1+INT(real(imaxbin)*ran1(idum))
              iphibin=1+INT(real(imaxbin)*ran1(idum))
              prob=ran1(idum)
              if (prob.le.dNdpt(i,ipTbin,irTbin,itbin,iphibin)) then
                 rap=rapgap*(ran1(idum)-.5d0)

                 pT=(dble(ipTbin)-ran1(idum)+.5)*dpt
                 xmT=sqrt(pT**2+MCLUST(i)**2)
                 psi=2d0*PI*ran1(idum)
                 px=pT*sin(psi)
                 py=pT*cos(psi)
                 pz=xmT*sinh(rap)
                 energy=xmT*cosh(rap)

                 rT=(dble(irTbin)-ran1(idum))*drtbin*RAu
                 phi=dble(iphibin-ran1(idum)+.5)*dphi
                 rx=rT*sin(psi-phi)
                 ry=rT*cos(psi-phi)

c define space-time rapidity
                 if (MCLUST(i).ge.1d0) then
                    strap=rap
                 else
                    strap=rap+(ran1(idum)-.5)
                 end if
                 time=(dble(itbin)-ran1(idum))*dtbin*RAu*cosh(strap)

                 rz=time*tanh(strap)

                 write(19,904) icount, id,
     &                px, py, pz, energy, MCLUST(i),     
     &                rx, ry, rz, time

              else
c momentum/coord. combination rejected, try another ...
                 goto 222
              end if
C ---------- end of loop over particles of given species & event --------
 77	     continue
           end do
C ++++++++++++++++ end of loop over particle species ++++++++++++++++++
        end do
C **************** end of loop over events **************
      end do
C
C That's it !
      STOP
 777  FORMAT(E10.3,15G13.5)
 904  format (i10,2x,i10,2x,9(e12.6,2x))
 99   write(6,*) 'ERROR in hypersurface file'
      STOP
      END


c ********************************************
c determine OSCAR-particle-id from internal id
c ********************************************
      function id_clust(i,idum)
      implicit real*8 (a-h,o-z)
      real ran1

      if (i.eq.1) then
	j=3.*ran1(idum)
        if (j.eq.0) id_clust=111		! pi0
        if (j.eq.1) id_clust=211		! pi+
        if (j.eq.2) id_clust=-211		! pi-
      else if (i.eq.63) then
	id_clust=60221				! sigma
      else if (i.eq.2) then
	id_clust=221				! eta
      else if (i.eq.6) then
	id_clust=331				! eta'
      else if (i.eq.3) then
	j=3.*ran1(idum)
        if (j.eq.0) id_clust=113		! rho0
        if (j.eq.1) id_clust=213		! rho+
        if (j.eq.2) id_clust=-213		! rho-
      else if (i.eq.77) then
	j=3.*ran1(idum)
        if (j.eq.0) id_clust=40113		! rho_1465
        if (j.eq.1) id_clust=40213		! 
        if (j.eq.2) id_clust=-40213		! 
      else if (i.eq.81) then
	j=3.*ran1(idum)
        if (j.eq.0) id_clust=30113		! rho_1700
        if (j.eq.1) id_clust=30213		! 
        if (j.eq.2) id_clust=-30213		! 
      else if (i.eq.4) then
	id_clust=223				! omega
      else if (i.eq.79) then
	id_clust=50223				! omega_1420
      else if (i.eq.83) then
	id_clust=60223				! omega_1600
      else if (i.eq.61) then
	j=3.*ran1(idum)
        if (j.eq.0) id_clust=10111		! a_0 (980)
        if (j.eq.1) id_clust=10211		! 
        if (j.eq.2) id_clust=-10211		! 
      else if (i.eq.65) then
	j=3.*ran1(idum)
        if (j.eq.0) id_clust=20113		! a_1^0 (1260)
        if (j.eq.1) id_clust=20213		! a_1^+
        if (j.eq.2) id_clust=-20213		! a_1^-
      else if (i.eq.69) then
	j=3.*ran1(idum)
        if (j.eq.0) id_clust=115		! a_2 (1320)
        if (j.eq.1) id_clust=215		! 
        if (j.eq.2) id_clust=-215		! 
      else if (i.eq.73) then
	j=3.*ran1(idum)
        if (j.eq.0) id_clust=10113		! b_1(1235)
        if (j.eq.1) id_clust=10213		! 
        if (j.eq.2) id_clust=-10213		! 
      else if (i.eq.75) then
	id_clust=10223                  	! h_1(1170)
      else if (i.eq.64) then
	id_clust=10221                  	! f_0(980)
      else if (i.eq.71) then
	id_clust=225				! f_2(1270)
      else if (i.eq.67) then
	id_clust=20223				! f_1(1285)
      else if (i.eq.68) then
	id_clust=40223				! f_1(1510)
      else if (i.eq.72) then
	id_clust=335				! f_2'(1525)
      else if (i.eq.39) then
	id_clust=311+10*NINT(ran1(idum))	! K0 -> 311, K+ -> 321
      else if (i.eq.38) then
	id_clust=-311-10*NINT(ran1(idum))	! A-K0 -> -311, K- -> -321
      else if (i.eq.66) then
        id_clust=10313+10*NINT(ran1(idum))      ! K_1(1270)
      else if (i.eq.86) then
        id_clust=-10313-10*NINT(ran1(idum))     ! A-K_1(1270)
      else if (i.eq.74) then
        id_clust=20313+10*NINT(ran1(idum))      ! K_1(1400)
      else if (i.eq.88) then
        id_clust=-20313-10*NINT(ran1(idum))     ! A-K_1(1400)
      else if (i.eq.62) then
	id_clust=10311+10*NINT(ran1(idum))	! K_0(1430)
      else if (i.eq.85) then
	id_clust=-10311-10*NINT(ran1(idum))	! A-K_0(1430)
      else if (i.eq.70) then
	id_clust=315+10*NINT(ran1(idum))	! K_2(1430)
      else if (i.eq.87) then
	id_clust=-315-10*NINT(ran1(idum))	! A-K_2(1430)
      else if (i.eq.49) then
        id_clust=313+10*NINT(ran1(idum))        ! K*(895)
      else if (i.eq.48) then
        id_clust=-313-10*NINT(ran1(idum))       ! A-K*(895)
      else if (i.eq.78) then
        id_clust=30313+10*NINT(ran1(idum))      ! K*(1410)
      else if (i.eq.89) then
        id_clust=-30313-10*NINT(ran1(idum))     ! A-K*(1410)
      else if (i.eq.82) then
        id_clust=40313+10*NINT(ran1(idum))      ! K*(1680)
      else if (i.eq.90) then
        id_clust=-40313-10*NINT(ran1(idum))     ! A-K*(1680)
      else if (i.eq.56) then
	id_clust=333				! phi
      else if (i.eq.80) then
	id_clust=10333				! phi_1680
      else if (i.eq.5) then
	id_clust=2112+100*NINT(ran1(idum))	! n -> 2112, p -> 2212
      else if (i.eq.22) then
	id_clust=-2112-100*NINT(ran1(idum))	! A-n -> -2112, A-p -> -2212
      else if (i.eq.8) then
	id_clust=12112+100*NINT(ran1(idum))	! N_1440
      else if (i.eq.24) then
	id_clust=-12112-100*NINT(ran1(idum))	! A-N_1440
      else if (i.eq.9) then
	id_clust=1214+910*NINT(ran1(idum))	! N_1520
      else if (i.eq.25) then
	id_clust=-1214-910*NINT(ran1(idum))	! A-N_1520
      else if (i.eq.10) then
	id_clust=22112+100*NINT(ran1(idum))	! N_1535
      else if (i.eq.26) then
	id_clust=-22112-100*NINT(ran1(idum))	! A-N_1535
      else if (i.eq.12) then
	id_clust=32112+100*NINT(ran1(idum))	! N_1650
      else if (i.eq.28) then
	id_clust=-32112-100*NINT(ran1(idum))	! A-N_1650
      else if (i.eq.13) then
	id_clust=2116+100*NINT(ran1(idum))	! N_1675
      else if (i.eq.29) then
	id_clust=-2116-100*NINT(ran1(idum))	! A-N_1675
      else if (i.eq.14) then
	id_clust=12116+100*NINT(ran1(idum))	! N_1680
      else if (i.eq.30) then
	id_clust=-12116-100*NINT(ran1(idum))	! A-N_1680
      else if (i.eq.15) then
	id_clust=21214+910*NINT(ran1(idum))	! N_1700
      else if (i.eq.31) then
	id_clust=-21214-910*NINT(ran1(idum))	! A-N_1700
      else if (i.eq.16) then
	id_clust=42112+100*NINT(ran1(idum))	! N_1710
      else if (i.eq.32) then
	id_clust=-42112-100*NINT(ran1(idum))	! A-N_1710
      else if (i.eq.17) then
	id_clust=31214+910*NINT(ran1(idum))	! N_1720
      else if (i.eq.33) then
	id_clust=-31214-910*NINT(ran1(idum))	! A-N_1720
      else if (i.eq.7) then
	j=4.*ran1(idum)
        if (j.eq.0) id_clust=1114		! Delta
        if (j.eq.1) id_clust=2114		! 
        if (j.eq.2) id_clust=2214		! 
        if (j.eq.3) id_clust=2224		! 
      else if (i.eq.23) then
	j=4.*ran1(idum)
        if (j.eq.0) id_clust=-1114		! A-Delta
        if (j.eq.1) id_clust=-2114		! 
        if (j.eq.2) id_clust=-2214		! 
        if (j.eq.3) id_clust=-2224		! 
      else if (i.eq.91) then
	j=4.*ran1(idum)
        if (j.eq.0) id_clust=31114		! Delta_1600
        if (j.eq.1) id_clust=32114		! 
        if (j.eq.2) id_clust=32214		! 
        if (j.eq.3) id_clust=32224		! 
      else if (i.eq.95) then
	j=4.*ran1(idum)
        if (j.eq.0) id_clust=-31114		! A-Delta_1600
        if (j.eq.1) id_clust=-32114		! 
        if (j.eq.2) id_clust=-32214		! 
        if (j.eq.3) id_clust=-32224		! 
      else if (i.eq.11) then
	j=4.*ran1(idum)
        if (j.eq.0) id_clust=1112		! Delta_1620
        if (j.eq.1) id_clust=1212		! 
        if (j.eq.2) id_clust=2122		! 
        if (j.eq.3) id_clust=2222		! 
      else if (i.eq.27) then
	j=4.*ran1(idum)
        if (j.eq.0) id_clust=-1112		! A-Delta_1620
        if (j.eq.1) id_clust=-1212		! 
        if (j.eq.2) id_clust=-2122		! 
        if (j.eq.3) id_clust=-2222		! 
      else if (i.eq.92) then
	j=4.*ran1(idum)
        if (j.eq.0) id_clust=11114		! Delta_1700
        if (j.eq.1) id_clust=12114		! 
        if (j.eq.2) id_clust=12214		! 
        if (j.eq.3) id_clust=12224		! 
      else if (i.eq.96) then
	j=4.*ran1(idum)
        if (j.eq.0) id_clust=-11114		! A-Delta_1700
        if (j.eq.1) id_clust=-12114		! 
        if (j.eq.2) id_clust=-12214		! 
        if (j.eq.3) id_clust=-12224		! 
      else if (i.eq.93) then
	j=4.*ran1(idum)
        if (j.eq.0) id_clust=11112		! Delta_1900
        if (j.eq.1) id_clust=11212		! 
        if (j.eq.2) id_clust=12122		! 
        if (j.eq.3) id_clust=12222		! 
      else if (i.eq.97) then
	j=4.*ran1(idum)
        if (j.eq.0) id_clust=-11112		! A-Delta_1900
        if (j.eq.1) id_clust=-11212		! 
        if (j.eq.2) id_clust=-12122		! 
        if (j.eq.3) id_clust=-12222		! 
      else if (i.eq.18) then
	j=4.*ran1(idum)
        if (j.eq.0) id_clust=1116		! Delta_1905
        if (j.eq.1) id_clust=1216		! 
        if (j.eq.2) id_clust=2126		! 
        if (j.eq.3) id_clust=2226		! 
      else if (i.eq.34) then
	j=4.*ran1(idum)
        if (j.eq.0) id_clust=-1116		! A-Delta_1905
        if (j.eq.1) id_clust=-1216		! 
        if (j.eq.2) id_clust=-2126		! 
        if (j.eq.3) id_clust=-2226		! 
      else if (i.eq.19) then
	j=4.*ran1(idum)
        if (j.eq.0) id_clust=21112		! Delta_1910
        if (j.eq.1) id_clust=21212              ! 
        if (j.eq.2) id_clust=22122		! 
        if (j.eq.3) id_clust=22222		! 
      else if (i.eq.35) then
	j=4.*ran1(idum)
        if (j.eq.0) id_clust=-21112		! A-Delta_1910
        if (j.eq.1) id_clust=-21212		! 
        if (j.eq.2) id_clust=-22122		! 
        if (j.eq.3) id_clust=-22222		! 
      else if (i.eq.94) then
	j=4.*ran1(idum)
        if (j.eq.0) id_clust=21114		! Delta_1920
        if (j.eq.1) id_clust=22114		! 
        if (j.eq.2) id_clust=22214		! 
        if (j.eq.3) id_clust=22224		! 
      else if (i.eq.98) then
	j=4.*ran1(idum)
        if (j.eq.0) id_clust=-21114		! A-Delta_1920
        if (j.eq.1) id_clust=-22114		! 
        if (j.eq.2) id_clust=-22214		! 
        if (j.eq.3) id_clust=-22224		! 
      else if (i.eq.20) then
	j=4.*ran1(idum)
        if (j.eq.0) id_clust=11116		! Delta_1930
        if (j.eq.1) id_clust=11216		! 
        if (j.eq.2) id_clust=12126		! 
        if (j.eq.3) id_clust=12226		! 
      else if (i.eq.36) then
	j=4.*ran1(idum)
        if (j.eq.0) id_clust=-11116		! A-Delta_1930
        if (j.eq.1) id_clust=-11216		! 
        if (j.eq.2) id_clust=-12126		! 
        if (j.eq.3) id_clust=-12226		! 
      else if (i.eq.21) then
	j=4.*ran1(idum)
        if (j.eq.0) id_clust=1118		! Delta_1950
        if (j.eq.1) id_clust=2118		! 
        if (j.eq.2) id_clust=2218		! 
        if (j.eq.3) id_clust=2228		! 
      else if (i.eq.37) then
	j=4.*ran1(idum)
        if (j.eq.0) id_clust=-1118		! A-Delta_1950
        if (j.eq.1) id_clust=-2118		! 
        if (j.eq.2) id_clust=-2218		! 
        if (j.eq.3) id_clust=-2228		! 
      else if (i.eq.40) then
	id_clust=3122				! L
      else if (i.eq.41) then
	id_clust=-3122				! A-L
      else if (i.eq.50) then
	id_clust=13122				! L_1405
      else if (i.eq.51) then
	id_clust=-13122				! A-L_1405
      else if (i.eq.42) then
	j=3.*ran1(idum)
        if (j.eq.0) id_clust=3222		! Sigma+,0,-
        if (j.eq.1) id_clust=3212		! 
        if (j.eq.2) id_clust=3112		! 
      else if (i.eq.43) then
	j=3.*ran1(idum)
        if (j.eq.0) id_clust=-3222		! A-Sigma+,0,-
        if (j.eq.1) id_clust=-3212		! 
        if (j.eq.2) id_clust=-3112		! 
      else if (i.eq.52) then
	j=3.*ran1(idum)
        if (j.eq.0) id_clust=3114		! Sigma_1385
        if (j.eq.1) id_clust=3214		! 
        if (j.eq.2) id_clust=3224		! 
      else if (i.eq.53) then
	j=3.*ran1(idum)
        if (j.eq.0) id_clust=-3114		! A-Sigma_1385
        if (j.eq.1) id_clust=-3214		! 
        if (j.eq.2) id_clust=-3224		! 
      else if (i.eq.44) then
        id_clust=3312+10*NINT(ran1(idum))       ! Casc.0,-
      else if (i.eq.45) then
        id_clust=-3312-10*NINT(ran1(idum))      ! A-Casc.0,-
      else if (i.eq.54) then
        id_clust=3314+10*NINT(ran1(idum))       ! Casc._1530
      else if (i.eq.55) then
        id_clust=-3314-10*NINT(ran1(idum))      ! A-Casc._1530
      else if (i.eq.46) then
	id_clust=3334				! Omega-
      else if (i.eq.47) then
	id_clust=-3334				! A-Omega-
      else
	id_clust=0				! uninteresting particle
      end if
      return
      END


