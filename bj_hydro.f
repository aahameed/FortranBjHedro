c**********************************************************************
c                                                                     *
c    This program calculates the 1+1d expansion of a blob in vacuum.  *
c                                                                     *
c    The correction due to cylindrical geometry and longitudinal      *
c    scaling is implemented via Sod's method.                         *
c                                                                     *
c    The Eos and the velocity of sound are to be specified in         *
c    appropriate subprograms.                                         *
c    However, there is no propagation of (mass) density, so p and     *
c    cs should not depend explicitly on (mass) density.               *
c                                                                     *
c    The signal velocities are estimated according to possibility 6,  *
c    i.e., a viscosity proportional to the velocity gradient is       *
c    added to the velocity of sound.                                  *
c                                                                     *
c    Original Code by D.H. Rischke                                    *
c    some modifications and the tab.-EoS by A. Dumitru                *
c**********************************************************************
      subroutine sinit(e,gamma,ks)
c
c  This subroutine-subprogram reads initialization data.
c  It is used in the main program.
c
c  Type declarations for variables in common-blocks
c
      real*8 gm,B
      real*8 dx,dt,t
      real*8 eta
      real*8 alp
      integer eos
c
c  Type declarations for variables in subroutine sinit
c
      real*8 e(1000),gamma(1000)
      real*8 r,epc,ra
      real*8 lam
      real*8 ks
      integer i
c
c  common-blocks
c
      common /adind/ gm,B
      common /gitter/ dx,dt,t
      common /visko/ eta
      common /flags/ eos 
      common /geom/ alp
c
c  Start reading data.
c  Initialize cell size.
c
      read(5,*) dx
c
c  Initialize ratio time step / cell size
c
      read(5,*) lam
      dt = lam * dx
c
c  Initialize eta in calculation of mean sound velocity
c
      read(5,*) eta
c
c  Initialize ratio of degrees of freedom in QGP to hadron matter,
c  adiabatic index of the EoS, epc, the initial energy density of 
c  the slab in units of the critical pressure (resp. in units of the
c  e0=147.7MeV/fm^3 for tab.-EoS), and the radius ra of the system.
c  Then initialize Bag constant B (in units of the critical pressure), 
c  and the fields e and gamma.
c
      r = 37d0/3d0
      gm = 1d0/3d0
      read(5,*) epc,ra
      B = r - 1d0
      do 113 i = 1,1000
        if (dabs(dfloat(i)-500.5d0)*dx.le.ra) then
          e(i) = epc
        else
          e(i) = 0d0
        end if
        gamma(i) = 1d0
 113  continue
c
c  Initialize slope approximation
c
      ks = 1d0
c
c  Initialize flags:
c     eos=1         equation of state is ultrarelativistic ideal gas
c     eos=0         ur-idgas + Bag Model QGP, 1st order phase transition
c     eos=2         all hadrons + Bag Model QGP, 1st order phase transition
      read(5,*) eos
c
c  Initialize geometry: alp = 0d0     longitudinal expansion
c                       alp = 1d0     cylindrical expansion
c                       alp = 2d0     spherical expansion
c                       alp = 3d0     cylindrical plus long. scaling
c
      read(5,*) alp
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine init(e,p,v,gamma,elab,m)
c
c  This subroutine-subprogram calculates the laboratory frame
c  quantities from the local rest frame quantities read in subroutine
c  sinit.
c  It is used in the main program.
c
c  Type declarations for variables in common-blocks
      real*8 total0(2),alp
      real*8 dx,dt,t
      integer kp(1000),km(1000)
c
c  Type declarations for variables in subroutine init
c
      real*8 e(1000),p(1000),v(1000),gamma(1000)
      real*8 elab(1000),m(1000),press
      integer i
c
c  common-blocks
c
      common /indpoi/ kp,km
      common /conser/ total0
      common /gitter/ dx,dt,t
      common /geom/ alp 
c
c  Start calculation. Null vector total0.
c
      do 1000 i = 1,2
        total0(i) = 0d0
 1000  continue
c
c  Start calculation. Each cell is treated separately.
c
      do 1001 i = 1,1000
c
c  Calculate index pointer (used in subroutine prop).
c
        kp(i) = i+1
        if (i.eq.1000) kp(i) = 1000
        km(i) = i-1
        if (i.eq.1) km(i) = 1
        if (e(i).le.0d0) then
c
c  Zero or negative temperatures correspond to vacuum.
c  The vacuum is supposed to have velocity -1 in (-x)-direction
c  and velocity +1 in (+x)-direction.
c
          e(i) = 0d0
          p(i) = 0d0
          if (i.le.500) v(i)=-1d0
          if (i.gt.500) v(i)=1d0
          elab(i) = 0d0
          m(i) = 0d0
          gamma(i) = 1d8
        else
          if (dabs(gamma(i)-1d0).lt.1d-16) then
c
c  If gamma is nearly 1.0, set v=0.
c
            v(i) = 0d0
          else if (gamma(i).lt.0d0) then
c
c  If gamma is negative, this corresponds to velocities in
c  (-x)-direction.
c
            v(i) = -dsqrt(1d0-1d0/(gamma(i)*gamma(i)))
          else
            v(i) = dsqrt(1d0-1d0/(gamma(i)*gamma(i)))
          end if
c
c  Afterwards, gamma gets its mathematical (i.e. positive) value.
c
          gamma(i) = dabs(gamma(i))
c
c  The pressure is calculated according to the equation of state.
c
          p(i) = press( e(i) )
          m(i) = gamma(i)*gamma(i) * (e(i)+p(i)) * v(i)
          elab(i) = gamma(i)*gamma(i) * (e(i)+p(i)) - p(i)
        end if
        if (alp.ne.0d0) then
          total0(1) = total0(1) + 
     1         elab(i)*(dx*dabs(dfloat(i)-500.5d0))**alp*dx
          total0(2) = total0(2) +
     1         m(i)*(dx*dabs(dfloat(i)-500.5d0))**alp*dx
        else
          total0(1) = total0(1) + 
     1         elab(i)*dx
          total0(2) = total0(2) +
     1         m(i)*dx
        end if
 1001 continue
      return
      end
c
c---------------------------------------------------------------------
c
      subroutine untang(elab,m,e,p,v,gamma)
c
c  This subroutine-subprogram calculates local rest frame quantities
c  from laboratory frame quantities.
c  It is used in subroutine prop and the main program.
c
c  Type declarations for variables in subroutine untang.
c
      real*8 e(1000),p(1000),v(1000),gamma(1000)
      real*8 elab(1000),m(1000)
      real*8 velo,press
      integer i
c
c  Start calculation. Each cell is separately treated.
c  Vacuum is assumed if the absolute amount of the 
c  lab-fields elab, m is smaller than 1e-16.
c
      do 100 i = 1,1000
        if ((dabs(m(i)).lt.1d-16).and.(dabs(elab(i)).lt.1d-16)) then
          e(i) = 0d0
          if (i.le.500) v(i) = -1d0
          if (i.gt.500) v(i) = 1d0
          gamma(i) = 1d8
        else if (dabs(m(i)).lt.1d-16) then
          e(i) = elab(i)
          v(i) = 0d0
          gamma(i) = 1d0
        else 
          v(i) = velo(elab(i),m(i))
          if (dabs(v(i)).lt.1d0-1d-16) then
            e(i) = elab(i) - v(i)*m(i)
            gamma(i) = 1d0/dsqrt(1d0-v(i)*v(i))
          else
            e(i) = 0d0
            if (i.le.500) v(i)=-1d0
            if (i.gt.500) v(i)=1d0
            gamma(i) = 1d8 
          end if
        end if
        p(i) = press( e(i) )
 100  continue
      return
      end
c
c----------------------------------------------------------------------
c
      function velo(el,ml)
c
c  This function-subprogram calculates the velocity as a fix
c  point of v =  M / ( Elab + p ) .
c  It is used in the subroutine untang.
c  el,ml are the values of Elab,M.
c
c  Type declarations for variables in function velo
c
      real*8 velo,el,ml
      real*8 v,f,root,press
      integer i
c
c  Fix-point iteration of the function v = f(v), with
c  f(x) = M / ( Elab + p( Elab-v*M ) )
c
c      if (el.lt.dabs(ml)) write(6,*) el,ml 
      i = 0
      f = 0d0
  10  v = f
      root = dsqrt( 1d0 - v*v )
      f = ml / ( el + press(el-v*ml) )
      i = i + 1
      if (i.gt.100) then
        write(6,*) ' No root found for Elab=',el,', M=',ml
        write(6,*) ' Program terminated.'
        stop
      end if
      if (dabs((v-f)/f).gt.1d-16) goto 10
      velo = f
      return
      end
c
c----------------------------------------------------------------------
c
      function press(a)
c
c  This function-subprogram determines the pressure of the
c  underlying EoS.
c  It is used in subroutine velo, untang.
c  a is the energy density in the local rest frame.
c
c  Type declarations for variables in common-blocks
c
      real*8 gm,B
      integer eos
c
c  Type declarations for variables used in function press
c
      real*8 a,press,press_table
c
c  common-blocks
c
      common /adind/ gm,B
      common /flags/ eos
c
c  Equation of state is ultrarelativistic ideal gas if flag eos=1,
c  otherwise, equation of state has phase transition
c
      if (eos.eq.1) then
        press = gm*a
      else if (eos.eq.0) then
        if (a.ge.(3d0+4d0*B)) then
          press = gm*(a-4d0*B)
        else if (a.ge.3d0) then
          press = 1d0
        else
          press = gm*a   
        end if
      else if (eos.eq.2) then
         press=press_table(a)
      end if
      return
      end
c
c----------------------------------------------------------------------
c This function determines the pressure for the tab.-EoS (same units
c  as energy density!)
c
      function press_table(e)
      implicit real*8 (a-h,o-z)
c
c  common-blocks
c
      common /pressuretab/ de,presstab(0:5000),
     &			temptab(0:5000)


      i=int(e/de)
      offe=e-DBLE(i)*de
      if (i.ge.0) then
         if (i.gt.4999) then
C *** 10.5=4B/e0 for B=380MeV/fm**3 ***
            press_table=(e-10.5d0)/3d0
            return
         end if
         p1 = presstab(i)
         p2 = presstab(i+1)
         press_table = (p2-p1)/de*offe + p1
      else
         press_table=0d0
      end if

      if (press.lt.0d0) press=0d0
      return
      end
c
c----------------------------------------------------------------------
c This function determines the temperature
c   - in units of T_C for eos=0 (ur-id.gas+QGP EoS)
c   - in units of p or e ^1/4 for eos=1 (ur-id.gas EoS)
c   - in units of MeV for eos=2 (tab.-EoS)
c
      function temp(e)
      implicit real*8 (a-h,o-z)
c
c  Type declaration for variables used in common-blocks
c
      real*8 gm,B
      real*8 dx,dt,t
      real*8 total0(2),alp
      integer eos
c
c  common-blocks
c
      common /adind/ gm,B
      common /gitter/ dx,dt,t
      common /conser/ total0
      common /flags/ eos
      common /geom/ alp

        if (eos.eq.1) then
          ttc = (e/3d0)**(gm/(1d0+gm))
        else if (eos.eq.0) then
          if (e.ge.(3d0+4d0*B)) then
            ttc = ((e-B)/3d0/(1d0+B))**(gm/(1d0+gm))
         else if (e.ge.3d0) then
            ttc = 1d0
         else
            ttc = (e/3d0)**(gm/(1d0+gm)) 
         end if
        else if (eos.eq.2) then
           ttc=temp_table(e)
        end if
        temp=ttc
        return
        end
c----------------------------------------------------------------------
c This function determines the temperature for the tab.-EoS in MeV!
c
      function temp_table(e)
      implicit real*8 (a-h,o-z)
c
c  common-blocks
c
      common /pressuretab/ de,presstab(0:5000),
     &			temptab(0:5000)

      i=int(e/de)
      offe=e-DBLE(i)*de
      if (i.ge.0) then
         if (i.gt.4999) then
            temp_table=1d3            ! out of bounds, return 1GeV
            return
         end if
         p1 = temptab(i)
         p2 = temptab(i+1)
         temp_table = (p2-p1)/de*offe + p1
      else
         temp_table=0d0
      end if
      return
      end
c----------------------------------------------------------------------
c
      function sound(a)
c
c  This function-subprogram determines the velocity of sound of the
c  underlying EoS.
c  It is used in subroutine prop, in the estimate for the signal
c  velocities.
c  a is the energy density in the local rest frame.
c
c  Type declarations for variables in common-blocks
c
      real*8 gm,B
      integer eos
c
c  Type declarations for variables used in function sound
c
      real*8 a,sound
c
c  common-blocks
c
      common /adind/ gm,B 
      common /flags/ eos
c
c  Velocity of sound is that of ultrarelativistic ideal gas 
c  if flag eos=1,
c  otherwise, velocity of sound is zero in mixed phase.
c     
      if (eos.eq.1) then
        sound = dsqrt(gm)
      else
        sound=sqrt(1d0/3d0)
      end if
c        if (a.ge.(3d0+4d0*B)) then
c          sound = dsqrt(gm)
c        else if (a.ge.3d0) then
c
c  one could put in physical sound velocity in the mixed phase, but this
c  causes instabilities for more or less stationary rarefaction shock waves.
c
c          sound = 1d-6
c          sound = dsqrt(gm)
c        else
c          sound = dsqrt(gm)
c        end if
c      end if
      return
      end
c
c-----------------------------------------------------------------------
c
      function minmod(a,b)
c
c  This function-subprogram calculates the minmod-function defined in
c  eq.(A.2).
c  It is used in the function-subprogram slope.
c
c  Type declarations for function minmod
c
      real*8 a,b,minmod
      if (a*b.le.0d0) then
        minmod = 0d0
      else
        if (a.lt.0d0) then
          minmod = dmax1(a,b)
        else
          minmod = dmin1(a,b)
        end if
      end if
      return
      end
c
c----------------------------------------------------------------------
c
      function slope(a,b,ks)
c
c  This function-subprogram calculates the slope approximation
c  according to appendix A.
c  It is used in subroutine prop.
c  a is the difference of the conserved quantity U at the surface
c  between cell i and i-1, divided by dx.
c  b is the corresponding expression between cell i+1 and i.
c  ks is the corresponding slope approximation factor.
c
c  Type-declarations for variables used in function-subprogram slope
c
      real*8 a,b,ks,slope
      real*8 minmod
      slope = minmod(a,b)*ks
      return
      end
c
c---------------------------------------------------------------------
c
      subroutine prop(elab,m,ks)
c
c  This subroutine-subprogram propagates the laboratory frame
c  quantities Elab,M. It is used in the main program.
c  Since a second order scheme is used, laboratory frame quantities
c  are computed at cell interfaces.
c
c  Type declarations for variables in common-blocks.
c
      real*8 dx,dt,t
      real*8 eta
      real*8 alp
      integer kp(1000),km(1000)
c
c  Type declarations for variables in subroutine prop.
c
      real*8 elab(1000),m(1000),elabt(1000),mt(1000),pos
      real*8 dd,sslope,bp,bm,flux,csmean,vmean,sqrl,sqrr,slope
      real*8 denomi,sound
      real*8 csp(1000),csm(1000)
      real*8 a1(1000),a2(1000),a3(1000)
      real*8 elm(1000),mm(1000),em(1000),
     1       pm(1000),vm(1000),gamm(1000)
      real*8 elp(1000),mp(1000),ep(1000),
     1       pp(1000),vp(1000),gamp(1000)
      real*8 dpv(1000),gv(1000)
      real*8 fem(1000),fmm(1000)
      real*8 fep(1000),fmp(1000)
      real*8 ks
      integer i
c
c  common-blocks
c
      common /indpoi/ kp,km
      common /gitter/ dx,dt,t
      common /visko/ eta
      common /geom/ alp
c
c  Start calculation. dd gets its half-step value
c
      dd = 0.5d0*dt/dx
c
c  Step 1: Slope and interface values are calculated.
c
c  For Elab:
c
      do 2000 i = 1,1000
        dpv(i) = ( elab(kp(i)) - elab(i) )/dx
2000  continue
      do 2001 i = 1,1000
        sslope = 0.5d0*dx*slope(dpv(km(i)),dpv(i),ks)
        elp(i) = elab(i) + sslope
        elm(i) = elab(i) - sslope
2001  continue
c
c  For M:
c
      do 2004 i = 1,1000
        dpv(i) = ( m(kp(i))-m(i) )/dx
2004  continue
      do 2005 i = 1,1000
        sslope = 0.5d0*dx*slope(dpv(km(i)),dpv(i),ks)
        if ( m(i)+sslope .gt. 0d0 ) then
          mp(i) = dmin1( m(i)+sslope , elp(i) )
        else
          mp(i) = dmax1( m(i)+sslope , - elp(i) )
        end if
        if ( 2d0*m(i) - mp(i) .gt. 0d0 ) then
          mm(i) = dmin1( 2d0*m(i) - mp(i), elm(i) )
        else
          mm(i) = dmax1( 2d0*m(i) - mp(i), - elm(i) )
        end if
2005  continue
c
c  Step 2: Calculate rest frame variables via untang.
c
      call untang(elp,mp,ep,pp,vp,gamp)
      call untang(elm,mm,em,pm,vm,gamm)
c
c  Step 3: Calculate interface values in half step.
c
      do 2006 i = 1,1000
        flux = dd*(mp(i)*vp(i) + pp(i) - mm(i)*vm(i) - pm(i))
        mp(i) = mp(i)-flux
        mm(i) = mm(i)-flux
        flux = dd*( (pp(i) + elp(i))*vp(i) - (pm(i) + elm(i))*vm(i))
        elp(i) = elp(i)-flux
        elm(i) = elm(i)-flux
c
c  Check physical consistence of half step values.
c
        if (elp(i)-dabs(mp(i)).lt.0d0) then
          if (dabs(mp(i)).gt.1d-16) then
            mp(i) = elp(i)*mp(i)/dabs(mp(i))
c           write(6,*) '  M+(',i,') changed to Elab+ at Step 3'
c
c  After correction, linear slope is enforced (but not if it
c  produces acausalities in elm(i), mm(i) ).
c
            if (mm(i).gt.0d0) then
              mm(i) = dmin1( 2d0*m(i) - mp(i) , elm(i) )
            else
              mm(i) = dmax1( 2d0*m(i) - mp(i) , -elm(i) )
            end if
          else
            mp(i) = 0d0
          end if
        end if
        if (elm(i)-dabs(mm(i)).lt.0d0) then
          elm(i) = dabs(mm(i))
c         write(6,*) '  Elab-(',i,') changed to M- at Step 3'
c
c  After correction, linear slope is enforced (but not if it
c  produces acausalities in elp(i), mp(i) ).
c
          elp(i) = dmax1( 2d0*elab(i) - elm(i) , dabs(mp(i)) )
        end if
 2006 continue
c
c  Step 4: Update local rest frame variables at cell interfaces
c  via untang.
c
      call untang(elp,mp,ep,pp,vp,gamp)
      call untang(elm,mm,em,pm,vm,gamm)
c
c  Step 5: Calculate sound velocity at the interface.
c
      do 2007 i=1,1000
        csp(i) = sound(ep(i))
        csm(i) = sound(em(kp(i)))
 2007 continue
c
c  Step 6: Estimate signal velocities.
c
      do 2008 i = 1,1000
c
c  Roe-mean for the velocities
c
        sqrl = dsqrt(elp(i))
        sqrr = dsqrt(elm(kp(i)))
        denomi = sqrl + sqrr
        if(denomi.gt.1d-16) then
         vmean = ( vp(i)*sqrl + vm(kp(i))*sqrr )/denomi
         csmean = dsqrt((csp(i)**2*sqrl + csm(kp(i))**2*sqrr)/denomi +
     1     eta*sqrl*sqrr/denomi/denomi*( vp(i) - vm(kp(i)) )**2 )
         bp = dmax1(0d0 , (vmean + csmean)/(1d0 + vmean*csmean) ,
     1          (vm(kp(i)) + csm(kp(i)))/(1d0 + vm(kp(i))*csm(kp(i))))
         bm = dmin1(0d0 , (vmean - csmean)/(1d0 - vmean*csmean) ,
     1          (vp(i) - csp(i))/(1d0 - vp(i)*csp(i)))
        else
         bp = 1d0
         bm = -1d0
        end if
c
c  Step 7, first part: calculate coefficients a.
c
        a1(i) = bp/(bp-bm)
        a2(i) = bm/(bp-bm)
        a3(i) = bp*bm/(bp-bm)
 2008 continue
c
c  Step 7, second part: calculate fluxes.
c
      do 2009 i=1,1000
        fep(i) = (pp(i) + elp(i))*vp(i)
        fmp(i) = mp(i)*vp(i) + pp(i)
        fem(i) = (pm(i) + elm(i))*vm(i)
        fmm(i) = mm(i)*vm(i) + pm(i)
 2009 continue
c
c Step 8 and 9 : Calculation of the total numerical fluxes and the
c                transported quantities at full time step.
c
      dd = dd+dd
      do 2012 i = 1,1000
        gv(i) = a1(i)*fep(i)-a2(i)*fem(kp(i))+a3(i)*(em(kp(i))-ep(i))
 2012 continue
      do 2013 i = 1,1000
        elabt(i) = elab(i) + dd*(gv(km(i))-gv(i))
 2013 continue
      do 2014 i = 1,1000
        gv(i) = a1(i)*fmp(i)-a2(i)*fmm(kp(i))+a3(i)*(mm(kp(i))-mp(i))
 2014 continue
      do 2015 i = 1,1000
        mt(i) = m(i) + dd*(gv(km(i))-gv(i))
 2015 continue
c
c  Step 9a: remove acausalities
c
      do 2116 i = 1,1000
        if (elabt(i)-dabs(mt(i)).lt.0d0) then
          if (dabs(mt(i)).gt.1d-16) then
            mt(i) = elabt(i)*mt(i)/dabs(mt(i))*(1d0-1d-10)
c           write(6,*) '  Mt(',i,') changed to Elabt at Step 9a'
          else
            mt(i) = 0d0
          end if
        end if
 2116 continue
c
c  Step 9b: Sod's method to correct for geometry
c           First call untangle, since the source terms
c           contain variables in the local rest frame.
c           These variables are temporarily stored in
c           the fields ep,pp,vp,gamp.
c      
      call untang(elabt,mt,ep,pp,vp,gamp)
      do 2115 i = 1,1000 
          pos = (dfloat(i) - 500.5d0)*dx
          if(alp.eq.3d0) then
            elab(i) = elabt(i) - dt*(vp(i)/pos+1d0/t)*(elabt(i)+pp(i))
            m(i) = mt(i) - dt*(vp(i)/pos+1d0/t)*mt(i)
          else 
            elab(i) = elabt(i) - dt*alp*vp(i)/pos*(elabt(i)+pp(i))
            m(i) = mt(i) - dt*alp*vp(i)/pos*mt(i)
          end if
 2115 continue
c
c  Step 10: The last consistency check is to remove acausalities
c           produced in the full step propagation.
c
      do 2016 i = 1,1000
        if (elab(i)-dabs(m(i)).lt.0d0) then
          if (dabs(m(i)).gt.1d-16) then
            m(i) = elab(i)*m(i)/dabs(m(i))*(1d0-1d-10)
c           write(6,*) '  M(',i,') changed to Elab at Step 10'
          else
            m(i) = 0d0
          end if
        end if
 2016 continue
c
c  Step 11: Cells with smaller lab energy density than a cut-off
c           are regarded as vacuum.
c
      do 2017 i = 1,1000
        if (elab(i).lt.1d-16) then
          elab(i) = 0d0
          m(i) = 0d0
        end if
 2017 continue
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine fileo(elab,m,e,p,v,gamma,it,ind)
c
c  This subroutine-subprogram writes data on File
c  It is used in the main program.
c
c  Type declaration for variables used in common-blocks
c
      real*8 gm,B
      real*8 dx,dt,t
      real*8 total0(2),alp
      integer eos
c
c  Type declaration for variables used in fileo
c
      real*8 e(1000),p(1000),v(1000),gamma(1000)
      real*8 elab(1000),m(1000)
      real*8 total(2),pos(1000)
      real*8 ttc,temp
      integer i,it,index,ind
c
c  common-blocks
c
      common /adind/ gm,B
      common /gitter/ dx,dt,t
      common /conser/ total0
      common /flags/ eos
      common /geom/ alp
c
c  Total is a vector which contains the sum of the primary variables
c  over all cells. It is set zero in the beginning.
c
      do 50 i = 1,2
        total(i) = 0d0
 50   continue
c
c  Pos is a vector whose element i contains the centre position of
c  cell no. i.
c
      pos(1) = 0.5d0*dx
      do 49 i = 2,1000
c        pos(i) = pos(i-1) + dx
         pos(i)= (dble(i)-.5d0)*dx
 49   continue
c
c  Start writing.
c
      index = 7 + it/ind
      write(6,*) ' time=',t,' file index=',index-(7-1)
c
c  Write quantities for all cells with pos(i)>0.
c
      do 51 i = 501,1000
        ttc=temp(e(i))
        write(index,7777) pos(i)-5d2*dx,ttc,elab(i),v(i),
     ?                      e(i),m(i) 

        if (alp.ne.0d0) then
          total(1) = total(1) + 
     1         elab(i)*(dx*dabs(dfloat(i)-500.5d0))**alp*dx
          total(2) = total(2) +
     1         m(i)*(dx*dabs(dfloat(i)-500.5d0))**alp*dx
        else
          total(1) = total(1) + 
     1         elab(i)*dx
          total(2) = total(2) +
     1         m(i)*dx
        end if
 51   continue
c      write(6,*) ' conservation law balance'
c      write(6,*) '      E             M'
c      write(6,7777) ((total(i)-total0(i))/total0(1),i=1,2)
      return
 7777 format(1x,6(e13.5,1x))
      end
c----------------------------------------------------------------------
c
c This subroutine reads the pressure/temp. table from file "eos.dat"
c
       subroutine PressInit
       implicit real*8 (a-h,o-z)
c
c Common Blocks
c
      common /pressuretab/ de,presstab(0:5000),
     &			temptab(0:5000)

	open(unit=81,file='eos.dat')
	read(81,*) de
	do i=0,5000
	  read(81,*) d1,presstab(i),d2,temptab(i)
	enddo
	close(81)
	return
	end
c
c----------------------------------------------------------------------
c  This subroutine writes hypersurface-data on File.
c  It is called by the main program.
c  (write out cells on h.s. in each timestep)
c
      subroutine wrhyps(e,v,tc)
      implicit real*8 (a-h,o-z)
c
c
c  Type declaration for variables used in common-blocks.
c
      real*8 gm,B
      real*8 dx,dt,t
c
c  Type declaration for variables used in writc.
c
      real*8 e(1000),v(1000),pos
      real*8 tc,ttc,ni
      integer i,index
c
c  Common-blocks.
c
      common /adind/ gm,B
      common /gitter/ dx,dt,t


      pos = -(1000d0*0.5d0+0.5d0)*dx
      index = 1
      pos = pos + dx
      do 51 i = 2,1000
        pos = pos + dx
c
c  Find points which have given energy density
	ttc=e(i)
	ttcm=e(i-1)

c  Search routine uses the fact that for given time T (or e)
c  first increases for increasing cell index.
c
        if ((ttc.ge.tc).and.(index.eq.1)) then
	  j=abs(i-500)
	  write(30,77) pos,t,e(i),v(i),0d0
          index = 2
        else if ((ttc.lt.tc).and.(index.eq.2)) then
	  j=abs(i-500)
	  del_x=(tc-ttc)/(ttc-ttcm)	! interpolate to position between
	  posi=pos+del_x*dx		! the cells
	  ei=e(i)+del_x*(e(i)-e(i-1))
	  vi=v(i)+del_x*(v(i)-v(i-1))
	  write(30,77) posi,t,ei,vi,0d0
          index = 1
        end if
 51   continue

      return
 77   FORMAT(5G14.6)
      end
c
c----------------------------------------------------------------------
c
c sort points from file fort.30, save hypersurface to hyper.dat
c called by main program
c
      subroutine sort(t0)
      implicit real*8 (a-h,o-z)
      real*8 dx,dt
      real*8 x(10000),t(10000),v(10000),abs(10000)
	real*8 e(10000),rho(10000)
	real*8 pos,time,vel,t0,eps
      real*8 xn(10000),tn(10000),vn(10000),min,xx
	real*8 en(10000),rhon(10000)
      real*8 xn2(10000),tn2(10000),vn2(10000)
	real*8 en2(10000),rhon2(10000)
      integer i,j,k,l,n,ind
      common /gitter/ dx,dt,dummy
c
c
c  Read in cells on hypersurface.
c
      OPEN(30,FILE='fort.30',STATUS='OLD')
      OPEN(40,FILE='hyper.dat')
      n = 0
      do 14 i=1,10000
        read(30,*,END=20) pos,time,eps,vel,rhob
        if (pos.ge.0d0) then
          n = n+1
          x(n) = pos
          t(n) = time
          v(n) = vel
	  e(n)=eps
	  rho(n)=rhob
        end if
 14   continue
 20   continue
      if (n.eq.0) then
        return
      end if
c
c  Sort (x,t) combinations.
c
      xn(1) = x(1)
      tn(1) = t(1)
      vn(1) = v(1)
      en(1) = e(1)
      rhon(1)= rho(1)

      do 30 i=1,n-1
c
c  Calculate vector of distances between remaining points
c  and point i.
c
        do 40 j=i+1,n
          abs(j) = (tn(i)-t(j))**2+(xn(i)-x(j))**2
 40     continue
        min = abs(i+1)
        l = i+1
        do 50 k=i+2,n
          if(abs(k).lt.min) then
            min = abs(k)
            l = k
          end if
 50     continue

        xn(i+1) = x(l)
        tn(i+1) = t(l)
        vn(i+1) = v(l)
        en(i+1) = e(l)
	rhon(i+1)= rho(l)

        x(l) = x(i+1)
        t(l) = t(i+1)
        v(l) = v(i+1)
        e(l) = e(i+1)
	rho(l)= rho(i+1)

 30   continue
c
c  If positions are the same for different times, take
c  only the arithmetic mean time at that position.
c
      write(40,17) 1d0,t0,0d0,0d0,0d0,0d0,0d0
      i = 1
      l = 1

 31   xn2(i) = xn(l)
      tn2(i) = tn(l)
      en2(i) = en(l)
      vn2(i) = vn(l)
      rhon2(i)= rhon(l)

      count = 1
c
c  Search next 100 entries if position is equal.
c
      do 32 j = 1,100
        if (xn(j+l).eq.xn(l)) then
          tn2(i) = tn2(i)+tn(j+l)
          en2(i) = en2(i)+en(j+l)
          vn2(i) = vn2(i)+vn(j+l)
          rhon2(i)= rhon2(i)+rhon(j+l)
          count = count + 1
        else
          goto 33
        end if
 32   continue
 33   continue
      tn2(i) = tn2(i)/count
      en2(i) = en2(i)/count
      vn2(i) = vn2(i)/count
      rhon2(i)= rhon2(i)/count

c      fo_temp=temp(en2(i),rhon2(i))
      fo_temp=temp(en2(i))
      fo_muq=0d0
      fo_mus=0d0

      write(40,17) xn2(i),tn2(i),en2(i),vn2(i),
     &			fo_temp,fo_muq,fo_mus

      i = i + 1
      l = l + count
      if (l.lt.n) goto 31

      if (xn2(i-1).gt.dx) then
        xx = xn2(i-1)-dx
 34     continue
        write(40,17) xx,tn2(i-1),en2(i-1),0d0,fo_temp,fo_muq
     &				,fo_mus
        if (xx.gt.dx) then
          xx = xx-dx
          goto 34
        end if
      end if
      write(40,17) 0d0,tn2(i-1),0d0,0d0,0d0,0d0,0d0
      CLOSE(30,STATUS='DELETE')
      CLOSE(40)

      return
 17   FORMAT(8G14.6)
      END
c
c----------------------------------------------------------------------
c
c     main program
c
c  Type declarations for variables in common-blocks
c
      real*8 gm,B
      real*8 dx,dt,t
      real*8 eta
      integer eos
c
c  Type declarations for variables used in main
c
      real*8 e(1000),p(1000),v(1000),gamma(1000)
      real*8 elab(1000),m(1000),t0
      real*8 ks
      integer itplot,it,ind
      real*8 eps_hyper
c
c  common-blocks
c
      common /adind/ gm,B
      common /gitter/ dx,dt,t
      common /visko/ eta
      common /flags/ eos
c
c  First part of initialization (continued in sinit)
c
c  Initialize t0 and number of time steps
c
      open(unit=5,file='par.dat')
      read(5,*) t0,itplot
      call sinit(e,gamma,ks)
      if (eos.eq.2) call PressInit
      call init(e,p,v,gamma,elab,m)
      read(5,*) eps_hyper
      close(5)

      t = t0
      open(unit=6,file='rhllew00.dat')
      open(unit=7,file='rhllew01.dat')
      open(unit=8,file='rhllew02.dat')
      open(unit=9,file='rhllew03.dat')
      open(unit=10,file='rhllew04.dat')
      open(unit=11,file='rhllew05.dat')
      open(unit=12,file='rhllew06.dat')
      open(unit=13,file='rhllew07.dat')
      open(unit=14,file='rhllew08.dat')
      open(unit=15,file='rhllew09.dat')
      open(unit=16,file='rhllew10.dat')
      open(unit=17,file='rhllew11.dat')
      open(unit=18,file='rhllew12.dat')
      open(unit=19,file='rhllew13.dat')
      open(unit=20,file='rhllew14.dat')
      open(unit=21,file='rhllew15.dat')
      open(unit=22,file='rhllew16.dat')
      open(unit=23,file='rhllew17.dat')
      open(unit=24,file='rhllew18.dat')
      open(unit=25,file='rhllew19.dat')
      open(unit=26,file='rhllew20.dat')
      open(unit=27,file='rhllew21.dat') 
c intermediate hypersurface file
      open(unit=30,file='fort.30')

      write(6,1212)  eta,gm
      write(6,1213)  itplot,dx,dt
      ind = itplot/50
      call fileo(elab,m,e,p,v,gamma,0,ind)

      do 10 it = 1,itplot    ! main loop (time-steps)
        t = t+dt
        call prop(elab,m,ks)
        call untang(elab,m,e,p,v,gamma)
        if ((it.eq.ind).or.(it.eq.ind*2).or.(it.eq.ind*3).or.
     1      (it.eq.ind*4).or.(it.eq.ind*5).or.(it.eq.ind*6).or.
     2      (it.eq.ind*7).or.(it.eq.ind*8).or.(it.eq.ind*9).or.
     3      (it.eq.ind*10).or.(it.eq.ind*11).or.(it.eq.ind*12).or.
     4      (it.eq.ind*13).or.(it.eq.ind*14).or.(it.eq.ind*15).or.
     5      (it.eq.ind*16).or.(it.eq.ind*17).or.(it.eq.ind*18).or.
     6      (it.eq.ind*19).or.(it.eq.ind*20)) then
          call fileo(elab,m,e,p,v,gamma,it,ind)
        end if

c write out cells on the hypersurface
        call wrhyps(e,v,eps_hyper)
 10   continue    ! main loop (time-steps)

      close(6)      
      close(7)
      close(8)
      close(9)
      close(10)
      close(11)
      close(12)
      close(13)
      close(14)
      close(15)
      close(16)
      close(17)
      close(18)
      close(19)
      close(20)
      close(21)
      close(22)
      close(23)
      close(24)
      close(25)
      close(26)
      close(27)
      close(30)

      call sort(t0)
      stop
 1212 format(1x,' et=',f8.3,'  Gamma-1=',e16.8)
 1213 format(1x,' max. no. of time steps=',i4,' dx=',f5.2,' dt=',f5.2)
      end
