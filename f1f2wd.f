      SUBROUTINE RESCSP(w2,q2,sigtp,siglp,dsigtp,dsiglp)
      IMPLICIT none

      real*8 w2,q2,sigtp,dsigtp,siglp,dsiglp
      real*8 xvalp(100),xval1(50),xvalL(50)
      Integer i
      real*8 sigtdis,sigLdis

      data xvalp / 

     & 0.12298E+01,0.15304E+01,0.15057E+01,0.16980E+01,0.16650E+01,
     & 0.14333E+01,0.13573E+00,0.22000E+00,0.82956E-01,0.95782E-01,
     & 0.10936E+00,0.37944E+00,0.77805E+01,0.42291E+01,0.12598E+01,
c     & 0.21242E+01,0.63351E+01,0.68232E+04,0.33521E+05,0.25686E+01, !Oleksandr Koshchii
     & 0.21242E+01,0.63351E+01,0.11400E+01,0.55900E+01,0.25686E+01,
     & 0.60347E+00,0.21240E+02,0.55746E-01,0.24886E+01,0.23305E+01,
     & -.28789E+00,0.18607E+00,0.63534E-01,0.19790E+01,-.56175E+00,
     & 0.38964E+00,0.54883E+00,0.22506E-01,0.46213E+03,0.19221E+00,
     & 0.19141E+01,0.24606E+03,0.67469E-01,0.13501E+01,0.12054E+00,
     & -.89360E+02,0.20977E+00,0.15715E+01,0.90736E-01,-.38495E-02,
     & 0.10362E-01,0.19341E+01,0.38000E+00,0.34187E+01,0.14462E+00,
     & 0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,
     & 0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,0.00000E+00,
     & 0.00000E+00,0.00000E+00,0.29414E+02,0.19910E+02,0.22587E+00,
     & 0.00000E+00,0.00000E+00,0.38565E+04,0.65717E+00,0.00000E+00,
     & 0.15792E+03,0.97046E+02,0.31042E+00,0.00000E+00,0.42160E+01,
     & 0.38200E-01,0.12182E+01,0.00000E+00,0.13764E+02,0.31393E+00,
     & 0.29997E+01,0.00000E+00,0.55124E+01,0.53743E-01,0.13091E+01,
     & 0.00000E+00,0.86746E+02,0.40864E-04,0.40294E+01,0.31285E+01,
     & 0.33403E+00,0.49623E+01,0.00000E+00,0.00000E+00,0.11000E+02,
     & 0.18951E+01,0.51376E+00,0.00000E+00,0.42802E+01,0.00000E+00 /

      do i=1,50
        xval1(i) = xvalp(i)
        xvalL(i) = xvalp(50+i) 
        if(i.LE.12) xvalL(i) = xval1(i)  !!! use same masses for L as T !!!
      enddo
      xvalL(43) = xval1(47)
      xvalL(44) = xval1(48)
      xvalL(50) = xval1(50)
   
      call resmodp(1,w2,q2,xval1,sigTp,dsigTp)
      call resmodp(2,w2,q2,xvalL,sigLp,dsigLp)
c      call disp(w2,q2,sigTdis,sigLdis)
      if(w2.GT.8.0) then 
c        write(6,*) "resmod: ",w2,q2,sigTp,sigLp
        call disp(w2,q2,sigTdis,sigLdis)
c        write(6,*) "dismod: ",w2,q2,sigtdis,sigLdis
        if(w2.LE.10.0) then
         sigTp = (10.0-w2)*sigTp+(w2-8.0)*sigTdis
         sigLp = (10.0-w2)*sigLp+(w2-8.0)*sigLdis
         sigTp = sigTp/2.0
         sigLp = sigLp/2.0
         else
          sigTp = sigTdis
          sigLp = sigLdis
        endif
      endif

      return
      end




      SUBROUTINE RESMODP(sf,w2,q2,xval,sig,dsig) 

      IMPLICIT NONE
      REAL*8 W,w2,q2,q22,mp,mp2,mpi2,xb,dxb,xth(4),sig,dsig,xval(50),mass(7),width(7)
      REAL*8 height(7),heightprime(7),sig_del,sig_21,sig_22,sig_31,sig_32,rescoef(6,4)
      REAL*8 nr_coef(3,4),sigr(7),dsigr(7),wdif(2),sig_nr,dsig_nr,sig_4
      REAL*8 mpi,meta,intwidth(7),k,kcm,kr(7),kcmr(7),ppicm,ppi2cm
      REAL*8 petacm,ppicmr(7),ppi2cmr(7),petacmr(7),epicmr(7),epi2cmr(7)
      REAL*8 eetacmr(7),epicm,epi2cm,eetacm,br(7,3),spin(7),ang(7)
      REAL*8 pgam(7),pwid(7,3),x0(7),dip,mon,q20,h_nr(3),dh_nr(3)
      REAL*8 sig_res,dsig_res,sig_4L,sigtemp,slope,t,dt,xpr(2),dxpr(2),m0
      real*8 sigrsv(7),sig_nrsv
      INTEGER i,j,l,num,sf
      real*8 sig_mec
      logical first/.true./
      common/tst2/sigrsv,sig_nrsv,sig_mec


      mp = 0.9382727
      mpi = 0.135
      mpi2 = mpi*mpi
      meta = 0.547
      mp2 = mp*mp
      W = sqrt(w2)
      wdif(1) = w - (mp + mpi)
      wdif(2) = w - (mp + 2.*mpi)

      m0 = 0.125
      if(sf.EQ.2) m0 = xval(49)

      if(sf.EQ.1) then
        q20 = 0.05
      else
        q20 = 0.125
      endif 
   

CCCC   single pion branching ratios  CCCC

      br(1,1) = 1.0       !!!  P33(1232)       
      br(2,1) = 0.45      !!!  S11(1535)   
      br(3,1) = 0.65      !!!  D13(1520)
      br(4,1) = 0.65      !!!  F15(1680)
      br(5,1) = 0.4       !!!  S11(1650)
      br(6,1) = 0.65      !!!  P11(1440) roper 
      br(7,1) = 0.50      !!!  F37(1950)

CCCC  eta branching ratios   CCCC

      br(1,3) = 0.0       !!!  P33(1232)
      br(2,3) = 0.45      !!!  S11(1535) 
      br(3,3) = 0.0       !!!  D13(1520)
      br(4,3) = 0.0       !!!  F15(1680)
      br(5,3) = 0.1       !!!  S11(1650)
      br(6,3) = 0.0       !!!  P11(1440) roper   
      br(7,3) = 0.0       !!!  F37(1950)

CCCC  2-pion branching ratios  CCCC

      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo


CCCC   Meson angular momentum   CCCC



      ang(1) = 1.       !!!  P33(1232)
      ang(2) = 0.       !!!  S11(1535)
      ang(3) = 2.       !!!  D13(1520)
      ang(4) = 3.       !!!  F15(1680)
      ang(5) = 0.       !!!  S15(1650)
      ang(6) = 1.       !!!  P11(1440) roper   
      ang(7) = 3.       !!!  F37(1950)

      do i=1,7     !!!  resonance damping parameter  !!!
        x0(i) = 0.215
c        x0(i) = xval(50)
      enddo

      x0(1) = 0.15
      x0(1) = xval(50)   

      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo
 
        q22=q2!0.0!q2!OK    
        
 
      dip = 1./(1.+q22/0.71)**2.             !!!  Dipole parameterization  !!!
      mon = 1./(1.+q22/0.71)**1.

  
      xb = q22/(q22+w2-mp2)
      dxb = (w2-mp2)/(q22+w2-mp2)**2
      xpr(1) = 1.+(w2-(mp+mpi)**2)/(q22+q20)
      xpr(1) = 1./xpr(1)
      dxpr(1) = xpr(1)**2*(w2-(mp+mpi)**2)/(q22+q20)**2 !OK 
      xpr(2) = 1.+(w2-(mp+mpi+mpi)**2)/(q22+q20)
      xpr(2) = 1./xpr(2)


      t = log(log((q22+m0)/0.330**2)/log(m0/0.330**2))
      dt = 1./(m0+q22)/log((q22+m0)/0.330**2)
CCC    Calculate kinematics needed for threshold Relativistic B-W  CCC

      k = (w2 - mp2)/2./mp
      kcm = (w2-mp2)/2./w

      epicm = (W2 + mpi**2 -mp2 )/2./w
      ppicm = SQRT(MAX(0.0,(epicm**2 - mpi**2)))
      epi2cm = (W2 + (2.*mpi)**2 -mp2 )/2./w
      ppi2cm = SQRT(MAX(0.0,(epi2cm**2 - (2.*mpi)**2)))
      eetacm = (W2 + meta*meta -mp2 )/2./w
      petacm =  SQRT(MAX(0.0,(eetacm**2 - meta**2)))

      num = 0

      do i=1,6              !!!  Read in resonance masses     !!!
        num = num + 1
        mass(i) = xval(i)
      enddo
      do i=1,6              !!!  Read in resonance widths     !!!
        num = num + 1
        intwidth(i) = xval(num)
        width(i) = intwidth(i)
      enddo

      if(sf.EQ.2) then      !!!  Put in 4th resonance region  !!!
        mass(7) = xval(43)
        intwidth(7) = xval(44)
        width(7) = intwidth(7)
      else
        mass(7) = xval(47)
        intwidth(7) = xval(48)
        width(7) = intwidth(7) 
      endif

      do i=1,7
        kr(i) = (mass(i)**2-mp2)/2./mp
        kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
        epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
        ppicmr(i) = SQRT(MAX(0.0,(epicmr(i)**2 - mpi**2)))
        epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
        ppi2cmr(i) = SQRT(MAX(0.0,(epi2cmr(i)**2 - (2.*mpi)**2)))
        eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
        petacmr(i) =  SQRT(MAX(0.0,(eetacmr(i)**2 - meta**2)))

CCC   Calculate partial widths   CCC

        pwid(i,1) = intwidth(i)*(ppicm/ppicmr(i))**(2.*ang(i)+1.)
     &           *((ppicmr(i)**2+x0(i)**2)/(ppicm**2+x0(i)**2))**ang(i)
c         !!!  1-pion decay mode


        pwid(i,2) = intwidth(i)*(ppi2cm/ppi2cmr(i))**(2.*ang(i)+4.)
     &         *((ppi2cmr(i)**2+x0(i)**2)/(ppi2cm**2+x0(i)**2))
     &         **(ang(i)+2)   !!!  2-pion decay mode

        pwid(i,2) = W/mass(i)*pwid(i,2)


        pwid(i,3) = 0.          !!!  eta decay mode


        if(i.EQ.2.OR.i.EQ.5) then
          pwid(i,3) =  intwidth(i)*(petacm/petacmr(i))**(2.*ang(i)+1.)
     &          *((petacmr(i)**2+x0(i)**2)/(petacm**2+x0(i)**2))**ang(i)
c         !!!  eta decay only for S11's 
        endif 



        pgam(i) = (kcm/kcmr(i))**2*
     &                   (kcmr(i)**2+x0(i)**2)/(kcm**2+x0(i)**2)

        pgam(i) = intwidth(i)*pgam(i)

        width(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)+br(i,3)*pwid(i,3)

      enddo

CCC    End resonance kinematics and Widths calculations   CCC


CCC    Begin resonance Q^2 dependence calculations   CCC
           
      do i=1,6
        do j=1,4
          num = num + 1
          rescoef(i,j)=xval(num)
        enddo

        if(sf.EQ.1) then
 
          height(i) = rescoef(i,1)*
     &                (1.+rescoef(i,2)*q22/(1.+rescoef(i,3)*q22))/
     &                (1.+q22/0.91)**rescoef(i,4)
CCC     Modification by Oleksandr Koshchii
          !A^i_T prime below
          heightprime(i) = rescoef(i,1)*
     &                     (1./(1.+q22/0.91)**rescoef(i,4)*
     &                     (rescoef(i,2)/(1.+rescoef(i,3)*q22)-
     &          rescoef(i,2)*rescoef(i,3)*q22/(1.+rescoef(i,3)*q22)**2)-
     &          rescoef(i,4)/0.91*(1.+q22/0.91)**(-1.-rescoef(i,4))*
     &          (1.+rescoef(i,2)*q22/(1.+rescoef(i,3)*q22))
     &         )
CCC     End of the modification
        else

          height(i) = rescoef(i,1)*q22/(1.+rescoef(i,2)*q22)
     &                             *exp(-1.*rescoef(i,3)*q22)
          heightprime(i) = rescoef(i,1)*exp(-1.*rescoef(i,3)*q22)*
     &                     (1./(1.+rescoef(i,2)*q22)**2
     &                      -rescoef(i,3)*q22/(1.+rescoef(i,2)*q22))
        endif

        heightprime(i)=2.*heightprime(i)*height(i)!OK dsig/dq2 ~ 2*Ai*Aipr 
        height(i) = height(i)*height(i) !sig ~ Ai*Ai
        
      enddo
                  
      if(sf.EQ.2) then      !!!  4th resonance region  !!!
        height(7) = xval(45)*q22/(1.+xval(46)*q22)*exp(-1.*xval(47)*q22)
        heightprime(7)=xval(45)*exp(-1.*xval(47)*q22)*
     &                     (1./(1.+xval(46)*q22)**2
     &                      -xval(47)*q22/(1.+xval(46)*q22))
      else
        height(7) = xval(49)/(1.+q22/0.91)**1. 
        heightprime(7) = -xval(49)/0.91/(1.+q22/0.91)**2
      endif
      heightprime(7)=2.*heightprime(7)*height(7)
      height(7) = height(7)*height(7)

      !if(sf.EQ.2) then
      !      write(*,*)'h1=',heightprime(1)!height(1)!
      !      write(*,*)'h2=',heightprime(2)!height(2)!
      !      write(*,*)'h3=',heightprime(3)!height(3)!
      !      write(*,*)'h4=',heightprime(4)!height(4)!
      !      write(*,*)'h5=',heightprime(5)!height(5)!
      !      write(*,*)'h6=',heightprime(6)!height(6)!
      !      write(*,*)'h7=',heightprime(7)!height(6)!
      !end if

        
CCC    End resonance Q^2 dependence calculations   CCC <--checked OK

     
      do i=1,3               !!!  Non-Res coefficients  !!!
        do j=1,4
          num = num + 1
          nr_coef(i,j)=xval(num)
        enddo
      enddo


CCC   Calculate Breit-Wigners for all resonances   CCC

      sig_res = 0.0
      dsig_res = 0.0
      
      do i=1,7
        sigr(i) = width(i)*pgam(i)/((W2 - mass(i)**2.)**2. 
     &              + (mass(i)*width(i))**2.)
        sigr(i) = height(i)*kr(i)/k*kcmr(i)/kcm*sigr(i)/intwidth(i)
CCC     OLeksandr Koshchii
        dsigr(i) = width(i)*pgam(i)/((W2 - mass(i)**2.)**2. 
     &              + (mass(i)*width(i))**2.)
        dsigr(i) = heightprime(i)*kr(i)/k*kcmr(i)/kcm*dsigr(i)/
     &             intwidth(i)
CCC     end Oleksandr Koshchii   
        if(sf.eq.1) sigrsv(i) = sigr(i)
        sig_res = sig_res + sigr(i)
        dsig_res = dsig_res + dsigr(i)!Oleksandr Koshchii
      enddo

      sig_res = sig_res*w
      dsig_res = dsig_res*w!OLeksandr Koshchii
        !write(*,*) 'w2=',w2,'q2=',q22,'sigR=',sig_res,'dsigR=',dsig_res
CCC    Finish resonances / start non-res background calculation   CCC checled OK

 
      sig_nr = 0.
      dsig_nr = 0.
      
      if(sf.EQ.1) then
        do i=1,2  
          h_nr(i) = nr_coef(i,1)/     
     &       (q22+nr_coef(i,2))**
     &       (nr_coef(i,3)+nr_coef(i,4)*q22+xval(44+i)*q22**2)
          dh_nr(i) = nr_coef(i,1)*     
     &       (q22+nr_coef(i,2))**
     &       (-nr_coef(i,3)-nr_coef(i,4)*q22-xval(44+i)*q22**2)*
     &       ((-nr_coef(i,3)-nr_coef(i,4)*q22-xval(44+i)*q22**2)/
     &        (q22+nr_coef(i,2))-(nr_coef(i,4)+2.*xval(44+i)*q22)*
     &        log(q22+nr_coef(i,2)))     
          sig_nr = sig_nr + h_nr(i)*(wdif(1))**(float(2*i+1)/2.)
          dsig_nr = dsig_nr + dh_nr(i)*(wdif(1))**(float(2*i+1)/2.)
        enddo

        dsig_nr = dxpr(1)*sig_nr+xpr(1)*dsig_nr!OK
        sig_nr = sig_nr*xpr(1)
        sig_nrsv = sig_nr

      elseif(sf.EQ.2) then

        do i=1,1
          sig_nr = sig_nr + nr_coef(i,1)*
     &      (1.-xpr(i))**(nr_coef(i,3)+nr_coef(i,2)*t)
     &               /(1.-xb)/(q22+q20)
     &        *(q22/(q22+q20))**nr_coef(i,4)*xpr(i)**(xval(41)+xval(42)*t)
           
           !a=nr_coef(i,2)   b=nr_coef(i,3)      c=nr_coef(i,4)     d=xval(41)  e=xval(42)
          dsig_nr = dsig_nr + nr_coef(i,1)*
     &      1./(xb-1.)**2*(q22)**nr_coef(i,4)*(q22+q20)**(-2.-nr_coef(i,4))*
     &      (1.-xpr(i))**(nr_coef(i,3)+nr_coef(i,2)*t)*xpr(i)**(xval(41)+xval(42)*t)*
     &      ( (1.+nr_coef(i,4))*(-1.+xb)-nr_coef(i,4)*(q22+q20)*(xb-1.)/q22+
     &        (q22+q20)*dxb+(q22+q20)*(1.-xb)*(nr_coef(i,2)*log(1.-xpr(i))*dt+
     &        (nr_coef(i,3)+nr_coef(i,2)*t)*dxpr(1)/(xpr(i)-1.))+
     &        (q22+q20)*(1.-xb)*(xval(42)*log(xpr(i))*dt+(xval(41)+xval(42)*t)*
     &        dxpr(1)/xpr(i))
     &       )
        !write(*,*) 'sig0=',nr_coef(i,1),'a=',nr_coef(i,2),'b=',nr_coef(i,3),'c=',nr_coef(i,4),'d=',xval(41),'e=',xval(42)
        enddo
      !write(*,*) 'w2=',w2,'q2=',q22,'sigNR=',sig_nr,'dsigNR=',dsig_nr

      endif


      sig = sig_res + sig_nr
      dsig = dsig_res + dsig_nr
      
      if(w2.LT.1.159) sig = 0.0
      if(w2.LT.1.159) dsig = 0.0
      
      !if (sf.eq.1) then
      !  write(*,*) 'w2=',w2,'q2=',q22,'sigT=',sig,'dsigT=',dsig
      !else if (sf.eq.2) then
      !  write(*,*) 'w2=',w2,'q2=',q22,'sigL=',sig,'dsigL=',dsig
      !endif
 1000  format(8f12.5)

      RETURN 
      END 
  

      SUBROUTINE DISP(w2,q2,sigt,sigl)
      IMPLICIT none

      real*8 w2,q2,q2t,x,sigt,sigl,f1,f2,fL,r,dr,r1,r2
      Integer i,modt
      character*1 targ
      real*8 mp2,pi,pi2,alpha,t1,t2
      logical goodfit


      targ = 'P'
      modt = 12
      mp2 = 0.938272
      mp2 = mp2*mp2
      pi = 3.14159
      pi2 = pi*pi
      alpha = 1/137.0359

      x = q2/(w2+q2-mp2)

      call f2allm(x,q2,f2)
      call f2glob(x,q2,targ,modt,f2)
      call r1998(x,q2,r,dr,goodfit)
      if(q2.LE.0.15) then
        q2t = 0.15
        call r1998(x,q2t,r1,dr,goodfit)
        r = r1/0.25*q2
c        write(6,*) r1,r
        call f2allm(x,q2,f2)
      endif


      f1 = f2/2./x/(r+1.0)*(1.0+4.0*mp2*x*x/q2)
      fL = 2.*x*r*f1    

      sigt = 0.3894e3*f1*pi2*alpha*8.0/abs(w2-mp2)
      sigL = r*sigt

      return
      end



      SUBROUTINE SFNU(nu,q2,F1p,FLp,F2p,sigTp,sigLp,dsigTp,dsigLp)
cccc    This subroutine is created by Oleksandr Koshchii on April 2020 in order to add 
cccc    proton photoabsorption XS into calculation of beam asymmetry
CCCC   Converts reduced cross sections to structure functions for protons and neutrons  CCCCC 

      IMPLICIT none

      real*8 nu,w2,q2,x,sigtp,dsigtp,siglp,dsiglp,f1p,f2p,fLp,der
      real*8 pi,pi2,alpha,mp,mp2
      Integer i

      mp = 0.938272
      mp2 = mp*mp
      pi = 3.14159
      pi2 = pi*pi
      alpha = 1/137.03599 
      !x = q2/(q2+w2-mp2)
      !write(*,*)'w2=',w2
      !return
      w2=2.0d0*mp*nu+mp2-q2
      x = q2/(q2+w2-mp2)
      !write(*,*)'w2=',w2
      call rescsp(w2,q2,sigTp,sigLp,dsigTp,dsigLp)
      !w2=2.2081956d0
      !call xsder(w2,der)
      !write(*,*)'sing=',der!SignF(-1.0d0)
      
      f1p = sigTp/0.3894e3/pi2/alpha/8.0*abs(w2-mp2)!/0.3894e3
      fLp = sigLp*2.0*x/0.3894e3/pi2/alpha/8.0*abs(w2-mp2)!/0.3894e3
      f2p = (2.*x*f1p+fLp)/(1.+4.*mp2*x*x/q2)

      !The line below converts transverse and longitudinal XS to inverse GeV^2. Inseted by OK.
      sigTp=sigTp/0.3894e3
      dsigTp=dsigTp/0.3894e3
      sigLp=sigLp/0.3894e3
      dsigLp=dsigLp/0.3894e3
      return
      end


