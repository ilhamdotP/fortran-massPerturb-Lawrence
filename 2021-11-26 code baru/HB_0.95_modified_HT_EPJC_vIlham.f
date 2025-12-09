      PROGRAM TOVSOLVER                   
C     *********************************************************
C     
C     RMF BSP parameter set
c      for p-anisotropic slow rotating isotropic NS, 24 Nov 2015
C     R. L. Bowers and E. P. T. Liang, Astrophys. J 88, 657 (1974)
C     LBL=-2 and 2             
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z) 
C      INTEGER I, J, IM, NS, LI, N, IL, IN    
C      DIMENSION YA(10), EK(5,10), Y(10)
       INTEGER IL
C--------------------------------------------------------
C result save in file X.dat 
C--------------------
      OPEN (unit=8,STATUS='unknown',FILE='HB 0.95 BSP_MIAIBLX_1.txt')
      OPEN (unit=9,STATUS='unknown',FILE='HB 0.95 BSP_CRMIAIBLX_1.txt')
      OPEN (unit=10,STATUS='unknown',FILE='HB 0.95 MASSCORR&DEFORM.txt')

      HC= 197.327D0
      PI= 3.14159265358979D0      
      GS=1.325D-12
      MSS=1.1155D15
  
C---------------------------------------------------------------------
C  PCC is pressure in center
C  TUCT is trail initial condition for matric component nu

c       DO 10 IL=2,450,5
C        DO 10 IL=2,500,1
       DO 10 IL=2,500,1
       PCC=1.0D0*IL
       TUCT=0.1D-8

c       PCC=20D0
C------------------------------------------------------------------------
c Initial steps to calculate the correct initial condition for 
C matric component nu
C   ROS is radius in km
C   GMOS is mass in solar mass
C   MNURT is matric component nu in R
C   Rot is the correct initial condition for matric component nu
C
      CALL TOVI(PCC,TUCT,ROS,GMOS,MNURT)


C----------------------------------------------------------------------
      RNS=ROS*1.D3
      DBS=1.D0-2.D0*GS*GMOS*MSS/RNS

      KC=0.5D0*DLOG(DBS)-MNURT
      ROT= TUCT + KC

C----------------------------------------
C  Check matric nu
C----------------------------------------- 
c      CALL TOVI(PCC,ROT,ROS,GMOS,MAMA)
c      CMNUR=0.5D0*DLOG(1.D0-2.D0*GS*GMOS*MSS/RNS)
c      WRITE(*,*)MAMA,MNURT,CMNUR,PCC,ROT,ROS,GMOS,MAMA
C--------------------------------------------------
c Moment of Inertia related properties
c----------------------------------------------------
C   ROS2 is radius in km calculate using initial correct nu
C   GMOS2 is mass in solar mass calculate using initial correct nu
C   MAMA is matric component nu in R calculate using initial correct nu
C   OMEGA is rotation frequency and KAPPA is nedded to calculate moment
C   of inertia
C   MOMIN is I/MR^2 dimesionless of NS
C   MI is moment of inertia in 10^45 g cm^2 of NS

C   Defining O, i.e. OmegaK
       Z=GS*MSS*GMOS/(((ROS*1.0D3)**3))
       QA=(24.0D-2*2.0D0*GS*GMOS/ROS)
       QB=2.0D0*(24.0D-2*24.0D-2)*(GS*GMOS/ROS)**2
       Q=(1.0D0/(SQRT(1.0D0+QA-QB)))
       O=Q*(SQRT(Z))
       
       PMIN=1.0D-9
       CALL TOVMI(PCC,ROT,PMIN,ROS2,GMOS2,MNURT,OMGA,KP,O,Po,Mo,V,W,Z,X)
       RNS2=ROS2*1.D3

C------------------------------------------------------
       KAPA=KP
       MOMIN=KAPA/(GS*GMOS*MSS*RNS2*RNS2)
       MI=MOMIN*1.98892D33*1.D10*GMOS*ROS2*ROS2/1.0D45
       
       EJMI=O*KAPA
       DELTAM=Mo/MSS+((EJMI**2)/(MSS*RNS**3))
      
c-----------------------------------------------------------
c     CRUST Properties
c----------------------------------------------------------
C     PT pressure at core-crust transition
C     RNST is core radius
c     MOMINC I/MR^2 dimesionless of the core of NS
c     MIC is moment of Inertia    in 10^45 g cm^2 of NS core
c     MICR,MGCR,RPMIC,RPMGC are crust related properties
c--------------------------------------------------------------       
       PT=2.863D-1
       CALL TOVMI(PCC,ROT,PT,RT,GMC,MNUC,OMGAC,KC,O,PoC,MoC,VC,WC,ZC,XC)
       RNST=RT*1.D3
       
       KAPAC=KC
       MOMINC=KAPAC/(GS*GMC*MSS*RNST*RNST)
       MIC=MOMINC*1.98892D33*1.D10*GMC*RT*RT/1.0D45
C------------------------------------------------------------
       MICR=MI-MIC
       MGCR=GMOS2-GMC
       RPMIC=MICR/MI*1D2
       RPMGC= MGCR/GMOS2*1D2
       
       OS=O*3*1.0D8
       OMS=O*OMGA*3*1.0D8

C      total mass, polar radius, equator radius, aveage radius, eccentricity, and ellipticity        
       MTOT=GMOS2+DELTAM
       RTOTP=ROS2+Z
       RTOTE=ROS2+X
       RAVE=(RTOTP+RTOTE)/2.0D0
       ECCEN=SQRT(1.0D0-((RTOTP/RTOTE)**2))
       ELLIPT=(1.0D0/2.0D0)*(ECCEN**2)
    

        WRITE(*,*)IL,ROS2,GMOS2,Po,DELTAM,V,W,Z,X,OMS,OS
       WRITE(8,*)IL,PCC,ROS2,RT,GMOS2,MOMIN,MI,Po,Mo,V,W,DELTAM,O,OMGA
       WRITE(9,*)IL,GMOS2,MGCR,MICR,RPMIC,RPMGC,PoC,MoC,VC,WC
       WRITE(10,*)IL,ROS2,RTOTP,RTOTE,RAVE,MTOT,OMS,OS,ECCEN,ELLIPT

 10   CONTINUE

       
      STOP
      END

C----------------------------------------------------------
C TOV Inertia Moment 
C--------------------------------------------
           
      SUBROUTINE TOVMI(PCC,TUCT,PMIN,ROS,GMOS,MNURT,D,K,O,Po,Mo,V,W,Z,X)
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z) 
      INTEGER I, J, IM, NS, LI, N, IN    
      DIMENSION YA(10), EK(5,10), Y(10)
      HC  = 197.327D0


C     IM = NUMBER OF EQUATIONS (10)
C     IN =  NUMBER OF FIRST ORDER DIFFERENTIAL EQUATIONS (9)


      IM=10
      IN=IM-1

C---------------------------------------------------
C     Y(1)=Pressure, 
C     Y(2)=Star Mass,
C     Y(3)=Matric Nu
C     Y(4)=Omega, 
C     Y(5)=Kappa,
C     Y(6)=Po
C     Y(7)=Mo
C     Y(8)=W
C     Y(9)=V
C----------------------------------------------------
C     Y(10)=Energy Density,

         
      Y(1)=PCC
      Y(3)=TUCT

      Y(2)=0.1D-8     
      Y(4)=0.1D-8
      Y(5)=0.1D-9
      Y(6)=0.1D-9
      Y(7)=0.1D-9
      Y(8)=0.1D-9
      Y(9)=0.1D-9
      
c--------------------------------------------------------------------------
c Compute initial Energy density
c---------------------------------------------------------------------      

      
      Y(10)=FED(PCC)


     
C     PU=INTERVAL OF radius FOR PRINTING, NS=NUMBER OF STEP IN ONE PRINT
C     INTERVAL OF T,XL=MAXIMUM radius TO STOP CALCULATION

      PU=1.0D-1
      NS=32
      XL=30.0D3

      H=PU/NS
      XP=1.0D-3
      HH=H/(2.0D0)

C     LINE NUMBER INITIALIZATION

      LI=0
 

 28   LI=LI+1

C     XB=OLD radius, XP=NEW radius, XM=MIDPOINT radius

      DO N=1,NS
         XB=XP
         XP=XP+H
         XM=XB+HH


C    COMPUTE K1, L1, M1, N1, O1

         J=1
         DO I=1,IM
            YA(I)=Y(I)
         END DO
         XA=XB

   
         CALL FUNCX(EK,J,YA,XA,H,O,PTWO,XIP,XIE)

C    COMPUTE K2, L2, M2, N2, O2

         J=2
         
         DO I=1,IN
            YA(I)=Y(I)+EK(1,I)/(2.D0)
         END DO
            P0=YA(1)
c--------------------------------------------------------------------------
c Call Presure vs energy and baryon densities parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
        YA(10)=ED
c---------------------------------------------------------------------
                     
         XA=XM

         CALL FUNCX(EK,J,YA,XA,H,O,PTWO,XIP,XIE)

C    COMPUTE K3, L3, M3, N3, O3

         J=3
         DO I=1,IN
            YA(I)=Y(I)+EK(2,I)/(2.D0)
         END DO
     

          P0=YA(1)

c--------------------------------------------------------------------------
c Call Presure vs energy and baryon densities parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
        YA(10)=ED    
c---------------------------------------------------------------------- 
                               
         XA=XM


         CALL FUNCX(EK,J,YA,XA,H,O,PTWO,XIP,XIE)

C    COMPUTE K4, L4, M4, N4, O4

         J=4
         DO I=1,IM
            YA(I)=Y(I)+EK(3,I)
         END DO

          P0=YA(1)

c--------------------------------------------------------------------------
c Call Presure vs energy and baryon densities parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
        YA(10)=ED

c---------------------------------------------------------------------
                 

         XA=XP

         CALL FUNCX(EK,J,YA,XA,H,O,PTWO,XIP,XIE)

C    4-TH ORDER RUNGE KUTTA SCHEME

         DO I=1,IN
            Y(I)=Y(I)+(EK(1,I)+2.D0*EK(2,I)+2.D0*EK(3,I)+EK(4,I))/6.D0
         END DO
          P0=Y(1)

c--------------------------------------------------------------------------
c Call Presure vs energy and baryon density parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
        Y(10)=ED
c---------------------------------------------------------------------  
                   
          
       END DO



       PS=Y(1)
         
c       PT=2.863D-1 
c       PMIN=1.0D-9 
      
       IF (PS .GT. PMIN  ) GOTO 28
     
       
      ROS=(XP/1.D3)
      GMOS=Y(2)
      MNURT=Y(3)
     
      OMET=Y(4)
      KAPAT=Y(5)
      Po=Y(6)
      Mo=Y(7)
      W=Y(8)
      V=Y(9)
      ETA=1.D0/(OMET+2.D0*KAPAT/(XP*XP*XP))

      D=OMET
      K=ETA*KAPAT
      
C      X represents deformation at the pole 
      Z=XIP
      X=XIE
   
      RETURN
      END
      
C--------------------------------------------------------------
C  TOV initial condition for nu
C----------------------------------------------------------

      SUBROUTINE TOVI(PCC,TUCT,ROS,GMOS,MNURT)
      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z) 
      INTEGER I, J, IM, NS, LI, N, IN    
      DIMENSION YA(10), EK(5,10), Y(10)

      HC  = 197.327D0

C     IM = NUMBER OF EQUATIONS
C------------------------------------------
C     IN =  NUMBER OF FIRST ORDER DIFFERENTIAL EQUATIONS
C     Y(1)=Pressure, 
C     Y(2)=Star Mass,
C     Y(3)=Matric Nu,  
C----------------------------------------------
C     Y(4)=Energy Density, 


      IM=4
      IN=IM-1  

      Y(1)=PCC
      Y(2)=0.1D-8
      Y(3)=TUCT

      
c--------------------------------------------------------------------------
c Call Presure vs energy and baryon densities parametrization (EOS)
c---------------------------------------------------------------------      
    
      Y(4)=FED(PCC)
     
     
C     PU=INTERVAL OF radius FOR PRINTING, NS=NUMBER OF STEP IN ONE PRINT
C     INTERVAL OF T,XL=MAXIMUM radius TO STOP CALCULATION

      PU=1.0D-1
      NS=16
      XL=20.0D3

      H=PU/NS
      XP=1.0D-3
      HH=H/(2.0D0)

C     LINE NUMBER INITIALIZATION

      LI=0
 

 28   LI=LI+1

C     XB=OLD radius, XP=NEW radius, XM=MIDPOINT radius

      DO N=1,NS
         XB=XP
         XP=XP+H
         XM=XB+HH

C    COMPUTE K1, L1, M1

         J=1
         DO I=1,IM
            YA(I)=Y(I)
         END DO
         XA=XB
   
         CALL FUNCT(EK,J,YA,XA,H)

C    COMPUTE K2, L2, M2

         J=2
         
         DO I=1,IN
            YA(I)=Y(I)+EK(1,I)/(2.D0)
         END DO
         
        P0=YA(1)
c--------------------------------------------------------------------------
c Call Presure vs energy and baryon densities parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
        YA(4)=ED
c---------------------------------------------------------------------

        
         XA=XM

         CALL FUNCT(EK,J,YA,XA,H)

C    COMPUTE K3, L3, M3

         J=3
         DO I=1,IN
            YA(I)=Y(I)+EK(2,I)/(2.D0)
         END DO
     

          P0=YA(1)

c--------------------------------------------------------------------------
c Call Presure vs energy and baryon densities parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
        YA(4)=ED 
                 
c---------------------------------------------------------------------- 
           
         XA=XM

         CALL FUNCT(EK,J,YA,XA,H)

C    COMPUTE K4, L4, M4

         J=4
         DO I=1,IM
            YA(I)=Y(I)+EK(3,I)
         END DO

          P0=YA(1)

c--------------------------------------------------------------------------
c Call Presure vs energy and baryon densities parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
        YA(4)=ED

c---------------------------------------------------------------------
                

         XA=XP

         CALL FUNCT(EK,J,YA,XA,H)

C    4-TH ORDER RUNGE KUTTA SCHEME

         DO I=1,IN
            Y(I)=Y(I)+(EK(1,I)+2.D0*EK(2,I)+2.D0*EK(3,I)+EK(4,I))/6.D0
         END DO
          P0=Y(1)

c--------------------------------------------------------------------------
c Call Presure vs energy and baryon density parametrization (EOS)
c---------------------------------------------------------------------      

        ED=FED(P0)
        Y(4)=ED
c---------------------------------------------------------------------  
                  
          
       END DO

  

       PS=Y(1)
       PMIN=1.0D-9
c       PMIN=2.0D-5
     
      IF (PS .GT. PMIN  ) GOTO 28

      ROS=(XP/1.D3)
      GMOS=Y(2)
      MNURT=Y(3)
     
      
      RETURN
      END
      
       SUBROUTINE FUNCX(EK,J,YA,XA,H,O,PTWO,XIP,XIE)
C     *********************************************************
C     DEFINES TOV EQUATIONS AIP BL 
C     P-unisotropic used to calculate Moment of inertia
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K-Z) 
      INTEGER J   
      DIMENSION EK(5,10), YA(10)
      PI  = 3.14159265358979D0      
      GS=1.325D-12
      MSS=1.1155D15
C ----------------------------------
c    model parameter !! Mmax=2.08 Ms
  
      Con=0.950D0

C--------------------------------------       
      PRESS=YA(1)
      MASST=YA(2)
      MNU=YA(3)
      OMGT=YA(4)
      KAPT=YA(5)
      Po=YA(6)
      Mo=YA(7)
      W=YA(8)
      V=YA(9)
      EDEN=YA(10) 
      EXPMNU=EXP(MNU)
      EXPMMNU=EXP(-MNU)
C      j = EJ
      EJ=EXPMMNU*SQRT(1.D0-(2.D0*GS*MASST*MSS/XA))
C      dj/dr = DJ
      DJ=EXPMMNU*EXPMMNU/EJ
C      (dj/dr)^2 = DJKUADRAT
      DJKUADRAT=-8.D0*PI*GS*XA*EXPMMNU*EXPMMNU*(EDEN+PRESS)
C      (dnu/dr) = MNUAKSEN
      MNUAKSEN=GS*(MASST*MSS+4.D0*PI*(XA**3)*PRESS)/(XA*(XA-2.D0
     &         *GS*MSS*MASST))
C      d(omega bar)/dr = DOMEGABAR
      DOMEGABAR=O*(6.D0/(EJ*XA*XA*XA*XA))*KAPT
C      omega bar = OMEGABAR
      OMEGABAR=O*OMGT

      AICORR=1.D0+(1.D0-Con)*GS*MASST*MSS*H/(2.D0*XA)
     &       *(1.D0+4.D0*PI*XA*XA*XA*PRESS/(MSS*MASST))
     &       /(1.D0-2.D0*GS*MASST*MSS/XA)
     
C     Define things related to derivative of AICORR with rescpect to P

      PAKSEN=-Con*GS*EDEN*MASST*MSS/(XA*XA)
     &        *(1.D0+PRESS/EDEN)
     &        *(1.D0+4.D0*PI*XA*XA*XA*PRESS/(MSS*MASST))
     &        /(1.D0-2.D0*GS*MASST*MSS/XA)
     
      SATU=-1.D0*(MSS*MASST+4.D0*PI*PRESS*(XA**3)
     &     *(EDEN+PRESS)*(1.D0/PAKSEN))/(2*Con*(1.0D0-2.D0*GS*MSS
     &     *MASST/XA)*XA*XA)
     
      DUA=(PRESS+EDEN)*(4.D0*PI*(XA**3)+12.D0*PRESS*PI*(1/PAKSEN)
     &    +(1/PAKSEN)*4.D0*PI*XA*XA*EDEN)/(2*Con*(1.0D0-2.D0*GS*MSS
     &    *MASST/XA)*XA)
     
      TIGA=-1.D0*(MSS*MASST+4.D0*PRESS*PI*XA**3)*(PRESS+EDEN)
     &     *((2.D0*GS*MSS*MASST*(1/PAKSEN)/(XA**2))
     &     -(2*GS*4.D0*PI*XA*XA*EDEN*(1/PAKSEN)/XA))
     &     /(2*Con*((1.0D0-2.D0*GS*MSS*MASST/XA)**2)*XA)
    
      EMPAT=(MASST*MSS+4.D0*PRESS*PI*XA**3)*(1.D0+DFED(PRESS))
     &      /(2*Con*(1.0D0-2.D0*GS*MSS*MASST/XA)*XA)
     
C      DAICORR is derivative of sigma with respect to pressure
     
      DAICORR=(GS*(1-Con))*(SATU+DUA+TIGA+EMPAT)

C      solving TOV eq (dp/dr)

      EK(J,1)=-Con*GS*EDEN*MASST*MSS*H/(XA*XA)
     &        *(1.D0+PRESS/EDEN)
     &       *(1.D0+4.D0*PI*XA*XA*XA*PRESS/(MSS*MASST))
     &       /(1.D0-2.D0*GS*MASST*MSS/XA)

C      solving dm/dr
 
      EK(J,2)=4.D0*PI*XA*XA*EDEN*H/MSS

C      solving d(nu)/dr
 
      EK(J,3)=GS*MASST*MSS*H/(XA*XA)
     &       *(1.D0+4.D0*PI*XA*XA*XA*PRESS/(MSS*MASST))
     &       /(1.D0-2.D0*GS*MASST*MSS/XA)

C      solving d(kappatilde)/dr 

      EK(J,4)=6.D0*EXPMNU*H/DSQRT(1.D0-2.D0*GS*MASST*MSS/XA)
     &     *KAPT/(XA*XA*XA*XA)

C      solving d(omegatilde)/dr

      EK(J,5)=GS*8.D0*PI*(XA*XA*XA*XA)*EXPMMNU*EDEN*H/3.D0
     &        *(1.D0+PRESS/EDEN)*AICORR*OMGT
     &        /DSQRT(1.D0-2.D0*GS*MASST*MSS/XA)

C      solving dp_0/dr; one, two, three, four are terms in the eq
     
      ONE=GS*Mo*(1.D0+8.D0*GS*PI*XA*XA*PRESS)/((XA-2.D0*GS*MSS
     &    *MASST)*(XA-2.D0*GS*MSS*MASST))
      TWO=-1.D0*GS*4.D0*PI*XA*XA*((EDEN+PRESS)/(XA-2.D0*GS*MSS*MASST))
     &    *AICORR*Po
      THREE=(1.D0/12.D0)*((XA**4)*(EJ**2)/(XA-2.D0*GS*MSS*MASST))
     &      *(DOMEGABAR**2)
     
      FOURONE=3.D0*XA*XA*(EJ**2)
     &     *(OMEGABAR**2)
      FOURTWO=2.D0*(XA**3)*(OMEGABAR**2)*EJ*DJ
      FOURTHREE=2.D0*(XA**3)
     &     *(EJ**2)*(OMEGABAR*DOMEGABAR)*(XA-2.D0*GS*MSS*MASST)
      FOURFOUR=((XA**3)*(EJ**2)*(OMEGABAR**2))
     &     *(1.D0-2.D0*GS*4.D0*PI*XA*XA
     &     *EDEN)
      FOUR=(1.D0/3.D0)
     &     *((FOURONE+FOURTWO+FOURTHREE)/(XA-2.D0*GS*MSS*MASST))
     &     -(FOURFOUR)/((XA-2.D0*GS*MSS*MASST)**2)

      EK(J,6)=H*(ONE+TWO+THREE+FOUR)
 
C      solving dm_0/dr; five, six, seven are terms in the eq
      
      FIVE=4.D0*GS*PI*DFED(PRESS)*(EDEN+PRESS)*Po*AICORR
     &     /(1.D0-DAICORR)
      SIX=(1.D0/12.D0)*(EJ**2)*(XA**4)*(DOMEGABAR**2)/GS
      SEVEN=(-1.D0/3.D0)*(XA**3)*(DJKUADRAT)*(OMEGABAR**2)*AICORR/GS
      EK(J,7)=H*(FIVE+SIX+SEVEN)

C      solving dh_2/dr; eight, nine, ten, eleven are terms in the eq      
c      W represents h_2
C      V represents v_2
      
      EIGHT=W*(-2.D0*MNUAKSEN+(XA/((XA-2.D0*GS*MASST*MSS)*(2.D0
     &      *2.D0*MNUAKSEN)))*(8.D0*GS*PI*(EDEN+PRESS)*AICORR
     &      /(1.D0-DAICORR)-4.D0*GS*MSS*MASST/(XA**3)))
      MNINE=-4.D0*V/(XA*XA*(1.D0-2.D0*GS*MSS*MASST/XA)*2.D0*MNUAKSEN)
      TENONE=(1.D0/2.D0)*2.D0*MNUAKSEN*XA*(XA**3)*(EJ**2)*(DOMEGABAR**2)
      TENTWO=(1.D0/((XA-2*GS
     &       *MASST*MSS)*2.D0*MNUAKSEN))*(XA**3)*(EJ**2)*(DOMEGABAR**2)
      TEN=(1.D0/6.D0)*(TENONE-TENTWO)
      ELEVEN=(-1.D0/3.D0)*((1.D0-(6.D0*AICORR)+6.D0)*(1.D0
     &       /2.D0)*2.D0*MNUAKSEN*XA+(1.D0/((XA-2.D0*GS*MSS*MASST)
     &       *2.D0*MNUAKSEN))*AICORR*(1.D0
     &       -DAICORR))*(XA**2)*DJKUADRAT*(OMEGABAR**2)
      EK(J,8)=H*(EIGHT+MNINE+TEN+ELEVEN)
      

C      solving dv_2/dr; twelve, thirteen, fourteen are terms in the eq
c      W represents h_2
C      V represents v_2
      
      TWELVE=-2.D0*MNUAKSEN*W
      THIRTEEN=(2.D0*AICORR)-2.D0
      FOURTEEN=(1.D0/6.D0)*(EJ**2)*(XA**4)*(DOMEGABAR**2)
      EK(J,9)=H*(TWELVE+((1.D0/XA)+(1.D0/2.D0)*2.D0
     &        *MNUAKSEN)*(-1.D0*((1.D0/3.D0)-THIRTEEN)
     &        *(XA**3)*DJKUADRAT*(OMEGABAR**2)+FOURTEEN))
     
c      WRITE(*,*) EK(J,9),EK(J,8),EK(j,7),EK(J,6)
C      declaring p_2, fifteen & sixteen are terms in p_2     
     
      FIFTEEN=-1.D0*W
      SIXTEEN=-(1.D0/3.D0)*(XA**2)*(EXPMMNU**2)*(OMEGABAR**2)
      
C      declaring P2
      
      PTWO=-FIFTEEN-SIXTEEN
      
C     KOMPXI is a mathematical thing that is used to find XI 
      
      KOMPXI=(-Con*GS*EDEN*MASST*MSS*H/(XA*XA)
     &        *(1.D0+PRESS/EDEN)
     &       *(1.D0+4.D0*PI*XA*XA*XA*PRESS/(MSS*MASST))
     &       /(1.D0-2.D0*GS*MASST*MSS/XA)-DAICORR)/((EDEN+PRESS)
     &       *AICORR)
     
C     declaring XIP for the deformation of the pole

      XIP=(-Po-PTWO)/(KOMPXI*1.D3) 

C     declaring XIE for the deformation of the equator
     
      XIE=(-Po+(1.D0/2.D0)*PTWO)/(KOMPXI*1.D3)
      
      RETURN
      END


      SUBROUTINE FUNCT(EK,J,YA,XA,H)
C     *********************************************************
C     DEFINES TOV EQUATIONS AIP BL 
C     P-anisotropic used to calculate initial condition for nu
C     *********************************************************
      IMPLICIT DOUBLE PRECISION (A-H,K-Z) 
      INTEGER J   
      DIMENSION EK(5,10), YA(10)
      PI  = 3.14159265358979D0     
      GS=1.325D-12
      MSS=1.1155D15
C ----------------------------------
c    model parameter !! Mmax=2.08 Ms
  
      Con=0.950D0 
c--------------------------------------------      
      PRESS=YA(1)
      MASST=YA(2)
      MNU=YA(3)
      EDEN=YA(4) 
    
      EK(J,1)=-Con*GS*EDEN*MASST*MSS*H/(XA*XA)
     &        *(1.D0+PRESS/EDEN)
     &       *(1.D0+4.D0*PI*XA*XA*XA*PRESS/(MSS*MASST))
     &       /(1.D0-2.D0*GS*MASST*MSS/XA)
     

      EK(J,2)=4.D0*PI*XA*XA*EDEN*H/MSS
 
      EK(J,3)=GS*MASST*MSS*H/(XA*XA)
     &       *(1.D0+4.D0*PI*XA*XA*XA*PRESS/(MSS*MASST))
     &       /(1.D0-2.D0*GS*MASST*MSS/XA)
   
      RETURN
      END

 
C---------------------------------------------------------------------
C  Energy density as a function of pressure
C------------------------------------------------------------------
      FUNCTION FED(P0)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)

        IF ( P0 .GT. 50.D0 ) THEN


        FED=752.1779439721782D0 - 40.108620203769D0*P0 + 
     -  1.4752000864096837D0*P0**2 - 
     -  0.027036299898738337D0*P0**3 + 
     -  0.0003125837415657684D0*P0**4 - 
     -  2.441653734444954D-6*P0**5 + 
     -  1.3308776414303D-8*P0**6 - 
     -  5.1246996017212506D-11*P0**7 + 
     -  1.3886907993633848D-13*P0**8 - 
     -  2.590695651603915D-16*P0**9 + 
     -  3.166994006343349D-19*P0**10 - 
     -  2.283247911494129D-22*P0**11 + 
     -  7.357464623695704D-26*P0**12



  

        ELSE IF ( P0 .GT.  2.863D-1 .AND. P0 .LE. 50.D0 ) THEN


        FED=65.94376092512762D0 + 57.952630512664356D0*P0 - 
     -  15.437797576854498D0*P0**2 + 
     -  2.7511762122977346D0*P0**3 - 
     -  0.299643408613215D0*P0**4 + 
     -  0.02089083577091044D0*P0**5 - 
     -  0.0009679766327308126D0*P0**6 + 
     -  0.000030437894581131215D0*P0**7 - 
     -  6.521271499180282D-7*P0**8 + 
     -  9.372308503275227D-9*P0**9 - 
     -  8.643173476154485D-11*P0**10 + 
     -  4.621195967220976D-13*P0**11 - 
     -  1.0889239514378402D-15*P0**12
     
        ELSE IF (P0 .GT. 4.99313436D-4 .AND. P0 .LE. 2.863D-1) THEN


        FED=0.05015663787134234D0 + 836.2363942486941D0*P0 - 
     -  9315.146969977652D0*P0**2 + 79689.1930322726D0*P0**3 - 
     -  412197.6475732246D0*P0**4 + 1.116366190255507D6*P0**5 - 
     -  1.1988188657021397D6*P0**6


        ELSE
        
        FED=0.00020663104786863406D0 + 985.7550962048012D0*P0 - 
     -  6.452649548410687D6*P0**2 + 4.045493650683396D10*P0**3 - 
     -  1.2422017897554384D14*P0**4 + 1.765493095354931D17*P0**5 - 
     -  9.15294170121406D19*P0**6

        
        END IF

   
      RETURN
      END
      
C-------------------------------------------------------------------------
C  Derivative of energy density with respect to pressure
C-------------------------------------------------------------------------


      FUNCTION DFED(P0)
      IMPLICIT DOUBLE PRECISION (A-H,K-Z)

        IF ( P0 .GT. 50.D0 ) THEN


        DFED= - 40.108620203769D0 + 
     -  2.0D0*1.4752000864096837D0*P0 - 
     -  3.0D0*0.027036299898738337D0*P0**2 + 
     -  4.0D0*0.0003125837415657684D0*P0**3 - 
     -  5.0D0*2.441653734444954D-6*P0**4 + 
     -  6.0D0*1.3308776414303D-8*P0**5 - 
     -  7.0D0*5.1246996017212506D-11*P0**6 + 
     -  8.0D0*1.3886907993633848D-13*P0**7 - 
     -  9.0D0*2.590695651603915D-16*P0**8 + 
     -  10.0D0*3.166994006343349D-19*P0**9 - 
     -  11.0D0*2.283247911494129D-22*P0**10 + 
     -  12.0D0*7.357464623695704D-26*P0**11



  

        ELSE IF ( P0 .GT.  2.863D-1 .AND. P0 .LE. 50.D0 ) THEN


        DFED= 57.952630512664356D0 - 
     -  2.0D0*15.437797576854498D0*P0 + 
     -  3.0D0*2.7511762122977346D0*P0**2 - 
     -  4.0D0*0.299643408613215D0*P0**3 + 
     -  5.0D0*0.02089083577091044D0*P0**4 - 
     -  6.0D0*0.0009679766327308126D0*P0**5 + 
     -  7.0D0*0.000030437894581131215D0*P0**6 - 
     -  8.0D0*6.521271499180282D-7*P0**7 + 
     -  9.0D0*9.372308503275227D-9*P0**8 - 
     -  10.0D0*8.643173476154485D-11*P0**9 + 
     -  11.0D0*4.621195967220976D-13*P0**10 - 
     -  12.0D0*1.0889239514378402D-15*P0**11
     
        ELSE IF (P0 .GT. 4.99313436D-4 .AND. P0 .LE. 2.863D-1) THEN


        DFED= 836.2363942486941D0 - 
     -  2.0D0*9315.146969977652D0*P0 + 
     -  3.0D0*79689.1930322726D0*P0**2 - 
     -  4.0D0*412197.6475732246D0*P0**3 + 
     -  5.0D0*1.116366190255507D6*P0**4 - 
     -  6.0D0*1.1988188657021397D6*P0**5


        ELSE
        
        DFED=985.7550962048012D0 - 
     -  2.0D0*6.452649548410687D6*P0 +
     -  3.0D0*4.045493650683396D10*P0**2 - 
     -  4.0D0*1.2422017897554384D14*P0**3 +
     -  5.0D0*1.765493095354931D17*P0**4 - 
     -  6.0D0*9.15294170121406D19*P0**5

        
        END IF

   
      RETURN
      END

 
c---------------------------------------------------------------------
C end of the code
