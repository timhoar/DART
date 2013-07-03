      SUBROUTINE IONDEN2(QTOT,QI,Z,RZUR,PHI,TMO,RLT,RLTM,RLGM,DIP)     
C
C	This is a modified version of the CHIU(1975) empirical ionospheric
C	model.  Changes in the FORTRAN code are contained between the
C  ****************************************************************** ' s
C	below.  Included are a midlatitude correction to the height of the
C	F2-layer, and equatorial-region (between +20 and -20 degrees
C	magnetic latitude) corrections to foF2 and hmF2, and modified
C	scale heights to be used in the Chapman-layer specification of 
C	the bottomside and topside F-region.  Added subroutines are
C	CORRECT and HERMW.  CORRECT computes the equatorial corrections
C	and HERMW generates Hermite polynomials used in CORRECT.
C
C	
      DIMENSION QI(3),AMP(3),ALF(3),VAR(3),ZMAX(3),HMAX(3)  
      COMMON/HEIGHT/ZMF2,FOF2 
      COMMON /MODEi/NZONAL,MOIS,NSS
      COMMON /MODEr/PERIOD,FREQ,FLUX
      DATA AMP/1.36,2.44,.66/ 
      DATA ALF/.5,.5,1./      
      DATA ZMAX(1),ZMAX(2)/110.,180./   
      DATA HMAX(1),HMAX(2)/10.,34./     
      
      PI=ACOS(-1.)  
      SDEC=.39795*SIN(PI*(TMO-3.167)/6.)
      DEC=ASIN(SDEC)
      DELP=ABS(ABS(RLT)-PI/2.)
      IF(DELP.GT.1.E-03) GOTO 5         
      SEASN=RLT*DEC 
      IF(SEASN.LT.0.) PHI=0.  
      IF(SEASN.GE.0.) PHI=PI
C     write *,"RZUR= ", RZUR
 5    R=RZUR/100.   
C     write *,"R= ", R
      CLTM=COS(RLTM)
      SLTM=SIN(RLTM)
      XLAM=1.+.5*ALOG(1.+30.*R)         
      VAR(1)=TVEF1(1.15,0.,.4,2.,RLT,R,PHI,DEC)   
      VAR(2)=TVEF1(1.24,.25,.25,XLAM,RLT,R,PHI,DEC)         
      VAR(3)=TVARF2(RLTM,RLGM,DIP,R,PHI,TMO,DEC,PI,CLTM,SLTM)         
      IF(Z.EQ.0.) GOTO 100    
      ZALF=-4.5*ABS(RLTM)-PI  
      ZBA=240.+10.*CLTM*COS(PI*(TMO/3.-1.5))      
      ZBAR=ZBA+R*(75.+83.*CLTM*ZETA(DEC,RLTM))    
      ZMAX(3)=ZBAR+30.*COS(PHI+ZALF)    
C***************************************************************      
C         
C     MIDLATITUDE CORRECTION  
C         
C         
      THM = RLTM*180./PI 
      DEGRAD=PI/180.      
      CLTM40N=COS(40.*DEGRAD) 
      CLTM40S=COS(-40.*DEGRAD)
      SLTM40N=SIN(40.*DEGRAD) 
      SLTM40S=SIN(-40.*DEGRAD)
      CLTM20N=COS(20.*DEGRAD) 
      CLTM20S=COS(-20.*DEGRAD)
      SLTM20N=SIN(20.*DEGRAD) 
      SLTM20S=SIN(-20.*DEGRAD)
      RLTM40N=40.*PI/180.
      RLTM40S=-40.*PI/180.         
      RLTM20N=20.*PI/180.
      RLTM20S=-20.*PI/180.         
      IF(THM.LT.40..AND.THM.GT.20.)GO TO 112      
      GO TO 113     
112   VAR40N=TVARF2(RLTM40N,RLGM,DIP,R,PHI,TMO,DEC,PI,      
     1CLTM40N,SLTM40N)        
      VAR20N=TVARF2(RLTM20N,RLGM,DIP,R,PHI,TMO,DEC,PI,      
     1CLTM20N,SLTM20N)        
      VAR3N=VAR40N-(VAR40N-VAR20N)*(40.-THM)/20.  
113   CONTINUE      
      IF(THM.GT.-40..AND.THM.LT.-20.) GO TO 114   
      GO TO 115     
114   VAR40S=TVARF2(RLTM40S,RLGM,DIP,R,PHI,TMO,DEC,PI,      
     1CLTM40S,SLTM40S)        
      VAR20S=TVARF2(RLTM20S,RLGM,DIP,R,PHI,TMO,DEC,PI,      
     1CLTM20S,SLTM20S)        
      VAR3S=VAR40S-(VAR40S-VAR20S)*(-40.-THM)/(-20.)        
115   CONTINUE      
      ZALF = -2.*PI + .75*PI*EXP(-THM*THM/400.)   
      ZBAR = ZBA + 100.*SQRT(R) - 10. + 83.*R*CLTM*ZETA(DEC,RLTM)     
      ZMAX(3) = ZBAR + (56.-12.*SQRT(R))*COS(PHI+ZALF)      
     1 + (9.+12.*SQRT(R))*COS(2.*PHI-PI/4.)       
     2 + 12.*COS(3.*PHI-.75*PI)         
C         
C     LOW-LATITUDE/EQUATORIAL CORRECTION
C
C	CHMF2 = correction to hmF2 (km)
C	CFOF2 = correction to foF2 (MHz)
C	SHUP  = topside plasma scale height (km)
C	SHLO  = bottomside plasma scale height (km)
C         
      FOF2=2.85*SQRT(AMP(3)*VAR(3))     
      IF(ABS(THM).GT.40.) GO TO 111     
      NSOL=(1.9*(R+.62))      
      NSEAS=1       
      IF(TMO.LE.8.5.AND.TMO.GE.4.5) NSEAS=2       
      IF(TMO.LT.2.5) NSEAS=3  
      IF(TMO.GT.10.5) NSEAS=3 
      NS=NSEAS+3*(NSOL-1)     
      IF(THM.GT.20.) GO TO 121
      GO TO 122     
121   CALL CORRECT(20.,PHI,CHMF2,CFOF2,NS,SHUP,SHLO)        
      FOF220N=2.85*SQRT(AMP(3)*VAR20N)+CFOF2      
      FOF240N=2.85*SQRT(AMP(3)*VAR40N)  
      DFDLN=(FOF240N-FOF220N)/20.       
122   CONTINUE      
      IF(THM.LT.-20.) GO TO 123         
      GO TO 124     
123   CALL CORRECT(-20.,PHI,CHMF2,CFOF2,NS,SHUP,SHLO)       
      FOF220S=2.85*SQRT(AMP(3)*VAR20S)+CFOF2      
      FOF240S=2.85*SQRT(AMP(3)*VAR40S)  
      DFDLS=(FOF240S-FOF220S)/(-20.)    
124   CONTINUE      
      CALL CORRECT(THM,PHI,CHMF2,CFOF2,NS,SHUP,SHLO)        
      ZMAX(3)=ZMAX(3)+CHMF2   
      FOF2=FOF2+CFOF2         
      IF(THM.LT.40..AND.THM.GT.20.) FOF2=FOF240N-(40.-THM)*DFDLN      
      IF(THM.GT.-40..AND.THM.LT.-20.) FOF2=FOF240S-(-40.-THM)*DFDLS   
120   CONTINUE      
      VAR(3)=((FOF2/2.85)**2.)/AMP(3)   
111   CONTINUE      
      ZMF2=ZMAX(3)  
C         
C***************************************************************      
      HMAX(3)=.2*ZMAX(3)+40.  
      IF(Z.LT.ZMAX(3)) HMAX(3)=.2*Z+40. 
      IF(ABS(THM).LE.40.)HMAX(3)=SHUP   
      IF(ABS(THM).LE.40..AND.Z.LT.ZMAX(3))HMAX(3)=SHLO      
      SUM=0.        
      DO 10 I=1,3   
      ZP=(Z-ZMAX(I))/HMAX(I)  
      CZ=EXP(ALF(I)*(1.-ZP-EXP(-ZP)))   
      QI(I)=AMP(I)*VAR(I)*CZ  
      SUM=SUM+QI(I) 
 10   CONTINUE      
      QTOT=SUM      
      RETURN        
 100  CONTINUE      
      DO 110 I=1,3  
      QI(I)=AMP(I)*VAR(I)     
 110  CONTINUE      
      QTOT=0.       
      RETURN        
      END 
      SUBROUTINE CORRECT(THM,PHI,CHMF2,CFOF2,NS,SHUP,SHLO)  
      DIMENSION HH(9,9,5),HN(9,9,5),A(9),B(9)     
      DIMENSION HU(9,19,4),HL(9,19,4),AA(13),BB(13)         
C
C	Hermite polynomial coefficients (k,i,j) (left to right, k = 0
C	to k = 4) for each local-time Fourier component (top to bottom
C	a0,a1,b1,.....a4,b4 corresponding to an*cos(nt)+bn*sin(nt) where
C	t = local time in radians and n = 0 to n = 4  (i = 1,9).
C
C	HH corresponds to hmF2 (km)
C	HN corresponds to foF2 (mHz)
C
C	NS = 1,2,3 implies equinox,summer,winter   (solar minimum)
C	NS = 4,5,6 implies equinox,summer,winter   (average solar activity)
C	NS = 7,8,9 implies equinox,summer,winter   (solar maximum)
C
C         
      DATA ((HH(1,I,J),J=1,5),I=1,9)/   
     1  25.24      ,38.24      ,-.06      ,.42      ,-.63,      
     2  -53.85     ,62.57      ,.24       ,20.02    ,2.04,    
     3  -58.61     ,148.85     ,.05       ,85.33    ,6.89,   
     4  2.44       ,17.81      ,-1.29     ,14.98    ,.73,     
     5  12.29      ,-56.44     ,-.10      ,-28.13   ,-2.35, 
     6  21.20      ,-52.28     ,.23       ,-20.48   ,-1.96,  
     7  -11.38     ,5.25       ,-.22      ,-4.53    ,.28,     
     8  -9.96      ,43.64      ,1.30      ,24.20    ,1.78,    
     9  11.39      ,-37.72     ,.46       ,-17.09   ,-1.04/  
C         
      DATA ((HH(2,I,J),J=1,5),I=1,9)/   
     1  46.74      ,-40.44     ,30.37     ,-40.18   ,-3.98,
     2  -87.05     ,149.87     ,-14.09    ,64.08    ,5.24,
     3  -58.08     ,180.13     ,-4.45     ,95.71    ,7.16, 
     4   42.23     ,-85.03     ,11.33     ,-46.37   ,-4.07,         
     5  49.88      ,-138.16    ,-1.01     ,-81.87   ,-7.23,         
     6  -9.27      ,33.72      ,-.02      ,27.33    ,2.41,    
     7  -24.09     ,40.18      ,9.81      ,12.36    ,1.13,   
     8  25.50      ,-48.42     ,-7.13     ,-23.90   ,-1.89,
     9  -5.94      ,-.58       ,3.05      ,1.67     ,-.10/      
C         
      DATA ((HH(3,I,J),J=1,5),I=1,9)/   
     1  7.58       ,73.0       ,-37.24    ,24.03    ,2.39,    
     2  -55.02     ,111.50     ,12.55     ,25.05    ,1.30, 
     3  -31.31     ,109.44     ,3.80      ,55.87    ,3.17,  
     4  -.20       ,6.33       ,-10.01    ,8.88     ,.49,      
     5  12.73      ,-61.40     ,-2.39     ,-34.44   ,-2.56,
     6  18.06      ,-35.37     ,1.32      ,-11.63   ,-1.37, 
     7  -22.68     ,39.39      ,-8.34     ,16.69    ,1.66,  
     8  -.99       ,11.00      ,5.98      ,12.66    ,1.12,     
     9  19.64      ,-51.94     ,-.81      ,-33.27   ,-3.22/ 
C         
      DATA ((HN(1,I,J),J=1,5),I=1,9)/   
     1  -2.02      ,2.21       ,-.08      ,-.59     ,-.19,      
     2  .84        ,-.11       ,.11       ,-1.20    ,-.13,        
     3  1.99       ,-3.81      ,-.04      ,-1.49    ,.01,      
     4  1.66       ,-2.16      ,.05       ,-1.77    ,-.20,      
     5  -.29       ,-.18       ,.00       ,.35      ,.01,
     6  .32        ,.41        ,-.02      ,.07      ,.01, 
     7  -.07       ,.44        ,.03       ,.27      ,.05, 
     8  .56        ,-1.77      ,-.02      ,-.99     ,-.09,       
     9  .02        ,.89        ,-.01      ,.53      ,.02/ 
C         
      DATA ((HN(2,I,J),J=1,5),I=1,9)/   
     1  -1.10      ,1.06       ,-.25      ,-.18     ,-.05,      
     2  -.90       ,1.58       ,.00       ,.60      ,.09,
     3  -.08       ,1.10       ,.24       ,.91      ,.13,
     4  .44        ,-.67       ,-.12      ,-.35     ,-.05,        
     5  -.82       ,1.29       ,.39       ,1.11     ,.07,         
     6  -.28       ,1.49       ,.01       ,.82      ,.08,
     7  .13        ,.12        ,-.18      ,.05      ,.02, 
     8  .26        ,-.90       ,.08       ,-.52     ,-.05,         
     9  -.23       ,1.20       ,.01       ,.61      ,.04/
C         
      DATA ((HN(3,I,J),J=1,5),I=1,9)/   
     1  -.98       ,-2.71      ,.46       ,-2.50    ,-.21,      
     2  .22        ,2.08       ,.52       ,.54      ,.02, 
     3  -.93       ,5.28       ,-.36      ,3.07     ,.29,        
     4  .18        ,.39        ,.31       ,-.02     ,-.02,
     5  .04        ,-.81       ,-.53      ,-.18     ,.00,         
     6  -.22       ,1.35       ,-.05      ,.70      ,.07,         
     7  -.32       ,1.23       ,.14       ,.76      ,.08,
     8  -.16       ,-.02       ,-.08      ,.01      ,.00,         
     9  -.12       ,.83        ,.02       ,.46      ,.03/ 
C         
      DATA ((HH(4,I,J),J=1,5),I=1,9)/   
     1  48.74      ,-27.19     ,-.08      ,-49.72   ,-4.64, 
     2  -63.36     ,97.79      ,.52       ,46.53    ,4.39,    
     3  -50.74     ,82.53      ,.33       ,59.56    ,4.84,    
     4  25.90      ,-55.34     ,-1.93     ,-10.15   ,-.34, 
     5  16.94      ,-96.03     ,.31       ,-45.75   ,-4.31,  
     6  -3.71      ,6.12       ,-.02      ,23.5     ,1.73,      
     7  -31.76     ,59.39      ,-.39      ,25.90    ,2.58,   
     8  24.88      ,-51.50     ,-.21      ,-31.16   ,-3.03, 
     9  29.27      ,-69.92     ,.37       ,-45.78   ,-4.36/  
C         
      DATA ((HH(5,I,J),J=1,5),I=1,9)/   
     1  22.46      ,-1.48      ,29.53     ,-22.56   ,-2.39, 
     2  -46.25     ,53.46      ,-14.14    ,21.14    ,1.65, 
     3  -38.36     ,100.46     ,-1.94     ,55.95    ,4.44, 
     4  1.62       ,-3.00      ,9.22      ,6.12     ,.33,       
     5  10.02      ,-59.75     ,-.90      ,-37.14   ,-3.77, 
     6  10.53      ,-13.32     ,1.12      ,3.50     ,.48,     
     7  -38.60     ,76.58      ,9.44      ,31.55    ,2.78,   
     8  7.05       ,-4.88      ,-6.99     ,-.19     ,.09,      
     9  2.13       ,-9.96      ,-.89      ,-7.50    ,-.67/     
C         
      DATA ((HH(6,I,J),J=1,5),I=1,9)/   
     1  34.22      ,-1.26      ,-35.17    ,-26.44   ,-1.94,
     2  -53.98     ,107.60     ,14.28     ,29.74    ,2.63, 
     3  -32.66     ,64.96      ,2.35      ,44.88    ,2.52,   
     4  -4.59      ,11.62      ,-12.35    ,13.73    ,1.11,  
     5  2.63       ,-72.47     ,.97       ,-27.90   ,-2.54,   
     6  -4.08      ,17.57      ,-2.29     ,20.39    ,1.87,   
     7  -11.82     ,2.19       ,-9.66     ,-4.44    ,-.67,   
     8  5.59       ,-3.17      ,6.63      ,1.06     ,.03,       
     9  -.06       ,-.82       ,.34       ,-.47     ,-.06/        
C         
      DATA ((HN(4,I,J),J=1,5),I=1,9)/   
     1  -3.02      ,6.69       ,-.18      ,1.98     ,-.03,      
     2  .61        ,3.42       ,.19       ,.41      ,-.12,
     3  3.86       ,-8.29      ,-.02      ,-4.87    ,-.25,     
     4  2.59       ,-3.48      ,.04       ,-2.98    ,-.34,      
     5  1.03       ,-3.48      ,-.01      ,-1.63    ,-.13,     
     6  .77        ,-.35       ,-.04      ,-.65     ,-.06,        
     7  .46        ,-.29       ,.03       ,-.11     ,.02,
     8  .80        ,-2.76      ,-.03      ,-1.61    ,-.13,      
     9  -.26       ,1.57       ,-.01      ,1.10     ,.08/        
C         
      DATA ((HN(5,I,J),J=1,5),I=1,9)/   
     1  -1.73      ,2.42       ,-.60      ,.27      ,-.04,       
     2  -.24       ,1.81       ,-.08      ,.36      ,.05,         
     3  1.88       ,-3.03      ,.55       ,-1.46    ,-.02,      
     4  1.20       ,-2.22      ,-.34      ,-1.50    ,-.15,     
     5  -1.11      ,1.50       ,.56       ,1.43     ,.12,        
     6  .19        ,.70        ,.07       ,.30      ,.04,  
     7  .34        ,-.14       ,-.21      ,.04      ,.02,
     8  .10        ,-.71       ,.08       ,-.37     ,-.04,         
     9  -.49       ,1.90       ,-.03      ,1.08     ,.08/        
C         
      DATA ((HN(6,I,J),J=1,5),I=1,9)/   
     1  -1.44      ,-1.02      ,.75       ,-2.06    ,-.26,     
     2  .45        ,3.67       ,.85       ,1.22     ,.03,
     3  .92        ,1.28       ,-.69      ,.80      ,.16,
     4  1.05       ,-.84       ,.71       ,-1.05    ,-.15,       
     5  -.09       ,-.77       ,-.82      ,-.27     ,.00,        
     6  .35        ,.55        ,-.01      ,.02      ,.00, 
     7  .27        ,.18        ,.12       ,.17      ,.05,  
     8  .16        ,-.95       ,-.02      ,-.58     ,-.05,        
     9  -.31       ,1.52       ,.07       ,.86      ,.06/
C         
C         
      DATA ((HH(7,I,J),J=1,5),I=1,9)/   
     1  46.72      ,3.31       ,1.13      ,-46.51   ,-4.20,   
     2  -29.51     ,30.58      ,.68       ,15.87    ,2.10,    
     3  -57.42     ,40.63      ,-.04      ,56.22    ,5.00,   
     4  13.01      ,-31.89     ,-1.16     ,6.37     ,1.47,   
     5  -.39       ,-80.12     ,.06       ,-31.92   ,-2.93,   
     6  7.55       ,-48.25     ,.15       ,4.38     ,.38,       
     7  -20.53     ,39.17      ,-1.52     ,9.86     ,.89,    
     8  -36.62     ,105.62     ,.61       ,58.5     ,4.87,    
     9  6.07       ,5.03       ,.21       ,-.46     ,-.92/        
C         
      DATA ((HH(8,I,J),J=1,5),I=1,9)/   
     1  46.12      ,-46.14     ,33.10     ,-55.52   ,-5.34,
     2  -52.76     ,70.52      ,-15.89    ,36.88    ,3.42, 
     3  -47.51     ,93.67      ,-4.71     ,61.90    ,4.67,  
     4  34.38      ,-96.27     ,12.93     ,-34.42   ,-2.67,
     5  -1.04      ,-48.90     ,-4.88     ,-27.94   ,-2.45,
     6  1.16       ,6.32       ,2.55      ,17.84    ,1.72,      
     7  -42.83     ,84.45      ,4.78      ,34.07    ,2.82,   
     8  19.57      ,-31.43     ,-7.02     ,-16.33   ,-1.44,
     9  -3.27      ,13.83      ,-2.66     ,-.23     ,-.34/    
C         
      DATA ((HH(9,I,J),J=1,5),I=1,9)/   
     1  27.48      ,37.94      ,-35.88    ,-19.02   ,-1.18,
     2  -57.31     ,148.25     ,19.27     ,51.14    ,5.15, 
     3  -17.04     ,-12.06     ,2.62      ,12.75    ,.05,   
     4  -4.20      ,23.25      ,-11.66    ,17.17    ,2.31,  
     5  -6.14      ,-91.00     ,2.50      ,-29.57   ,-2.57, 
     6  1.75       ,3.05       ,-5.93     ,11.32    ,1.39,     
     7  -27.80     ,30.15      ,-9.60     ,15.15    ,1.29,  
     8  8.30       ,-6.78      ,6.16      ,-1.95    ,-.33,     
     9  23.54      ,-55.27     ,-2.76     ,-30.58   ,-2.96/
C         
C         
      DATA ((HN(7,I,J),J=1,5),I=1,9)/   
     1  -3.71      ,11.27      ,-.26      ,5.01     ,.19,      
     2  -.65       ,9.03       ,.24       ,4.14     ,.11,         
     3  6.15       ,-13.44     ,.01       ,-8.92    ,-.62,     
     4  2.60       ,-2.54      ,.02       ,-2.78    ,-.37,      
     5  2.42       ,-6.97      ,-.03      ,-4.13    ,-.35,     
     6  1.26       ,-1.57      ,-.05      ,-1.78    ,-.18,     
     7  1.23       ,-1.63      ,.03       ,-.94     ,-.04,       
     8  1.23       ,-4.30      ,.03       ,-2.75    ,-.23,      
     9  -.20       ,1.37       ,-.02      ,1.04     ,.08/        
C         
      DATA ((HN(8,I,J),J=1,5),I=1,9)/   
     1  -1.98      ,4.70       ,-.91      ,1.52     ,.04,       
     2  -.30       ,3.50       ,-.39      ,1.11     ,.07,        
     3  2.96       ,-5.59      ,.81       ,-3.34    ,-.17,      
     4  1.70       ,-3.09      ,-.50      ,-2.45    ,-.25,     
     5  -.72       ,.07        ,.76       ,.61      ,.07, 
     6  .49        ,.16        ,.05       ,-.14     ,-.01,
     7  .62        ,-.52       ,-.24      ,-.03     ,.02,         
     8  .08        ,-.92       ,.12       ,-.51     ,-.04,         
     9  -.56       ,2.13       ,-.05      ,1.35     ,.10/        
C         
      DATA ((HN(9,I,J),J=1,5),I=1,9)/   
     1  -1.87      ,1.97       ,.86       ,-.47     ,-.17,       
     2  .02        ,6.79       ,1.28      ,3.35     ,.14,         
     3  1.77       ,-.34       ,-.81      ,-.45     ,.05,        
     4  .60        ,.94        ,1.00      ,-.13     ,-.13,         
     5  .35        ,-1.88      ,-1.00     ,-1.22    ,-.09,     
     6  -.03       ,1.90       ,.04       ,.74      ,.03,
     7  1.11       ,-1.26      ,.20       ,-.81     ,-.04,       
     8  .08        ,-1.09      ,-.01      ,-.69     ,-.07,       
     9  -.19       ,1.46       ,.15       ,.77      ,.05/
C         
C         
C	Hermite polynomial coefficients (k,i,j) (left to right, k = 0
C	to k = 3) for each local-time Fourier component (top to bottom
C	a0,a1,b1,.....a6,b6 corresponding to an*cos(nt)+bn*sin(nt) where
C	t = local time in radians and n = 0 to n = 6  (i = 1,13).
C
C	HU corresponds to plasma scale height for topside F-layer
C	HL corresponds to plasma scale height for bottomside F-layer
C
C	NS = 1,2,3 implies equinox,summer,winter   (solar minimum)
C	NS = 4,5,6 implies equinox,summer,winter   (average solar activity)
C	NS = 7,8,9 implies equinox,summer,winter   (solar maximum)
C
C         
C         
      DATA ((HU(1,I,J),J=1,4),I=1,13)/   
     1  96.44      ,-34.94     ,0.17      ,-13.41,         
     2  2.56       ,-4.82      ,-0.28     ,-0.81, 
     3  -6.21      ,8.79       ,0.92      ,10.98,  
     4  5.62       ,-14.08     ,-0.77     ,-8.13,
     5  7.45       ,-3.20      ,-0.07     ,-4.05, 
     6  -2.40      ,-1.12      ,-0.12     ,0.24, 
     7  -0.41      ,2.33       ,-0.58     ,2.10,  
     8  -0.27      ,-1.45      ,-0.07     ,-2.34,
     9  -2.25      ,-2.46      ,0.04      ,-2.53, 
     1  -1.31      ,10.63      ,-0.04     ,4.84, 
     1  -3.65      ,4.32       ,0.41      ,1.14,   
     2  2.10       ,-4.66      ,0.01      ,-1.87,  
     3  -1.62      ,9.16       ,0.06      ,2.79/   
C
      DATA ((HU(2,I,J),J=1,4),I=1,13)/   
     1  94.93      ,-16.62     ,-2.34     ,-4.84,         
     2  -0.31      ,4.25       ,1.30      ,4.21,   
     3  -10.18     ,22.56      ,0.72      ,10.42,
     4  -1.30      ,3.12       ,0.60      ,-0.47,  
     5  3.42       ,7.53       ,-2.69     ,-0.80,  
     6  -2.39      ,2.59       ,2.87      ,2.49,   
     7  -0.36      ,2.70       ,-0.63     ,0.70,  
     8  1.52       ,-4.84      ,-1.70     ,-2.02, 
     9  -2.71      ,1.05       ,2.39      ,0.23,   
     1  -0.55      ,4.49       ,-0.48     ,1.86,  
     2  -0.38      ,-2.12      ,-2.93     ,-1.07,
     3  2.49       ,-4.91      ,2.64      ,-1.47,  
     4  0.01       ,3.04       ,-1.57     ,0.63/   
C      
      DATA ((HU(3,I,J),J=1,4),I=1,13)/   
     1  96.02      ,-25.26     ,0.62      ,-3.77,
     2  3.00       ,-13.19     ,-4.53     ,-3.06,
     3  -9.83      ,18.14      ,-0.28     ,9.55, 
     4  1.48       ,-5.10      ,-2.96     ,-2.82, 
     5  3.63       ,2.56       ,4.52      ,0.76,    
     6  -0.37      ,-1.69      ,-2.96     ,-0.96,
     7  -3.79      ,7.83       ,0.43      ,3.31,   
     8  -1.21      ,2.85       ,1.30      ,0.14,   
     9  2.59       ,-10.62     ,-2.03     ,-4.24,
     1  -1.51      ,6.07       ,2.16      ,2.07,   
     2  0.65       ,-2.40      ,0.47      ,-0.86,  
     3  -1.21      ,3.20       ,-2.14     ,1.55,  
     4  -0.12      ,3.90       ,0.59      ,0.59/   
C      
      DATA ((HL(1,I,J),J=1,4),I=1,13)/   
     1  75.55      ,-3.50      ,0.14      ,-1.86, 
     2  -5.25      ,-28.27     ,-0.08     ,-1.26,         
     3  1.45       ,-12.56     ,0.48      ,-0.39, 
     4  -1.77      ,11.70      ,-0.80     ,2.03, 
     5  11.10      ,-12.84     ,0.12      ,-2.15,
     6  -6.58      ,1.27       ,0.04      ,3.05,   
     7  4.79       ,-10.51     ,-0.38     ,-5.59,
     8  3.16       ,-1.93      ,0.74      ,-0.26,  
     9  0.18       ,-13.30     ,0.08      ,-5.49, 
     1  1.65       ,0.13       ,0.04      ,-0.48,   
     1  3.75       ,-4.96      ,0.68      ,-3.18,  
     2  -0.70      ,1.23       ,-0.87     ,1.22,  
     3  3.15       ,-2.48      ,0.00      ,-1.13/  
C
      DATA ((HL(2,I,J),J=1,4),I=1,13)/  
     1  73.44      ,1.58       ,7.56      ,-0.89,  
     2  8.46       ,-38.46     ,-10.33    ,-3.53,         
     3  -5.09      ,9.52       ,-2.21     ,5.23,  
     4  3.22       ,-3.28      ,3.58      ,-2.73,  
     5  -1.38      ,16.78      ,-5.50     ,5.50, 
     6  -2.20      ,-4.52      ,2.53      ,0.26,  
     7  -3.30      ,10.67      ,2.76      ,1.61,  
     8  0.98       ,0.72       ,-3.04     ,0.24,   
     9  -4.30      ,-0.72      ,0.13      ,1.54,  
     1  0.74       ,0.52       ,0.56      ,-1.71,   
     2  0.25       ,-0.83      ,0.00      ,-0.52,  
     3  -3.55      ,6.49       ,2.19      ,1.89,   
     4  4.28       ,-4.46      ,3.15      ,-1.40/  
C
      DATA ((HL(3,I,J),J=1,4),I=1,13)/  
     1  74.26      ,-4.66      ,-10.28    ,-1.90,         
     2  3.59       ,-27.30     ,8.46      ,-2.45, 
     3  -5.59      ,6.72       ,3.71      ,5.47,   
     4  0.20       ,1.53       ,-1.91     ,0.43,   
     5  0.19       ,5.53       ,1.39      ,1.72,    
     6  -4.93      ,4.57       ,-2.52     ,2.51,  
     7  4.16       ,-6.23      ,-5.06     ,-3.63, 
     8  3.26       ,-4.93      ,2.63      ,0.11,   
     9  -4.93      ,6.16       ,-3.03     ,1.40,  
     1  2.14       ,-3.12      ,-1.09     ,-0.44, 
     2  -1.08      ,4.72       ,-1.13     ,0.03,  
     3  -3.71      ,3.79       ,-0.45     ,0.99,  
     4  1.47       ,1.61       ,1.62      ,1.90/    
C
      DATA ((HU(4,I,J),J=1,4),I=1,13)/   
     1  98.50      ,-16.10     ,0.80      ,-11.91,         
     2  2.09       ,-4.86      ,0.02      ,0.46,   
     3  -6.53      ,-13.89     ,0.41      ,5.99, 
     4  3.53       ,-0.45      ,-0.19     ,-3.96, 
     5  8.78       ,-2.50      ,-0.18     ,-7.33, 
     6  -7.08      ,4.89       ,0.16      ,1.66,   
     7  -3.34      ,5.70       ,-0.69     ,3.53,  
     8  4.74       ,-7.28      ,0.79      ,-2.13,  
     9  1.94       ,-6.60      ,-0.22     ,-3.08, 
     1  2.05       ,0.19       ,-0.07     ,1.59,   
     2  -1.88      ,4.77       ,0.41      ,3.06,   
     3  2.63       ,-5.12      ,-0.61     ,-1.39, 
     4  5.67       ,-7.94      ,-0.32     ,-3.03/ 
C
      DATA ((HU(5,I,J),J=1,4),I=1,13)/  
     1  102.70     ,-16.13     ,-0.37     ,-8.20,        
     2  -3.03      ,3.13       ,-0.28     ,5.27,  
     3  -7.13      ,2.01       ,-1.15     ,8.29,  
     4  5.55       ,-6.55      ,2.49      ,-5.14,  
     5  4.61       ,3.09       ,-2.79     ,-2.96,  
     6  -2.24      ,-0.43      ,1.15      ,0.87,  
     7  -3.82      ,6.84       ,0.01      ,3.70,   
     8  1.01       ,-2.06      ,-0.81     ,-0.94, 
     9  0.03       ,-3.02      ,1.51      ,-2.06,  
     1  -0.72      ,5.64       ,-0.64     ,1.93,  
     2  -1.81      ,3.03       ,-0.12     ,1.20,  
     3  -0.95      ,2.61       ,1.14      ,1.18,   
     4  1.29       ,-1.38      ,-0.17     ,-2.05/ 
C
      DATA ((HU(6,I,J),J=1,4),I=1,13)/   
     1  103.62     ,-21.26     ,1.33      ,-9.01,         
     2  -3.52      ,3.47       ,-0.85     ,2.78,  
     3  -6.28      ,-0.26      ,1.82      ,5.55,  
     4  5.24       ,-5.57      ,-3.89     ,-6.17, 
     5  6.57       ,-4.63      ,2.87      ,-3.21,  
     6  -2.14      ,-1.55      ,-2.30     ,0.25, 
     7  -0.36      ,-0.05      ,0.40      ,0.12,  
     8  2.07       ,-2.87      ,0.41      ,-0.71,  
     9  -1.53      ,1.01       ,-1.92     ,-0.89, 
     1  0.76       ,1.12       ,2.47      ,0.78,    
     2  -1.26      ,2.57       ,-0.56     ,1.75,  
     3  0.35       ,-0.51      ,-1.36     ,-0.44, 
     4  1.08       ,-0.37      ,1.42      ,-0.68/  
C
      DATA ((HL(4,I,J),J=1,4),I=1,13)/   
     1  79.30      ,8.65       ,-0.47     ,-4.50, 
     2  -1.68      ,-36.25     ,0.15      ,-3.62,
     3  -3.96      ,-12.56     ,0.34      ,-1.26,
     4  12.86      ,-13.72     ,-1.25     ,-6.33,         
     5  3.91       ,-3.56      ,-0.13     ,-0.67, 
     6  1.75       ,-6.66      ,-0.20     ,-0.29, 
     7  4.00       ,-10.87     ,-0.33     ,-3.29,
     8  -1.49      ,3.19       ,-0.57     ,-0.27, 
     9  -7.43      ,10.38      ,0.21      ,3.44,  
     1  -0.22      ,0.47       ,0.01      ,0.82,   
     2  -0.76      ,3.63       ,0.47      ,0.06,   
     3  -4.37      ,5.69       ,-0.07     ,2.02,  
     4  -2.82      ,6.51       ,0.72      ,2.59/   
C
      DATA ((HL(5,I,J),J=1,4),I=1,13)/  
     1  82.80      ,0.07       ,7.58      ,-3.71,  
     2  0.39       ,-33.17     ,-10.09    ,-1.37,         
     3  -4.12      ,-4.15      ,0.26      ,0.56,  
     4  5.59       ,-2.82      ,3.40      ,-3.53,  
     5  1.25       ,5.06       ,-3.44     ,3.20,   
     6  1.20       ,-9.46      ,4.97      ,-2.55,  
     7  1.05       ,0.49       ,4.13      ,0.00,    
     8  1.57       ,-3.37      ,-2.71     ,-1.19, 
     9  -5.53      ,6.74       ,0.80      ,1.78,   
     1  0.87       ,0.21       ,-1.52     ,-0.14,  
     2  -0.60      ,3.90       ,-3.52     ,1.49,  
     3  1.43       ,-4.78      ,2.51      ,-1.88,  
     4  -1.47      ,5.52       ,-0.87     ,2.98/  
C
      DATA ((HL(6,I,J),J=1,4),I=1,13)/  
     1  82.90      ,-1.18      ,-11.67    ,-4.28,         
     2  -5.12      ,-19.33     ,8.77      ,-1.76,
     3  -6.78      ,-0.78      ,-0.94     ,2.43, 
     4  2.00       ,4.87       ,-5.20     ,-0.93,  
     5  4.80       ,-9.43      ,3.14      ,-1.12,  
     6  -0.94      ,-0.39      ,-2.70     ,-0.86,
     7  0.97       ,-1.69      ,-3.06     ,-0.35, 
     8  2.15       ,-1.96      ,2.90      ,-1.38,  
     9  -1.72      ,-2.96      ,0.17      ,-0.86, 
     1  -2.06      ,7.65       ,-0.33     ,3.17,  
     2  0.12       ,1.53       ,3.42      ,-0.62,   
     3  0.62       ,-2.73      ,-2.16     ,-0.44, 
     4  -3.76      ,9.78       ,-0.35     ,3.73/  
C
      DATA ((HU(7,I,J),J=1,4),I=1,13)/  
     1  103.61     ,1.74       ,0.66      ,-13.03,         
     2  3.32       ,0.05       ,-0.44     ,-1.27,
     3  -8.35      ,-30.64     ,0.29      ,5.49,         
     4  5.68       ,-11.42     ,-0.33     ,-6.69,        
     5  12.28      ,-4.85      ,0.21      ,-9.87,         
     6  -4.23      ,6.73       ,0.05      ,-0.09,
     7  -5.44      ,7.13       ,0.20      ,6.05, 
     8  1.29       ,-2.40      ,0.40      ,0.21, 
     9  4.21       ,-12.49     ,-0.23     ,-6.55,        
     1  1.65       ,-0.72      ,-0.24     ,-0.44,         
     2  -1.52      ,0.24       ,0.21      ,2.90, 
     3  -1.33      ,0.99       ,-0.97     ,3.21,
     4  5.58       ,-7.61      ,0.26      ,-3.14/
C
      DATA ((HU(8,I,J),J=1,4),I=1,13)/  
     1  104.96     ,-2.39      ,-0.16     ,-9.49,         
     2  -5.83      ,5.63       ,-3.12     ,9.31,  
     3  -6.70      ,-18.05     ,-0.45     ,4.61,
     4  6.72       ,-11.75     ,0.57      ,-5.61, 
     5  3.30       ,9.50       ,-0.02     ,-1.39,  
     6  -4.76      ,4.56       ,1.36      ,0.97,   
     7  -2.15      ,3.78       ,1.65      ,3.20,   
     8  1.60       ,-3.15      ,-0.62     ,-0.98, 
     9  0.29       ,-3.12      ,1.11      ,-1.70,  
     1  0.70       ,3.12       ,-0.80     ,1.67,   
     2  -0.68      ,1.25       ,-1.65     ,1.34,  
     3  1.04       ,-2.25      ,0.66      ,-0.79,  
     4  0.27       ,-0.55      ,0.29      ,0.08/   
C         
      DATA ((HU(9,I,J),J=1,4),I=1,13)/  
     1  107.72     ,-12.13     ,1.51      ,-10.25,        
     2  -6.33      ,9.64       ,2.33      ,5.76,   
     3  -8.99      ,-8.63      ,1.90      ,4.51,  
     4  3.13       ,-0.91      ,-3.47     ,-3.77, 
     5  5.25       ,1.05       ,0.68      ,-3.33,   
     6  -4.02      ,1.97       ,-0.34     ,0.33,  
     7  -3.48      ,5.68       ,-0.45     ,3.86,  
     8  3.96       ,-6.81      ,0.89      ,-3.28,  
     9  1.92       ,-4.38      ,-0.50     ,-3.09, 
     1  -0.50      ,4.69       ,0.83      ,2.27,   
     2  -0.03      ,-0.12      ,1.18      ,1.38,  
     3  1.54       ,-4.13      ,-0.28     ,-1.48, 
     4  -1.39      ,6.26       ,-0.42     ,1.56/  
C
       DATA ((HL(7,I,J),J=1,4),I=1,13)/  
     1  81.32      ,28.28      ,-0.12     ,-3.92,        
     2  -3.78      ,-38.35     ,0.10      ,-1.08,        
     3  -4.23      ,-30.42     ,0.52      ,-1.05,        
     4  12.35      ,-3.91      ,-0.55     ,-3.58,        
     5  -0.14      ,4.47       ,0.09      ,-0.58,
     6  -0.30      ,-14.61     ,0.39      ,0.21,         
     7  3.83       ,-4.68      ,-1.40     ,-3.75,         
     8  2.44       ,3.74       ,0.20      ,-1.84, 
     9  -2.34      ,4.69       ,-0.07     ,5.01,
     1  -5.56      ,13.41      ,0.17      ,6.91,
     2  0.70       ,5.32       ,0.36      ,0.24,  
     3  -0.90      ,-0.50      ,0.20      ,-3.80,         
     4  -7.85      ,15.04      ,-0.35     ,5.68/         
C
      DATA ((HL(8,I,J),J=1,4),I=1,13)/  
     1  86.46      ,10.72      ,12.32     ,-2.27,
     2  -2.34      ,-36.70     ,-11.46    ,-3.74,        
     3  -4.22      ,-8.85      ,-0.84     ,1.90, 
     4  8.30       ,-7.48      ,4.06      ,-4.03,  
     5  1.00       ,4.20       ,-5.48     ,0.93,   
     6  0.48       ,-3.59      ,5.76      ,-2.68,  
     7  -0.32      ,2.15       ,2.73      ,0.04,   
     8  2.02       ,-2.55      ,-1.82     ,-0.31, 
     9  -5.04      ,9.37       ,-2.28     ,3.20,  
     1  0.11       ,-0.72      ,-0.28     ,1.03,  
     2  0.12       ,2.35       ,-3.09     ,0.99,   
     3  -2.26      ,0.66       ,1.85      ,0.39,   
     4  -2.31      ,6.63       ,-1.24     ,2.25/  
C         
      DATA ((HL(9,I,J),J=1,4),I=1,13)/  
     1  84.85      ,11.73      ,-13.20    ,-2.68,         
     2  -6.96      ,-18.38     ,14.36     ,0.47,
     3  -5.50      ,-6.58      ,-0.64     ,0.99, 
     4  4.11       ,4.22       ,-3.38     ,-1.09,  
     5  3.06       ,-6.78      ,2.65      ,-2.43,  
     6  -3.33      ,4.17       ,-4.62     ,0.25,  
     7  2.25       ,-2.92      ,-4.86     ,-3.37, 
     8  1.47       ,1.12       ,3.85      ,1.03,    
     9  -2.12      ,1.39       ,-1.02     ,0.77,  
     1  4.10       ,-8.08      ,1.19      ,-2.12,  
     2  3.18       ,-3.38      ,1.90      ,-1.50,  
     3  -3.32      ,6.01       ,-1.51     ,1.77,  
     4  5.48       ,-11.98     ,0.26      ,-3.94/ 
C
      FC=1.         
      ATHM=ABS(THM) 
      IF(ATHM.GT.20.) FC=EXP(-(ATHM-20.)/4.)      
      CALL HERMW(THM,H0,H1,H2,H4)       
      DO 1 I=1,9    
      A(I)=HH(NS,I,1)+HH(NS,I,2)*H0+HH(NS,I,3)*H1+HH(NS,I,4)*H2       
     1  +HH(NS,I,5)*H4        
      B(I)=HN(NS,I,1)+HN(NS,I,2)*H0+HN(NS,I,3)*H1+HN(NS,I,4)*H2       
     1  +HN(NS,I,5)*H4        
1     CONTINUE      
      DO 10 I=1,13  
      AA(I)=HU(NS,I,1)+HU(NS,I,2)*H0+HU(NS,I,3)*H1+HU(NS,I,4)*H2      
      BB(I)=HL(NS,I,1)+HL(NS,I,2)*H0+HL(NS,I,3)*H1+HL(NS,I,4)*H2      
10    CONTINUE      
C         
      CHMF2=A(1)    
      CFOF2=B(1)    
      K=0 
      DO 2 N=2,8,2  
      K=K+1         
      W=K*PHI       
      CSF=COS(W)    
      SNF=SIN(W)    
C         
      CHMF2=CHMF2+A(N)*CSF+A(N+1)*SNF   
      CFOF2=CFOF2+B(N)*CSF+B(N+1)*SNF   
2     CONTINUE      
      K=0 
      SHUP=AA(1)    
      SHLO=BB(1)    
      DO 20 N=2,12,2
      K=K+1         
      W=K*PHI       
      CSF=COS(W)    
      SNF=SIN(W)    
      SHUP=SHUP+AA(N)*CSF+AA(N+1)*SNF   
      SHLO=SHLO+BB(N)*CSF+BB(N+1)*SNF   
20    CONTINUE      
      CHMF2=CHMF2*FC
      CFOF2=CFOF2*FC
      RETURN        
      END 
      SUBROUTINE HERMW(DEGLAT,H0,H1,H2,H4)        
C
C	Generates Hermite polynomials for use in CORRECT
C
      XSI=DEGLAT/14.
      XSI2=XSI*XSI  
      WT=EXP(-XSI2) 
      H0=WT         
      H1=2.*XSI*WT  
      H2=(-2.+4.*XSI2)*WT     
      H4=(12.-48.*XSI2+16.*XSI2*XSI2)*WT
      RETURN        
      END 