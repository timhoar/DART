
      SUBROUTINE VARX(Z,F1,F2)

C  DEFINES COEFF. TO INCORPORATE STRETCHING IN VERTICAL COORDINATE (FOR
C  DETAILS SEE FORBES, JGR, 1982, 5222-5240, ATMOSPHERIC TIDES, PART I)
      
      REAL CS(3),DS(2),ZS(2)
      DATA CS(1),CS(2),CS(3),DS(1),ZS(1),DS(2),ZS(2)
     1     /4.,.15,.025,.5,1.5,40.,120./
      CVT=1.0E-03
      F1=CS(1)
      F2=0.0
      DO 10 I=1,2
      ARG=(Z-ZS(I))/DS(I)
      ARG2=ZS(I)/DS(I)
      FACT1=.5*DS(I)*(CS(I+1)-CS(I))
      IF(ARG.GT.10.) ARG=10.
      FACT3=COSH(ARG)
      FACT4=TANH(ARG)
      F1=F1+FACT1*(1.+FACT4)/DS(I)
10    F2=F2-(FACT1/(DS(I)*DS(I)))/(FACT3*FACT3)
      F1=F1*CVT
      F2=F2*CVT*CVT
      RETURN
      END
