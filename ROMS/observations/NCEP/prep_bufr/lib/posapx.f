      SUBROUTINE POSAPX(LUNIT)

C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    POSAPX
C   PRGMMR: WOOLLEN          ORG: NP20       DATE: 1994-01-06
C
C ABSTRACT: THIS SUBROUTINE READS TO THE END OF THE FILE AND BACKSPACES
C   IN ORDER TO POSITION FOR APPEND.
C
C PROGRAM HISTORY LOG:
C 1994-01-06  J. WOOLLEN -- ORIGINAL AUTHOR
C 1998-07-08  J. WOOLLEN -- REPLACED CALL TO CRAY LIBRARY ROUTINE
C                           "ABORT" WITH CALL TO NEW INTERNAL BUFRLIB
C                           ROUTINE "BORT"
C 2000-09-19  J. WOOLLEN -- MAXIMUM MESSAGE LENGTH INCREASED FROM
C                           10,000 TO 20,000 BYTES
C 2003-11-04  S. BENDER  -- ADDED REMARKS/BUFRLIB ROUTINE
C                           INTERDEPENDENCIES
C 2003-11-04  D. KEYSER  -- UNIFIED/PORTABLE FOR WRF; ADDED
C                           DOCUMENTATION (INCLUDING HISTORY); OUTPUTS
C                           MORE COMPLETE DIAGNOSTIC INFO WHEN ROUTINE
C                           TERMINATES ABNORMALLY
C 2004-08-09  J. ATOR    -- MAXIMUM MESSAGE LENGTH INCREASED FROM
C                           20,000 TO 50,000 BYTES
C DART $Id$
C
C USAGE:    CALL POSAPX (LUNIT)
C   INPUT ARGUMENT LIST:
C     LUNIT    - INTEGER: FORTRAN LOGICAL UNIT NUMBER FOR BUFR FILE
C
C   INPUT FILES:
C     UNIT "LUNIT" - BUFR FILE
C
C REMARKS:
C    THIS ROUTINE CALLS:        BORT     RDMSGW
C    THIS ROUTINE IS CALLED BY: OPENBF
C                               Normally not called by any application
C                               programs.
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 77
C   MACHINE:  PORTABLE TO ALL PLATFORMS
C
C$$$

      INCLUDE 'bufrlib.prm'

      DIMENSION   MBAY(MXMSGLD4)

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      REWIND LUNIT
      IREC = 0

C  TRY TO READ TO THE END OF THE FILE
C  ----------------------------------

1     CALL RDMSGW(LUNIT,MBAY,IER)
      IF(IER.EQ.-1) GOTO 2
      IF(IER.EQ.-2) GOTO 3
      IREC = IREC+1
      GOTO 1

C  IF SUCCESSFUL, BACKSPACE FOR APPENDING AND RETURN
C  -------------------------------------------------

2     BACKSPACE LUNIT
      GOTO 100

C  IF AN I/O ERROR IS ENCOUNTERED, THEN REREAD THE GOOD RECORDS,
C  BACKSPACE FOR APPENDING AND RETURN
C  -----------------------------------------------------------------

3     REWIND LUNIT
      DO J=1,IREC
        CALL RDMSGW(LUNIT,MBAY,IER)
        IF(IER.EQ.-1) GOTO 2
        IF(IER.EQ.-2) GOTO 900
      ENDDO

C  EXITS
C  -----

100   RETURN
900   CALL BORT('BUFRLIB: POSAPX - ERROR READING A BUFR MESSAGE')
      END
