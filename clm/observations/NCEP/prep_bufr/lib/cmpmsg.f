      SUBROUTINE CMPMSG(CF)

C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    CMPMSG
C   PRGMMR: ATOR            ORG: NP12       DATE: 2005-03-09
C
C ABSTRACT: THIS SUBROUTINE IS USED TO SPECIFY WHETHER OR NOT BUFR
C   MESSAGES CREATED BY FUTURE CALLS TO EITHER OF THE BUFR ARCHIVE
C   LIBRARY SUBROUTINES WRITSB OR WRITSA ARE TO BE COMPRESSED.
C   THIS SUBROUTINE CAN BE CALLED AT ANY TIME AFTER THE FIRST CALL
C   TO BUFR ARCHIVE LIBRARY SUBROUTINE OPENBF, AND THE POSSIBLE VALUES
C   FOR CF ARE 'N' (= 'NO', WHICH IS THE DEFAULT) AND 'Y' (= 'YES').
C
C PROGRAM HISTORY LOG:
C 2005-03-09  J. ATOR    -- ORIGINAL AUTHOR
C DART $Id$
C
C USAGE:    CALL CMPMSG (CF)
C   INPUT ARGUMENT LIST:
C     CF       - CHARACTER*1: FLAG INDICATING WHETHER BUFR MESSAGES
C                OUTPUT BY FUTURE CALLS TO WRITSB OR WRITSA ARE TO
C                BE COMPRESSED:
C                       'N' = 'NO' (THE DEFAULT)
C                       'Y' = 'YES'
C
C REMARKS:
C    THIS ROUTINE CALLS:        BORT     CAPIT
C    THIS ROUTINE IS CALLED BY: COPYSB   WRITCA   WRITCP
C                               Also called by application programs.
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 77
C   MACHINE:  PORTABLE TO ALL PLATFORMS
C
C$$$

      COMMON /MSGCMP/ CCMF

      CHARACTER*128 BORT_STR
      CHARACTER*1   CCMF, CF

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

      CALL CAPIT(CF)
      IF(CF.NE.'Y'.AND. CF.NE.'N') GOTO 900
      CCMF = CF 

C  EXITS
C  -----

      RETURN
900   WRITE(BORT_STR,'("BUFRLIB: CMPMSG - INPUT ARGUMENT IS ",A1,'//
     . '", IT MUST BE EITHER Y OR N")') CF
      CALL BORT(BORT_STR)
      END
