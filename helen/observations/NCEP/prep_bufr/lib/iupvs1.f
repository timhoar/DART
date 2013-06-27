      FUNCTION IUPVS1(LUNIT,IL)

C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    IUPVS1
C   PRGMMR: ATOR             ORG: NP12       DATE: 2004-08-18
C
C ABSTRACT: THIS FUNCTION UNPACKS AND RETURNS THE BINARY INTEGER WORD
C   CONTAINED WITHIN BYTE IL OF SECTION 1 (OR BYTE 8 OF SECTION 0, IF
C   IL = 0) OF THE LAST BUFR MESSAGE THAT WAS READ FROM LOGICAL UNIT
C   NUMBER LUNIT VIA BUFR ARCHIVE LIBRARY SUBROUTINE READMG, READERME,
C   READIBM OR EQUIVALENT.  THIS FUNCTION IS SIMILAR TO BUFR ARCHIVE
C   LIBRARY FUNCTION IUPBS1 EXCEPT THAT IT OPERATES ON A BUFR MESSAGE
C   THAT HAS ALREADY BEEN READ INTO THE INTERNAL BUFR ARCHIVE LIBRARY
C   ARRAYS (VIA A PREVIOUS CALL TO READMG, READERME, READIBM, ETC.)
C   RATHER THAN ON A BUFR MESSAGE PASSED DIRECTLY INTO THE FUNCTION
C   VIA A MEMORY ARRAY.  NOTE THAT THIS FUNCTION IS CONSIDERED
C   OBSOLETE AND MAY BE REMOVED FROM THE BUFR ARCHIVE LIBRARY IN A
C   FUTURE VERSION; USERS SHOULD INSTEAD MIGRATE TO THE USE OF BUFR
C   ARCHIVE LIBRARY FUNCTION IUPVS01.
C
C PROGRAM HISTORY LOG:
C 2004-08-18  J. ATOR    -- ORIGINAL AUTHOR
C 2005-11-29  J. ATOR    -- MARKED AS OBSOLETE
C DART $Id$
C
C USAGE:    IUPVS1 (LUNIT, IL)
C   INPUT ARGUMENT LIST:
C     LUNIT    - INTEGER: FORTRAN LOGICAL UNIT NUMBER FOR BUFR FILE
C     IL       - INTEGER: BYTE TO UNPACK WITHIN SECTION 1 OF BUFR MSG
C                       0 = UNPACK BYTE 8 OF SECTION 0
C
C   OUTPUT ARGUMENT LIST:
C     IUPVS1   - INTEGER: UNPACKED INTEGER WORD
C
C REMARKS:
C    THIS ROUTINE CALLS:        BORT     IUPBS1   STATUS
C    THIS ROUTINE IS CALLED BY: None
C                               Normally called only by application
C                               programs.
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 77
C   MACHINE:  PORTABLE TO ALL PLATFORMS
C
C$$$

      INCLUDE 'bufrlib.prm'

      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(MXMSGLD4),MBYT(NFILES),
     .                MBAY(MXMSGLD4,NFILES)

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C  CHECK THE FILE STATUS
C  ---------------------

      CALL STATUS(LUNIT,LUN,ILST,IMST)
      IF(ILST.EQ.0) GOTO 900
      IF(ILST.GT.0) GOTO 901
      IF(IMST.EQ.0) GOTO 902

C  UNPACK THE REQUESTED BYTE
C  -------------------------

      IUPVS1 = IUPBS1(MBAY(1,LUN),IL)

C  EXITS
C  -----

      RETURN
900   CALL BORT('BUFRLIB: IUPVS1 - INPUT BUFR FILE IS CLOSED, IT '//
     . 'MUST BE OPEN FOR INPUT')
901   CALL BORT('BUFRLIB: IUPVS1 - INPUT BUFR FILE IS OPEN FOR '//
     . 'OUTPUT, IT MUST BE OPEN FOR INPUT')
902   CALL BORT('BUFRLIB: IUPVS1 - A MESSAGE MUST BE OPEN IN INPUT '//
     . 'BUFR FILE, NONE ARE')
      END
