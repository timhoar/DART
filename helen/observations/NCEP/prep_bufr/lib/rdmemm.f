      SUBROUTINE RDMEMM(IMSG,SUBSET,JDATE,IRET)

C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    RDMEMM
C   PRGMMR: WOOLLEN          ORG: NP20       DATE: 1994-01-06
C
C ABSTRACT: THIS SUBROUTINE READS A PARTICULAR BUFR MESSAGE FROM
C   INTERNAL MEMORY (ARRAY MSGS IN COMMON BLOCK /MSGMEM/) INTO A
C   MESSAGE BUFFER (ARRAY MBAY IN COMMON BLOCK /BITBUF/).  IT IS
C   IDENTICAL TO BUFR ARCHIVE LIBRARY SUBROUTINE READMM EXCEPT IT DOES
C   NOT ADVANCE THE VALUE OF IMSG PRIOR TO RETURNING TO CALLING
C   PROGRAM.
C
C PROGRAM HISTORY LOG:
C 1994-01-06  J. WOOLLEN -- ORIGINAL AUTHOR
C 1998-07-08  J. WOOLLEN -- REPLACED CALL TO CRAY LIBRARY ROUTINE
C                           "ABORT" WITH CALL TO NEW INTERNAL BUFRLIB
C                           ROUTINE "BORT"; MODIFIED TO MAKE Y2K
C                           COMPLIANT
C 1999-11-18  J. WOOLLEN -- THE NUMBER OF BUFR FILES WHICH CAN BE
C                           OPENED AT ONE TIME INCREASED FROM 10 TO 32
C                           (NECESSARY IN ORDER TO PROCESS MULTIPLE
C                           BUFR FILES UNDER THE MPI); THE MAXIMUM
C                           NUMBER OF BYTES REQUIRED TO STORE ALL
C                           MESSAGES INTERNALLY WAS INCREASED FROM 4
C                           MBYTES TO 8 MBYTES
C 2000-09-19  J. WOOLLEN -- REMOVED MESSAGE DECODING LOGIC THAT HAD
C                           BEEN REPLICATED IN THIS AND OTHER READ
C                           ROUTINES AND CONSOLIDATED IT INTO A NEW
C                           ROUTINE CKTABA, CALLED HERE, WHICH IS
C                           ENHANCED TO ALLOW COMPRESSED AND STANDARD
C                           BUFR MESSAGES TO BE READ; MAXIMUM MESSAGE
C                           LENGTH INCREASED FROM 10,000 TO 20,000
C                           BYTES
C 2001-08-15  D. KEYSER  -- PARAMETER MAXMEM (THE MAXIMUM NUMBER OF
C                           BYTES REQUIRED TO STORE ALL MESSAGES
C                           INTERNALLY) WAS INCREASED FROM 8 MBYTES TO
C                           16 MBYTES
C 2003-11-04  S. BENDER  -- ADDED REMARKS/BUFRLIB ROUTINE
C                           INTERDEPENDENCIES
C 2003-11-04  D. KEYSER  -- PARAMETER MAXMSG (THE MAXIMUM NUMBER OF
C                           BUFR MESSAGES WHICH CAN BE STORED
C                           INTERNALLY) INCREASED FROM 50000 TO 200000;
C                           UNIFIED/PORTABLE FOR WRF; ADDED
C                           DOCUMENTATION (INCLUDING HISTORY); OUTPUTS
C                           MORE COMPLETE DIAGNOSTIC INFO WHEN ROUTINE
C                           TERMINATES ABNORMALLY OR UNUSUAL THINGS
C                           HAPPEN
C 2004-08-09  J. ATOR    -- MAXIMUM MESSAGE LENGTH INCREASED FROM
C                           20,000 TO 50,000 BYTES
C 2004-11-15  D. KEYSER  -- PARAMETER MAXMEM (THE MAXIMUM NUMBER OF
C                           BYTES REQUIRED TO STORE ALL MESSAGES
C                           INTERNALLY) WAS INCREASED FROM 16 MBYTES TO
C                           50 MBYTES
C DART $Id$
C
C USAGE:    CALL RDMEMM (IMSG, SUBSET, JDATE, IRET)
C   INPUT ARGUMENT LIST:
C     IMSG     - INTEGER: POINTER TO BUFR MESSAGE NUMBER (RECORD) IN
C                         STORAGE
C                       0 = resets the memory file
C
C   OUTPUT ARGUMENT LIST:
C     SUBSET   - CHARACTER*8: TABLE A MNEMONIC FOR TYPE OF BUFR MESSAGE
C                BEING READ
C     JDATE    - INTEGER: DATE-TIME STORED WITHIN SECTION 1 OF BUFR
C                MESSAGE BEING READ, IN FORMAT OF EITHER YYMMDDHH OR
C                YYYYMMDDHH, DEPENDING ON DATELEN() VALUE
C     IRET     - INTEGER: RETURN CODE:
C                       0 = normal return
C                      -1 = IMSG is either zero or greater than the
C                           number of messages in memory
C
C   OUTPUT FILES:
C     UNIT 06  - STANDARD OUTPUT PRINT
C
C REMARKS:
C    NOTE THAT UFBMEM IS CALLED PRIOR TO THIS TO STORE THE BUFR
C    MESSAGES INTO INTERNAL MEMORY.
C
C    THIS ROUTINE CALLS:        BORT     CKTABA   STATUS   WTSTAT
C    THIS ROUTINE IS CALLED BY: UFBMMS   UFBMNS   UFBRMS   UFBTAM
C                               Also called by application programs.
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 77
C   MACHINE:  PORTABLE TO ALL PLATFORMS
C
C$$$

      INCLUDE 'bufrlib.prm'

      COMMON /MSGCWD/ NMSG(NFILES),NSUB(NFILES),MSUB(NFILES),
     .                INODE(NFILES),IDATE(NFILES)
      COMMON /BITBUF/ MAXBYT,IBIT,IBAY(MXMSGLD4),MBYT(NFILES),
     .                MBAY(MXMSGLD4,NFILES)
      COMMON /MSGMEM/ MUNIT,MLAST,MSGP(0:MAXMSG),MSGS(MAXMEM)
      COMMON /QUIET / IPRT

      CHARACTER*8 SUBSET

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C  CHECK THE MESSAGE REQUEST AND FILE STATUS
C  -----------------------------------------

      CALL STATUS(MUNIT,LUN,IL,IM)
      CALL WTSTAT(MUNIT,LUN,IL, 1)
      IF(IL.EQ.0) GOTO 900
      IF(IL.GT.0) GOTO 901
      IRET = 0

      IF(IMSG.EQ.0 .OR.IMSG.GT.MSGP(0)) THEN
         CALL WTSTAT(MUNIT,LUN,IL,0)
         IF(IPRT.GE.1) THEN
      PRINT*
      PRINT*,'+++++++++++++++++++++++WARNING+++++++++++++++++++++++++'
            IF(IMSG.EQ.0)  THEN
               PRINT*, 'BUFRLIB: RDMEMM - REQUESTED MEMORY MESSAGE ',
     .          'NUMBER {FIRST (INPUT) ARGUMENT} IS 0, RETURN WITH ',
     .          'IRET = -1'
            ELSE
               PRINT*, 'BUFRLIB: RDMEMM - REQ. MEMORY MESSAGE NO. {',
     .          IMSG,' - {1ST (INPUT) ARG.} > NO. OF MESSAGES IN ',
     .          'MEMORY (',MSGP(0),'), RETURN WITH IRET = -1'
            ENDIF
      PRINT*,'+++++++++++++++++++++++WARNING+++++++++++++++++++++++++'
      PRINT*
         ENDIF
         IRET = -1
         GOTO 100
      ENDIF

C  READ MEMORY MESSAGE NUMBER IMSG INTO A MESSAGE BUFFER
C  -----------------------------------------------------

      IPTR = MSGP(IMSG)
      IF(IMSG.LT.MSGP(0)) LPTR = MSGP(IMSG+1)-IPTR
      IF(IMSG.EQ.MSGP(0)) LPTR = MLAST-IPTR+1
      IPTR = IPTR-1

      DO I=1,LPTR
      MBAY(I,LUN) = MSGS(IPTR+I)
      ENDDO

C  PARSE THE MESSAGE SECTION CONTENTS
C  ----------------------------------

      CALL CKTABA(LUN,SUBSET,JDATE,JRET)
      NMSG(LUN) = IMSG

C  EXITS
C  -----

100   RETURN
900   CALL BORT('BUFRLIB: RDMEMM - INPUT BUFR FILE IS CLOSED, IT '//
     . 'MUST BE OPEN FOR INPUT')
901   CALL BORT('BUFRLIB: RDMEMM - INPUT BUFR FILE IS OPEN FOR '//
     . 'OUTPUT, IT MUST BE OPEN FOR INPUT')
      END
