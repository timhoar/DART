      SUBROUTINE ELEMDX(CARD,LUN)

C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    ELEMDX
C   PRGMMR: WOOLLEN          ORG: NP20       DATE: 1994-01-06
C
C ABSTRACT: THIS SUBROUTINE DECODES THE SCALE FACTOR, REFERENCE VALUE,
C   BIT WIDTH AND UNITS (I.E., THE "ELEMENTS") FROM A TABLE B MNEMONIC
C   DEFINITION CARD THAT WAS PREVIOUSLY READ FROM A USER-SUPPLIED BUFR
C   DICTIONARY TABLE FILE IN CHARACTER FORMAT BY BUFR ARCHIVE LIBRARY
C   SUBROUTINE RDUSDX.  THESE DECODED VALUES ARE THEN ADDED TO THE
C   ALREADY-EXISTING ENTRY FOR THAT MNEMONIC (BUILT IN RDUSDX) WITHIN
C   THE INTERNAL BUFR TABLE B ARRAY TABB(*,LUN) IN COMMON BLOCK
C   /TABABD/.
C
C PROGRAM HISTORY LOG:
C 1994-01-06  J. WOOLLEN -- ORIGINAL AUTHOR
C 1995-06-28  J. WOOLLEN -- INCREASED THE SIZE OF INTERNAL BUFR TABLE
C                           ARRAYS IN ORDER TO HANDLE BIGGER FILES
C 1998-07-08  J. WOOLLEN -- REPLACED CALL TO CRAY LIBRARY ROUTINE
C                           "ABORT" WITH CALL TO NEW INTERNAL BUFRLIB
C                           ROUTINE "BORT"
C 1999-11-18  J. WOOLLEN -- THE NUMBER OF BUFR FILES WHICH CAN BE
C                           OPENED AT ONE TIME INCREASED FROM 10 TO 32
C                           (NECESSARY IN ORDER TO PROCESS MULTIPLE
C                           BUFR FILES UNDER THE MPI)
C 2003-11-04  J. ATOR    -- ADDED DOCUMENTATION
C 2003-11-04  S. BENDER  -- ADDED REMARKS/BUFRLIB ROUTINE
C                           INTERDEPENDENCIES
C 2003-11-04  D. KEYSER  -- UNIFIED/PORTABLE FOR WRF; ADDED HISTORY
C                           DOCUMENTATION; OUTPUTS MORE COMPLETE
C                           DIAGNOSTIC INFO WHEN ROUTINE TERMINATES
C                           ABNORMALLY; CHANGED CALL FROM BORT TO BORT2
C DART $Id$
C
C USAGE:    CALL ELEMDX (CARD, LUN)
C   INPUT ARGUMENT LIST:
C     CARD     - CHARACTER*80: MNEMONIC DEFINITION CARD THAT WAS READ
C                FROM A USER-SUPPLIED BUFR DICTIONARY TABLE
C     LUN      - INTEGER: I/O STREAM INDEX INTO INTERNAL MEMORY ARRAYS 
C
C REMARKS:
C    THIS ROUTINE CALLS:        BORT2    CAPIT    JSTCHR   JSTNUM
C                               NEMTAB
C    THIS ROUTINE IS CALLED BY: RDUSDX
C                               Normally not called by any application
C                               programs.
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 77
C   MACHINE:  PORTABLE TO ALL PLATFORMS
C
C$$$

      INCLUDE 'bufrlib.prm'

      COMMON /TABABD/ NTBA(0:NFILES),NTBB(0:NFILES),NTBD(0:NFILES),
     .                MTAB(MAXTBA,NFILES),IDNA(MAXTBA,NFILES,2),
     .                IDNB(MAXTBB,NFILES),IDND(MAXTBD,NFILES),
     .                TABA(MAXTBA,NFILES),TABB(MAXTBB,NFILES),
     .                TABD(MAXTBD,NFILES)

      CHARACTER*600 TABD
      CHARACTER*128 BORT_STR1,BORT_STR2
      CHARACTER*128 TABB
      CHARACTER*128 TABA
      CHARACTER*80  CARD
      CHARACTER*24  UNIT
      CHARACTER*11  REFR,REFR_ORIG
      CHARACTER*8   NEMO
      CHARACTER*4   SCAL,SCAL_ORIG
      CHARACTER*3   BITW,BITW_ORIG
      CHARACTER*1   SIGN,TAB

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C  CAPTURE THE VARIOUS ELEMENTS CHARACTERISTICS
C  --------------------------------------------

      NEMO = CARD( 3:10)
      SCAL = CARD(14:17)
      REFR = CARD(21:31)
      BITW = CARD(35:37)
      UNIT = CARD(41:64)
c  .... Make sure the units are all capitalized
      CALL CAPIT(UNIT)

C  FIND THE ELEMENT TAG IN TABLE B
C  -------------------------------

C     Note that an entry for this mnemonic should already exist within
C     the internal BUFR Table B array TABB(*,LUN); this entry should
C     have been created by subroutine RDUSDX when the mnemonic and its
C     associated FXY value and description were initially defined within
C     a card read from the "Descriptor Definition" section at the top of
C     the user-supplied BUFR dictionary table in character format.  Now,
C     we need to retrieve the positional index for that entry within
C     TABB(*,LUN) so that we can access the entry and then add the scale
C     factor, reference value, bit width, and units to it.

      CALL NEMTAB(LUN,NEMO,IDSN,TAB,IELE)
      IF(TAB.NE.'B') GOTO 900

C  LEFT JUSTIFY AND STORE CHARACTERISTICS
C  --------------------------------------

      CALL JSTCHR(UNIT)
      TABB(IELE,LUN)(71:94) = UNIT

      SCAL_ORIG=SCAL
      CALL JSTNUM(SCAL,SIGN,IRET)
      IF(IRET.NE.0) GOTO 901
      TABB(IELE,LUN)(95:95) = SIGN
      TABB(IELE,LUN)(96:98) = SCAL

      REFR_ORIG=REFR
      CALL JSTNUM(REFR,SIGN,IRET)
      IF(IRET.NE.0) GOTO 902
      TABB(IELE,LUN)( 99: 99) = SIGN
      TABB(IELE,LUN)(100:109) = REFR

      BITW_ORIG=BITW
      CALL JSTNUM(BITW,SIGN,IRET)
      IF(IRET.NE.0  ) GOTO 903
      IF(SIGN.EQ.'-') GOTO 903
      TABB(IELE,LUN)(110:112) = BITW

C  EXITS
C  -----

      RETURN
900   WRITE(BORT_STR1,'("BUFRLIB: ELEMDX - CARD READ IN IS: ",A)') CARD
      WRITE(BORT_STR2,'(18X,"MNEMONIC ",A," IS NOT A TABLE B ENTRY '//
     . '(UNDEFINED, TAB=",A,")")') NEMO,TAB
      CALL BORT2(BORT_STR1,BORT_STR2)
901   WRITE(BORT_STR1,'("BUFRLIB: ELEMDX - CARD READ IN IS: ",A)') CARD
      WRITE(BORT_STR2,'(18X,"PARSED SCALE VALUE (=",A,") IS NOT '//
     . 'NUMERIC")') SCAL_ORIG
      CALL BORT2(BORT_STR1,BORT_STR2)
902   WRITE(BORT_STR1,'("BUFRLIB: ELEMDX - CARD READ IN IS: ",A)') CARD
      WRITE(BORT_STR2,'(18X,"PARSED REFERENCE VALUE (=",A,") IS NOT '//
     . 'NUMERIC")') REFR_ORIG
      CALL BORT2(BORT_STR1,BORT_STR2)
903   WRITE(BORT_STR1,'("BUFRLIB: ELEMDX - CARD READ IN IS: ",A)') CARD
      WRITE(BORT_STR2,'(18X,"PARSED BIT WIDTH VALUE (=",A,") IS NOT '//
     . 'NUMERIC")') BITW_ORIG
      CALL BORT2(BORT_STR1,BORT_STR2)
      END