      SUBROUTINE DXDUMP(LUNIT,LDXOT)

C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    DXDUMP
C   PRGMMR: ATOR             ORG: NP12       DATE: 2004-08-18
C
C ABSTRACT: THIS SUBROUTINE WRITES, TO LOGICAL UNIT LDXOT, AN ASCII
C   COPY OF THE BUFR DICTIONARY TABLE INFORMATION ASSOCIATED WITH
C   THE BUFR FILE DEFINED BY LOGICAL UNIT LUNIT.  IT IS ESPECIALLY
C   USEFUL FOR DETERMINING THE CONTENTS OF ARCHIVE BUFR FILES WHICH
C   MAY HAVE SUCH INFORMATION EMBEDDED AS DX MESSAGES AT THE FRONT
C   OF THE FILE.  THE OUTPUT FILE WILL BE IN A FORMAT SUITABLE FOR
C   SUBSEQUENT INPUT AS A USER-DEFINED DICTIONARY TABLES FILE TO
C   BUFR ARCHIVE LIBRARY SUBROUTINE OPENBF, AND IN THAT SENSE THIS
C   SUBROUTINE CAN BE VIEWED AS THE LOGICAL INVERSE OF BUFR ARCHIVE
C   LIBRARY SUBROUTINE RDUSDX.  NOTE THAT THE BUFR FILE ASSOCIATED
C   WITH LOGICAL UNIT LUNIT MUST HAVE ALREADY BEEN IDENTIFIED TO
C   THE BUFR ARCHIVE LIBRARY SOFTWARE VIA A PRIOR CALL TO OPENBF.
C
C PROGRAM HISTORY LOG:
C 2004-08-18  J. ATOR    -- ORIGINAL AUTHOR
C DART $Id$
C
C USAGE:    CALL DXDUMP (LUNIT, LDXOT)
C   INPUT ARGUMENT LIST:
C     LUNIT    - INTEGER: FORTRAN LOGICAL UNIT NUMBER FOR BUFR FILE
C     LDXOT    - INTEGER: FORTRAN LOGICAL UNIT NUMBER FOR OUTPUT FILE
C
C   INPUT FILES:
C     UNIT "LUNIT" - BUFR FILE WITH EMBEDDED DX DICTIONARY MESSAGES
C
C   OUTPUT FILES:
C     UNIT "LDXOT" - ASCII VERSION OF DX DICTIONARY INFORMATION, IN
C                    FORMAT SUITABLE FOR SUBSEQUENT INPUT TO OPENBF
C
C REMARKS:
C    THIS ROUTINE CALLS:        BORT     NEMTBD   STATUS   STRSUC
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

      COMMON /REPTAB/ IDNR(5,2),TYPS(5,2),REPS(5,2),LENS(5)
      COMMON /TABABD/ NTBA(0:NFILES),NTBB(0:NFILES),NTBD(0:NFILES),
     .                MTAB(MAXTBA,NFILES),IDNA(MAXTBA,NFILES,2),
     .                IDNB(MAXTBB,NFILES),IDND(MAXTBD,NFILES),
     .                TABA(MAXTBA,NFILES),TABB(MAXTBB,NFILES),
     .                TABD(MAXTBD,NFILES)

      CHARACTER*600 TABD
      CHARACTER*128 TABB
      CHARACTER*128 TABA
      CHARACTER*80  CARD,CARDI1,CARDI2,CARDI3,CARDI4
      CHARACTER*20  CMSTR
      CHARACTER*8   NEMS(MAXCD),WRK1,WRK2
      CHARACTER*6   ADN
      CHARACTER*3   TYPS
      CHARACTER*1   REPS

      DIMENSION     IRPS(MAXCD),KNTS(MAXCD)

      LOGICAL       TBSKIP, TDSKIP, XTRCI1

      DATA          CARDI1( 1:40)
     .              /'|          |        |                   '/
      DATA          CARDI1(41:80)
     .              /'                                       |'/
      DATA          CARDI2( 1:40)
     .              /'|          |                            '/
      DATA          CARDI2(41:80)
     .              /'                                       |'/
      DATA          CARDI3( 1:40)
     .              /'|          |      |             |     | '/
      DATA          CARDI3(41:80)
     .              /'                         |-------------|'/
      DATA          CARDI4( 1:40)
     .              /'|---------------------------------------'/
      DATA          CARDI4(41:80)
     .              /'---------------------------------------|'/

C-----------------------------------------------------------------------
      TBSKIP(ADN) = ((ADN.EQ.'063000').OR.(ADN.EQ.'063255').OR.
     .               (ADN.EQ.'031000').OR.(ADN.EQ.'031001').OR.
     .               (ADN.EQ.'031002'))
      TDSKIP(ADN) = ((ADN.EQ.'360001').OR.(ADN.EQ.'360002').OR.
     .               (ADN.EQ.'360003').OR.(ADN.EQ.'360004'))
C-----------------------------------------------------------------------

C     DETERMINE LUN FROM LUNIT.

      CALL STATUS(LUNIT,LUN,IL,IM)
      IF(IL.EQ.0) GOTO 900

C     CREATE AND WRITE OUT (TO LDXOT) THE HEADER CARDS FOR THE
C     DESCRIPTOR DEFINITION SECTION.

      CARD=CARDI4
      CARD( 1: 1)='.'
      CARD(80:80)='.'
      WRITE (LDXOT,'(A)') CARD

      CARD=CARDI4
      CARD( 2: 2)=' '
      CARD(79:79)=' '
      CARD(15:64)='   USER DEFINITIONS FOR TABLE-A TABLE-B TABLE D   '
      WRITE (LDXOT,'(A)') CARD

      WRITE (LDXOT,'(A)') CARDI4

      CARD=CARDI1
      CARD( 3:10)='MNEMONIC' 
      CARD(14:19)='NUMBER' 
      CARD(23:33)='DESCRIPTION'
      WRITE (LDXOT,'(A)') CARD

      CARD=CARDI4
      CARD(12:12)='|'
      CARD(21:21)='|'
      WRITE (LDXOT,'(A)') CARD

C     CREATE AND WRITE OUT (TO LDXOT) THE TABLE D DESCRIPTOR
C     DEFINITION CARDS.

      WRITE (LDXOT,'(A)') CARDI1

      XTRCI1=.FALSE.
      DO N=1,NTBD(LUN)
        IF(.NOT.TDSKIP(TABD(N,LUN)(1:6))) THEN
          CARD=CARDI1
          CARD( 3:10)=TABD(N,LUN)( 7:14) 
          CARD(14:19)=TABD(N,LUN)( 1: 6) 
          CARD(23:77)=TABD(N,LUN)(16:70)

C         CHECK IF THIS TABLE D MNEMONIC IS ALSO A TABLE A MNEMONIC.
C         IF SO, THEN LABEL IT AS SUCH AND ALSO CHECK IF IT IS THE
C         LAST OF THE TABLE A MNEMONICS, IN WHICH CASE AN EXTRA
C         CARDI1 LINE WILL BE WRITTEN TO LDXOT IN ORDER TO SEPARATE
C         THE TABLE A MNEMONICS FROM THE OTHER TABLE D MNEMONICS.

          DO NA=1,NTBA(LUN)
            IF(TABA(NA,LUN)(4:11).EQ.TABD(N,LUN)(7:14)) THEN
              CARD(14:14)='A'
              IF(NA.EQ.NTBA(LUN)) XTRCI1=.TRUE.
              GOTO 10
            END IF
          END DO
   10     WRITE (LDXOT,'(A)') CARD
          IF(XTRCI1) THEN
            WRITE (LDXOT,'(A)') CARDI1
            XTRCI1=.FALSE.
          END IF
        END IF
      END DO

C     CREATE AND WRITE OUT (TO LDXOT) THE TABLE B DESCRIPTOR
C     DEFINITION CARDS.

      WRITE (LDXOT,'(A)') CARDI1

      DO N=1,NTBB(LUN)
        IF(.NOT.TBSKIP(TABB(N,LUN)(1:6))) THEN
          CARD=CARDI1
          CARD( 3:10)=TABB(N,LUN)( 7:14) 
          CARD(14:19)=TABB(N,LUN)( 1: 6) 
          CARD(23:77)=TABB(N,LUN)(16:70)
          WRITE (LDXOT,'(A)') CARD
        END IF
      END DO

      WRITE (LDXOT,'(A)') CARDI1

C     CREATE AND WRITE OUT (TO LDXOT) THE HEADER CARDS FOR THE
C     SEQUENCE DEFINITION SECTION.

      WRITE (LDXOT,'(A)') CARDI4

      CARD=CARDI2
      CARD( 3:10)='MNEMONIC' 
      CARD(14:21)='SEQUENCE' 
      WRITE (LDXOT,'(A)') CARD

      CARD=CARDI4
      CARD(12:12)='|'
      WRITE (LDXOT,'(A)') CARD

C     CREATE AND WRITE OUT (TO LDXOT) THE TABLE D SEQUENCE
C     DEFINITION CARDS.

      WRITE (LDXOT,'(A)') CARDI2

      DO N=1,NTBD(LUN)
        IF(.NOT.TDSKIP(TABD(N,LUN)(1:6))) THEN
          CARD=CARDI2
          CARD( 3:10)=TABD(N,LUN)( 7:14) 
          IC = 14

C         GET THE LIST OF CHILD MNEMONICS FOR THIS TABLE D DESCRIPTOR,
C         AND THEN ADD EACH ONE (INCLUDING ANY REPLICATION TAGS) TO
C         THE SEQUENCE DEFINITION CARD FOR THIS TABLE D DESCRIPTOR.
 
          CALL NEMTBD(LUN,N,NSEQ,NEMS,IRPS,KNTS)
          IF(NSEQ.GT.0) THEN
            DO NC=1,NSEQ
              CMSTR=' '
              ICMS=0
              CALL STRSUC(NEMS(NC),WRK2,NCH)
              IF(IRPS(NC).NE.0) THEN

C               ADD THE OPENING REPLICATION TAG.

                ICMS=ICMS+1
                CMSTR(ICMS:ICMS)=REPS(IRPS(NC),1)
              END IF
              CMSTR(ICMS+1:ICMS+NCH)=WRK2(1:NCH)
              ICMS=ICMS+NCH
              IF(IRPS(NC).NE.0) THEN

C               ADD THE CLOSING REPLICATION TAG.

                ICMS=ICMS+1
                CMSTR(ICMS:ICMS)=REPS(IRPS(NC),2)
              END IF
              IF(KNTS(NC).NE.0) THEN

C               ADD THE FIXED REPLICATION COUNT.

                WRK1=' '
                WRITE (WRK1,'(I3)') KNTS(NC)
                CALL STRSUC(WRK1,WRK2,NCH)
                CMSTR(ICMS+1:ICMS+NCH)=WRK2(1:NCH)
                ICMS=ICMS+NCH
              END IF

C             WILL THIS CHILD (AND ITS REPLICATION TAGS, IF ANY) FIT
C             INTO THE CURRENT SEQUENCE DEFINITION CARD?  IF NOT, THEN
C             WRITE OUT (TO LDXOT) THE CURRENT CARD AND INITIALIZE A
C             NEW ONE TO HOLD THIS CHILD.

              IF(IC.GT.(79-ICMS)) THEN
                WRITE (LDXOT,'(A)') CARD
                CARD=CARDI2
                CARD( 3:10)=TABD(N,LUN)( 7:14) 
                IC = 14
              END IF
              CARD(IC:IC+ICMS-1)=CMSTR(1:ICMS)

C             NOTE THAT WE WANT TO LEAVE 2 BLANK SPACES BETWEEN EACH
C             CHILD WITHIN THE SEQUENCE DEFINITION CARD (TO IMPROVE
C             READABILITY).

              IC=IC+ICMS+2
            END DO
            WRITE (LDXOT,'(A)') CARD
            WRITE (LDXOT,'(A)') CARDI2
          END IF
        END IF
      END DO

C     CREATE AND WRITE OUT (TO LDXOT) THE HEADER CARDS FOR THE
C     ELEMENT DEFINITION SECTION.

      WRITE (LDXOT,'(A)') CARDI4

      CARD=CARDI3
      CARD( 3:10)='MNEMONIC' 
      CARD(14:17)='SCAL' 
      CARD(21:29)='REFERENCE' 
      CARD(35:37)='BIT' 
      CARD(41:45)='UNITS' 
      WRITE (LDXOT,'(A)') CARD

      CARD=CARDI4
      CARD(12:12)='|'
      CARD(19:19)='|'
      CARD(33:33)='|'
      CARD(39:39)='|'
      CARD(66:66)='|'
      WRITE (LDXOT,'(A)') CARD

C     CREATE AND WRITE OUT (TO LDXOT) THE TABLE B ELEMENT
C     DEFINITION CARDS.

      WRITE (LDXOT,'(A)') CARDI3

      DO N=1,NTBB(LUN)
        IF(.NOT.TBSKIP(TABB(N,LUN)(1:6))) THEN
          CARD=CARDI3
          CARD( 3:10)=TABB(N,LUN)( 7:14) 
          CARD(41:64)=TABB(N,LUN)(71:94) 

C         ADD THE SCALE FACTOR.

          CALL STRSUC(TABB(N,LUN)(96:98),WRK2,NCH)
          CARD(17-NCH+1:17)=WRK2 
          IF(TABB(N,LUN)(95:95).EQ.'-') CARD(17-NCH:17-NCH)='-'

C         ADD THE REFERENCE VALUE.

          CALL STRSUC(TABB(N,LUN)(100:109),WRK2,NCH)
          CARD(31-NCH+1:31)=WRK2 
          IF(TABB(N,LUN)(99:99).EQ.'-') CARD(31-NCH:31-NCH)='-'

C         ADD THE BIT WIDTH.

          CALL STRSUC(TABB(N,LUN)(110:112),WRK2,NCH)
          CARD(37-NCH+1:37)=WRK2 
          WRITE (LDXOT,'(A)') CARD
        END IF
      END DO

      WRITE (LDXOT,'(A)') CARDI3

C     CREATE AND WRITE OUT (TO LDXOT) THE CLOSING CARD.

      CARD=CARDI4
      CARD( 1: 1)='`'
      CARD(80:80)=''''
      WRITE (LDXOT,'(A)') CARD

      RETURN
900   CALL BORT('BUFRLIB: DXDUMP - BUFR FILE IS CLOSED, IT MUST BE'//
     . ' OPEN')

      END
