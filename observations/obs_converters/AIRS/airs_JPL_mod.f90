! This code is not protected by the DART copyright agreement.
! DART $Id$

! adapted from original JPL code, example AIRS readers
!
! updated for version 6 of the AIRS data formats
! added fields needed to support radiances
! removed unused items to streamline the code.
!
! april 2019, nsc

module airs_JPL_mod

! the contents of this file are an amalgam of:
!  airs_ret_typ.inc
!  airs_ret_struct.inc
!  airs_ret_rdr.f
! although in several cases they were modified to use
! fortran 90 derived types and free format continuation
! lines, and to not use the unknown syntax 'double precision'.

implicit none
public

      integer   AIRS_RET_GEOXTRACK                      
      parameter(AIRS_RET_GEOXTRACK                      =    30)
      integer   AIRS_RET_GEOTRACK                       
      parameter(AIRS_RET_GEOTRACK                       =    45)
      integer   AIRS_RET_STDPRESSURELEV                 
      parameter(AIRS_RET_STDPRESSURELEV                 =    28)
      integer   AIRS_RET_STDPRESSURELAY                 
      parameter(AIRS_RET_STDPRESSURELAY                 =    28)
      integer   AIRS_RET_AIRSXTRACK                     
      parameter(AIRS_RET_AIRSXTRACK                     =     3)
      integer   AIRS_RET_AIRSTRACK                      
      parameter(AIRS_RET_AIRSTRACK                      =     3)
      integer   AIRS_RET_CLOUD                          
      parameter(AIRS_RET_CLOUD                          =     2)
      integer   AIRS_RET_H2OPRESSURELEV                 
      parameter(AIRS_RET_H2OPRESSURELEV                 =    15)
      integer   AIRS_RET_H2OPRESSURELAY                 
      parameter(AIRS_RET_H2OPRESSURELAY                 =    14)

! Record holds an entire granule of airs_ret
type airs_granule_type

! Attributes
        integer      NumLandSurface
        integer      NumOceanSurface
        integer      start_year
        integer      start_month
        integer      start_day
        integer      start_hour
        integer      start_minute
        real         start_sec
        integer      granule_number

! Geolocation fields
        real*8 Latitude(AIRS_RET_GEOXTRACK, &
                                  AIRS_RET_GEOTRACK)
        real*8 Longitude(AIRS_RET_GEOXTRACK, & 
                                   AIRS_RET_GEOTRACK)
        real*8 Time(AIRS_RET_GEOXTRACK,  &
                              AIRS_RET_GEOTRACK)

! Data Fields
        real*8         sat_lat( 45)
        real*8         sat_lon( 45)
        byte           scan_node_type( 45)
        real           satzen( 30, 45)
        real           satazi( 30, 45)
        real           solzen( 30, 45)
        real           solazi( 30, 45)
        real           glintlat( 45)
        real           glintlon( 45)
        integer*2      sun_glint_distance( 30, 45)
        real           topog( 30, 45)
        real           topog_err( 30, 45)
        real           landFrac( 30, 45)
        real           landFrac_err( 30, 45)
        real           pressStd( 28)
        real           pressH2O( 15)
        real           latAIRS( 3, 3, 30, 45)
        real           lonAIRS( 3, 3, 30, 45)
        integer*2      Qual_Guess_PSurf( 30, 45)
        real           PSurfStd( 30, 45)
        integer        nSurfStd( 30, 45)
        real           PBest( 30, 45)
        real           PGood( 30, 45)
        integer*2      nBestStd( 30, 45)
        integer*2      nGoodStd( 30, 45)
        real           TAirStd( 28, 30, 45)
        real           TAirStdErr( 28, 30, 45)
        integer*2      TAirStd_QC( 28, 30, 45)
        real           TSurfAir( 30, 45)
        real           TSurfAirErr( 30, 45)
        integer*2      TSurfAir_QC( 30, 45)
        integer*2      Qual_Surf( 30, 45)
        real           TSurfStd( 30, 45)
        real           TSurfStdErr( 30, 45)
        integer*2      Qual_H2O( 30, 45)
        real           H2OMMRStd( 14, 30, 45)
        real           H2OMMRStdErr( 14, 30, 45)
        integer*2      H2OMMRStd_QC( 14, 30, 45)
        real           totH2OStd( 30, 45)
        real           totH2OStdErr( 30, 45)
        integer*2      totH2OStd_QC( 30, 45)
        integer        numCloud( 30, 45)
        real           TCldTopStd( 2, 30, 45)
        real           TCldTopStdErr( 2, 30, 45)
        real           PCldTopStd( 2, 30, 45)
        real           PCldTopStdErr( 2, 30, 45)
        real           CldFrcStd( 2, 3, 3, 30, 45)
        real           CldFrcStdErr( 2, 3, 3, 30, 45)
        real           totCldH2OStd( 30, 45)
        real           totCldH2OStdErr( 30, 45)
        byte           retrieval_type( 30, 45)
        byte           Startup( 30, 45)
END type airs_granule_type


contains

! This function is autogenerated by the mkezio program to read
! an AIRS swath of type "L2_Standard_atmospheric&surface_product" from file given by the
! file_name argument into a buffer pointed to by the airs_ret_gran
! argument.  The caller owns the buffer.  The entire granule
! is read -- every attribute and field, the whole lat/lon/time
! extent.
!
! Errors opening the file, etc. are fatal and cause STOP.
! Problems reading individual attributes or fields are reported to
! the console but do not interrupt program flow.

      subroutine airs_ret_rdr(file_name, airs_ret_gran, ver)
      character(len=*),        intent(in) :: file_name
      type(airs_granule_type), intent(out) :: airs_ret_gran
      integer,                 intent(in) :: ver

      integer :: statn           ! HDF-EOS status. 0 for success
      integer :: fid             ! HDF-EOS file ID
      integer :: swid            ! HDF-EOS swath ID
      integer :: nchar           ! Number of characters
      character*256 :: swathname   ! Name of swath
      integer :: nswath          ! Number of swaths
      integer :: start(10) = (/0,0,0,0,0, 0,0,0,0,0/)
                                ! start of each dimensions for Swath I/O
                                ! 0 => start with first element
      integer :: stride(10) = (/1,1,1,1,1, 1,1,1,1,1/)
                                ! stride of each dimensions for Swath I/O
                                ! 1 => use every element
      integer :: edge(10)       ! size of each dimension for swath I/O
                                ! will be set for each individual read
      integer :: swopen, swinqswath, swattach
      integer :: swrdfld, swrdattr
      integer :: swdetach, swclose

      fid = swopen(file_name, 1)
      if (fid .eq. -1) then
        print *, "Error ", fid, " opening file ", file_name
        stop
      end if

      ! Get name of swath(s)
      nswath = swinqswath(file_name, swathname, nchar)
      if (nswath .ne. 1) then
        print *, "swinqswath found ", nswath, " swaths for file ", &
                 file_name, " Need exactly 1"
        stop
      end if

      ! There's exactly one swath.  Make sure it is the right one.
      if (swathname .ne. &
         'L2_Standard_atmospheric&surface_product') then
        print *, "Error: bad swath name ", swathname, " in file ", &
                 file_name
        print *, "Expected L2_Standard_atmospheric&surface_product"
        stop
      end if

      ! Attach to (open) the one swath.
      swid = swattach(fid, swathname)
      if (swid .eq. -1) then
        print *, "Failed to attach to swath ", swathname, &
                 " in file ", file_name
        stop
      end if

! Attributes
      statn = swrdattr(swid, "NumLandSurface", &
                   airs_ret_gran%NumLandSurface)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "NumLandSurface"

      statn = swrdattr(swid, "NumOceanSurface", &
                   airs_ret_gran%NumOceanSurface)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "NumOceanSurface"

      statn = swrdattr(swid, "start_year", &
                   airs_ret_gran%start_year)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "start_year"

      statn = swrdattr(swid, "start_month", &
                   airs_ret_gran%start_month)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "start_month"

      statn = swrdattr(swid, "start_day", &
                   airs_ret_gran%start_day)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "start_day"

      statn = swrdattr(swid, "start_hour", &
                   airs_ret_gran%start_hour)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "start_hour"

      statn = swrdattr(swid, "start_minute", &
                   airs_ret_gran%start_minute)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "start_minute"

      statn = swrdattr(swid, "start_sec", &
                   airs_ret_gran%start_sec)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "start_sec"

      statn = swrdattr(swid, "granule_number", &
                   airs_ret_gran%granule_number)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "granule_number"

! Geolocation fields
      edge(1) = AIRS_RET_GEOXTRACK
      edge(2) = AIRS_RET_GEOTRACK
      statn = swrdfld(swid, "Latitude", start, stride, edge, &
                      airs_ret_gran%Latitude)
      if (statn .ne. 0)  &
        print *, "Error ", statn, " reading field Latitude"

      statn = swrdfld(swid, "Longitude", start, stride, edge, &
                      airs_ret_gran%Longitude)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field Longitude"

      statn = swrdfld(swid, "Time", start, stride, edge, &
                      airs_ret_gran%Time)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field Time"


! Data Fields
      edge(1) = 45
      statn = SWrdfld(swid, "sat_lat", &
                   start, stride, edge, &
                   airs_ret_gran%sat_lat)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "sat_lat"

      edge(1) = 45
      statn = SWrdfld(swid, "sat_lon", &
                   start, stride, edge, &
                   airs_ret_gran%sat_lon)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "sat_lon"

      edge(1) = 45
      statn = SWrdfld(swid, "scan_node_type", &
                   start, stride, edge, &
                   airs_ret_gran%scan_node_type)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "scan_node_type"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "satzen", &
                   start, stride, edge, &
                   airs_ret_gran%satzen)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "satzen"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "satazi", &
                   start, stride, edge, &
                   airs_ret_gran%satazi)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "satazi"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "solzen", &
                   start, stride, edge, &
                   airs_ret_gran%solzen)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "solzen"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "solazi", &
                   start, stride, edge, &
                   airs_ret_gran%solazi)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "solazi"

      edge(1) = 45
      statn = SWrdfld(swid, "glintlat", &
                   start, stride, edge, &
                   airs_ret_gran%glintlat)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "glintlat"

      edge(1) = 45
      statn = SWrdfld(swid, "glintlon", &
                   start, stride, edge, &
                   airs_ret_gran%glintlon)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "glintlon"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "sun_glint_distance", &
                   start, stride, edge, &
                   airs_ret_gran%sun_glint_distance)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "sun_glint_distance"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "topog", &
                   start, stride, edge, &
                   airs_ret_gran%topog)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "topog"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "topog_err", &
                   start, stride, edge, &
                   airs_ret_gran%topog_err)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "topog_err"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "landFrac", &
                   start, stride, edge, &
                   airs_ret_gran%landFrac)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "landFrac"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "landFrac_err", &
                   start, stride, edge, &
                   airs_ret_gran%landFrac_err)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "landFrac_err"

      edge(1) = 28
      statn = SWrdfld(swid, "pressStd", &
                   start, stride, edge, &
                   airs_ret_gran%pressStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "pressStd"

      edge(1) = 15
      statn = SWrdfld(swid, "pressH2O", &
                   start, stride, edge, &
                   airs_ret_gran%pressH2O)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "pressH2O"

      edge(4) = 45
      edge(3) = 30
      edge(2) = 3
      edge(1) = 3
      statn = SWrdfld(swid, "latAIRS", &
                   start, stride, edge, &
                   airs_ret_gran%latAIRS)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "latAIRS"

      edge(4) = 45
      edge(3) = 30
      edge(2) = 3
      edge(1) = 3
      statn = SWrdfld(swid, "lonAIRS", &
                   start, stride, edge, &
                   airs_ret_gran%lonAIRS)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "lonAIRS"

      if (ver == 5) then
         edge(2) = 45
         edge(1) = 30
         statn = SWrdfld(swid, "Qual_Guess_PSurf", &
                      start, stride, edge, &
                      airs_ret_gran%Qual_Guess_PSurf)
         if (statn .ne. 0) &
           print *, "Error ", statn, " reading field ", &
                     "Qual_Guess_PSurf"
      else
         airs_ret_gran%Qual_Guess_PSurf = -1
      endif  

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "PSurfStd", &
                   start, stride, edge, &
                   airs_ret_gran%PSurfStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "PSurfStd"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "nSurfStd", &
                   start, stride, edge, &
                   airs_ret_gran%nSurfStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "nSurfStd"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "PBest", &
                   start, stride, edge, &
                   airs_ret_gran%PBest)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "PBest"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "PGood", &
                   start, stride, edge, &
                   airs_ret_gran%PGood)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "PGood"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "nBestStd", &
                   start, stride, edge, &
                   airs_ret_gran%nBestStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "nBestStd"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "nGoodStd", &
                   start, stride, edge, &
                   airs_ret_gran%nGoodStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "nGoodStd"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 28
      statn = SWrdfld(swid, "TAirStd", &
                   start, stride, edge, &
                   airs_ret_gran%TAirStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "TAirStd"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 28
      statn = SWrdfld(swid, "TAirStdErr", &
                   start, stride, edge, &
                   airs_ret_gran%TAirStdErr)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "TAirStdErr"

      if (ver == 6) then
         edge(3) = 45
         edge(2) = 30
         edge(1) = 28
         statn = SWrdfld(swid, "TAirStd_QC", &
                      start, stride, edge, &
                      airs_ret_gran%TAirStd_QC)
         if (statn .ne. 0) &
           print *, "Error ", statn, " reading field ", &
                     "TAirStd_QC"
      else
         airs_ret_gran%TAirStd_QC = -1
      endif  

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "TSurfAir", &
                   start, stride, edge, &
                   airs_ret_gran%TSurfAir)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "TSurfAir"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "TSurfAirErr", &
                   start, stride, edge, &
                   airs_ret_gran%TSurfAirErr)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "TSurfAirErr"

      if (ver == 6) then
         edge(2) = 45
         edge(1) = 30
         statn = SWrdfld(swid, "TSurfAir_QC", &
                      start, stride, edge, &
                      airs_ret_gran%TSurfAir_QC)
         if (statn .ne. 0) &
           print *, "Error ", statn, " reading field ", &
                     "TSurfAir_QC"
      else
         airs_ret_gran%TSurfAir_QC = -1
      endif  

      if (ver == 5) then
         edge(2) = 45
         edge(1) = 30
         statn = SWrdfld(swid, "Qual_Surf", &
                      start, stride, edge, &
                      airs_ret_gran%Qual_Surf)
         if (statn .ne. 0) &
           print *, "Error ", statn, " reading field ", &
                     "Qual_Surf"
      else
         airs_ret_gran%Qual_Surf = -1
      endif  

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "TSurfStd", &
                   start, stride, edge, &
                   airs_ret_gran%TSurfStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "TSurfStd"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "TSurfStdErr", &
                   start, stride, edge, &
                   airs_ret_gran%TSurfStdErr)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "TSurfStdErr"


      if (ver == 5) then
         edge(2) = 45
         edge(1) = 30
         statn = SWrdfld(swid, "Qual_H2O", &
                      start, stride, edge, &
                      airs_ret_gran%Qual_H2O)
         if (statn .ne. 0) &
           print *, "Error ", statn, " reading field ", &
                     "Qual_H2O"
      else
         airs_ret_gran%Qual_H2O = -1
      endif  

      edge(3) = 45
      edge(2) = 30
      edge(1) = 14
      statn = SWrdfld(swid, "H2OMMRStd", &
                   start, stride, edge, &
                   airs_ret_gran%H2OMMRStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "H2OMMRStd"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 14
      statn = SWrdfld(swid, "H2OMMRStdErr", &
                   start, stride, edge, &
                   airs_ret_gran%H2OMMRStdErr)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "H2OMMRStdErr"

      if (ver == 6) then
         edge(3) = 45
         edge(2) = 30
         edge(1) = 14
         statn = SWrdfld(swid, "H2OMMRStd_QC", &
                      start, stride, edge, &
                      airs_ret_gran%H2OMMRStd_QC)
         if (statn .ne. 0) &
           print *, "Error ", statn, " reading field ", &
                     "H2OMMRStd_QC"
      else
         airs_ret_gran%H2OMMRStd_QC = -1
      endif  

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "totH2OStd", &
                   start, stride, edge, &
                   airs_ret_gran%totH2OStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "totH2OStd"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "totH2OStdErr", &
                   start, stride, edge, &
                   airs_ret_gran%totH2OStdErr)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "totH2OStdErr"

      if (ver == 6) then
         edge(2) = 45
         edge(1) = 30
         statn = SWrdfld(swid, "totH2OStd_QC", &
                      start, stride, edge, &
                      airs_ret_gran%totH2OStd_QC)
         if (statn .ne. 0) &
           print *, "Error ", statn, " reading field ", &
                     "totH2OStd_QC"
      else
         airs_ret_gran%totH2OStd_QC = -1
      endif  

      if (ver == 5) then
         edge(2) = 45
         edge(1) = 30
         statn = SWrdfld(swid, "numCloud", &
                      start, stride, edge, &
                      airs_ret_gran%numCloud)
         if (statn .ne. 0) &
           print *, "Error ", statn, " reading field ", &
                     "numCloud"
      else
         airs_ret_gran%numCloud = -1
      endif  

      if (ver == 5) then
         edge(3) = 45
         edge(2) = 30
         edge(1) = 2
         statn = SWrdfld(swid, "TCldTopStd", &
                      start, stride, edge, &
                      airs_ret_gran%TCldTopStd)
         if (statn .ne. 0) &
           print *, "Error ", statn, " reading field ", &
                     "TCldTopStd"

         edge(3) = 45
         edge(2) = 30
         edge(1) = 2
         statn = SWrdfld(swid, "TCldTopStdErr", &
                      start, stride, edge, &
                      airs_ret_gran%TCldTopStdErr)
         if (statn .ne. 0) &
           print *, "Error ", statn, " reading field ", &
                     "TCldTopStdErr"
   
         edge(3) = 45
         edge(2) = 30
         edge(1) = 2
         statn = SWrdfld(swid, "PCldTopStd", &
                      start, stride, edge, &
                      airs_ret_gran%PCldTopStd)
         if (statn .ne. 0) &
           print *, "Error ", statn, " reading field ", &
                     "PCldTopStd"
   
         edge(3) = 45
         edge(2) = 30
         edge(1) = 2
         statn = SWrdfld(swid, "PCldTopStdErr", &
                      start, stride, edge, &
                      airs_ret_gran%PCldTopStdErr)
         if (statn .ne. 0) &
           print *, "Error ", statn, " reading field ", &
                  "PCldTopStdErr"
      else
         airs_ret_gran%TCldTopStd = -1
         airs_ret_gran%TCldTopStdErr = -1
         airs_ret_gran%PCldTopStd = -1
         airs_ret_gran%PCldTopStdErr = -1
      endif  

      edge(5) = 45
      edge(4) = 30
      edge(3) = 3
      edge(2) = 3
      edge(1) = 2
      statn = SWrdfld(swid, "CldFrcStd", &
                   start, stride, edge, &
                   airs_ret_gran%CldFrcStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "CldFrcStd"

      edge(5) = 45
      edge(4) = 30
      edge(3) = 3
      edge(2) = 3
      edge(1) = 2
      statn = SWrdfld(swid, "CldFrcStdErr", &
                   start, stride, edge, &
                   airs_ret_gran%CldFrcStdErr)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "CldFrcStdErr"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "totCldH2OStd", &
                   start, stride, edge, &
                   airs_ret_gran%totCldH2OStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "totCldH2OStd"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "totCldH2OStdErr", &
                   start, stride, edge, &
                   airs_ret_gran%totCldH2OStdErr)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "totCldH2OStdErr"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "retrieval_type", &
                   start, stride, edge, &
                   airs_ret_gran%retrieval_type)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "retrieval_type"

      if (ver == 5) then
         edge(2) = 45
         edge(1) = 30
         statn = SWrdfld(swid, "Startup", &
                      start, stride, edge, &
                      airs_ret_gran%Startup)
         if (statn .ne. 0) &
           print *, "Error ", statn, " reading field ", &
                     "Startup"
      else
         airs_ret_gran%Startup = -1
      endif  


      ! Final clean-up
      statn = swdetach(swid)
      if (statn .ne. 0) &
        print *, "Error detaching from input file ", file_name
      statn = swclose(fid)
      if (statn .ne. 0) &
        print *, "Error closing input file ", file_name

      return
end subroutine

end module

