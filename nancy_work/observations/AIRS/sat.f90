      implicit none
      integer*4  iqc, ipc , itype, iprof
      character  subset*6, mday*2, mday1*2
!
      real*8 lon, lat, lev, obs, aug, prof, time, obserr, pi
      integer iday, ii

!--------------------------------------------------------
      pi = 3.1415926
!
      do 2000 iday =13, 16
      write(mday, '(i2.2)') iday
      write(mday1, '(i2.2)') iday+1

      open(unit = 75, file = 'satemp_obs.200708'//mday)
      open(unit = 76, file = 'airs_obs.200708'//mday)     

!----------------------------------------------------------
!   read in GPS RO data locations first
!
      aug = 0.0
      itype = 173
      iqc = 0
      ipc = 1
      subset = 'AIRX2R'

obsloop:  do ii=1, 100000
      read(75, *, end=102) lon, lat, lev, obs, prof, time

      iprof = prof

      if (iprof < 5000000 ) then
!  temp errors
!   100       -99.000        0
!   150         1.445     3287
!   200         1.454     3287
!   250         1.056     3285
!   300         1.090     3286
!   400         1.061     3287
!   500         0.890     3287
!   700         0.990     6572
!   850         1.384     3284
!   925         1.380     3269
!  1000         1.026     3214

      obserr = 1.0                     ! for temperature (K)
      if(lev == 250.0 ) obserr = 1.2
      if(lev  < 250.0 ) obserr = 1.5
      if(lev  >= 850.0 ) obserr = 1.3

      else                             ! for moisture
!     use an approxiamte estimae of the error from (O-F)
!   300         0.207     3287
!   400         0.545     3287
!   500         0.980     3287
!   700         1.714     6574
!   850         3.270     3286
!   925         1.622     3286
!  1000         2.017     3252

      if(lev < 300.00) cycle obsloop   !  skip moisture observatons above 300mb

      if(lev >=750.0) obserr = 3.5
      if(lev < 750.0 .and. lev >= 600.0 ) obserr = 2.0
      if(lev < 600.0 .and. lev >= 450.0 ) obserr = 1.0
      if(lev < 450.0 .and. lev >= 350.0 ) obserr = 0.5
      if(lev <=350.0) obserr = 0.2

!     obs = obs/(1.0 + 0.001*obs)      ! convert to specific humidity (g/kg)
      endif

!----------------------------------------------
!   select the obs for the specific domain
!
      if( obs >= 0.0 .and. obs < 400.0) then
!  for monsoon domain
!     if (lon > -95.0 .and. lon < -35.0) then
!     if (lat > 5.0 .and. lat < 27.0) then
!  for monsoon domain
!     if (lon > 85.0 .and. lon < 165.0) then
!     if (lat > -20.0 .and. lat < 40.0) then
!
!  for Dean domain
      if (lon > -90.0 .and. lon < -10.0) then
      if (lat > -8.0 .and. lat < 30.0) then
!----------------------------------------------

      if(lon < 0.0 ) lon = lon + 360.0
      lon = lon * pi/180.0
      lat = lat * pi/180.0
      write(76, 880) obserr, lon, lat, lev, obs, aug,  &
                     iprof, time, itype, iqc, subset, ipc
      endif
      endif
      endif

  enddo obsloop

  102 continue
       close(75)
       close(76)
 2000 continue

  880  format(f6.2, 2f8.3, e14.5, 2f9.2, i9, f7.2, 2i4,1x,a6,i2)

      stop           
      end
