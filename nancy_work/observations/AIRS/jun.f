!The data is consisted of the AIRS single field retrievals created by CIMSS sounding team. It is a direct access binary file (big-endian). The following information should help user access file content.
!
!
!Total record number: 11008
!record length: 311
!data type: float32
!missing value: -9999.
!
!For each record:
!1        latitude
!2        longitude
!3        year
!4        month
!5        day
!6        hour
!7        minute
!8        second
!9:109    pressure level (101 levels)      (unit: mb)
!110:210  temperature profile (101 levels) (unit: K)
!211:311  moisture profile (101 levels)    (unit: g/kg, specific humidity)
!

!__________________________________________________________
!Sample of Fortran to read data file
!__________________________________________________________

!C This is a test to readin data from a binary file created for NCAR
!--------------------------------------
!   ifort -convert big_endian rrr.f
!   plot 'airs_obs.200706' u ($2)*180/3.14:($3)*180/3.14
!--------------------------------------

      implicit none
      integer leng, nrec, nlev
      parameter (leng=311, nrec=11008, nlev= 101)

      real*4 lon, lat, lev, obs, aug, time, obserr, pi
      real*4 prof(leng), pressure(nlev), temp(nlev), qv(nlev)

      integer*4  iqc, ipc , itype, num, numq, j, k, nr, iuf
      character file*100,   subset*6
!

      aug = 0.0
      itype = 173
      iqc = 0
      ipc = 1
      pi  = 3.1415926
      subset = 'AIRSLI'
      num = 1000000 
      numq= 5000000 

        file = 'AIRS_rtv_081507.bin'
        iuf = 90
        open(iuf,file=file,recl=leng,access='direct',
     $        status='old')

       open(unit = 76, file = 'airs_jun_obs.20070815')

!     take one obs from the nine AMSU pixels of the one AIRS pixel
      do 1000 nr=1,nrec, 6
!  2008/0405/
!     do 1000 nr=1,nrec

          read(iuf,rec=nr) prof

      lat   = prof(1)
      lon   = prof(2)
      time  = prof(6) + prof(7)/60.0 + prof(8)/3600.0

      if(lon < 0.0 ) lon = lon + 360.0
      lon = lon * pi/180.0
      lat = lat * pi/180.0

      do 100 j=9, 109
      k = j-8
      pressure(k) = prof(j)
  100 continue

      do 120 j=110, 210
      k = j-8-nlev
      temp(k) = prof(j)
  120 continue

      do 130 j=211, 311
      k = j-8-nlev*2
      qv(k) = prof(j)
  130 continue

!     choose the same (11) levels below 150mb as from the AIRS standard data.
      do 200 k=nlev-6, 1, -8

      lev = pressure(k)

!   for temperature
      if(lev >= 150.0) then
      num = num + 1
      obserr = 1.0                     ! for temperature (K)
      obs = temp(k)
      
      if(lev == 250.0 ) obserr = 1.2
      if(lev  < 250.0 ) obserr = 1.5
      if(lev  >= 850.0 ) obserr = 1.3

      if( obs .ne. -9999.0) then
      write(76, 880) obserr, lon, lat, lev, obs, aug,  
     &               num, time, itype, iqc, subset, ipc
      endif
      endif

!   for moisture
      if(lev >= 300.0) then
!     obserr = 2.0                     ! for specific humidity (g/kg/)
      if(lev >=750.0) obserr = 3.5
      if(lev < 750.0 .and. lev >= 600.0 ) obserr = 2.0
      if(lev < 600.0 .and. lev >= 450.0 ) obserr = 1.0
      if(lev < 450.0 .and. lev >= 350.0 ) obserr = 0.5
      if(lev <=350.0) obserr = 0.2

!  test
      obserr = 0.5 * obserr
!  test

      obs = qv(k)
      numq = numq + 1

      if( obs .ne. -9999.0) then
      write(76, 880) obserr, lon, lat, lev, obs, aug,  
     &               numq, time, itype, iqc, subset, ipc
      endif
      endif

  200   continue
 1000   continue

        stop 
  880  format(f6.2, 2f8.3, e14.5, 2f9.2, i9, f7.2, 2i4,1x,a6,i2)

        end

