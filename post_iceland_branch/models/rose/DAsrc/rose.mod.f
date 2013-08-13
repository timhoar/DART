c----------------------------------------------------------------------
      module time_gcm 
c----------------------------------------------------------------------
      
      use params

      implicit none
      
      save

c... N(2D) global mean from TIME-GCM

      real ::  n2d(nz) = (/
     $ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     $ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     $ 0.00,      3.607e-17, 9.562e-17, 2.204e-16, 4.369e-16, 8.689e-16,
     $ 3.655e-15, 3.154e-14, 8.326e-14, 8.010e-13, 4.401e-12, 1.424e-11,
     $ 4.013e-11, 1.453e-10
     $ /)

c      real ::  n2d(nz) = (/
c     $ 2.062e-20, 2.062e-20, 2.062e-20, 2.062e-20, 2.062e-20, 2.062e-20,
c     $ 2.062e-20, 2.062e-20, 2.062e-20, 2.062e-20, 2.062e-20, 2.062e-20,
c     $ 2.062e-20, 2.079e-20, 2.180e-20, 2.522e-20, 3.749e-20, 6.325e-20,
c     $ 1.004e-19, 2.268e-19, 5.393e-19, 1.183e-18, 2.649e-18, 7.671e-18,
c     $ 1.910e-17, 3.607e-17, 9.562e-17, 2.204e-16, 4.369e-16, 8.689e-16,
c     $ 3.655e-15, 3.154e-14, 8.326e-14, 8.010e-13, 4.401e-12, 1.424e-11,
c     $ 4.013e-11, 1.453e-10
c     $ /)

      end module time_gcm 
c----------------------------------------------------------------------