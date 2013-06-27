! get_name_info
! -------------
! Given the value of module-specific constant values and a name of 
! a variable in the COAMPS restart file, returns an integer 
! representing the dimension type and a two-dimensional integer 
! array containing the position of that variable in that particular 
! dimension set in the COAMPS restart file written by either 
! single-processor or multi-processor I/O. 
!  PARAMETERS
!   IN  DIM_TYPE_2D         Constant for 2-D dimension
!   IN  DIM_TYPE_3D         Constant for 3-D dimension
!   IN  DIM_TYPE_3DW        Constant for 3-D (w level) dimension
!   IN  SINGLEIO            Constant for single-processor I/O
!   IN  MULTIIO             Constant for multi-processor I/O
!   IN  var_name            the name of the variable to look up
!   OUT var_dim_type        the dimension type (2d/3d/3dw)
!                           1: 2D  2: 3D  3: 3DW
!   OUT var_record_num      the position of the variable in its
!                           particular dimension
subroutine get_name_info(DIM_TYPE_2D, DIM_TYPE_3D, DIM_TYPE_3DW,   &
                         SINGLEIO, MULTIIO, var_name, var_dim_type,&
                         var_record_num)
  integer, intent(in)                :: DIM_TYPE_2D
  integer, intent(in)                :: DIM_TYPE_3D
  integer, intent(in)                :: DIM_TYPE_3DW
  integer, intent(in)                :: SINGLEIO
  integer, intent(in)                :: MULTIIO
  character(len=*), intent(in)       :: var_name
  integer, intent(out)               :: var_dim_type
  integer, dimension(2), intent(out) :: var_record_num

  select case (trim(var_name))
  case('aalhs')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 1
    var_record_num(MULTIIO)  = 1
  case('acx')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 103
  case('acy')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 104
  case('akhm')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 2
    var_record_num(MULTIIO)  = 2
  case('akhu')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 3
    var_record_num(MULTIIO)  = 3
  case('akhv')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 4
    var_record_num(MULTIIO)  = 4
  case('albed')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 5
    var_record_num(MULTIIO)  = 5
  case('aln')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 6
    var_record_num(MULTIIO)  = 6
  case('amoob')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 7
    var_record_num(MULTIIO)  = 7
  case('anca')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 1
    var_record_num(MULTIIO)  = 1
  case('aoz')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 8
    var_record_num(MULTIIO)  = 8
  case('auz')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 9
    var_record_num(MULTIIO)  = 9
  case('avlhf')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 11
    var_record_num(MULTIIO)  = 11
  case('avshf')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 12
    var_record_num(MULTIIO)  = 12
  case('avstx')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 13
    var_record_num(MULTIIO)  = 13
  case('avsty')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 14
    var_record_num(MULTIIO)  = 14
  case('avz')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 10
    var_record_num(MULTIIO)  = 10
  case('bblhs')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 2
    var_record_num(MULTIIO)  = 2
  case('bcr')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 105
  case('bcx')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 106
  case('bcy')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 107
  case('bflux')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 15
    var_record_num(MULTIIO)  = 15
  case('blht')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 16
    var_record_num(MULTIIO)  = 16
  case('cclhs')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 3
    var_record_num(MULTIIO)  = 3
  case('ccx')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 108
  case('ccy')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 109
  case('cond')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 1
    var_record_num(MULTIIO)  = 1
  case('diab')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 2
    var_record_num(MULTIIO)  = 2
  case('dimpl')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 4
    var_record_num(MULTIIO)  = 4
  case('distx')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 17
    var_record_num(MULTIIO)  = 17
  case('disty')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 18
    var_record_num(MULTIIO)  = 18
  case('div3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 3
    var_record_num(MULTIIO)  = 3
  case('dqcdt')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 29
    var_record_num(MULTIIO)  = 29
  case('dqrdt')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 30
    var_record_num(MULTIIO)  = 30
  case('dqvdt')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 31
    var_record_num(MULTIIO)  = 31
  case('dtdt')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 32
    var_record_num(MULTIIO)  = 32
  case('dtdts')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 19
    var_record_num(MULTIIO)  = 19
  case('dxav1')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 20
    var_record_num(MULTIIO)  = 20
  case('dxav2')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 21
    var_record_num(MULTIIO)  = 21
  case('dxm')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 22
    var_record_num(MULTIIO)  = 22
  case('dxu')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 23
    var_record_num(MULTIIO)  = 23
  case('dxv')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 24
    var_record_num(MULTIIO)  = 24
  case('dyav1')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 25
    var_record_num(MULTIIO)  = 25
  case('dyav2')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 26
    var_record_num(MULTIIO)  = 26
  case('dym')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 27
    var_record_num(MULTIIO)  = 27
  case('dyu')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 28
    var_record_num(MULTIIO)  = 28
  case('dyv')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 29
    var_record_num(MULTIIO)  = 29
  case('dzsdx')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 30
    var_record_num(MULTIIO)  = 30
  case('dzsdxm')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 31
  case('dzsdy')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 31
    var_record_num(MULTIIO)  = 32
  case('dzsdym')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 33
  case('e1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 8
    var_record_num(MULTIIO)  = 8
  case('e2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 9
    var_record_num(MULTIIO)  = 9
  case('e3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 10
    var_record_num(MULTIIO)  = 10
  case('ebi')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 48
    var_record_num(MULTIIO)  = -1
  case('ebs')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 49
    var_record_num(MULTIIO)  = -1
  case('edis')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 32
    var_record_num(MULTIIO)  = 34
  case('eimpl')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 5
    var_record_num(MULTIIO)  = 5
  case('ekh')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 4
    var_record_num(MULTIIO)  = 4
  case('ekm')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 5
    var_record_num(MULTIIO)  = 5
  case('eks')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 6
    var_record_num(MULTIIO)  = 6
  case('exbm')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 7
    var_record_num(MULTIIO)  = 7
  case('exbw')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 6
    var_record_num(MULTIIO)  = 6
  case('fcd')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 110
  case('fcx')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 111
  case('fcy')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 112
  case('fimpl')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 7
    var_record_num(MULTIIO)  = 7
  case('fm')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 33
    var_record_num(MULTIIO)  = 35
  case('fu')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 34
    var_record_num(MULTIIO)  = 36
  case('fv')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 35
    var_record_num(MULTIIO)  = 37
  case('gimpl')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 8
    var_record_num(MULTIIO)  = 8
  case('grdi')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 36
    var_record_num(MULTIIO)  = 38
  case('grdj')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 37
    var_record_num(MULTIIO)  = 39
  case('grdrot')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 38
    var_record_num(MULTIIO)  = 40
  case('gwet')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 39
    var_record_num(MULTIIO)  = 41
  case('hflxl')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 40
    var_record_num(MULTIIO)  = 42
  case('hflxs')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 41
    var_record_num(MULTIIO)  = 43
  case('himpl')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 9
    var_record_num(MULTIIO)  = 9
  case('hlong')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 42
    var_record_num(MULTIIO)  = 44
  case('hxm')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 43
    var_record_num(MULTIIO)  = 45
  case('hxu')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 46
    var_record_num(MULTIIO)  = 48
  case('hxv')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 44
    var_record_num(MULTIIO)  = 46
  case('hym')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 45
    var_record_num(MULTIIO)  = 47
  case('hyu')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 47
    var_record_num(MULTIIO)  = 49
  case('hyv')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 48
    var_record_num(MULTIIO)  = 50
  case('nc1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 48
  case('nc2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 49
  case('nc3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 50
  case('ncn1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 54
  case('ncn2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 55
  case('ncn3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 56
  case('nr1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 51
  case('nr2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 52
  case('nr3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 53
  case('p1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 11
    var_record_num(MULTIIO)  = 11
  case('p2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 12
    var_record_num(MULTIIO)  = 12
  case('p3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 13
    var_record_num(MULTIIO)  = 13
  case('phi')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 49
    var_record_num(MULTIIO)  = 51
  case('phim')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 50
    var_record_num(MULTIIO)  = 52
  case('phis')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 51
    var_record_num(MULTIIO)  = 53
  case('ppbi')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 50
    var_record_num(MULTIIO)  = -1
  case('ppbs')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 51
    var_record_num(MULTIIO)  = -1
  case('precp')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 52
    var_record_num(MULTIIO)  = 54
  case('prlcl')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 53
    var_record_num(MULTIIO)  = 55
  case('prtop')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 54
    var_record_num(MULTIIO)  = 56
  case('psfc')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 55
    var_record_num(MULTIIO)  = 57
  case('psih')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 56
    var_record_num(MULTIIO)  = 58
  case('psim')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 57
    var_record_num(MULTIIO)  = 59
  case('psr')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 58
    var_record_num(MULTIIO)  = 60
  case('qc1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 33
    var_record_num(MULTIIO)  = 33
  case('qc2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 34
    var_record_num(MULTIIO)  = 34
  case('qc3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 35
    var_record_num(MULTIIO)  = 35
  case('qcbi')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 52
    var_record_num(MULTIIO)  = -1
  case('qcbs')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 53
    var_record_num(MULTIIO)  = -1
  case('qg1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 45
    var_record_num(MULTIIO)  = 45
  case('qg2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 46
    var_record_num(MULTIIO)  = 46
  case('qg3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 47
    var_record_num(MULTIIO)  = 47
  case('qgbi')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 60
    var_record_num(MULTIIO)  = -1
  case('qgbs')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 61
    var_record_num(MULTIIO)  = -1
  case('qi1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 36
    var_record_num(MULTIIO)  = 36
  case('qi2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 37
    var_record_num(MULTIIO)  = 37
  case('qi3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 38
    var_record_num(MULTIIO)  = 38
  case('qibi')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 54
    var_record_num(MULTIIO)  = -1
  case('qibs')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 55
    var_record_num(MULTIIO)  = -1
  case('qr1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 39
    var_record_num(MULTIIO)  = 39
  case('qr2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 40
    var_record_num(MULTIIO)  = 40
  case('qr3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 41
    var_record_num(MULTIIO)  = 41
  case('qrbi')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 56
    var_record_num(MULTIIO)  = -1
  case('qrbs')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 57
    var_record_num(MULTIIO)  = -1
  case('qs1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 42
    var_record_num(MULTIIO)  = 42
  case('qs2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 43
    var_record_num(MULTIIO)  = 43
  case('qs3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 44
    var_record_num(MULTIIO)  = 44
  case('qsbi')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 58
    var_record_num(MULTIIO)  = -1
  case('qsbs')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 59
    var_record_num(MULTIIO)  = -1
  case('qstar')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 59
    var_record_num(MULTIIO)  = 61
  case('qv1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 14
    var_record_num(MULTIIO)  = 14
  case('qv2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 15
    var_record_num(MULTIIO)  = 15
  case('qv3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 16
    var_record_num(MULTIIO)  = 16
  case('qvbi')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 62
    var_record_num(MULTIIO)  = -1
  case('qvbs')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 63
    var_record_num(MULTIIO)  = -1
  case('qvsea')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 60
    var_record_num(MULTIIO)  = 62
  case('rainc')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 61
    var_record_num(MULTIIO)  = 63
  case('raink')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 62
    var_record_num(MULTIIO)  = 64
  case('rainp')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 64
    var_record_num(MULTIIO)  = 66
  case('rains')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 63
    var_record_num(MULTIIO)  = 65
  case('raint')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 65
    var_record_num(MULTIIO)  = 67
  case('rbm')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 17
    var_record_num(MULTIIO)  = 17
  case('rbw')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 10
    var_record_num(MULTIIO)  = 10
  case('rib')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 66
    var_record_num(MULTIIO)  = 68
  case('sigdt')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 11
    var_record_num(MULTIIO)  = 11
  case('slpr')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 67
    var_record_num(MULTIIO)  = 69
  case('slw')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 68
    var_record_num(MULTIIO)  = 70
  case('snow')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 69
    var_record_num(MULTIIO)  = 71
  case('snowd')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 70
    var_record_num(MULTIIO)  = 72
  case('snowx')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 71
    var_record_num(MULTIIO)  = 73
  case('spdmax')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 101
    var_record_num(MULTIIO)  = -1
  case('ssw')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 72
    var_record_num(MULTIIO)  = 74
  case('stres')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 73
    var_record_num(MULTIIO)  = 75
  case('strx')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 74
    var_record_num(MULTIIO)  = 76
  case('stry')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 75
    var_record_num(MULTIIO)  = 77
  case('t1000')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 80
    var_record_num(MULTIIO)  = 82
  case('th1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 19
    var_record_num(MULTIIO)  = 19
  case('th2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 20
    var_record_num(MULTIIO)  = 20
  case('th3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 21
    var_record_num(MULTIIO)  = 21
  case('thbi')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 64
    var_record_num(MULTIIO)  = -1
  case('thbm')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 18
    var_record_num(MULTIIO)  = 18
  case('thbs')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 65
    var_record_num(MULTIIO)  = -1
  case('thbw')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 12
    var_record_num(MULTIIO)  = 12
  case('trad')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 22
    var_record_num(MULTIIO)  = 22
  case('tsea')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 76
    var_record_num(MULTIIO)  = 78
  case('tsoil')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 77
    var_record_num(MULTIIO)  = 79
  case('tstar')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 78
    var_record_num(MULTIIO)  = 80
  case('ttau0')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 79
    var_record_num(MULTIIO)  = 81
  case('u1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 23
    var_record_num(MULTIIO)  = 23
  case('u10m')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 97
    var_record_num(MULTIIO)  = 99
  case('u2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 24
    var_record_num(MULTIIO)  = 24
  case('u3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 25
    var_record_num(MULTIIO)  = 25
  case('ubi')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 66
    var_record_num(MULTIIO)  = -1
  case('ubs')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 67
    var_record_num(MULTIIO)  = -1
  case('ustar')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 81
    var_record_num(MULTIIO)  = 83
  case('uzgeu')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 16
  case('uzgev')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 17
  case('v1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 26
    var_record_num(MULTIIO)  = 26
  case('v10m')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 98
    var_record_num(MULTIIO)  = 100
  case('v2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 27
    var_record_num(MULTIIO)  = 27
  case('v3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 28
    var_record_num(MULTIIO)  = 28
  case('vbi')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 68
    var_record_num(MULTIIO)  = -1
  case('vbs')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 69
    var_record_num(MULTIIO)  = -1
  case('vzgeu')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 18
  case('vzgev')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 19
  case('w0avg')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 16
    var_record_num(MULTIIO)  = 20
  case('w1')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 13
    var_record_num(MULTIIO)  = 13
  case('w2')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 14
    var_record_num(MULTIIO)  = 14
  case('w3')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 15
    var_record_num(MULTIIO)  = 15
  case('wbi')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 17
    var_record_num(MULTIIO)  = -1
  case('wbs')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 18
    var_record_num(MULTIIO)  = -1
  case('wtm')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 82
    var_record_num(MULTIIO)  = 84
  case('wtu')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 83
    var_record_num(MULTIIO)  = 85
  case('wtv')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 84
    var_record_num(MULTIIO)  = 86
  case('xland')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 85
    var_record_num(MULTIIO)  = 87
  case('xnr')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 99
    var_record_num(MULTIIO)  = 101
  case('xpos')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 86
    var_record_num(MULTIIO)  = 88
  case('ynr')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 100
    var_record_num(MULTIIO)  = 102
  case('ypos')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 87
    var_record_num(MULTIIO)  = 89
  case('z0')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 96
    var_record_num(MULTIIO)  = 98
  case('zaol')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 88
    var_record_num(MULTIIO)  = 90
  case('zerom')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 89
    var_record_num(MULTIIO)  = 91
  case('zerou')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 90
    var_record_num(MULTIIO)  = 92
  case('zerov')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 91
    var_record_num(MULTIIO)  = 93
  case('zsfc')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 92
    var_record_num(MULTIIO)  = 94
  case('zsfcu')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 93
    var_record_num(MULTIIO)  = 95
  case('zsfcv')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 94
    var_record_num(MULTIIO)  = 96
  case('zvkrm')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 95
    var_record_num(MULTIIO)  = 97
  case default
    write (*,*) "Can't match name " // var_name // "in restart"
  end select
end subroutine get_name_info
