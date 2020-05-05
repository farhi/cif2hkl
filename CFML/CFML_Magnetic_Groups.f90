!!-------------------------------------------------------
!!---- Crystallographic Fortran Modules Library (CrysFML)
!!-------------------------------------------------------
!!---- The CrysFML project is distributed under LGPL. In agreement with the
!!---- Intergovernmental Convention of the ILL, this software cannot be used
!!---- in military applications.
!!----
!!---- Copyright (C) 1999-2012  Institut Laue-Langevin (ILL), Grenoble, FRANCE
!!----                          Universidad de La Laguna (ULL), Tenerife, SPAIN
!!----                          Laboratoire Leon Brillouin(LLB), Saclay, FRANCE
!!----
!!----                          Laboratoire Leon Brillouin(LLB), Saclay, FRANCE
!!----
!!---- Authors: Juan Rodriguez-Carvajal (ILL)
!!----          Javier Gonzalez-Platas  (ULL)
!!----          Nebil Ayape Katcho      (ILL)
!!----
!!---- Contributors: Laurent Chapon     (ILL)
!!----               Marc Janoschek     (Los Alamos National Laboratory, USA)
!!----               Oksana Zaharko     (Paul Scherrer Institute, Switzerland)
!!----               Tierry Roisnel     (CDIFX,Rennes France)
!!----               Eric Pellegrini    (ILL)
!!----
!!---- This library is free software; you can redistribute it and/or
!!---- modify it under the terms of the GNU Lesser General Public
!!---- License as published by the Free Software Foundation; either
!!---- version 3.0 of the License, or (at your option) any later version.
!!----
!!---- This library is distributed in the hope that it will be useful,
!!---- but WITHOUT ANY WARRANTY; without even the implied warranty of
!!---- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!!---- Lesser General Public License for more details.
!!----
!!---- You should have received a copy of the GNU Lesser General Public
!!---- License along with this library; if not, see <http://www.gnu.org/licenses/>.
!!----
!!---- Special module Created by JRC in November 2012.
!!---- Code based in the program "read_magnetic_data.f" and the data base "magnetic_data.txt"
!!---- provided by Harold T. Stokes and Branton J. Campbell
!!---- Brigham Young University, Provo, Utah, USA
!!----
 Module CFML_Magnetic_Groups
   Use CFML_GlobalDeps
   Use CFML_Crystallographic_Symmetry, only: Get_Seitz_Symbol,Get_Trasfm_Symbol, &
                                             Get_Shubnikov_Operator_Symbol, Get_Transl_Symbol
   Use CFML_String_Utilities,          only: Get_LogUnit
   Use CFML_Math_3D,                   only: matrix_inverse
   Use CFML_Math_General,              only: modulo_lat

   !
   ! Declaration of all the variables that are defined in the file "magnetic_data.txt"
   !
   Implicit None
   public

   logical             :: err_magg=.false.,database_allocated=.false.
   character(len=256)  :: err_magg_mess=" "
   character(len=2)    :: end_line
   logical             :: mcif=.false.
   integer             :: current_group=1
   ! For the ith nonhexagonal point operator:
   Character(Len=8),  dimension(:),allocatable :: point_op_label   ! point_op_label(i): point operator symbol (from Litvin)
   Character(Len=10), dimension(:),allocatable :: point_op_xyz
   Integer,       dimension(:,:,:),allocatable :: point_op_matrix  ! point_op_matrix(i): point operator matrix
   ! For the ith hexagonal point operator:
   Character(Len=8),  dimension(:),allocatable :: point_op_hex_label  ! point_op_hex_label(i): point operator symbol (from Litvin)
   Character(Len=10), dimension(:),allocatable :: point_op_hex_xyz    ! point_op_hex_xyz(i): point operator in x,y,z notation
   Integer,       dimension(:,:,:),allocatable :: point_op_hex_matrix ! point_op_hex_matrix(i): point operator matrix

   ! Number of magnetic space groups
   Integer, Parameter :: magcount=1651
   ! For the ith magnetic space group
   Character(Len=12), dimension(  :),allocatable :: nlabel_bns           ! nlabel_bns(i): numerical label in BNS setting
   Integer,           dimension(:,:),allocatable :: nlabelparts_bns      ! nlabel_parts_bns(j,i): jth part of nlabel_bns
   Character(Len=14), dimension(  :),allocatable :: spacegroup_label_bns ! label_bns(i): group symbol
   Character(Len=12), dimension(  :),allocatable :: nlabel_og            ! nlabel_og(i): numerical label in OG setting
   Integer,           dimension(:,:),allocatable :: nlabelparts_og       ! nlabel_parts_og(j,i): jth part of nlabel_og
   Character(Len=14), dimension(  :),allocatable :: spacegroup_label_og  ! label_og(i): group symbol
   Integer,           dimension(  :),allocatable :: magtype              ! magtype(i): type of magnetic space group (1-4)
   ! BNS-OG transformation (if type-4)
   Integer,         dimension(:,:,:),allocatable :: bnsog_point_op     ! bnsog_point_op(j,k,i): 3x3 point operator part of transformation
   Integer,         dimension(:,  :),allocatable :: bnsog_origin       ! bnsog_origin(j,i): translation part of transformation
   Integer,         dimension(    :),allocatable :: bnsog_origin_denom ! bnsog_point_origin(i): common denominator
   Integer,         dimension(    :),allocatable :: ops_count          ! iops_count(i): number of point operators
   Integer,         dimension(    :),allocatable :: wyckoff_site_count ! wyckoff_count(i): number of wyckoff sites
   Integer,         dimension(: , :),allocatable :: wyckoff_pos_count  ! wyckoff_pos_count(j,i): number of positions in jth wyckoff site
   Integer,         dimension(: , :),allocatable :: wyckoff_mult       ! wyckoff_mult(j,i): multiplicity for jth wyckoff site
   Character(Len=1),dimension(: , :),allocatable :: wyckoff_label      ! wyckoff_label(j,i): symbol (a,b,c,...,z,alpha) for jth wyckoff site

   ! For BNS setting
   Integer, dimension(    :),allocatable :: lattice_bns_vectors_count  ! number of lattice vectors defining the lattice
   Integer, dimension(:,:,:),allocatable :: lattice_bns_vectors        ! (k,j,i): kth component of the jth lattice vector
   Integer, dimension(:,  :),allocatable :: lattice_bns_vectors_denom  !(j,i): common denominator
   ! For jth operator
   Integer, dimension(  :,:),allocatable :: ops_bns_point_op    ! ops_bns_point_op(j,i): point operator part
   Integer, dimension(:,:,:),allocatable :: ops_bns_trans       ! ops_bns_trans(k,j,i): kth component of translation part
   Integer, dimension(  :,:),allocatable :: ops_bns_trans_denom ! ops_bns_trans_denom(j,i): common denominator
   Integer, dimension(  :,:),allocatable :: ops_bns_timeinv     ! ops_bns_timeinv(j,i): 1=no time inversion, -1=time inversion
   ! For jth wyckoff site
   Integer, dimension(:,  :,:,:),allocatable :: wyckoff_bns_fract       ! wyckoff_bns_fract(k,j,i): kth component of fractional part of wyckoff position
   Integer, dimension(    :,:,:),allocatable :: wyckoff_bns_fract_denom ! wyckoff_bns_fract_denom(j,i): common denominator
   Integer, dimension(:,:,:,:,:),allocatable :: wyckoff_bns_xyz         ! wyckoff_bns_xyz(m,k,j,i): mth component to coeffcient of kth parameter (x,y,z)
   Integer, dimension(:,:,:,:,:),allocatable :: wyckoff_bns_mag  ! wyckoff_bns_mag(m,k,j,i): mth component to coeffcient of kth magnetic parameter (mx,my,mz)
   ! For OG setting (for type-4 groups)
   Integer, dimension(    :),allocatable :: lattice_og_vectors_count  ! lattice_og_vectors_count(i): number of lattice vectors defining the lattice
   Integer, dimension(:,:,:),allocatable :: lattice_og_vectors   ! lattice_og_vectors(k,j,i): kth component of the jth lattice vector
   Integer, dimension(:,  :),allocatable :: lattice_og_vectors_denom  ! lattice_og_vectors_denom(j,i): common denominator
   ! For jth operator
   Integer, dimension(  :,:),allocatable :: ops_og_point_op    ! ops_og_point_op(j,i): point operator part
   Integer, dimension(:,:,:),allocatable :: ops_og_trans       ! ops_og_trans(k,j,i): kth component of translation part
   Integer, dimension(  :,:),allocatable :: ops_og_trans_denom ! ops_og_trans_denom(j,i): common denominator
   Integer, dimension(  :,:),allocatable :: ops_og_timeinv     ! ops_og_timeinv(j,i): 1=no time inversion, -1=time inversion
   ! For jth wyckoff site
   Integer, dimension(:,  :,:,:),allocatable :: wyckoff_og_fract        ! wyckoff_og_fract(k,j,i): kth component of fractional part of wyckoff position
   Integer, dimension(    :,:,:),allocatable :: wyckoff_og_fract_denom  ! wyckoff_og_fract_denom(j,i): common denominator
   Integer, dimension(:,:,:,:,:),allocatable :: wyckoff_og_xyz          ! wyckoff_og_xyz(m,k,j,i): mth component to coefficient of kth parameter (x,y,z)
   Integer, dimension(:,:,:,:,:),allocatable :: wyckoff_og_mag          ! wyckoff_og_mag(m,k,j,i): mth component to coefficient of kth magnetic parameter (mx,my,mz)

  contains

  Subroutine Allocate_DataBase()
    if(database_allocated) return
    if(.not. allocated(point_op_label))             Allocate(point_op_label(48))
    if(.not. allocated(point_op_xyz))               Allocate(point_op_xyz(48))
    if(.not. allocated(point_op_matrix))            Allocate(point_op_matrix(3,3,48))
    if(.not. allocated(point_op_hex_label))         Allocate(point_op_hex_label(24))
    if(.not. allocated(point_op_hex_xyz))           Allocate(point_op_hex_xyz(24))
    if(.not. allocated(point_op_hex_matrix))        Allocate(point_op_hex_matrix(3,3,24))
    if(.not. allocated(nlabel_bns))                 Allocate(nlabel_bns(magcount))
    if(.not. allocated(nlabelparts_bns))            Allocate(nlabelparts_bns(2,magcount))
    if(.not. allocated(spacegroup_label_bns))       Allocate(spacegroup_label_bns(magcount))
    if(.not. allocated(nlabel_og))                  Allocate(nlabel_og(magcount))
    if(.not. allocated(nlabelparts_og))             Allocate(nlabelparts_og(3,magcount))
    if(.not. allocated(spacegroup_label_og))        Allocate(spacegroup_label_og(magcount))
    if(.not. allocated(magtype))                    Allocate(magtype(magcount))
    if(.not. allocated(bnsog_point_op))             Allocate(bnsog_point_op(3,3,magcount))
    if(.not. allocated(bnsog_origin))               Allocate(bnsog_origin(3,magcount))
    if(.not. allocated(bnsog_origin_denom))         Allocate(bnsog_origin_denom(magcount))
    if(.not. allocated(ops_count))                  Allocate(ops_count(magcount))
    if(.not. allocated(wyckoff_site_count))         Allocate(wyckoff_site_count(magcount))
    if(.not. allocated(wyckoff_pos_count))          Allocate(wyckoff_pos_count (27,magcount))
    if(.not. allocated(wyckoff_mult))               Allocate(wyckoff_mult(27,magcount))
    if(.not. allocated(wyckoff_label))              Allocate(wyckoff_label(27,magcount))
    if(.not. allocated(lattice_bns_vectors_count))  Allocate(lattice_bns_vectors_count (magcount))
    if(.not. allocated(lattice_bns_vectors))        Allocate(lattice_bns_vectors(3,6,magcount))
    if(.not. allocated(lattice_bns_vectors_denom))  Allocate(lattice_bns_vectors_denom(6,magcount))
    if(.not. allocated(ops_bns_point_op))           Allocate(ops_bns_point_op(96,magcount))
    if(.not. allocated(ops_bns_trans))              Allocate(ops_bns_trans(3,96,magcount))
    if(.not. allocated(ops_bns_trans_denom))        Allocate(ops_bns_trans_denom(96,magcount))
    if(.not. allocated(ops_bns_timeinv))            Allocate(ops_bns_timeinv(96,magcount))
    if(.not. allocated(wyckoff_bns_fract))          Allocate(wyckoff_bns_fract(3,96,27,magcount))
    if(.not. allocated(wyckoff_bns_fract_denom ))   Allocate(wyckoff_bns_fract_denom(96,27,magcount))
    if(.not. allocated(wyckoff_bns_xyz))            Allocate(wyckoff_bns_xyz(3,3,96,27,magcount))
    if(.not. allocated(wyckoff_bns_mag))            Allocate(wyckoff_bns_mag(3,3,96,27,magcount))
    if(.not. allocated(lattice_og_vectors_count))   Allocate(lattice_og_vectors_count(magcount))
    if(.not. allocated(lattice_og_vectors))         Allocate(lattice_og_vectors(3,6,magcount))
    if(.not. allocated(lattice_og_vectors_denom))   Allocate(lattice_og_vectors_denom(6,magcount))
    if(.not. allocated(ops_og_point_op))            Allocate(ops_og_point_op(96,magcount))
    if(.not. allocated(ops_og_trans))               Allocate(ops_og_trans(3,96,magcount))
    if(.not. allocated(ops_og_trans_denom))         Allocate(ops_og_trans_denom(96,magcount))
    if(.not. allocated(ops_og_timeinv))             Allocate(ops_og_timeinv(96,magcount))
    if(.not. allocated(wyckoff_og_fract))           Allocate(wyckoff_og_fract(3,96,27,magcount))
    if(.not. allocated(wyckoff_og_fract_denom))     Allocate(wyckoff_og_fract_denom (96,27,magcount))
    if(.not. allocated(wyckoff_og_xyz))             Allocate(wyckoff_og_xyz(3,3,96,27,magcount))
    if(.not. allocated(wyckoff_og_mag))             Allocate(wyckoff_og_mag(3,3,96,27,magcount))
    database_allocated=.true.
  End Subroutine Allocate_DataBase

  Subroutine deAllocate_DataBase()
    if(.not. database_allocated) return
    if(allocated(point_op_label))            deAllocate(point_op_label)
    if(allocated(point_op_xyz))              deAllocate(point_op_xyz)
    if(allocated(point_op_matrix))           deAllocate(point_op_matrix)
    if(allocated(point_op_hex_label))        deAllocate(point_op_hex_label)
    if(allocated(point_op_hex_xyz))          deAllocate(point_op_hex_xyz)
    if(allocated(point_op_hex_matrix))       deAllocate(point_op_hex_matrix)
    if(allocated(nlabel_bns))                deAllocate(nlabel_bns)
    if(allocated(nlabelparts_bns))           deAllocate(nlabelparts_bns)
    if(allocated(spacegroup_label_bns))      deAllocate(spacegroup_label_bns)
    if(allocated(nlabel_og))                 deAllocate(nlabel_og)
    if(allocated(nlabelparts_og))            deAllocate(nlabelparts_og)
    if(allocated(spacegroup_label_og ))      deAllocate(spacegroup_label_og)
    if(allocated(magtype))                   deAllocate(magtype)
    if(allocated(bnsog_point_op))            deAllocate(bnsog_point_op)
    if(allocated(bnsog_origin))              deAllocate(bnsog_origin)
    if(allocated(bnsog_origin_denom))        deAllocate(bnsog_origin_denom )
    if(allocated(ops_count))                 deAllocate(ops_count)
    if(allocated(wyckoff_site_count))        deAllocate(wyckoff_site_count)
    if(allocated(wyckoff_pos_count))         deAllocate(wyckoff_pos_count)
    if(allocated(wyckoff_mult))              deAllocate(wyckoff_mult)
    if(allocated(wyckoff_label))             deAllocate(wyckoff_label)
    if(allocated(lattice_bns_vectors_count)) deAllocate(lattice_bns_vectors_count)
    if(allocated(lattice_bns_vectors))       deAllocate(lattice_bns_vectors)
    if(allocated(lattice_bns_vectors_denom)) deAllocate(lattice_bns_vectors_denom)
    if(allocated(ops_bns_point_op))          deAllocate(ops_bns_point_op)
    if(allocated(ops_bns_trans))             deAllocate(ops_bns_trans)
    if(allocated(ops_bns_trans_denom))       deAllocate(ops_bns_trans_denom)
    if(allocated(ops_bns_timeinv))           deAllocate(ops_bns_timeinv)
    if(allocated(wyckoff_bns_fract))         deAllocate(wyckoff_bns_fract)
    if(allocated(wyckoff_bns_fract_denom))   deAllocate(wyckoff_bns_fract_denom )
    if(allocated(wyckoff_bns_xyz))           deAllocate(wyckoff_bns_xyz)
    if(allocated(wyckoff_bns_mag))           deAllocate(wyckoff_bns_mag)
    if(allocated(lattice_og_vectors_count))  deAllocate(lattice_og_vectors_count)
    if(allocated(lattice_og_vectors))        deAllocate(lattice_og_vectors)
    if(allocated(lattice_og_vectors_denom))  deAllocate(lattice_og_vectors_denom)
    if(allocated(ops_og_point_op))           deAllocate(ops_og_point_op)
    if(allocated(ops_og_trans))              deAllocate(ops_og_trans)
    if(allocated(ops_og_trans_denom))        deAllocate(ops_og_trans_denom)
    if(allocated(ops_og_timeinv))            deAllocate(ops_og_timeinv)
    if(allocated(wyckoff_og_fract))          deAllocate(wyckoff_og_fract)
    if(allocated(wyckoff_og_fract_denom))    deAllocate(wyckoff_og_fract_denom)
    if(allocated(wyckoff_og_xyz))            deAllocate(wyckoff_og_xyz)
    if(allocated(wyckoff_og_mag))            deAllocate(wyckoff_og_mag)
    database_allocated=.false.
  End Subroutine deAllocate_DataBase


  Subroutine read_magnetic_data()
    ! read data about magnetic space groups
    ! input data from magnetic_table.dat
    integer :: i,j,k,n,m,i_mag,ier
    character(len=512) :: fullprof_suite, database
    !*****************************************************************************
    ! open data file
    err_magg=.false.
    err_magg_mess=" "
    call Get_LogUnit(i_mag)
    call GET_ENVIRONMENT_VARIABLE("FULLPROF",fullprof_suite)
    n=len_trim(fullprof_suite)
    !write(*,*) trim(fullprof_suite)
    if(n == 0) then
      err_magg=.true.
      Select Case (OPS)
        Case(1)
           end_line=char(13)//char(10)
           write(unit=err_magg_mess,fmt="(a)") " => The FULLPROF environment variable is not defined! "//end_line// &
                                               "    This is needed for localizing the data base: magnetic_data.txt"//end_line// &
                                               "    that should be within the %FULLPROF%/Databases directory"
        Case Default
           end_line=char(10)
           write(unit=err_magg_mess,fmt="(a)") " => The FULLPROF environment variable is not defined! "//trim(end_line)// &
                                               "    This is needed for localizing the data base: magnetic_data.txt"//trim(end_line)// &
                                               "    that should be within the $FULLPROF/Databases directory"
      End Select
      return
    else
       if(fullprof_suite(n:n) /= OPS_SEP) then
         database=trim(fullprof_suite)//OPS_SEP//"Databases"//OPS_SEP//'magnetic_data.txt'
         !write(*,*) trim(database)
       else
         database=trim(fullprof_suite)//"Databases"//OPS_SEP//'magnetic_data.txt'
         !write(*,*) trim(database)
       end if
    end if
    Open(unit=i_mag,File=Trim(database),status="old",action="read",position="rewind",iostat=ier)
    if( ier /= 0) then
      err_magg=.true.
      err_magg_mess="     Error opening the data base: "//trim(database)
      return
    end if
    if(.not. database_allocated) call Allocate_DataBase()
    ! read nonhexagonal point operators
    err_magg=.false.
    Do i=1,48
      Read(i_mag,*)n,point_op_label(i),point_op_xyz(i),  &
          ((point_op_matrix(k,j,i),j=1,3),k=1,3)
      If(n /= i) then
        err_magg=.true.
        err_magg_mess= 'Error in numbering of nonhexagonal point operators'
        return
      End If
    End Do
    ! read hexagonal point operators
    Do i=1,24
      Read(i_mag,*) n,point_op_hex_label(i), point_op_hex_xyz(i),  &
          ((point_op_hex_matrix(k,j,i),j=1,3),k=1,3)
      If(n /= i)then
        err_magg=.true.
        err_magg_mess= 'Error in numbering of hexagonal point operators'
        return
      End If
    End Do
    ! read data for each magnetic space group
    Do i=1,1651
      Read(i_mag,*) (nlabelparts_bns(j,i),j=1,2),nlabel_bns(i),  &
          spacegroup_label_bns(i),(nlabelparts_og(j,i),j=1,3),  &
          nlabel_og(i),spacegroup_label_og(i)
      Read(i_mag,*) magtype(i)
      If(magtype(i) == 4) Then
        Read(i_mag,*) ((bnsog_point_op(j,k,i),j=1,3),k=1,3),  &
            (bnsog_origin(j,i),j=1,3),bnsog_origin_denom(i)
      End If
      Read(i_mag,*) ops_count(i)
      Read(i_mag,*) (ops_bns_point_op(j,i),(ops_bns_trans(k,j,i),k=1,3),  &
          ops_bns_trans_denom(j,i),ops_bns_timeinv(j,i), j=1,ops_count(i))
      Read(i_mag,*) lattice_bns_vectors_count(i)
      Read(i_mag,*) ((lattice_bns_vectors(k,j,i),k=1,3),  &
          lattice_bns_vectors_denom(j,i), j=1,lattice_bns_vectors_count(i))
      Read(i_mag,*) wyckoff_site_count(i)
      Do j=1,wyckoff_site_count(i)
        Read(i_mag,*) wyckoff_pos_count(j,i),wyckoff_mult(j,i), wyckoff_label(j,i)
        Do k=1,wyckoff_pos_count(j,i)
          Read(i_mag,*) (wyckoff_bns_fract(m,k,j,i),m=1,3),  &
              wyckoff_bns_fract_denom(k,j,i),  &
              ((wyckoff_bns_xyz(m,n,k,j,i),m=1,3),n=1,3),  &
              ((wyckoff_bns_mag(m,n,k,j,i),m=1,3),n=1,3)
        End Do
      End Do
      If(magtype(i) == 4) Then
        Read(i_mag,*) ops_count(i)
        Read(i_mag,*) (ops_og_point_op(j,i),(ops_og_trans(k,j,i),k=1,3),  &
            ops_og_trans_denom(j,i),ops_og_timeinv(j,i), j=1,ops_count(i))
        Read(i_mag,*) lattice_og_vectors_count(i)
        Read(i_mag,*) ((lattice_og_vectors(k,j,i),k=1,3),  &
            lattice_og_vectors_denom(j,i), j=1,lattice_og_vectors_count(i))
        Read(i_mag,*) wyckoff_site_count(i)
        Do j=1,wyckoff_site_count(i)
          Read(i_mag,*) wyckoff_pos_count(j,i),wyckoff_mult(j,i), wyckoff_label(j,i)
          Do k=1,wyckoff_pos_count(j,i)
            Read(i_mag,*) (wyckoff_og_fract(m,k,j,i),m=1,3),  &
                  wyckoff_og_fract_denom(k,j,i),              &
                ((wyckoff_og_xyz(m,n,k,j,i),m=1,3),n=1,3),    &
                ((wyckoff_og_mag(m,n,k,j,i),m=1,3),n=1,3)
          End Do
        End Do
      End If
    End Do
    ! close data file
    Close(i_mag)
    err_magg=.false.
  End Subroutine read_magnetic_data

  Subroutine read_magnetic_binary()
    ! read data about magnetic space groups
    ! input data from magnetic_table.dat
    integer :: i,j,k,n,m,i_mag
    !*****************************************************************************
    ! open data file
    Call Get_LogUnit(i_mag)
    Open(unit=i_mag,File='magnetic_data.bin',status="old",action="read",form="unformatted",access="stream")
    !For the old Lahey compiler use this
    !Open(unit=i_mag,File='magnetic_data.bin',status="old",action="read",form="unformatted",access="transparent") ! For Lahey
    ! read nonhexagonal point operators
    Do i=1,48
      Read(i_mag)n,point_op_label(i),point_op_xyz(i),  &
          ((point_op_matrix(k,j,i),j=1,3),k=1,3)
    End Do
    ! read hexagonal point operators
    Do i=1,24
      Read(i_mag) n,point_op_hex_label(i), point_op_hex_xyz(i),  &
          ((point_op_hex_matrix(k,j,i),j=1,3),k=1,3)
    End Do
    ! read data for each magnetic space group
    Do i=1,1651
      Read(i_mag) (nlabelparts_bns(j,i),j=1,2),nlabel_bns(i),  &
          spacegroup_label_bns(i),(nlabelparts_og(j,i),j=1,3),  &
          nlabel_og(i),spacegroup_label_og(i)
      Read(i_mag) magtype(i)
      If(magtype(i) == 4) Then
        Read(i_mag) ((bnsog_point_op(j,k,i),j=1,3),k=1,3),  &
            (bnsog_origin(j,i),j=1,3),bnsog_origin_denom(i)
      End If
      Read(i_mag) ops_count(i)
      Read(i_mag) (ops_bns_point_op(j,i),(ops_bns_trans(k,j,i),k=1,3),  &
          ops_bns_trans_denom(j,i),ops_bns_timeinv(j,i), j=1,ops_count(i))
      Read(i_mag) lattice_bns_vectors_count(i)
      Read(i_mag) ((lattice_bns_vectors(k,j,i),k=1,3),  &
          lattice_bns_vectors_denom(j,i), j=1,lattice_bns_vectors_count(i))
      Read(i_mag) wyckoff_site_count(i)
      Do j=1,wyckoff_site_count(i)
        Read(i_mag) wyckoff_pos_count(j,i),wyckoff_mult(j,i), wyckoff_label(j,i)
        Do k=1,wyckoff_pos_count(j,i)
          Read(i_mag) (wyckoff_bns_fract(m,k,j,i),m=1,3),  &
              wyckoff_bns_fract_denom(k,j,i),  &
              ((wyckoff_bns_xyz(m,n,k,j,i),m=1,3),n=1,3),  &
              ((wyckoff_bns_mag(m,n,k,j,i),m=1,3),n=1,3)
        End Do
      End Do
      If(magtype(i) == 4) Then
        Read(i_mag) ops_count(i)
        Read(i_mag) (ops_og_point_op(j,i),(ops_og_trans(k,j,i),k=1,3),  &
            ops_og_trans_denom(j,i),ops_og_timeinv(j,i), j=1,ops_count(i))
        Read(i_mag) lattice_og_vectors_count(i)
        Read(i_mag) ((lattice_og_vectors(k,j,i),k=1,3),  &
            lattice_og_vectors_denom(j,i), j=1,lattice_og_vectors_count(i))
        Read(i_mag) wyckoff_site_count(i)
        Do j=1,wyckoff_site_count(i)
          Read(i_mag) wyckoff_pos_count(j,i),wyckoff_mult(j,i), wyckoff_label(j,i)
          Do k=1,wyckoff_pos_count(j,i)
            Read(i_mag) (wyckoff_og_fract(m,k,j,i),m=1,3),  &
                  wyckoff_og_fract_denom(k,j,i),              &
                ((wyckoff_og_xyz(m,n,k,j,i),m=1,3),n=1,3),    &
                ((wyckoff_og_mag(m,n,k,j,i),m=1,3),n=1,3)
          End Do
        End Do
      End If
    End Do
    ! close data file
    Close(i_mag)
  End Subroutine read_magnetic_binary

  Subroutine write_magnetic_binary()
    ! read data about magnetic space groups
    ! input data from magnetic_table.dat
    integer :: i,j,k,n,m,i_mag
    !*****************************************************************************
    ! open data file
    Call Get_LogUnit(i_mag)
    Open(unit=i_mag,File='magnetic_data.bin',status="replace",action="write",access="stream",form="unformatted")
    !For the old Lahey compiler use this
    !Open(unit=i_mag,File='magnetic_data.bin',status="replace",action="write",access="transparent",form="unformatted")  !For Lahey
    ! read nonhexangonal point operators
    Do i=1,48
      Write(i_mag) i,point_op_label(i),point_op_xyz(i),  &
                  ((point_op_matrix(k,j,i),j=1,3),k=1,3)
    End Do
    ! read hexagonal point operators
    Do i=1,24
      Write(i_mag) i,point_op_hex_label(i), point_op_hex_xyz(i),  &
          ((point_op_hex_matrix(k,j,i),j=1,3),k=1,3)
    End Do
    ! read data for each magnetic space group
    Do i=1,1651
      Write(i_mag) (nlabelparts_bns(j,i),j=1,2),nlabel_bns(i),  &
          spacegroup_label_bns(i),(nlabelparts_og(j,i),j=1,3),  &
          nlabel_og(i),spacegroup_label_og(i)
      Write(i_mag) magtype(i)
      If(magtype(i) == 4) Then
        Write(i_mag) ((bnsog_point_op(j,k,i),j=1,3),k=1,3),  &
            (bnsog_origin(j,i),j=1,3),bnsog_origin_denom(i)
      End If
      Write(i_mag) ops_count(i)
      Write(i_mag) (ops_bns_point_op(j,i),(ops_bns_trans(k,j,i),k=1,3),  &
          ops_bns_trans_denom(j,i),ops_bns_timeinv(j,i), j=1,ops_count(i))
      Write(i_mag) lattice_bns_vectors_count(i)
      Write(i_mag) ((lattice_bns_vectors(k,j,i),k=1,3),  &
          lattice_bns_vectors_denom(j,i), j=1,lattice_bns_vectors_count(i))
      Write(i_mag) wyckoff_site_count(i)
      Do j=1,wyckoff_site_count(i)
        Write(i_mag) wyckoff_pos_count(j,i),wyckoff_mult(j,i), wyckoff_label(j,i)
        Do k=1,wyckoff_pos_count(j,i)
          Write(i_mag) (wyckoff_bns_fract(m,k,j,i),m=1,3),  &
              wyckoff_bns_fract_denom(k,j,i),  &
              ((wyckoff_bns_xyz(m,n,k,j,i),m=1,3),n=1,3),  &
              ((wyckoff_bns_mag(m,n,k,j,i),m=1,3),n=1,3)
        End Do
      End Do
      If(magtype(i) == 4) Then
        Write(i_mag) ops_count(i)
        Write(i_mag) (ops_og_point_op(j,i),(ops_og_trans(k,j,i),k=1,3),  &
            ops_og_trans_denom(j,i),ops_og_timeinv(j,i), j=1,ops_count(i))
        Write(i_mag) lattice_og_vectors_count(i)
        Write(i_mag) ((lattice_og_vectors(k,j,i),k=1,3),  &
            lattice_og_vectors_denom(j,i), j=1,lattice_og_vectors_count(i))
        Write(i_mag) wyckoff_site_count(i)
        Do j=1,wyckoff_site_count(i)
          Write(i_mag) wyckoff_pos_count(j,i),wyckoff_mult(j,i), wyckoff_label(j,i)
          Do k=1,wyckoff_pos_count(j,i)
            Write(i_mag) (wyckoff_og_fract(m,k,j,i),m=1,3),   &
                  wyckoff_og_fract_denom(k,j,i),              &
                ((wyckoff_og_xyz(m,n,k,j,i),m=1,3),n=1,3),    &
                ((wyckoff_og_mag(m,n,k,j,i),m=1,3),n=1,3)
          End Do
        End Do
      End If
    End Do
    ! close data file
    Close(i_mag)
  End Subroutine write_magnetic_binary

  Subroutine write_magnetic_data(num,lun,mat,orig,dir,OG_BNS)
    integer,                          intent(in) :: num
    integer,optional,                 intent(in) :: lun  !Logical unit to write
    real,   optional, dimension(3,3), intent(in) :: mat
    real,   optional, dimension(3),   intent(in) :: orig
    logical,optional,                 intent(in) :: dir
    character(len=*),optional,        intent(in) :: OG_BNS

    real, dimension (3,3), parameter :: e = reshape ((/1.0,0.0,0.0,  &
                                                       0.0,1.0,0.0,  &
                                                       0.0,0.0,1.0/),(/3,3/))
    Character(len=34)      :: OG_Symb, BNS_Symb
    Character(len=30)      :: abc_op, Strsym, abc_symb
    Character(len=35)      :: ShOp_symb
    Character(len=3)       :: symb_OGBNS
    integer                :: idem,inv_time, ipr,j,k, ifail
    real,   dimension(3)   :: tr
    integer,dimension(3,3) :: op,rot,trf
    real,dimension(3,3)    :: S,Sinv !,Sog,Soginv

    logical                :: change_setting,tr_dir

    change_setting=.false.
    if(present(mat) .and. present(orig)) then
      change_setting=.true.
      trf=nint(mat)
      S=transpose(Mat)
      call matrix_inverse(S,Sinv,ifail)
      if (ifail /= 0) then
         write(unit=*,fmt="(a)") " => Inversion Matrix Failed in: write_magnetic_data!"
         return
      end if
      op=bnsog_point_op(:,:,num) !Transformation from OG to BNS
      tr=real(bnsog_origin(:,num))/real(bnsog_origin_denom(num))
    end if

    tr_dir=.true.
    if(present(dir)) tr_dir=dir
    symb_OGBNS="BNS"
    if(present(OG_BNS)) then
      if(magtype(num) == 4) then
         symb_OGBNS=OG_BNS
      end if
    end if
    ipr=6
    if(present(lun)) ipr=lun
    BNS_Symb="BNS:"//nlabel_bns(num)//" "//trim(spacegroup_label_bns(num))
    OG_Symb= " OG:"//nlabel_og(num)//" "//trim(spacegroup_label_og(num))
    Write(unit=ipr,fmt="(a,/)")  " "
    Write(unit=ipr,fmt="(a,i4)") " Group Ordering Number (BNS): ",num
    Write(unit=ipr,fmt="(a,i4)") "         Magnetic Group type: ",magtype(num)
    Write(unit=ipr,fmt="(a)")    "      Magnetic Group Symbols: "//BNS_Symb//OG_Symb

    if(.not. change_setting) then

       If(magtype(num) == 4) Then  !OG-BNS transformation:
         op=bnsog_point_op(:,:,num)
         tr=real(bnsog_origin(:,num))/real(bnsog_origin_denom(num))
         call Get_Trasfm_Symbol(Op,Tr,abc_op)
         Write(unit=ipr,fmt="(a)") "       OG-BNS transformation: "//trim(abc_op)
       End If
       If(magtype(num) == 4) Then
         Write(unit=ipr,fmt="(/,a)")   "   ----------------------------------- "
         Write(unit=ipr,fmt="(a)")     "   Belov-Neronova-Smirnova Description "
         Write(unit=ipr,fmt="(a,/)")   "   ----------------------------------- "
         Write(unit=ipr,fmt="(a,i3)")  "   Number of BNS operators: ",ops_count(num)
       else
         Write(unit=ipr,fmt="(/,a)")   "   ------------------------- "
         Write(unit=ipr,fmt="(a)")     "   Common BNS-OG Description "
         Write(unit=ipr,fmt="(a,/)")   "   ------------------------- "
         Write(unit=ipr,fmt="(a,i3)")  "   Number of operators: ",ops_count(num)
       end if

       do j=1,ops_count(num)
          idem=ops_bns_trans_denom(j,num)
          inv_time=ops_bns_timeinv(j,num)
          tr= real(ops_bns_trans(:,j,num))/real(idem)
          Call Get_Seitz_symbol(ops_bns_point_op(j,num),inv_time,Tr,Strsym)
          Write(unit=ipr,fmt="(a,i3,a,a25)",advance="no")  "   #",j,": ",trim(Strsym)
          if(mod(j,4) == 0) write(unit=ipr,fmt="(a)") "  "
       end do
       write(unit=ipr,fmt="(a)") "  "
       If(magtype(num) == 4) Then
         Write(unit=ipr,fmt="(a,i3)") "  Number of BNS lattice translations: ",lattice_bns_vectors_count(num)
       else
         Write(unit=ipr,fmt="(a,i3)") "  Number of lattice translations: ",lattice_bns_vectors_count(num)
       end if
       do j=1,lattice_bns_vectors_count(num)
          tr= real(lattice_bns_vectors(:,j,num))/real(lattice_bns_vectors_denom(j,num))
          Call Get_Transl_Symbol(tr,Strsym)
          if(present(lun)) then
            Write(unit=ipr,fmt="(a,i3,a,a20,tr4,3f8.4)")  "   #",j,": ",trim(Strsym), tr
          else
            Write(unit=ipr,fmt="(a,i3,a,a20)",advance="no")  "   #",j,": ",trim(Strsym)
            if(mod(j,4) == 0) write(unit=ipr,fmt="(a)") "  "
          end if
       end do

    else  !Change setting

       if(symb_OGBNS == "BNS") then
          If(magtype(num) == 4) Then
             Write(unit=ipr,fmt="(/,a)")   "   ----------------------------------------------------------- "
             Write(unit=ipr,fmt="(a)")     "   Belov-Neronova-Smirnova Description in non-standard setting "
             Write(unit=ipr,fmt="(a,/)")   "   ----------------------------------------------------------- "
          else
             Write(unit=ipr,fmt="(/,a)")   "   ----------------------------------------------------------- "
             Write(unit=ipr,fmt="(a)")     "   Common BNS/OG Description in non-standard setting "
             Write(unit=ipr,fmt="(a,/)")   "   ----------------------------------------------------------- "
          end if
          if(tr_dir) then
            Call Get_Trasfm_Symbol(trf,orig,abc_symb,.true.)
          else
            trf=transpose(nint(Sinv))
            Call Get_Trasfm_Symbol(trf,-orig,abc_symb,.true.)
          end if
          Write(unit=ipr,fmt="(a)")     "   The operators are transformed from the standard setting to: "//trim(abc_symb)
       end if

    end if

    If(symb_OGBNS == "BNS") Then
       write(unit=ipr,fmt="(a)") "  "
       If(magtype(num) == 4) Then
         Write(unit=ipr,fmt="(a,i3)") "  Number of BNS Wyckoff positions: ",wyckoff_site_count(num)
       else
         Write(unit=ipr,fmt="(a,i3)") "  Number of Wyckoff positions: ",wyckoff_site_count(num)
       end if
       Do j=1,wyckoff_site_count(num)
         write(unit=ipr,fmt="(a)") "  "
         Write(unit=ipr,fmt="(a,i3,a,2i4,a)") "   Wyckoff site #",j,": ",wyckoff_pos_count(j,num),wyckoff_mult(j,num),&
                                             wyckoff_label(j,num)
         Do k=1,wyckoff_pos_count(j,num)
           idem=wyckoff_bns_fract_denom(k,j,num)
           tr=real(wyckoff_bns_fract(:,k,j,num))/real(idem)
           op = wyckoff_bns_xyz(:,:,k,j,num)
           Rot = wyckoff_bns_mag(:,:,k,j,num)
           if(change_setting) then
             if(tr_dir) then     !Here the modulo_lat function is applied because in BNS setting we use the magnetic cell
               op=nint(matmul(matmul(Sinv,op),S))
               Rot=nint(matmul(matmul(Sinv,Rot),S))
               tr=Modulo_Lat(matmul(Sinv,tr-matmul(e-op,orig)))
             else
               op=nint(matmul(matmul(S,op),Sinv))
               Rot=nint(matmul(matmul(S,Rot),Sinv))
               tr=Modulo_Lat(matmul(S,tr+matmul(e-op,orig)))
             end if
           end if
           if(mcif) then
              Call Get_Shubnikov_Operator_Symbol(Op,Rot,tr,ShOp_symb,mcif)
           else
              Call Get_Shubnikov_Operator_Symbol(Op,Rot,tr,ShOp_symb)
           end if
           if(present(lun)) then
              write(unit=ipr,fmt="(a,i3,a)")  "   #",k,": "//trim(ShOp_symb)
           else
              write(unit=ipr,fmt="(a,i3,a,a35)",advance="no")  "   #",k,": ",trim(ShOp_symb)
              if(mod(k,3) == 0) write(unit=ipr,fmt="(a)") "  "
           end if
         End Do
       End Do

    End If  ! If(symb_OGBNS="BNS") Then

    If(magtype(num) == 4) Then

      if(.not. change_setting) then
          write(unit=ipr,fmt="(a)") "  "
          Write(unit=ipr,fmt="(/,a)")   "   ------------------------------- "
          Write(unit=ipr,fmt="(a)")     "   Opechowski-Guccione Description "
          Write(unit=ipr,fmt="(a,/)")   "   ------------------------------- "
          Write(unit=ipr,fmt="(a,i3)") "   Number of OG operators: ",ops_count(num)
          do j=1,ops_count(num)
             idem=ops_og_trans_denom(j,num)
             inv_time=ops_og_timeinv(j,num)
             tr= real(ops_og_trans(:,j,num))/real(idem)
             Call Get_Seitz_symbol(ops_og_point_op(j,num),inv_time,Tr,Strsym)
             Write(unit=ipr,fmt="(a,i3,a,a25)",advance="no")  "   #",j,": ",trim(Strsym)
             if(mod(j,4) == 0) write(unit=ipr,fmt="(a)") "  "
          end do
          write(unit=ipr,fmt="(a)") "  "
          Write(unit=ipr,fmt="(a,i3)") "  Number of OG lattice translations: ",lattice_og_vectors_count(num)
          do j=1,lattice_og_vectors_count(num)
             tr= real(lattice_og_vectors(:,j,num))/real(lattice_og_vectors_denom(j,num))
             Call Get_Transl_Symbol(tr,Strsym)
             if(present(lun)) then
               Write(unit=ipr,fmt="(a,i3,a,a20,tr4,3f8.4)")  "   #",j,": ",trim(Strsym), tr
             else
               Write(unit=ipr,fmt="(a,i3,a,a20)",advance="no")  "   #",j,": ",trim(Strsym)
               if(mod(j,4) == 0) write(unit=ipr,fmt="(a)") "  "
             end if
          end do
      else
         If(symb_OGBNS == "OG") Then
           Write(unit=ipr,fmt="(/,a)")   "   --------------------------------------------------------- "
           Write(unit=ipr,fmt="(a)")     "   Opechowski-Guccione Description in a non-standard setting "
           Write(unit=ipr,fmt="(a,/)")   "   --------------------------------------------------------- "
           Write(unit=ipr,fmt="(a)")     "   The operators are transformed from the standard setting to: "//trim(abc_symb)
         End If
      end if

      If(symb_OGBNS == "OG") Then
         write(unit=ipr,fmt="(a)") "  "
         Write(unit=ipr,fmt="(a,i3)") "  Number of OG Wyckoff positions: ",wyckoff_site_count(num)
         Do j=1,wyckoff_site_count(num)
           write(unit=ipr,fmt="(a)") "  "
           Write(unit=ipr,fmt="(a,i3,a,2i4,a)") "   Wyckoff site #",j,": ",wyckoff_pos_count(j,num),wyckoff_mult(j,num),&
                                                  wyckoff_label(j,num)
           Do k=1,wyckoff_pos_count(j,num)
             idem=wyckoff_og_fract_denom(k,j,num)
             tr=real(wyckoff_og_fract(:,k,j,num))/real(idem)
             op = wyckoff_bns_xyz(:,:,k,j,num)
             Rot = wyckoff_bns_mag(:,:,k,j,num)
             if(change_setting) then
               if(tr_dir) then  !Here we cannot use modulo_lat function
                 op=nint(matmul(matmul(Sinv,op),S))
                 Rot=nint(matmul(matmul(Sinv,Rot),S))
                 tr=matmul(Sinv,tr-matmul(e-op,orig))
               else
                 op=nint(matmul(matmul(S,op),Sinv))
                 Rot=nint(matmul(matmul(S,Rot),Sinv))
                 tr=matmul(S,tr+matmul(e-op,orig))
               end if
             end if
             if(mcif) then
                Call Get_Shubnikov_Operator_Symbol(Op,Rot,tr,ShOp_symb,mcif)
             else
                Call Get_Shubnikov_Operator_Symbol(Op,Rot,tr,ShOp_symb)
             end if
             if(present(lun)) then
                write(unit=ipr,fmt="(a,i3,a,a35)")  "   #",k,": ",trim(ShOp_symb)
             else
                write(unit=ipr,fmt="(a,i3,a,a35)",advance="no")  "   #",k,": ",trim(ShOp_symb)
                if(mod(k,3) == 0) write(unit=ipr,fmt="(a)") "  "
             end if
           End Do
         End Do
      end if
      write(unit=ipr,fmt="(a)") "  "
    End If !If(magtype(num) == 4)

    return
  End Subroutine write_magnetic_data

  Subroutine write_magnetic_data_binary(num,ipr)
    integer, intent(in)  :: num
    integer, intent(in)  :: ipr

    Character(len=34)      :: OG_Symb, BNS_Symb
    Character(len=30)      :: abc_op, Strsym
    Character(len=35)      :: ShOp_symb
    integer                :: idem,inv_time,j,k
    real,   dimension(3)   :: tr
    integer,dimension(3,3) :: op,rot

    BNS_Symb="BNS:"//nlabel_bns(num)//" "//trim(spacegroup_label_bns(num))
    OG_Symb= " OG:"//nlabel_og(num)//" "//trim(spacegroup_label_og(num))
    Write(unit=ipr) num          ! Group Ordering Number (BNS)
    Write(unit=ipr) magtype(num) ! Magnetic Group type
    Write(unit=ipr) BNS_Symb     ! Magnetic Group BNS Symbol
    Write(unit=ipr) OG_Symb      ! Magnetic Group OG Symbol

    If(magtype(num) == 4) Then  !OG-BNS transformation:
      op=bnsog_point_op(:,:,num)
      tr=real(bnsog_origin(:,num))/real(bnsog_origin_denom(num))
      call Get_Trasfm_Symbol(Op,Tr,abc_op)
      Write(unit=ipr) trim(abc_op) !OG-BNS transformation
    End If
    Write(unit=ipr)  ops_count(num)  ! Number of operators
    do j=1,ops_count(num)
       idem=ops_bns_trans_denom(j,num)
       inv_time=ops_bns_timeinv(j,num)
       tr= real(ops_bns_trans(:,j,num))/real(idem)
       Call Get_Seitz_symbol(ops_bns_point_op(j,num),inv_time,Tr,Strsym)
       Write(unit=ipr) Strsym
    end do
    Write(unit=ipr) lattice_bns_vectors_count(num) ! Number of lattice translations
    do j=1,lattice_bns_vectors_count(num)
       tr= real(lattice_bns_vectors(:,j,num))/real(lattice_bns_vectors_denom(j,num))
       Call Get_Transl_Symbol(tr,Strsym)
       Write(unit=ipr) Strsym
    end do
    Write(unit=ipr) wyckoff_site_count(num) !Number of Wyckoff positions
    Do j=1,wyckoff_site_count(num)
      Write(unit=ipr) wyckoff_pos_count(j,num),wyckoff_mult(j,num),&
                      wyckoff_label(j,num)
      Do k=1,wyckoff_pos_count(j,num)
        idem=wyckoff_bns_fract_denom(k,j,num)
        tr=real(wyckoff_bns_fract(:,k,j,num))/real(idem)
        op = wyckoff_bns_xyz(:,:,k,j,num)
        Rot = wyckoff_bns_mag(:,:,k,j,num)
        Call Get_Shubnikov_Operator_Symbol(Op,Rot,tr,ShOp_symb)
        write(unit=ipr) ShOp_symb
      End Do
    End Do

    If(magtype(num) == 4) Then
      Write(unit=ipr) ops_count(num) ! Number of OG operators
      do j=1,ops_count(num)
         idem=ops_og_trans_denom(j,num)
         inv_time=ops_og_timeinv(j,num)
         tr= real(ops_og_trans(:,j,num))/real(idem)
         Call Get_Seitz_symbol(ops_og_point_op(j,num),inv_time,Tr,Strsym)
         Write(unit=ipr)  Strsym
      end do
      Write(unit=ipr) lattice_og_vectors_count(num) ! Number of OG lattice translations
      do j=1,lattice_og_vectors_count(num)
         tr= real(lattice_og_vectors(:,j,num))/real(lattice_og_vectors_denom(j,num))
         Call Get_Transl_Symbol(tr,Strsym)
         Write(unit=ipr) Strsym
      end do

      Write(unit=ipr) wyckoff_site_count(num) ! Number of OG Wyckoff positions
      Do j=1,wyckoff_site_count(num)
        Write(unit=ipr) wyckoff_pos_count(j,num),wyckoff_mult(j,num),&
                                            wyckoff_label(j,num)
        Do k=1,wyckoff_pos_count(j,num)
          idem=wyckoff_og_fract_denom(k,j,num)
          tr=real(wyckoff_og_fract(:,k,j,num))/real(idem)
          op = wyckoff_bns_xyz(:,:,k,j,num)
          Rot = wyckoff_bns_mag(:,:,k,j,num)
          Call Get_Shubnikov_Operator_Symbol(op,Rot,tr,ShOp_symb)
          write(unit=ipr) ShOp_symb
        End Do
      End Do
    End If
    return
  End Subroutine write_magnetic_data_binary

  Subroutine write_magnetic_data_ASCII(num,ipr)
    integer, intent(in)  :: num
    integer, intent(in)  :: ipr

    Character(len=34)      :: OG_Symb, BNS_Symb
    Character(len=30)      :: abc_op,Strsym
    Character(len=35)      :: ShOp_symb
    integer                :: idem,inv_time,j,k
    real,   dimension(3)   :: tr
    integer,dimension(3,3) :: op,rot

    BNS_Symb="BNS:"//trim(nlabel_bns(num))//" "//trim(spacegroup_label_bns(num))
    OG_Symb= " OG:"//trim(nlabel_og(num))//" "//trim(spacegroup_label_og(num))
    Write(unit=ipr,fmt="(a,i5,a)") "Group Number:",num, " "//trim(BNS_Symb)//" "//trim(OG_Symb)          ! Group Ordering Number (BNS)
    Write(unit=ipr,fmt="(a,i5)")   "Group Type:",magtype(num) ! Magnetic Group type

    If(magtype(num) == 4) Then  !OG-BNS transformation:
      op=bnsog_point_op(:,:,num)
      tr=real(bnsog_origin(:,num))/real(bnsog_origin_denom(num))
      call Get_Trasfm_Symbol(Op,Tr,abc_op)
      Write(unit=ipr,fmt="(a)") trim(abc_op) !OG-BNS transformation
    End If
    Write(unit=ipr,fmt="(a,i4)")"# Operators:",  ops_count(num)  ! Number of operators
    do j=1,ops_count(num)
       idem=ops_bns_trans_denom(j,num)
       inv_time=ops_bns_timeinv(j,num)
       tr= real(ops_bns_trans(:,j,num))/real(idem)
       Call Get_Seitz_symbol(ops_bns_point_op(j,num),inv_time,Tr,Strsym)
       Write(unit=ipr,fmt="(a)") trim(Strsym)
    end do
    Write(unit=ipr,fmt="(a,i4)")"# Lattice vectors:", lattice_bns_vectors_count(num) ! Number of lattice translations
    do j=1,lattice_bns_vectors_count(num)
       tr= real(lattice_bns_vectors(:,j,num))/real(lattice_bns_vectors_denom(j,num))
       Call Get_Transl_Symbol(tr,Strsym)
       Write(unit=ipr,fmt="(a)") trim(Strsym)
    end do
    Write(unit=ipr,fmt="(a,i4)")"# Wyckoff sites:", wyckoff_site_count(num) !Number of Wyckoff positions
    Do j=1,wyckoff_site_count(num)
      Write(unit=ipr,fmt="(a,i3,a,2i4,a)") "Wyckoff site #",j,": ",wyckoff_pos_count(j,num),wyckoff_mult(j,num),&
                                            wyckoff_label(j,num)
      Do k=1,wyckoff_pos_count(j,num)
        idem=wyckoff_bns_fract_denom(k,j,num)
        tr=real(wyckoff_bns_fract(:,k,j,num))/real(idem)
        op = wyckoff_bns_xyz(:,:,k,j,num)
        Rot = wyckoff_bns_mag(:,:,k,j,num)
        Call Get_Shubnikov_Operator_Symbol(Op,Rot,tr,ShOp_symb)
        write(unit=ipr,fmt="(a)") trim(ShOp_symb)
      End Do
    End Do

    If(magtype(num) == 4) Then
      Write(unit=ipr,fmt="(a,i4)")"# Operators (OG) :",  ops_count(num) ! Number of OG operators
      do j=1,ops_count(num)
         idem=ops_og_trans_denom(j,num)
         inv_time=ops_og_timeinv(j,num)
         tr= real(ops_og_trans(:,j,num))/real(idem)
         Call Get_Seitz_symbol(ops_og_point_op(j,num),inv_time,Tr,Strsym)
         Write(unit=ipr,fmt="(a)") trim(Strsym)
      end do
      Write(unit=ipr,fmt="(a,i4)")"# Lattice vectors (OG):", lattice_og_vectors_count(num) ! Number of OG lattice translations
      do j=1,lattice_og_vectors_count(num)
         tr= real(lattice_og_vectors(:,j,num))/real(lattice_og_vectors_denom(j,num))
         Call Get_Transl_Symbol(tr,Strsym)
         Write(unit=ipr,fmt="(a)") trim(Strsym)
      end do

      Write(unit=ipr,fmt="(a,i4)")"# Wyckoff sites (OG):", wyckoff_site_count(num) ! Number of OG Wyckoff positions
      Do j=1,wyckoff_site_count(num)
        Write(unit=ipr,fmt="(a,i3,a,2i4,a)") "Wyckoff site #",j,": ",wyckoff_pos_count(j,num),wyckoff_mult(j,num),&
                                            wyckoff_label(j,num)
        Do k=1,wyckoff_pos_count(j,num)
          idem=wyckoff_og_fract_denom(k,j,num)
          tr=real(wyckoff_og_fract(:,k,j,num))/real(idem)
          op = wyckoff_bns_xyz(:,:,k,j,num)
          Rot = wyckoff_bns_mag(:,:,k,j,num)
          Call Get_Shubnikov_Operator_Symbol(op,Rot,tr,ShOp_symb)
          write(unit=ipr,fmt="(a)") trim(ShOp_symb)
        End Do
      End Do
    End If
    return
  End Subroutine write_magnetic_data_ASCII

 End Module CFML_Magnetic_Groups