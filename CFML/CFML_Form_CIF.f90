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
!!----
!!---- MODULE: CFML_IO_Formats
!!----   INFO: Creation/Conversion for several formats
!!----
!!---- HISTORY
!!----    Update: 07/03/2011
!!----
!!----
!!---- DEPENDENCIES
!!----
!!---- VARIABLES
!!----    ERR_FORM
!!----    ERR_FORM_MESS
!!----    INTERVAL_TYPE
!!----    JOB_INFO_TYPE
!!----
!!---- PROCEDURES
!!----    Functions:
!!----
!!----    Subroutines:
!!----       FILE_TO_FILELIST
!!----       GET_JOB_INFO
!!----       GET_PHASES_FILE
!!--++       GET_NPHASES_CIFFILE
!!--++       GET_NPHASES_PCRFILE
!!----       INIT_ERR_FORM
!!----       READ_ATOM
!!----       READ_CELL
!!----       READ_CIF_ATOM
!!----       READ_CIF_CELL
!!----       READ_CIF_CHEMICALNAME
!!----       READ_CIF_CONT
!!----       READ_CIF_HALL
!!----       READ_CIF_HM
!!----       READ_CIF_LAMBDA
!!----       READ_CIF_PRESSURE
!!----       READ_CIF_TEMP
!!----       READ_CIF_SYMM
!!----       READ_CIF_TITLE
!!----       READ_CIF_Z
!!----       READ_FILE_ATOM
!!--++       READ_FILE_ATOMLIST             [Overloaded]
!!--++       READ_FILE_POINTLIST            [Overloaded]
!!----       READ_FILE_CELL
!!--++       READ_FILE_CELLc                [Overloaded]
!!--++       READ_FILE_CELLt                [Overloaded]
!!----       READ_FILE_LAMBDA
!!----       READ_FILE_RNGSINTL
!!----       READ_FILE_SPG
!!----       READ_FILE_TRANSF
!!----       READ_SHX_ATOM
!!----       READ_SHX_CELL
!!----       READ_SHX_CONT
!!----       READ_SHX_FVAR
!!----       READ_SHX_LATT
!!----       READ_SHX_SYMM
!!----       READ_SHX_TITL
!!----       READ_UVALS
!!----       READN_SET_MAGNETIC_SPACE_GROUP
!!--++       READN_SET_XTAL_CFL             [Private]
!!--++       READN_SET_XTAL_CFL_MOLEC       [Private]
!!--++       READN_SET_XTAL_CFL_SHUB        [Private]
!!--++       READN_SET_XTAL_CIF             [Private]
!!--++       READN_SET_XTAL_PCR             [Private]
!!--++       READN_SET_XTAL_SHX             [Private]
!!----       READN_SET_XTAL_STRUCTURE
!!--++       READN_SET_XTAL_STRUCTURE_MOLCR [Overloaded]
!!--++       READN_SET_XTAL_STRUCTURE_SPLIT [Overloaded]
!!----       SET_MAGNETIC_SPACE_GROUP
!!----       WRITE_CIF_POWDER_PROFILE
!!----       WRITE_CIF_TEMPLATE
!!----       WRITE_MCIF
!!----       WRITE_SHX_TEMPLATE
!!----
!!
 Module CFML_IO_Formats

    !---- Use modules ----!
    Use CFML_GlobalDeps,                only: cp,sp,dp,pi,eps,Write_Date_Time
    Use CFML_Math_General,              only: sind,equal_matrix,Sort
    Use CFML_Math_3D,                   only: determ_a
    Use CFML_String_Utilities
    Use CFML_Crystal_Metrics,           only: Crystal_Cell_Type, Set_Crystal_Cell, Convert_U_Betas, &
                                              Convert_B_Betas, U_Equiv, Convert_Betas_U
    Use CFML_Crystallographic_Symmetry, only: Space_Group_Type, Magnetic_Space_Group_Type,Set_SpaceGroup,     &
                                              Init_Magnetic_Space_Group_Type,Get_Multip_Pos,Get_MagMatSymb,   &
                                              Read_Xsym,Read_Msymm, Setting_Change, get_symsymb,Sym_Oper_type,&
                                              Msym_Oper_Type,is_Lattice_vec,Get_Stabilizer,Get_mOrbit
    Use CFML_Atom_TypeDef,              only: Atom_Type, Init_Atom_Type,atom_list_type,         &
                                              Allocate_atom_list, Deallocate_atom_list
    Use CFML_Molecular_Crystals,        only: Err_Molec, Err_Molec_Mess,Molecular_Crystal_Type, &
                                              Read_Molecule, Set_Euler_Matrix, Write_Molecule
    Use CFML_Geometry_Calc,             only: Point_List_Type, Get_Euler_from_Fract
    Use CFML_Diffraction_Patterns,      only: Diffraction_Pattern_type
    Use CFML_Scattering_Chemical_Tables,only: Set_Magnetic_Form, Remove_Magnetic_Form, num_mag_form, &
                                              Magnetic_Form, get_magnetic_form_factor
    Use CFML_Magnetic_Groups
    Use CFML_EisPack,                   only: rg_ort

    !---- Variables ----!
    implicit none

    private


    !---- List of public functions ----!

    !---- List of public subroutines ----!
    public :: Init_Err_Form, Read_Atom, Read_Cell, Read_Cif_Atom, Read_Cif_Cell,                 &
              Read_Cif_Cont, Read_Cif_Hall, Read_Cif_Hm, Read_Cif_Lambda, Read_Cif_Symm,         &
              Read_Cif_Title, Read_Cif_Z, Read_File_Atom, Read_File_Spg, Read_Cif_ChemicalName,  &
              Read_File_Transf, Read_Shx_Atom, Read_Shx_Cell, Read_Shx_Cont, Read_Shx_Fvar,      &
              Read_Shx_Latt, Read_Shx_Symm, Read_Shx_Titl, Read_Uvals, Write_Cif_Powder_Profile, &
              Write_Cif_Template, Write_Shx_Template, Read_File_rngSINTL, Read_File_Lambda,      &
              Get_job_info, File_To_FileList, Get_Phases_File, Read_Cif_Pressure,                &
              Read_Cif_Temp,Readn_Set_Magnetic_Space_Group, Set_Magnetic_Space_Group,            &
              Cleanup_Symmetry_Operators,Write_mCIF, Get_Refinement_Codes, Get_moment_ctr,       &
              Readn_Set_Magnetic_Structure_MCIF

    !---- List of public overloaded procedures: subroutines ----!
    public :: Read_File_Cell, Readn_Set_Xtal_Structure, Write_Atoms_CFL, Write_CFL

    !---- List of private functions ----!

    !---- List of private subroutines ----!
    private:: Read_File_Cellc, Read_File_Cellt, Read_File_Atomlist,Read_File_Pointlist,               &
              Readn_Set_Xtal_CFL, Readn_Set_Xtal_CIF, Readn_Set_Xtal_PCR,Readn_Set_Xtal_SHX,          &
              Readn_Set_Xtal_CFL_Molec, Readn_Set_Xtal_Structure_Split,                               &
              Readn_Set_Xtal_Structure_Molcr, Get_NPhases_CIFFile,Get_NPHases_PCRFile,                &
              Write_CFL_Molcrys, Write_CFL_Atom_List_Type, Write_Atoms_CFL_ATM, Write_Atoms_CFL_MOLX, &
              Write_Atoms_CFL_MOLX_orig, Readn_Set_XTal_CFL_Shub

    !---- Definitions ----!


    character (len=1), dimension(26),parameter, private   :: &
    cdd=(/'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r', &
        's','t','u','v','w','x','y','z'/)
    real(kind=dp), parameter, private :: epps=0.000001_dp
    !!----
    !!---- ERR_FORM
    !!----    logical, public :: err_form
    !!----
    !!----    Logical Variable indicating an error in CFML_IO_Formats
    !!----
    !!---- Update: February - 2005
    !!
    logical, public :: err_form

    !!----
    !!---- ERR_FORM_MESS
    !!----    character(len=150), public :: ERR_Form_Mess
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: February - 2005
    !!
    character(len=150),       public  :: ERR_Form_Mess

    !!----
    !!---- EPSV
    !!----    real(kind=cp), parameter, private :: epsv=1.0e-5_cp
    !!----
    !!----    Private small real number for floating point comparisons
    !!----
    !!---- Update: February - 2011
    !!
    real(kind=cp), parameter, private :: epsv=1.0e-5_cp

    !!----
    !!---- TYPE :: INTERVAL_TYPE
    !!--..
    !!---- Type, public :: interval_type
    !!----    real(kind=cp) :: mina  !low limit
    !!----    real(kind=cp) :: maxb  !high limit
    !!---- End Type interval_type
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: interval_type
       real(kind=cp) :: mina  !low limit
       real(kind=cp) :: maxb  !high limit
    End Type interval_type

    !!----
    !!---- TYPE :: JOB_INFO_TYPE
    !!--..
    !!---- Type, public :: Job_Info_type
    !!----    character(len=120)                            :: Title          ! Title
    !!----    integer                                       :: Num_Phases     ! Number of phases
    !!----    integer                                       :: Num_Patterns   ! Number of patterns
    !!----    integer                                       :: Num_cmd        ! Number of command lines
    !!----    character(len=16),  dimension(:), allocatable :: Patt_typ       ! Type of Pattern
    !!----    character(len=128), dimension(:), allocatable :: Phas_nam       ! Name of phases
    !!----    character(len=128), dimension(:), allocatable :: cmd            ! Command lines: text for actions
    !!----    type(interval_type),dimension(:), allocatable :: range_stl      ! Range in sinTheta/Lambda
    !!----    type(interval_type),dimension(:), allocatable :: range_q        ! Range in 4pi*sinTheta/Lambda
    !!----    type(interval_type),dimension(:), allocatable :: range_d        ! Range in d-spacing
    !!----    type(interval_type),dimension(:), allocatable :: range_2theta   ! Range in 2theta-spacing
    !!----    type(interval_type),dimension(:), allocatable :: range_Energy   ! Range in Energy
    !!----    type(interval_type),dimension(:), allocatable :: range_tof      ! Range in Time of Flight
    !!----    type(interval_type),dimension(:), allocatable :: Lambda         ! Lambda
    !!----    real(kind=cp)      ,dimension(:), allocatable :: ratio          ! ratio lambda2/lambda1
    !!----    real(kind=cp)      ,dimension(:), allocatable :: dtt1,dtt2      ! d-to-TOF coefficients
    !!---- End Type Job_Info_type
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: Job_Info_type
       character(len=120)                            :: Title
       integer                                       :: Num_Phases
       integer                                       :: Num_Patterns
       integer                                       :: Num_cmd
       character(len=16),  dimension(:), allocatable :: Patt_typ
       character(len=128), dimension(:), allocatable :: Phas_nam
       character(len=128), dimension(:), allocatable :: cmd
       type(interval_type),dimension(:), allocatable :: range_stl
       type(interval_type),dimension(:), allocatable :: range_q
       type(interval_type),dimension(:), allocatable :: range_d
       type(interval_type),dimension(:), allocatable :: range_2theta
       type(interval_type),dimension(:), allocatable :: range_Energy
       type(interval_type),dimension(:), allocatable :: range_tof
       type(interval_type),dimension(:), allocatable :: Lambda
       real(kind=cp)      ,dimension(:), allocatable :: ratio
       real(kind=cp)      ,dimension(:), allocatable :: dtt1,dtt2
    End Type Job_Info_type

    !!----
    !!---- TYPE :: FILE_LIST_TYPE
    !!--..
    !!---- Type,public :: File_List_Type
    !!----    integer                                       :: nlines ! Number of lines in the file
    !!----    character(len=256), allocatable, dimension(:) :: line   ! Content of the lines
    !!---- End Type file_list_type
    !!----
    !!---- Updated: February - 2005, November 2012
    !!
    Type,public :: File_List_Type
       integer                                       :: nlines
       character(len=256), allocatable, dimension(:) :: line
    End Type File_List_Type


    !---- Interfaces - Overloaded procedures--!
    Interface  Read_File_Cell
       Module Procedure Read_File_Cellc  !Last Output Argument Vector Of Six Component With The Cell Parameters
       Module Procedure Read_File_Cellt  !Last output argument object of type Crystal_cell_type
    End interface

    Interface Read_File_Atom
       Module Procedure Read_File_Atomlist   !Last Output Argument of type Atom_list_type
       Module Procedure Read_File_Pointlist  !Last output argument of type Point_list_type
    End Interface

    Interface Readn_Set_Xtal_Structure
       Module Procedure Readn_Set_Xtal_Structure_Molcr ! For Molecular Crystal Type
       Module Procedure Readn_Set_Xtal_Structure_Split ! For Cell, Spg, A types
       Module Procedure Readn_Set_Xtal_Structure_Magn  ! Use Shubnikov groups
    End Interface

    Interface Write_CFL
       Module Procedure Write_CFL_Molcrys        ! For Molecular Crystal Type
       Module Procedure Write_CFL_Atom_List_Type ! For Cell, Spg, A Types
    End Interface

    Interface Write_Atoms_CFL
       Module Procedure Write_Atoms_CFL_MOLX ! For Molecular Crystal Type
       Module Procedure Write_Atoms_CFL_ATM  ! For Cell, Spg, A Types
    End Interface

 Contains

    !---- Functions ----!

    !---- Subroutines ----!

    !!---- Subroutine Cleanup_Symmetry_Operators(MSpG)
    !!----   Type(Magnetic_Space_Group_Type), intent(in out) :: MSpG
    !!----
    !!----  Subroutine to re-organize symmetry operators extracting lattice translations
    !!----  and anti-translations and reordering the whole set of operators.
    !!----  (Still under development). It is supposed that the identity symmetry operator
    !!----  is provided in the input MSpG object, otherwise ok=.false. and
    !!----  no re-order is done.
    !!----
    !!----  Created: February 2014 (JRC), July 2016 (JRC), January 2020 (moved from CFML_Magsymm.f90)
    !!----
    Subroutine Cleanup_Symmetry_Operators(MSpG)
      Type(Magnetic_Space_Group_Type), intent(in out) :: MSpG
      !--- Local variables ---!
      integer,      dimension(    MSpG%Multip) :: ip,inp
      real,         dimension(    MSpG%Multip) :: tr
      logical,      dimension(    MSpG%Multip) :: nul
      real(kind=cp),dimension(3,192)           :: Lat_tr
      real(kind=cp),dimension(3,192)           :: aLat_tr
      integer :: i,j,k,kp,L,m, Ng,num_lat, num_alat,invt,nl,i_centre,centri
      integer,    dimension(3,3) :: identity, nulo, inver,mat,imat
      real(kind=cp),dimension(3) :: v
      character(len=80)          :: ShOp_symb !
      character(len=4)           :: ired !
      logical                    :: centrosymm
      Type(MSym_Oper_Type),dimension(MSpG%Multip) :: MSymOp
      Type(Sym_Oper_Type), dimension(MSpG%Multip) :: SymOp
      character(len=80),   dimension(MSpG%Multip) :: SymbSymOp,SymbMSymOp
      character (len=*),dimension(0:2), parameter  :: Centro = &
                                         (/"Centric (-1 not at origin)", &
                                           "Acentric                  ", &
                                           "Centric (-1 at origin)    "/)

      identity=0; nulo=0
      do i=1,3
        identity(i,i)=1
      end do
      inver=-identity
      num_lat=0; num_alat=0
      centrosymm=.false.
      nul=.false.
      MSpG%MagType=1
      centri=1 !Default value for non-centrosymmetric groups or for those having
               !the centre of symmetry out of the origin.

      !The code below is to re-order the symmetry operators provided in the input MSpG object
      !----Start re-ordering
      do i=1,MSpG%Multip
        tr(i)=sum(abs(MSpG%SymOp(i)%tr))
      end do
      ip=0
      call sort(tr,MSpG%Multip,ip)
      do i=1,MSpG%Multip
        j=ip(i)
        SymOp(i) = MSpG%SymOp(j)
        MSymOp(i)= MSpG%MSymOp(j)
        SymbSymOp(i)=MSpG%SymOpSymb(j)
        SymbMSymOp(i)=MSpG%MSymOpSymb(j)
      end do
      MSpG%SymOp(:)=SymOp(:)
      MSpG%MSymOp(:)=MSymOp(:)
      MSpG%SymOpSymb(:) = SymbSymOp(:)
      MSpG%MSymOpSymb(:)= SymbMSymOp(:)

      !Reorder again the operators in case the identity is not given as the first operator
      j=0
      imat=MSpG%SymOp(1)%Rot(:,:)
      if(.not. ( equal_matrix(imat,identity,3) .and.  sum(abs(MSpG%SymOp(1)%tr))  < 0.0001)) then
        do i=2,MSpG%Multip
          imat=MSpG%SymOp(i)%Rot(:,:)
          if(equal_matrix(imat,identity,3) .and. sum(abs(MSpG%SymOp(i)%tr))  < 0.0001) then
            j=i
            exit
          end if
        end do
        if(j == 0) then
          Err_Form=.true.
          Err_Form_Mess="The identity operator is not provided in the mCIF file"
          return
        end if
        MSpG%SymOp(j)=MSpG%SymOp(1)
        MSpG%MSymOp(j)=MSpG%MSymOp(1)
        MSpG%SymOpSymb(j)=MSpG%SymOpSymb(1)
        MSpG%MSymOpSymb(j)=MSpG%MSymOpSymb(1)
        MSpG%SymOp(1)%Rot=identity
        MSpG%SymOp(1)%tr=0.0
        MSpG%MSymOp(1)%Rot=identity
        MSpG%MSymOp(1)%phas=1.0
        MSpG%SymOpSymb(1)="x,y,z"
        MSpG%MSymOpSymb(1)="mx,my,mz"
      end if

      !Now look for centre of symmetry associated with time inversion and promote
      !it to the second position
      j=0
      do i=2,MSpG%Multip
        imat=MSpG%SymOp(i)%Rot(:,:)
        if(equal_matrix(imat,inver,3) .and. sum(abs(MSpG%SymOp(i)%tr))  < 0.0001 .and. MSpG%MSymOp(i)%phas < 0.0) then
          j=i
          exit
        end if
      end do
      if(j /= 0) then
        MSpG%SymOp(j)=MSpG%SymOp(2)
        MSpG%MSymOp(j)=MSpG%MSymOp(2)
        MSpG%SymOpSymb(j)=MSpG%SymOpSymb(2)
        MSpG%MSymOpSymb(j)=MSpG%MSymOpSymb(2)
        MSpG%SymOp(2)%Rot=inver
        MSpG%SymOp(2)%tr=0.0
        MSpG%MSymOp(2)%Rot=inver
        MSpG%MSymOp(2)%phas=-1.0
        MSpG%SymOpSymb(2)="-x,-y,-z"
        MSpG%MSymOpSymb(2)="-mx,-my,-mz"
      end if

      !----End re-ordering

      Err_Form=.false.
      ip=0

      !Determine the lattice translations and anti-translations
      !Eliminate lattice translations and anti-translations
      do j=2,MSpG%Multip
        if(nul(j)) cycle
        invt= nint(MSpG%MSymOp(j)%phas)
        if(invt < 0) MSpG%MagType=3
        if(equal_matrix(identity,MSpG%SymOp(j)%Rot(:,:),3)) then
           if(invt == 1) then
              num_lat=num_lat+1
              Lat_tr(:,num_lat)=MSpG%SymOp(j)%tr(:)
              nul(j)=.true.   !Nullify centring translations
           else
              num_alat=num_alat+1
              aLat_tr(:,num_alat)=MSpG%SymOp(j)%tr(:)
              nul(j)=.true.  !Nullify anti-centring translations
           end if
        end if
      end do  !j=2,MSpG%Multip

      if(allocated(MSpG%Latt_trans)) deallocate(MSpG%Latt_trans)
      allocate(MSpG%Latt_trans(3,num_lat+1))
      MSpG%Latt_trans=0.0
      m=1
      do j=1,num_lat
        m=m+1
        MSpG%Latt_trans(:,m) = Lat_tr(:,j)
      end do
      MSpG%Num_Lat=num_lat+1

      if(num_alat > 0) then
        MSpG%MagType=4
        if(allocated(MSpG%aLatt_trans)) deallocate(MSpG%aLatt_trans)
        allocate(MSpG%aLatt_trans(3,num_alat))
        MSpG%aLatt_trans = aLat_tr(:,1:num_alat)
        MSpG%Num_aLat=num_alat
      end if

      !Eliminate centre of symmetry {-1|t} and select that having
      !t=0 if it exist
      k=0; kp=0
      do j=2,MSpG%Multip
          invt= nint(MSpG%MSymOp(j)%phas)
          imat=MSpG%SymOp(j)%Rot(:,:)
          if(equal_matrix(imat,inver,3)) then
            if(invt == 1) then
              kp=kp+1
              ip(kp)=j
            else
              k=k+1
              inp(k)=j
            end if
          end if
      end do

      i_centre=0
      if(kp > 0) then  !Centre of symmetry exist!, select that without translations
         i_centre=ip(1)
         do j=1,kp
           i=ip(j)
           if(sum(abs(MSpG%SymOp(i)%tr))  < 0.0001) then
             i_centre=i    !localization of the -x,-y,-z,+1 operator within the list
             centri=2      !Now this value indicates that the operor -x,-y,-z,+1 exists
             centrosymm=.true.
             nul(i)=.true.
             exit
           end if
         end do
      end if

      !Nullify the operators of inversion centres associated with time inversion
      !and have a translation corresponding to a centring or anticentring vector

      do i=1,k
         j=inp(i)
         v=MSpG%SymOp(j)%tr(:)
         if(sum(abs(v)) < 0.0001) cycle !Maintain the operaror -x,-y,-z,-1
         if(is_Lattice_vec(V,Lat_tr,num_lat,nl)) then
            nul(j)=.true.
            cycle
         end if

         if(is_Lattice_vec(V,aLat_tr,num_alat,nl)) then
            nul(j)=.true.
            cycle
         end if
      end do

      !Nullify the operators that can be deduced from others by applying translations,
      !anti-translations and centre of symmetry

      ip=0
      do j=2,MSpG%Multip-1
         do i=j+1,MSpG%Multip
           if(nul(i)) cycle
           mat=MSpG%SymOp(i)%Rot(:,:)-MSpG%SymOp(j)%Rot(:,:)
           if(equal_matrix(mat,nulo,3) ) then  !Pure lattice translation or antitranslation
              v=MSpG%SymOp(i)%tr(:)-MSpG%SymOp(j)%tr(:)

              if(is_Lattice_vec(V,Lat_tr,num_lat,nl)) then
                 nul(i)=.true.
                 cycle
              end if

              if(is_Lattice_vec(V,aLat_tr,num_alat,nl)) then
                 nul(i)=.true.
                 cycle
              end if

           end if

           if(centrosymm) then
              imat=MSpG%SymOp(i)%Rot(:,:)+MSpG%SymOp(j)%Rot(:,:)
              k=nint(MSpG%MSymOp(i)%phas)
              invt=nint(MSpG%MSymOp(j)%phas)

              if(equal_matrix(imat,nulo,3) .and. k == invt) then
                 v=MSpG%SymOp(i_centre)%tr(:)-MSpG%SymOp(i)%tr(:)-MSpG%SymOp(j)%tr(:)
                 if(is_Lattice_vec(V,Lat_tr,num_lat,nl)) then
                    nul(i)=.true.
                    cycle
                 end if
              end if

              if(equal_matrix(imat,nulo,3) .and. k /= invt) then
                 if(is_Lattice_vec(V,aLat_tr,num_alat,nl)) then
                    nul(i)=.true.
                    cycle
                 end if
              end if

           end if
         end do
      end do
      j=0

      ! => This is the reduced set of symmetry operators"
      do i=1,MSpG%Multip
        !write(*,"(a,i4,2x,L)") "  "//trim(MSpG%SymOpSymb(i))//"   "//trim(MSpG%MSymOpSymb(i)), nint(MSpG%MSymOp(i)%phas), nul(i)
        if(nul(i)) cycle
        j=j+1
        SymOp(j) = MSpG%SymOp(i)
        MSymOp(j)= MSpG%MSymOp(i)
      end do

      m=j*centri*(num_alat+num_lat+1)
      if( m /= MSpG%Multip) then
        write(unit=Err_Form_Mess,fmt="(2(a,i4))") " Warning! Multip=",MSpG%Multip, " Calculated Multip: ",m
        Err_Form=.true.
        return
      end if
      !Promote the reduced set of symmetry operators to the top of the list
      MSpG%SymOp(1:j)=SymOp(1:j)
      MSpG%MSymOp(1:j)=MSymOp(1:j)
      MSpG%Numops=j

      !Re-Construct, in an ordered way, all the symmetry operators in MSpG
      !starting with the reduced set
      m=MSpG%Numops

      if(centrosymm) then   !First apply the centre of symmetry
        MSpG%Centred=2
        MSpG%centre= centro(MSpG%Centred)
        do i=1,MSpG%Numops
          m=m+1
          MSpG%SymOp(m)%Rot  = -MSpG%SymOp(i)%Rot
          MSpG%SymOp(m)%tr   =  modulo_lat(-MSpG%SymOp(i)%tr)
          MSpG%MSymOp(m)= MSpG%MSymOp(i)
        end do
      else
        if(i_centre /= 0) then
          MSpG%Centred      = 0
          MSpG%centre       = centro(MSpG%Centred)
          MSpG%Centre_coord = MSpG%SymOp(i_centre)%tr(:)/2.0
        else
          MSpG%Centred=1
          MSpG%centre= centro(MSpG%Centred)
        end if
      end if

      ng=m

      if(MSpG%Num_aLat > 0) then   !Second apply the lattice centring anti-translations
        do L=1,MSpG%Num_aLat
           do i=1,ng
             m=m+1
             v=MSpG%SymOp(i)%tr(:) + MSpG%aLatt_trans(:,L)
             MSpG%SymOp(m)%Rot  = MSpG%SymOp(i)%Rot
             MSpG%SymOp(m)%tr   = modulo_lat(v)
             MSpG%MSymOp(m)%Rot = -MSpG%MSymOp(i)%Rot
             MSpG%MSymOp(m)%phas= -MSpG%MSymOp(i)%phas
           end do
        end do
      end if

      if(MSpG%Num_Lat > 1) then  !Third apply the lattice centring translations
        do L=2,MSpG%Num_Lat
           do i=1,ng
             m=m+1
             v=MSpG%SymOp(i)%tr(:) + MSpG%Latt_trans(:,L)
             MSpG%SymOp(m)%Rot  = MSpG%SymOp(i)%Rot
             MSpG%SymOp(m)%tr   = modulo_lat(v)
             MSpG%MSymOp(m)     = MSpG%MSymOp(i)
           end do
        end do
      end if


      !Normally here the number of operators should be equal to multiplicity
      !Test that everything is OK
      ng=m
      if(ng /= MSpG%Multip) then
        write(unit=Err_Form_Mess,fmt="(2(a,i3))") " => Problem! the multiplicity ",MSpG%Multip," has not been recovered, value of ng=",ng
        Err_Form=.true.
        return
      end if
      !Now re-generate all symbols from the symmetry operators and magnetic matrices
      ired=" => "
      do i=1,MSpG%Multip
         if(i > MSpG%Numops) ired="    "
         call Get_Shubnikov_Operator_Symbol(MSpG%SymOp(i)%Rot,MSpG%MSymOp(i)%Rot,MSpG%SymOp(i)%tr,ShOp_symb,.true.)
         j=index(ShOp_symb," ")
         MSpG%SymOpSymb(i)=ShOp_symb(1:j-1)
         ShOp_symb=adjustl(ShOp_symb(j:))
         j=index(ShOp_symb," ")
         MSpG%MSymOpSymb(i)=ShOp_symb(1:j-1)
         !write(*,"(a,i6,a,i4)") ired,i, "  "//trim(MSpG%SymOpSymb(i))//"   "//trim(MSpG%MSymOpSymb(i)), nint(MSpG%MSymOp(i)%phas)
      end do
      return
    End Subroutine Cleanup_Symmetry_Operators

    !!----
    !!---- Subroutine File_To_FileList(File_dat,File_list)
    !!----   character(len=*),     intent( in) :: file_dat  !Input data file
    !!----   type(file_list_type), intent(out) :: file_list !File list structure
    !!----
    !!----    Charge an external file to an object of File_List_Type.
    !!----
    !!---- Update: August - 2008
    !!
    Subroutine File_To_FileList(File_dat,File_list)
       !---- Arguments ----!
       character(len=*),      intent( in) :: file_dat
       type(file_list_type),  intent(out) :: file_list

       !---- Local Variables ----!
       integer                           :: nlines

       !---- Number of Lines in the input file ----!
       call Number_Lines(trim(File_dat), nlines)

       if (nlines==0) then
          err_form=.true.
          ERR_Form_Mess="The file "//trim(File_dat)//" contains nothing"
          return
       else
          file_list%nlines=nlines
          if (allocated(file_list%line)) deallocate(file_list%line)
          allocate(file_list%line(nlines))
          call reading_Lines(trim(File_dat),nlines,file_list%line)
       end if

       return
    End Subroutine File_To_FileList

    !!----
    !!---- Subroutine Get_Job_Info(file_dat,i_ini,i_end,Job_info)
    !!----   character(len=*), dimension(:), intent( in) :: file_dat     !Lines of text (content of a file)
    !!----   integer,                        intent( in) :: i_ini,i_end  !Lines to explore
    !!----   type(job_info_type),            intent(out) :: Job_info     !Object to be constructed here
    !!----
    !!----
    !!----    Constructor of the object Job_info. The arrary of strings file_dat
    !!----    have to be provided as input. It contains lines corresponding to the
    !!----    input control file. The analysis of the command lines is not given here.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_Job_Info(file_dat,i_ini,i_end,Job_info)
       !---- Arguments ----!
       character(len=*), dimension(:), intent( in) :: file_dat
       integer,                        intent( in) :: i_ini,i_end
       type(job_info_type),            intent(out) :: Job_info

       !---- Local Variables ----!
       integer                           :: i,nphas, ncmd,n_pat,ier, j
       integer, dimension(i_end-i_ini+1) :: ip,ic,ipt
       real(kind=sp)                     :: a1,a2,a3,a4,a5
       character(len=120)                :: line, fmtfields, fmtformat

       !--- Initialize FindFMT
       call Init_FindFMT(i_ini)
       nphas=0
       ncmd=0
       n_pat=0
       ip=i_end
       ic=0
       ipt=0
       Job_info%title=" General Job: CrysFML"
       Job_info%Num_Patterns=1

       do i=i_ini,i_end
          line=u_case(adjustl(file_dat(i)))
          if (line(1:5) == "TITLE") Job_info%title=line(7:)
          if (line(1:5) == "NPATT") then
             read(unit=line(7:), fmt=*,iostat=ier) Job_info%Num_Patterns
             if (ier /= 0) Job_info%Num_Patterns=1
          end if
          if (line(1:6) == "PHASE_") then
             nphas=nphas+1
             ip(nphas)=i
          end if
          if (line(1:4) == "CMDL") then
             ncmd=ncmd+1
             ic(ncmd)=i
          end if
          if (line(1:5) == "PATT_") then
             n_pat=n_pat+1
             ipt(n_pat)=i
          end if
       end do

       if (nphas == 0) then
          nphas=1
          ip(nphas)=0
       end if
       if (n_pat == 0) then
          n_pat=1
          ipt(n_pat) = 0
       end if

       if (Job_info%Num_Patterns /= n_pat) Job_info%Num_Patterns = n_pat
       Job_info%Num_Phases=nphas
       Job_info%Num_Cmd=ncmd

       if (allocated(Job_Info%Patt_typ)) deallocate(Job_Info%Patt_typ)
       allocate(Job_Info%Patt_typ(n_pat))

       if (allocated(Job_Info%Phas_nam)) deallocate(Job_Info%Phas_nam)
       allocate(Job_Info%Phas_nam(nphas))

       if (allocated(Job_Info%range_stl)) deallocate(Job_Info%range_stl)
       allocate(Job_Info%range_stl(n_pat))

       if (allocated(Job_Info%range_q)) deallocate(Job_Info%range_q)
       allocate(Job_Info%range_q(n_pat))

       if (allocated(Job_Info%range_d)) deallocate(Job_Info%range_d)
       allocate(Job_Info%range_d(n_pat))

       if (allocated(Job_Info%range_2theta)) deallocate(Job_Info%range_2theta)
       allocate(Job_Info%range_2theta(n_pat))

       if (allocated(Job_Info%range_energy)) deallocate(Job_Info%range_energy)
       allocate(Job_Info%range_energy(n_pat))

       if (allocated(Job_Info%range_tof)) deallocate(Job_Info%range_tof)
       allocate(Job_Info%range_tof(n_pat))

       if (allocated(Job_Info%lambda)) deallocate(Job_Info%lambda)
       allocate(Job_Info%lambda(n_pat))

       if (allocated(Job_Info%ratio)) deallocate(Job_Info%ratio)
       allocate(Job_Info%ratio(n_pat))

       if (allocated(Job_Info%dtt1)) deallocate(Job_Info%dtt1)
       allocate(Job_Info%dtt1(n_pat))

       if (allocated(Job_Info%dtt2)) deallocate(Job_Info%dtt2)
       allocate(Job_Info%dtt2(n_pat))

       !---- Initialize all variables
       Job_Info%Patt_typ    =" "
       Job_Info%Phas_nam    =" "
       Job_Info%range_stl%mina=0.0
       Job_Info%range_stl%maxb=0.0
       Job_Info%range_q%mina=0.0
       Job_Info%range_q%maxb=0.0
       Job_Info%range_d%mina=0.0
       Job_Info%range_d%maxb=0.0
       Job_Info%range_2theta%mina=0.0
       Job_Info%range_2theta%maxb=0.0
       Job_Info%range_Energy%mina=0.0
       Job_Info%range_Energy%maxb=0.0
       Job_Info%range_tof%mina=0.0
       Job_Info%range_tof%maxb=0.0
       Job_Info%Lambda%mina=0.0
       Job_Info%Lambda%maxb=0.0
       Job_Info%ratio = 0.0
       Job_Info%dtt1 = 0.0
       Job_Info%dtt2 = 0.0
       if (ncmd > 0) then
          if (allocated(Job_Info%cmd)) deallocate(Job_Info%cmd)
          allocate(Job_Info%cmd(ncmd))
          Job_Info%cmd=" "
       end if

       !---- Fill the different fields of Job_Info
       !---- Start with patterns
       fmtfields = "9fffff"

       !---- First asks if there is a PATT_ card, if not a standard is taken
       if (ipt(1) /= 0) then
          do n_pat=1, Job_info%Num_Patterns
             i=ipt(n_pat)
             line=u_case(adjustl(file_dat(i)))
             line=line(8:)
             call findfmt(0,line,fmtfields,fmtformat)
             read(unit=line,fmt=fmtformat) Job_Info%Patt_typ(n_pat), a1,a2,a3,a4,a5
             if (ierr_fmt /= 0) return
             line=u_case(Job_Info%Patt_typ(n_pat))

             select case(line(1:9))
                case("XRAY_2THE","NEUT_2THE","XRAY_SXTA","NEUT_SXTA")
                   if ( a1 <= 0.000001) a1=1.5405
                   if ( a2 <= 0.000001) then
                      a2=a1
                      a3=0.0
                   end if
                   if (a5 <= a4) a5=120.0
                   Job_Info%Lambda(n_pat)%mina=a1
                   Job_Info%Lambda(n_pat)%maxb=a2
                   Job_Info%ratio(n_pat)=a3
                   Job_Info%range_2theta(n_pat)%mina=a4
                   Job_Info%range_2theta(n_pat)%maxb=a5
                   a4=sind(0.5*a4)/a1
                   a5=sind(0.5*a5)/a2
                   Job_Info%range_stl(n_pat)%mina=a4
                   Job_Info%range_stl(n_pat)%maxb=a5
                   Job_Info%range_q(n_pat)%mina=a4*4.0*pi
                   Job_Info%range_q(n_pat)%maxb=a5*4.0*pi
                   Job_Info%range_d(n_pat)%mina=0.5/a5
                   Job_Info%range_d(n_pat)%maxb=0.5/a4

                case("NEUT_TOF ")
                   if (a1 <= 0.000001) a1=1000.0
                   if (a4 <= a3) a4=2.0*abs(a3)
                   Job_Info%dtt1(n_pat)=a1
                   Job_Info%dtt2(n_pat)=a2
                   Job_Info%range_tof(n_pat)%mina=a3
                   Job_Info%range_tof(n_pat)%maxb=a4
                   Job_Info%range_d(n_pat)%mina=0.5*(-1.0+sqrt(1.0+4.0*a2*a3/a1/a1))
                   Job_Info%range_d(n_pat)%maxb=0.5*(-1.0+sqrt(1.0+4.0*a2*a4/a1/a1))
                   Job_Info%range_stl(n_pat)%mina=0.5/Job_Info%range_d(n_pat)%maxb
                   Job_Info%range_stl(n_pat)%maxb=0.5/Job_Info%range_d(n_pat)%mina
                   Job_Info%range_q(n_pat)%mina=Job_Info%range_stl(n_pat)%mina*4.0*pi
                   Job_Info%range_q(n_pat)%maxb=Job_Info%range_stl(n_pat)%maxb*4.0*pi

                case("XRAY_ENER")
                   if (a1 <= 0.000001) a1=12.4 !(=hc(keV.Angstr.)
                   Job_Info%dtt1(n_pat)=a1
                   Job_Info%dtt2(n_pat)=0.0
                   Job_Info%range_energy(n_pat)%mina=a3
                   Job_Info%range_energy(n_pat)%maxb=a4
                   if (a3 <= 0.00001) a3=0.01
                   if (a4 <= 0.00001) a4=2.00
                   Job_Info%range_d(n_pat)%mina=a1/a4
                   Job_Info%range_d(n_pat)%maxb=a1/a3
                   Job_Info%range_stl(n_pat)%mina=0.5/Job_Info%range_d(n_pat)%maxb
                   Job_Info%range_stl(n_pat)%maxb=0.5/Job_Info%range_d(n_pat)%mina
                   Job_Info%range_q(n_pat)%mina=Job_Info%range_stl(n_pat)%mina*4.0*pi
                   Job_Info%range_q(n_pat)%maxb=Job_Info%range_stl(n_pat)%maxb*4.0*pi

             end select
          end do

       else
          n_pat=1
          a1=1.5405
          a2=a1
          a3=0.0
          a4=0.0
          a5=120.0
          Job_Info%Patt_typ(n_pat)="XRAY_2THE"
          Job_Info%Lambda(n_pat)%mina=a1
          Job_Info%Lambda(n_pat)%maxb=a2
          Job_Info%ratio(n_pat)=a3
          Job_Info%range_2theta(n_pat)%mina=a4
          Job_Info%range_2theta(n_pat)%maxb=a5
          a4=sind(0.5*a4)/a1
          a5=sind(0.5*a5)/a2
          Job_Info%range_stl(n_pat)%mina=a4
          Job_Info%range_stl(n_pat)%maxb=a5
          Job_Info%range_q(n_pat)%mina=a4*4.0*pi
          Job_Info%range_q(n_pat)%maxb=a5*4.0*pi
          Job_Info%range_d(n_pat)%mina=0.5/a5
          Job_Info%range_d(n_pat)%maxb=0.5/a4
       end if

       !---- Phase names
       if (ip(1) /= 0) then
          do i=1,nphas
             j=ip(i)
             line=adjustl(file_dat(j))
             Job_Info%Phas_nam(i)=line(8:)
          end do
       else
          Job_Info%Phas_nam(1)= Job_info%title
       end if

       !---- Command Lines, stored but not analysed here
       do i=1,ncmd
          j=ic(i)
          line=adjustl(file_dat(j))
          Job_Info%cmd(i)=line(8:)
       end do

       return
    End Subroutine Get_Job_Info

    Subroutine Get_moment_ctr(xnr,moment,Spg,codini,codes,ord,ss,att,Ipr)
       real(kind=cp), dimension(3),            intent(in)     :: xnr
       real(kind=cp), dimension(:),            intent(in out) :: moment
       type(Magnetic_Space_Group_type),        intent(in)     :: Spg
       Integer,                                intent(in out) :: codini
       real(kind=cp), dimension(:),            intent(in out) :: codes
       integer,                       optional,intent(in)     :: ord
       integer, dimension(:),         optional,intent(in)     :: ss
       real(kind=cp), dimension(:,:), optional,intent(in)     :: att
       integer,                       optional,intent(in)     :: Ipr

       ! Local variables
       integer,           dimension(3,3) :: magm   !g, magm= delta * det(g) * g
       character(len=1),  dimension(3)   :: codd
       integer                           :: i,j,order,n,ig,is
       real(kind=cp)                     :: suma
       integer,           dimension(48)  :: ss_ptr
       real(kind=cp),     dimension(3,48):: atr
       real(kind=cp),     dimension(3)   :: cod,multi
       real(kind=cp),     dimension(3)   :: x
       real(kind=dp),     dimension(3,3) :: sCtr
       real(kind=cp),     dimension(3)   :: momentL,TotMom


       !Test if all codes are given ... in such a case the user constraints
       !are prevalent

       suma=0.0_cp
       !iq=0
       n=3 !Real moments -> three components
       do j=1,3
          suma=suma+abs(codes(j))
       end do

       if(suma < epps ) return  !No refinement is required
       if(present(Ipr)) then
         write(Ipr,"(/,a)")         " => Calculation of symmetry constraints for magnetic moments "
       end if
       x=xnr
       !where(x < 0.0) x=x+1.0
       !where(x > 1.0) x=x-1.0

       if(present(ord) .and. present(ss) .and. present(att)) then
         order=ord
         ss_ptr(1:order) = ss(1:ord)
         atr(:,1:order)  = att(:,1:ord)
       else
         call get_stabilizer(x,SpG,order,ss_ptr,atr)
         if(present(ipr)) Write(unit=ipr,fmt="(a,i3)") " => Stabilizer without identity, order:",order
       end if

       momentL=moment
       sCtr=0.0_cp
       if(order > 1) then
         do ig=1,order
           magm(:,:) = Spg%MSymOp(ss_ptr(ig))%Rot
           sCtr=sCtr+magm !Adding constraint matrices for each operator of stabilizer
           if(present(ipr)) then
             write(unit=ipr,fmt='(a,i2,a,t20,a,t55,a,t75,9f8.4)') '     Operator ',ig,": ",trim(Spg%SymopSymb(ss_ptr(ig))), &
              trim(Spg%MSymopSymb(ss_ptr(ig))), sCtr
           end if
         end do  !ig operators
         sCtr=sCtr/order
         suma=sum(abs(sCtr))
         !write(*,"(a,f10.4,a,i3)") " suma:",suma, "Mag_Type:", spg%mag_type
         if(suma < epps .or. spg%magtype == 2) then !This corresponds to a grey point group
            moment=0.0_cp
            codes=0.0_cp
            if(present(Ipr)) then
              write(Ipr,"(a)")         " Grey point group or symmetry paramagnetic site: the attached moment is zero "
              write(Ipr,"(a,24f14.6)") " Final codes: ",codes(1:n)
              write(Ipr,"(a,24f14.6)") " Constrained moment: ",moment
            end if
            return
         end if
         TotMom=matmul(sCtr,momentL)
         call Get_Refinement_Codes(n,TotMom,sCtr,is,multi,codd,momentL)
         cod=0.0
         do j=1,n
           if(codd(j) /= "0") then
             do i=1,is
               if(codd(j) == cdd(i)) then
                 cod(j)=codini+i
                 exit
               end if
             end do
           end if
         end do
         moment=momentL
         codes=0.0
         do j=1,n
           if(abs(multi(j)) > epps)  codes(j) = sign(1.0_cp, multi(j))*(abs(cod(j))*10.0_cp + abs(multi(j)) )
         end do
         codini=codini+is
         if(present(Ipr)) then
           Write(unit=Ipr,fmt="(a,i4)")       " Number of free parameters: ",is
           write(unit=Ipr,fmt="(a,3f14.6)")   " Multipliers: ",(multi(j), j=1,n)
           write(unit=Ipr,fmt="(28a)")        " String with free parameters: ( ",(codd(j)//", ",j=1,n-1),codd(n)//" )"
           write(unit=Ipr,fmt="(a,3i6)")      " Resulting integer codes: ", nint(cod(1:n))
           write(unit=Ipr,fmt="(a,3f14.6)")   " Final codes: ",codes(1:n)
           write(unit=Ipr,fmt="(a,3f14.6)")   " Constrained Moment: ",moment
         end if

       else !No restrictions

         codd(1:n)=cdd(1:n)
         multi(1:n)=1.0_cp
         do j=1,n
           cod(j)=codini+j
           codes(j) = (abs(cod(j))*10.0_cp + abs(multi(j)))
         end do
         codini=codini+n
         if(present(Ipr)) then
           write(unit=Ipr,fmt="(a)")         " General position, no constraints in moment "
           write(unit=Ipr,fmt="(28a)")       " String with free parameters: ( ",(codd(j)//", ",j=1,n-1),codd(n)//" )"
           write(unit=Ipr,fmt="(a,24i6)")    " Resulting integer codes: ", nint(cod(1:n))
           write(unit=Ipr,fmt="(a,24f14.6)") " Final codes: ",codes(1:n)
           write(unit=Ipr,fmt="(a,24f14.6)") " Constrained moment: ",moment
         end if

       end if
       return
    End Subroutine Get_moment_ctr

    Subroutine Get_Refinement_Codes(n,vect_val,Ctr,is,multi,codd,vect_out)
      integer,                       intent(in)    :: n !dimension of the vector and the matrix
      real(kind=cp), dimension(:),   intent(in)    :: vect_val
      real(kind=dp), dimension(:,:), intent(in out):: Ctr
      integer,                       intent(out)   :: is
      real(kind=cp), dimension(:),   intent(out)   :: multi
      character(len=*), dimension(:),intent(out)   :: codd
      real(kind=cp), dimension(:),   intent(out)   :: vect_out
      !--- Local variables ---!
      real(kind=cp), dimension(n)   :: val
      integer,       dimension(n)   :: pti
      real(kind=dp), dimension(n,n) :: zv
      integer                       :: i,j,k,kval,ier,ip
      real(kind=dp)                 :: zmi
      real(kind=dp), dimension(n)   :: Wr, Wi
      logical,       dimension(n)   :: done

      !Diagonalize the matrix and pickup the lambda=1 eigenvalues
      !The corresponding eigenvector contains the constraints of all moment components
      !Calling the general diagonalization subroutine from EisPack
      call rg_ort(n,Ctr,wr,wi,.true.,zv,ier)
      is=0
      pti=0
      kval=0
      do i=1,n
        if(abs(wr(i)-1.0_dp) < epps .and. abs(wi(i)) < epps) then
          is=is+1   !Number of eigenvalues = 1 => number of free parameters
          pti(is)=i !This points to the eigenvectors with eigenvalue equal to 1.
          zmi=1.0e6 !normalize the eigenvectors so that the minimum (non-zero value) is 1.
          j=1
          do k=1,n
            if(abs(zv(k,i)) < epps) cycle
            if(abs(zv(k,i)) < zmi) then
              zmi=abs(zv(k,i))
              kval=k  !This is the basis value
              j=nint(sign(1.0_dp,zv(k,i)))
            end if
          end do
          zv(:,i)=j*zv(:,i)/zmi  !This provides directly the multipliers for a single lambda=1 eigenvalue
          val(is)=vect_val(kval) !This is the basis value to construct the new Moment
        end if
      end do
      codd="0"
      vect_out=0.0
      multi=0.0
      done=.false.
      where(abs(vect_val) < epps) done=.true.
      Select Case(is)
        case(1)
          vect_out(:)=val(1)*zv(:,pti(1))
          where(abs(vect_out) > epps)  codd(:)=cdd(1)
          multi(:)=zv(:,pti(1))
        case(2:)
          ip=0
          do i=1,n
            if(.not. done(i)) then
              if(abs(vect_val(i)) > epps) then
                ip=ip+1
                codd(i)=cdd(ip)
                multi(i)=1.0
                vect_out(i)=vect_val(i)
                done(i)=.true.
                do j=i+1,n
                  if(.not. done(j)) then
                    if(abs(vect_val(i)-vect_val(j)) < epps) then
                      codd(j)=cdd(ip)
                      multi(j)=1.0
                      vect_out(j)=vect_val(i)
                      done(j)=.true.
                    else if(abs(vect_val(i)+vect_val(j)) < epps) then
                      codd(j)=cdd(ip)
                      multi(j)=-1.0
                      vect_out(j)=-vect_val(i)
                      done(j)=.true.
                    end if
                  end if
                end do
              end if
            end if
          end do
      End Select
    End Subroutine Get_Refinement_Codes

    !!----
    !!---- Subroutine Init_Err_Form()
    !!----
    !!----    Initialize Errors Variable for this module
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Init_Err_Form()

       err_form=.false.
       ERR_Form_Mess=" "

       return
    End Subroutine Init_Err_Form

    !!----
    !!---- Subroutine Read_Atom(Line,Atomo)
    !!----    character(len=*), intent(in out ) :: line    !  In -> Input String with ATOM directive
    !!----    Type (Atom_Type), intent(out)     :: Atomo   ! Out -> Parameters on variable Atomo
    !!----
    !!----    Subroutine to read the atom parameters from a given "line"
    !!----    it construct the object Atomo of type Atom.
    !!----    Control of error is present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Atom(line,Atomo)
       !---- Arguments ----!
       character(len=*), intent(in out ) :: line
       Type (Atom_Type), intent(out)     :: Atomo

       !---- Local variables -----!
       integer                           :: iv, nlong1,n,ier,q
       real(kind=cp), dimension (10)     :: vet1
       real(kind=cp), dimension (10)     :: vet2
       character(len=4)                  :: dire
       character(len=40)                 :: magmom
       character(len=5)                  :: label
       character(len=132), dimension(1)  :: filevar
       character(len=*), parameter       :: digpm="0123456789+-"

       !---- Init ----!
       call init_err_form()
       call init_atom_type(Atomo)
       q=0
       iv=index(line,"#")
       if(iv /= 0) then
        atomo%AtmInfo=line(iv+1:)
        line=line(1:iv-1)
       end if
       iv=index(line,"Moment:") !Attemp to read magnetic moment components
       magmom=" "
       if(iv /= 0) then
         magmom=line(iv+7:)  !magmon should contain magnetic moment
         line=line(1:iv-1)   !Line after removing "Moment:" and infor
       end if
       call cutst(line,nlong1,dire)
       if (u_case(dire) /= "ATOM") then
          err_form=.true.
          ERR_Form_Mess=" Error reading the ATOM keyword"
          return
       end if

       !---- Atom Label ----!
       call cutst(line,nlong1,label)
       atomo%lab=label(1:5)

       !---- Atom Type (Chemical symbol & Scattering Factor) ----!
       call cutst(line,nlong1,label)

       if(len_trim(magmom) == 0) then
          n=index(digpm,label(2:2))
          if (n /=0) then
            atomo%chemsymb=U_case(label(1:1))
          else
            atomo%chemsymb=U_case(label(1:1))//L_case(label(2:2))
          end if
       else
          n=index(digpm,label(4:4))
          if(U_case(label(1:1)) /= "M" .and. U_case(label(1:1)) /= "J") then
             err_form=.true.
             ERR_Form_Mess=" Error reading the magnetic form factor of ATOM: "//trim(atomo%lab)
             return
          end if
          atomo%chemsymb=U_case(label(2:2))//L_case(label(3:3))
       end if
       atomo%SfacSymb=label(1:4)

       !---- Parameters ----!
       filevar(1)="atm "//trim(line)

       n=1
       call Read_Key_ValueSTD(filevar,n,n,"atm",vet1,vet2,iv)
      ! call getnum(line,vet,ivet,iv)
       if (iv <= 0) then
          err_form=.true.
          ERR_Form_Mess= "Error reading parameters of atom:"//atomo%lab
          return
       end if

       !---- Coordinates  ----!
       if (iv < 3) then
          err_form=.true.
          ERR_Form_Mess= "Error reading Coordinates of atom:"//atomo%lab
          return
       end if

       atomo%x(:)=vet1(1:3)
       atomo%x_std(:)=vet2(1:3)

       !---- Biso ----!
       if (iv > 3) then
         atomo%biso=vet1(4)
         atomo%biso_std=vet2(4)
       end if

       !---- Occ ----!
       if (iv > 4) then
          atomo%occ=vet1(5)
          atomo%occ_std=vet2(5)
       end if

       !---- Moment ----!
       if (iv > 5) atomo%moment=vet1(6)

       !---- Charge ----!
       if (iv > 6) atomo%charge=vet1(7)

       !Attempt to get the oxidation state from "Label"
       if(abs(atomo%charge) < eps) then
         iv=index(label,"+")
         Select Case(iv)
           Case(0) !No + sign
             n=index(label,"-")
             Select Case(n)
               Case(2) !Element with a single character symbol F-1
                  read(unit=label(3:),fmt="(i1)",iostat=ier)  q
                  if (ier /= 0) q=0
               Case(3) !Element in the form: F1- or Br-1
                  read(unit=label(2:2),fmt="(i1)",iostat=ier)  q
                  if (ier /= 0) then
                        read(unit=label(4:4),fmt="(i1)",iostat=ier)  q
                        if (ier /= 0) q=0
                  end if
               Case(4) !Element in the form: Br1-
                  read(unit=label(3:3),fmt="(i1)",iostat=ier)  q
                  if (ier /= 0) q=0
             End Select
             q=-q   !anions
           Case(2) !Element with a single character symbol C+4
                  read(unit=label(3:),fmt="(i1)",iostat=ier)  q
                  if (ier /= 0) q=0
           Case(3) !Element in the form: C4+ or Fe+3
                  read(unit=label(2:2),fmt="(i1)",iostat=ier)  q
                  if (ier /= 0) then
                        read(unit=label(4:4),fmt="(i1)",iostat=ier)  q
                        if (ier /= 0) q=0
                  end if
           Case(4) !Element in the form: Fe3+
                  read(unit=label(3:3),fmt="(i1)",iostat=ier)  q
                  if (ier /= 0) q=0
         End Select
         atomo%charge=real(q)
       end if

       !Now read the components of the magnetic moment
       if(len_trim(magmom) /= 0) then
         filevar(1)="mom "//trim(magmom)
         call Read_Key_ValueSTD(filevar,n,n,"mom",vet1,vet2,iv)

         !---- Moment components  ----!
         if (iv < 3) then
            err_form=.true.
            ERR_Form_Mess= "Error reading magnetic moment components of atom:"//atomo%lab
            return
         end if
         atomo%m_xyz(:)=vet1(1:3)
         atomo%sm_xyz(:)=vet2(1:3)
         if(atomo%moment < 0.001) atomo%moment=maxval(abs(atomo%m_xyz(:)))
       end if

       return
    End Subroutine Read_Atom

    !!----
    !!---- Subroutine Read_Cell(Line,Celda)
    !!----    character(len=*),          intent(in out ) :: line   !  In -> Input String with CELL Directive
    !!----    real(kind=cp),dimension(6),intent(out)     :: Celda  !  In -> Parameters on Celda Variable
    !!----
    !!----    Subroutine to read the cell parameters from a given "line"
    !!----    it construct the object Celda of type Crystal_Cell.
    !!----    Assumes the string "line" has been read from a file and
    !!----    starts with the word "cell", that is removed before reading
    !!----    the values of the parameters.
    !!----    Control of error is present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Cell(line,Celda)
       !---- Arguments ----!
       character(len=*),          intent(in out ) :: line
       real(kind=cp),dimension(6),intent(out)     :: Celda

       !---- Local variables -----!
       integer, dimension (6)               :: ivet
       real(kind=cp), dimension (6)         :: vet
       integer                              :: nlong1,iv
       character(len=4)                     :: dire

       call init_err_form()

       call cutst(line,nlong1,dire)
       if (u_case(dire) /= "CELL") then
          err_form=.true.
          ERR_Form_Mess=" Error reading the CELL keyword"
          return
       end if

       call getnum(line,vet,ivet,iv)
       if (iv /= 6 ) then
          err_form=.true.
          ERR_Form_Mess=" Error reading the Cell Parameters"
          return
       else
          celda=vet
       end if

       return
    End Subroutine Read_Cell

    !!----
    !!---- Subroutine Read_Cif_Atom(Filevar,Nline_Ini,Nline_End,N_Atom,Atm_List)
    !!----    character(len=*),dimension(:), intent(in)     :: filevar    !  In -> Input strings information
    !!----    integer,                       intent(in out) :: nline_ini  !  In -> Line to beginning search
    !!----                                                                   Out -> Current line on Filevar
    !!----    integer,                       intent(in)     :: nline_end  !  In -> Line to the End search
    !!----    integer,                       intent(out)    :: n_atom     ! Out -> Actual number of atom
    !!----    type (atom_list_type),        intent(out)    :: Atm_List   ! Out -> Atom list
    !!----
    !!----    Obtaining Atoms parameters from Cif file. A control error is present.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Cif_Atom(filevar,nline_ini,nline_end,n_atom,Atm_List)
       !---- Arguments ----!
       character(len=*), dimension(:),   intent(in)      :: filevar
       integer,                          intent(in out)  :: nline_ini
       integer,                          intent(in)      :: nline_end
       integer,                          intent(out)     :: n_atom
       type (atom_list_type),            intent(out)     :: Atm_List

       !---- Local Variables ----!
       character(len=len(filevar(1)))      :: string
       character(len=20),dimension(15)     :: label
       type(Atom_Type)                     :: aux_atm
       integer                             :: i, j, nc, nct, nline, iv, First, nline_big,num_ini,mm
       integer, dimension( 8)              :: lugar   !   1 -> label
                                                      !   2 -> Symbol
                                                      ! 3-5 -> coordinates
                                                      !   6 -> occupancy
                                                      !   7 -> Uequi
                                                      !   8 -> Biso
       real(kind=cp), dimension(1)     :: vet1,vet2
       type(atom_list_type)            :: Atm

       !---- Estimacion Inicial ----!
       lugar=0
       call allocate_atom_list(nline_end-nline_ini+1,Atm)

       num_ini=nline_ini
       n_atom=0
       !Change of _atom_site by _atom_site_label in order to be able the reading
       !the atoms positions even when the anisotropic parameters are given before
       call Read_Key_StrVal(filevar,nline_ini,nline_end,"_atom_site_label",string)
       !Look for the possibility that _atom_site_label is not the first item in the loop
       do i=nline_ini,num_ini,-1
         string=adjustl(filevar(i))
         if(string(1:) == "loop_") then
          nline_ini=i+1
          exit
         end if
       end do
       j=0
       do i=nline_ini,nline_end
          string=adjustl(filevar(i))
          if ("_atom_site_label" == string(1:16)) then
             j=j+1
             lugar(1)=j
             cycle
          end if
          if ("_atom_site_type_symbol" == string(1:22)) then
             j=j+1
             lugar(2)=j
             cycle
          end if
          if ("_atom_site_fract_x" == string(1:18)) then
             j=j+1
             lugar(3)=j
             cycle
          end if
          if ("_atom_site_fract_y" == string(1:18)) then
             j=j+1
             lugar(4)=j
             cycle
          end if
          if ("_atom_site_fract_z" == string(1:18)) then
             j=j+1
             lugar(5)=j
             cycle
          end if
          if ("_atom_site_occupancy" == string(1:20)) then
             j=j+1
             lugar(6)=j
             cycle
          end if
          if ("_atom_site_U_iso_or_equiv" == string(1:25)) then
             j=j+1
             lugar(7)=j
             cycle
          end if
          if ("_atom_site_B_iso_or_equiv" == string(1:25)) then
             j=j+1
             lugar(8)=j
             cycle
          end if
          if ("_atom_site_" == string(1:11)) then
             j=j+1
             cycle
          end if

          if ("_oxford_atom_site_" == string(1:18)) then
             j=j+1
             cycle
          end if

          nline=i
          exit
       end do

       if (any(lugar(3:5) == 0)) then
          err_form=.true.
          ERR_Form_Mess=" Error reading atoms"
          return
       end if
       nct=count(lugar > 0)
       nline_big=nline
       nline_ini=nline
       do i=nline_ini,nline_end
          string=adjustl(filevar(i))
          if (string(1:1) == "#" .or. string(1:1) == "?") cycle
          if (len_trim(string) == 0) exit
          if (string(1:1) == "_" .or. string(1:5) == "loop_") exit
          call getword(string,label,nc)
          if (nc < nct) then
            nline=i
            exit
          end if
          n_atom=n_atom+1

          ! _atom_site_label
          atm%atom(n_atom)%lab=label(lugar(1))

          ! _atom_site_type_symbol
          if (lugar(2) /= 0) then
             atm%atom(n_atom)%SfacSymb=label(lugar(2))(1:4)
             if(index("1234567890+-",label(lugar(2))(2:2)) /= 0 ) then
                atm%atom(n_atom)%chemSymb=U_case(label(lugar(2))(1:1))
             else
                atm%atom(n_atom)%chemSymb=U_case(label(lugar(2))(1:1))//L_case(label(lugar(2))(2:2))
             end if
          else
             if(index("1234567890+-",label(lugar(1))(2:2)) /= 0 ) then
                atm%atom(n_atom)%chemSymb=U_case(label(lugar(1))(1:1))
             else
                atm%atom(n_atom)%chemSymb=U_case(label(lugar(1))(1:1))//L_case(label(lugar(1))(2:2))
             end if
             atm%atom(n_atom)%SfacSymb=atm%atom(n_atom)%chemSymb
          end if

          call getnum_std(label(lugar(3)),vet1,vet2,iv)    ! _atom_site_fract_x
          atm%atom(n_atom)%x(1)=vet1(1)
          atm%atom(n_atom)%x_std(1)=vet2(1)
          call getnum_std(label(lugar(4)),vet1,vet2,iv)    ! _atom_site_fract_y
          atm%atom(n_atom)%x(2)=vet1(1)
          atm%atom(n_atom)%x_std(2)=vet2(1)
          call getnum_std(label(lugar(5)),vet1,vet2,iv)    ! _atom_site_fract_z
          atm%atom(n_atom)%x(3)=vet1(1)
          atm%atom(n_atom)%x_std(3)=vet2(1)

          ! _atom_site_occupancy
          if (lugar(6) /= 0) then
             call getnum_std(label(lugar(6)),vet1,vet2,iv)
          else
             vet1=1.0
          end if
          atm%atom(n_atom)%occ=vet1(1)
          atm%atom(n_atom)%occ_std=vet2(1)

          if (lugar(7) /= 0) then
             call getnum_std(label(lugar(7)),vet1,vet2,iv)    ! _atom_site_U_iso_or_equiv
             atm%atom(n_atom)%ueq=vet1(1)
             atm%atom(n_atom)%Biso=vet1(1)*78.95683521     !If anisotropic they
             atm%atom(n_atom)%Biso_std=vet2(1)*78.95683521 !will be put to zero
          else if (lugar(8) /= 0) then
             call getnum_std(label(lugar(8)),vet1,vet2,iv)    ! _atom_site_B_iso_or_equiv
             atm%atom(n_atom)%ueq=vet1(1)/78.95683521
             atm%atom(n_atom)%Biso=vet1(1)     !If anisotropic they
             atm%atom(n_atom)%Biso_std=vet2(1) !will be put to zero
          else
             atm%atom(n_atom)%ueq=0.0
             atm%atom(n_atom)%Biso=0.0     !If anisotropic they
             atm%atom(n_atom)%Biso_std=0.0 !will be put to zero
          end if
          atm%atom(n_atom)%utype="u_ij"
       end do

       if(nline >= nline_big) nline_big=nline
       !---- Anisotropic parameters ----!
       nline_ini=num_ini !Changed to be able the reading of anisotropic parameters
                         !even if given before the coordinates
       lugar=0
       call Read_Key_StrVal(filevar,nline_ini,nline_end,"_atom_site_aniso_",string)

       j=0
       do i=nline_ini,nline_end
          string=adjustl(filevar(i))
          !write(*,"(i6,a)") i,"  "//trim(string)
          if ("_atom_site_aniso_label" == string(1:22)) then
             j=j+1
             lugar(1)=j
             cycle
          end if
          if ("_atom_site_aniso_type_symbol" == string(1:28)) then
             j=j+1
             lugar(8)=j
             cycle
          end if
          if ("_atom_site_aniso_U_11" == string(1:21)) then
             j=j+1
             lugar(2)=j
             cycle
          end if
          if ("_atom_site_aniso_U_22" == string(1:21)) then
             j=j+1
             lugar(3)=j
             cycle
          end if
          if ("_atom_site_aniso_U_33" == string(1:21)) then
             j=j+1
             lugar(4)=j
             cycle
          end if
          if ("_atom_site_aniso_U_12" == string(1:21)) then
             j=j+1
             lugar(5)=j
             cycle
          end if
          if ("_atom_site_aniso_U_13" == string(1:21)) then
             j=j+1
             lugar(6)=j
             cycle
          end if
          if ("_atom_site_aniso_U_23" == string(1:21)) then
             j=j+1
             lugar(7)=j
             cycle
          end if

          if ("_atom_site_aniso" == string(1:16) ) then
             j=j+1
             cycle
          end if

          nline=i
          exit
       end do
       if(nline >= nline_big) nline_big=nline
       !if (all(lugar > 0)) then
       if (all(lugar(1:7) > 0)) then        ! T.R. June 2017
          nct=count(lugar > 0)
          nline_ini=nline
          mm=0
          do i=nline_ini,nline_end
             string=adjustl(filevar(i))
             if (string(1:1) == "#" .or. string(1:1) =="?") cycle
             if (len_trim(string) == 0) exit
             call getword(string,label,nc)
             if (nc < nct) then
                nline=i
                exit
             end if
             do j=1,n_atom
                if (atm%atom(j)%thtype == "aniso") cycle ! already assigned
                if (trim(atm%atom(j)%lab) /= trim(label(lugar(1))) ) cycle
                call getnum_std(label(lugar(2)),vet1,vet2,iv)    ! _atom_site_aniso_U_11
                atm%atom(j)%u(1)=vet1(1)
                atm%atom(j)%u_std(1)=vet2(1)
                call getnum_std(label(lugar(3)),vet1,vet2,iv)    ! _atom_site_aniso_U_22
                atm%atom(j)%u(2)=vet1(1)
                atm%atom(j)%u_std(2)=vet2(1)
                call getnum_std(label(lugar(4)),vet1,vet2,iv)    ! _atom_site_aniso_U_33
                atm%atom(j)%u(3)=vet1(1)
                atm%atom(j)%u_std(3)=vet2(1)
                call getnum_std(label(lugar(5)),vet1,vet2,iv)    ! _atom_site_aniso_U_12
                atm%atom(j)%u(4)=vet1(1)
                atm%atom(j)%u_std(4)=vet2(1)
                call getnum_std(label(lugar(6)),vet1,vet2,iv)    ! _atom_site_aniso_U_13
                atm%atom(j)%u(5)=vet1(1)
                atm%atom(j)%u_std(5)=vet2(1)
                call getnum_std(label(lugar(7)),vet1,vet2,iv)    ! _atom_site_aniso_U_23
                atm%atom(j)%u(6)=vet1(1)
                atm%atom(j)%u_std(6)=vet2(1)

                atm%atom(j)%thtype="aniso"
                atm%atom(j)%Biso=0.0
                atm%atom(j)%Biso_std=0.0
                mm=mm+1
                exit
             end do
          end do

       end if
       if(nline >= nline_big) nline_big=nline
       nline_ini=nline_big

       !Look for the first atoms fully occupying the site and put it in first position
       !This is needed for properly calculating the occupation factors
       !after normalization in subroutine Readn_Set_XTal_CIF
       vet1(1)=maxval(atm%atom(:)%occ)  !Normalize occupancies
       atm%atom(:)%occ=atm%atom(:)%occ/vet1(1)
       First=1
       do i=1,n_atom
        if(abs(atm%atom(i)%occ-1.0) < 0.0001) then
          First=i
          exit
        end if
       end do
       !Swapping the orinal atom at the first position with the first having full occupation
       if(First /= 1) Then
         aux_atm=atm%atom(1)
         atm%atom(1)=atm%atom(First)
         atm%atom(First)=aux_atm
       end if

       !Put the first atom the first having a full occupation factor 1.0
       !---- Adjusting ... ----!
       if (n_atom > 0) then
          call allocate_atom_list(n_atom,Atm_list)
          atm_list%natoms=n_atom
          do i=1,n_atom
             atm_list%atom(i)=atm%atom(i)
          end do
       end if
       call Deallocate_atom_list(atm)

       return
    End Subroutine Read_Cif_Atom

    !!----
    !!---- Subroutine Read_Cif_Cell(Filevar,Nline_Ini,Nline_End,Celda,Stdcelda)
    !!----    character(len=*), dimension(:), intent(in)     :: filevar      !  In -> String vector input
    !!----    integer,                        intent(in out) :: nline_ini    !  In -> Line to start the search
    !!----                                                                      Out -> Current line on Filevar
    !!----    integer,                        intent(in)     :: nline_end    !  In -> Line to finish the search
    !!----    real(kind=cp),dimension(6),     intent (out)   :: Celda        ! Out -> Cell variable
    !!----    real(kind=cp),dimension(6),     intent (out)   :: StdCelda     ! Out -> Cell variable
    !!----
    !!----    Read Cell Parameters from Cif file
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Cif_Cell(Filevar,Nline_Ini,Nline_End,Celda,StdCelda)
       !---- Arguments ----!
       character(len=*),  dimension(:),     intent(in)     :: filevar
       integer,                             intent(in out) :: nline_ini
       integer,                             intent(in)     :: nline_end
       real(kind=cp),dimension(6),          intent(out)    :: Celda
       real(kind=cp),dimension(6),optional, intent(out)    :: StdCelda

       !---- Local Variables ----!
       integer                     :: iv,initl
       real(kind=cp), dimension(1) :: vet1,vet2
       real(kind=cp), dimension(6) :: a

       !---- Valores iniciales ----!
       celda=(/1.0,1.0,1.0,90.0,90.0,90.0/)
       a=0.0
       if (present(stdcelda)) stdcelda=0.0

       !---- Celda ----!
       initl=nline_ini  !Preserve initial line => some CIF files have random order for cell parameters
       call read_key_valueSTD(filevar,nline_ini,nline_end,"_cell_length_a",vet1,vet2,iv)
       if (iv == 1) then
          Celda(1)   =vet1(1)
          a(1)=vet2(1)
       end if

       nline_ini=initl
       call read_key_valueSTD(filevar,nline_ini,nline_end,"_cell_length_b",vet1,vet2,iv)
       if (iv == 1) then
          Celda(2)   =vet1(1)
          a(2)=vet2(1)
       end if

       nline_ini=initl
       call read_key_valueSTD(filevar,nline_ini,nline_end,"_cell_length_c",vet1,vet2,iv)
       if (iv == 1) then
          Celda(3)   =vet1(1)
         a(3)=vet2(1)
       end if

       nline_ini=initl
       call read_key_valueSTD(filevar,nline_ini,nline_end,"_cell_angle_alpha",vet1,vet2,iv)
       if (iv == 1) then
          Celda(4)   =vet1(1)
          a(4)=vet2(1)
       end if

       nline_ini=initl
       call read_key_valueSTD(filevar,nline_ini,nline_end,"_cell_angle_beta",vet1,vet2,iv)
       if (iv == 1) then
          Celda(5)   =vet1(1)
          a(5)=vet2(1)
       end if

       nline_ini=initl
       call read_key_valueSTD(filevar,nline_ini,nline_end,"_cell_angle_gamma",vet1,vet2,iv)
       if (iv == 1) then
          Celda(6)   =vet1(1)
          a(6)=vet2(1)
       end if
       if (present(stdcelda)) stdcelda=a

       return
    End Subroutine Read_Cif_Cell

    !!----
    !!---- Subroutine Read_Cif_ChemicalName(Filevar,Nline_Ini,Nline_End,ChemName)
    !!----    character(len=*),  dimension(:), intent(in) :: filevar      !  In -> String vector
    !!----    integer,           intent(in out)           :: nline_ini    !  In -> Line to start the search
    !!----                                                                  Out -> Actual line on Filevar
    !!----    integer,           intent(in)               :: nline_end    !  In -> Line to finish the search
    !!----    character(len=*),  intent(out)              :: ChemName     ! Out -> Title string
    !!----
    !!----    Obtaining Chemical Name from Cif file
    !!----
    !!---- Update: March - 2009
    !!
    Subroutine Read_Cif_ChemicalName(Filevar,Nline_Ini,Nline_End,ChemName)
       !---- Arguments ----!
       character(len=*),  dimension(:), intent(in) :: filevar
       integer,           intent(in out)           :: nline_ini
       integer,           intent(in)               :: nline_end
       character(len=*),  intent(out)              :: ChemName

       !---- Local variables ----!
       integer :: np1, np2

       ChemName=" "
       call Read_Key_StrVal(filevar,nline_ini,nline_end, &
                            "_chemical_name_common",ChemName)

       if (len_trim(chemname) == 0) then
          call Read_Key_StrVal(filevar,nline_ini,nline_end, &
                            "_chemical_name_systematic",ChemName)
       end if

       if (len_trim(chemname) > 0) then
          if (trim(chemname) =="; ?" .or. trim(chemname)=="#") chemname=" "
          np1=index(chemname,"'")
          np2=index(chemname,"'",back=.true.)
          if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
             chemname=chemname(np1+1:np2-1)
          else
             np1=index(chemname,'"')
             np2=index(chemname,'"',back=.true.)
             if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
                chemname=chemname(np1+1:np2-1)
             end if
          end if
       end if

       return
    End Subroutine Read_Cif_ChemicalName

    !!----
    !!---- Subroutine Read_Cif_Cont(Filevar,Nline_Ini,Nline_End,N_Elem_Type,Elem_Type,N_Elem)
    !!----    character(len=*), dimension(:),      intent(in)      :: filevar       !  In -> String vector input
    !!----    integer,                             intent(in out)  :: nline_ini     !  In -> Line to start the search
    !!----                                                                             Out -> Actual line on Filevar
    !!----    integer,                             intent(in)      :: nline_end     !  In -> Line to finish the search
    !!----    integer,                             intent(out)     :: n_elem_type   ! Out -> N. of different elements
    !!----    character(len=*), dimension(:),      intent(out)     :: elem_type     ! Out -> String for Element type
    !!----    real(kind=cp), dimension(:),optional,intent(out)     :: n_elem        ! Out -> Number of elements
    !!----
    !!----    Obtaining the chemical contents from Cif file
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Cif_Cont(Filevar,Nline_Ini,Nline_End,N_Elem_Type,Elem_Type,N_Elem)
       !---- Arguments ----!
       character(len=*), dimension(:),      intent(in)      :: filevar
       integer,                             intent(in out)  :: nline_ini
       integer,                             intent(in)      :: nline_end
       integer,                             intent(out)     :: n_elem_type
       character(len=*), dimension(:),      intent(out)     :: elem_type
       real(kind=cp), dimension(:),optional,intent(out)     :: n_elem

       !---- Local  variables ----!
       character(len=len(filevar(1)))      :: string
       character(len=10),dimension(15)     :: label

       integer                    :: iv
       integer                    :: i,np1,np2,nlabel,nlong
       integer, dimension(1)      :: ivet

       real(kind=cp),dimension(1) :: vet

       n_elem_type = 0
       elem_type   = " "
       if (present(n_elem)) n_elem = 0.0

       call Read_Key_StrVal(filevar,nline_ini,nline_end, &
                            "_chemical_formula_sum",string)
       if (len_trim(string) ==0) string=filevar(nline_ini+1)
       string=adjustl(string)
       if (string(1:1) == "?") return
       np1=index(string,"'")
       np2=index(string,"'",back=.true.)
       nlabel=0
       if (np1 /= 0 .and. np2 /= 0 .and. np2 > np1) then
          call getword(string(np1+1:np2-1),label,nlabel)
       end if
       if (nlabel /=0) then
          n_elem_type = nlabel
          do i=1,nlabel
             nlong=len_trim(label(i))
             select case (nlong)
                 case (1)
                    elem_type(i)=label(i)(1:1)
                    if (present(n_elem)) n_elem(i)   = 1.0

                 case (2)
                    call getnum(label(i)(2:),vet,ivet,iv)
                    if (iv == 1) then
                       elem_type(i)=label(i)(1:1)
                       if (present(n_elem)) n_elem(i)   =vet(1)
                    else
                       elem_type(i)=label(i)(1:2)
                       if (present(n_elem)) n_elem(i)   = 1.0
                    end if

                 case (3:)
                    call getnum(label(i)(2:),vet,ivet,iv)
                    if (iv == 1) then
                       elem_type(i)=label(i)(1:1)
                       if (present(n_elem)) n_elem(i)   =vet(1)
                    else
                       call getnum(label(i)(3:),vet,ivet,iv)
                       if (iv == 1) then
                          elem_type(i)=label(i)(1:2)
                          if (present(n_elem)) n_elem(i)   =vet(1)
                       else
                          elem_type(i)=label(i)(1:2)
                          if (present(n_elem)) n_elem(i)   = 1.0
                       end if

                    end if

             end select
          end do
       end if

       return
    End Subroutine Read_Cif_Cont

    !!----
    !!---- Subroutine Read_Cif_Hall(Filevar,Nline_Ini,Nline_End,Spgr_Ha)
    !!----    character(len=*), dimension(:), intent(in) :: filevar      !  In -> String vector input
    !!----    integer,          intent(in out)           :: nline_ini    !  In -> Line to start the search
    !!----                                                                 Out -> Actual line on Filevar
    !!----    integer,          intent(in)               :: nline_end    !  In -> Line to finish the search
    !!----    character(len=*), intent(out)              :: spgr_ha      ! Out -> Hall symbol
    !!----
    !!----    Obtaining the Hall symbol of the Space Group
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Cif_Hall(Filevar,Nline_Ini,Nline_End,Spgr_Ha)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in) :: filevar
       integer,          intent(in out)           :: nline_ini
       integer,          intent(in)               :: nline_end
       character(len=*), intent(out)              :: spgr_ha

       !---- Local variables ----!
       integer :: np1, np2

       spgr_ha=" "
       call Read_Key_StrVal(filevar,nline_ini,nline_end, &
                            "_symmetry_space_group_name_Hall",spgr_ha)
       if (len_trim(spgr_ha)==0) spgr_ha=adjustl(filevar(nline_ini+1))
       ! TR  feb. 2015 .(re-reading the same item with another name)
       if(len_trim(spgr_ha) == 0) then
        spgr_ha=" "
        call Read_Key_StrVal(filevar,nline_ini,nline_end, "_space_group_name_Hall",spgr_ha)
        if (len_trim(spgr_ha)==0) spgr_ha=adjustl(filevar(nline_ini+1))
       end if

       if (spgr_ha =="?" .or. spgr_ha=="#") then
          spgr_ha=" "
       else
          np1=index(spgr_ha,"'")
          np2=index(spgr_ha,"'",back=.true.)
          if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
             spgr_ha=spgr_ha(np1+1:np2-1)
          else
             np1=index(spgr_ha,'"')
             np2=index(spgr_ha,'"',back=.true.)
             if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
                spgr_ha=spgr_ha(np1+1:np2-1)
             else
                spgr_ha=" "
             end if
          end if
       end if

       return
    End Subroutine Read_Cif_Hall

    !!----
    !!---- Subroutine Read_Cif_Hm(Filevar,Nline_Ini,Nline_End,Spgr_Hm)
    !!----    character(len=*),  dimension(:), intent(in) :: filevar     !  In -> String vector
    !!----    integer,           intent(in out)           :: nline_ini   !  In -> Line to start the search
    !!----                                                                 Out -> Actual Line on Filevar
    !!----    integer,           intent(in)               :: nline_end   !  In -> Line to finish the search
    !!----    character(len=*),  intent(out)              :: spgr_hm     ! Out -> Hermann-Mauguin symbol
    !!----
    !!----    Obtaining the Herman-Mauguin symbol of Space Group
    !!----
    !!---- Update: March - 2010
    !!
    Subroutine Read_Cif_Hm(Filevar,Nline_Ini,Nline_End,Spgr_Hm)
       !---- Arguments ----!
       character(len=*),  dimension(:), intent(in) :: filevar
       integer,           intent(in out)           :: nline_ini
       integer,           intent(in)               :: nline_end
       character(len=*),  intent(out)              :: spgr_hm

       !---- Local variables ----!
       character(len=1) :: csym, csym2
       integer          :: np1, np2

       spgr_hm=" "
       np1=nline_ini
       call Read_Key_Str(filevar,nline_ini,nline_end, &
                            "_symmetry_space_group_name_H-M",spgr_hm)
       !if (len_trim(spgr_hm) ==0 ) spgr_hm=adjustl(filevar(nline_ini+1))
       !nline_ini=np1
       ! TR  feb. 2015 .(re-reading the same item with another name)
       if(len_trim(spgr_hm) == 0) then
        nline_ini=np1
        spgr_hm = " "
        call Read_Key_Str(filevar,nline_ini,nline_end, "_space_group_name_H-M_alt",spgr_hm)
        if (len_trim(spgr_hm) ==0 ) spgr_hm=adjustl(filevar(nline_ini+1))
       end if

       if (spgr_hm =="?" .or. spgr_hm=="#") then
          spgr_hm=" "
       else
          np1=index(spgr_hm,"'")
          np2=index(spgr_hm,"'",back=.true.)
          if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
             spgr_hm=spgr_hm(np1+1:np2-1)
          else
             np1=index(spgr_hm,'"')
             np2=index(spgr_hm,'"',back=.true.)
             if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
                spgr_hm=spgr_hm(np1+1:np2-1)
             else
                spgr_hm=" "
             end if
          end if
       end if

       !---- Adapting Nomenclature from ICSD to our model ----!
       np1=len_trim(spgr_hm)
       if (np1 > 0) then
          csym=u_case(spgr_hm(np1:np1))
          select case (csym)
             case("1")
                csym2=u_case(spgr_hm(np1-1:np1-1))
                if (csym2 == "Z" .or. csym2 =="S") then
                   spgr_hm=spgr_hm(:np1-2)//":1"
                end if

             case("S","Z")
                csym2=u_case(spgr_hm(np1-1:np1-1))
                select case (csym2)
                   case ("H")
                      spgr_hm=spgr_hm(:np1-2)
                   case ("R")
                      spgr_hm=spgr_hm(:np1-2)//":R"
                   case default
                      spgr_hm=spgr_hm(:np1-1)
                end select

             case("R")
                csym2=u_case(spgr_hm(np1-1:np1-1))
                if (csym2 == "H" ) then
                   spgr_hm=spgr_hm(:np1-2)
                else
                   spgr_hm=spgr_hm(:np1-1)//":R"
                end if
             case("H")
                spgr_hm=spgr_hm(:np1-1)
                csym2=u_case(spgr_hm(np1-1:np1-1))
                if(csym2 == ":") spgr_hm=spgr_hm(:np1-2)

          end select
       end if

       return
    End Subroutine Read_Cif_Hm

    !!----
    !!---- Subroutine Read_Cif_Lambda(Filevar,Nline_Ini,Nline_End,Lambda)
    !!----    character(len=*), dimension(:), intent(in) :: filevar      !  In -> String vector
    !!----    integer,           intent(in out)          :: nline_ini    !  In -> Line to start of search
    !!----                                                                  Out -> Actual line on Filevar
    !!----    integer,           intent(in)              :: nline_end    !  In -> Line to finish the search
    !!----    real(kind=cp),     intent(out)             :: lambda       !  Out -> lamda value
    !!----
    !!----    Radiation length
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Cif_Lambda(Filevar,Nline_Ini,Nline_End,Lambda)
       !---- Arguments ----!
       character(len=*),  dimension(:), intent(in) :: filevar
       integer,           intent(in out)           :: nline_ini
       integer,           intent(in)               :: nline_end
       real(kind=cp),     intent(out)              :: lambda

       !---- Local Variables ----!
       integer                    :: iv
       integer,dimension(1)       :: ivet
       real(kind=cp), dimension(1):: vet

       lambda=0.71073    ! Mo

       call read_key_value(filevar,nline_ini,nline_end, &
                           "_diffrn_radiation_wavelength",vet,ivet,iv)
       if (iv == 1) then
          lambda=vet(1)
       end if

       return
    End Subroutine Read_Cif_Lambda

    !!----
    !!---- Subroutine Read_Cif_Pressure(Filevar,Nline_Ini,Nline_End,P,SigP)
    !!----    character(len=*), dimension(:), intent(in) :: filevar      !  In -> String vector
    !!----    integer,           intent(in out)          :: nline_ini    !  In -> Line to start of search
    !!----                                                                  Out -> Actual line on Filevar
    !!----    integer,           intent(in)              :: nline_end    !  In -> Line to finish the search
    !!----    real(kind=cp),     intent(out)             :: P            !  Out -> Pressure (GPa)value
    !!----    real(kind=cp),     intent(out)             :: SigP         !  Out -> Sigma Pressure
    !!----
    !!----    Pressure and Sigma
    !!----
    !!---- Update: October - 2016
    !!
    Subroutine Read_Cif_Pressure(Filevar,Nline_Ini,Nline_End,P,SigP)
       !---- Arguments ----!
       character(len=*),  dimension(:), intent(in) :: filevar
       integer,           intent(in out)           :: nline_ini
       integer,           intent(in)               :: nline_end
       real(kind=cp),     intent(out)              :: p
       real(kind=cp),     intent(out)              :: sigp

       !---- Local Variables ----!
       integer                    :: iv
       real(kind=cp),dimension(1) :: vet1,vet2

       !> Default values
       p=0.0
       sigp=1.0e-5

       call read_key_valuestd(filevar,nline_ini,nline_end, &
                           "_diffrn_ambient_pressure",vet1,vet2,iv)
       if (iv == 1) then
          p=vet1(1)*1.0e-6
          sigp=vet2(1)*1.0e-6
       end if

       return
    End Subroutine Read_Cif_Pressure

    !!----
    !!---- Subroutine Read_Cif_Temp(Filevar,Nline_Ini,Nline_End,T,SigT)
    !!----    character(len=*), dimension(:), intent(in) :: filevar      !  In -> String vector
    !!----    integer,           intent(in out)          :: nline_ini    !  In -> Line to start of search
    !!----                                                                  Out -> Actual line on Filevar
    !!----    integer,           intent(in)              :: nline_end    !  In -> Line to finish the search
    !!----    real(kind=cp),     intent(out)             :: T            !  Out -> Temp (K) value
    !!----    real(kind=cp),     intent(out)             :: SigT         !  Out -> Sigma Temp
    !!----
    !!----    Temperature and Sigma
    !!----
    !!---- Update: October - 2016
    !!
    Subroutine Read_Cif_Temp(Filevar,Nline_Ini,Nline_End,T,SigT)
       !---- Arguments ----!
       character(len=*),  dimension(:), intent(in) :: filevar
       integer,           intent(in out)           :: nline_ini
       integer,           intent(in)               :: nline_end
       real(kind=cp),     intent(out)              :: T
       real(kind=cp),     intent(out)              :: sigT

       !---- Local Variables ----!
       integer                    :: iv
       real(kind=cp),dimension(1) :: vet1,vet2

       !> Default values
       t=298.0
       sigt=1.0

       call read_key_valuestd(filevar,nline_ini,nline_end, &
                           "_diffrn_ambient_temperature",vet1,vet2,iv)
       if (iv == 1) then
          t=vet1(1)
          sigt=vet2(1)
       end if

       return
    End Subroutine Read_Cif_Temp
    !!----
    !!---- Subroutine Read_Cif_Symm(Filevar,Nline_Ini,Nline_End,N_Oper,Oper_Symm)
    !!----    character(len=*), dimension(:), intent(in) :: filevar       !  In -> String vector
    !!----    integer,          intent(in out)           :: nline_ini     !  In -> Line to start the search
    !!----                                                                  Out -> Actual line on Filevar
    !!----    integer,          intent(in)               :: nline_end     !  In -> Line to finish the search
    !!----    integer,          intent(out)              :: n_oper        ! Out -> Number of Operators
    !!----    character(len=*), dimension(:),intent(out) :: oper_symm     ! Out -> Vector with Symmetry Operators
    !!----
    !!----    Obtaining Symmetry Operators from Cif file
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Cif_Symm(Filevar,Nline_Ini,Nline_End,N_Oper,Oper_Symm)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in) :: filevar
       integer,          intent(in out)           :: nline_ini
       integer,          intent(in)               :: nline_end
       integer,          intent(out)              :: n_oper
       character(len=*), dimension(:),intent(out) :: oper_symm

       !---- Local variables ----!
       character(len=len(filevar(1))) :: string
       integer                        :: i,np1,np2

       n_oper=0
       oper_symm=" "
       np1=nline_ini
       call Read_Key_StrVal(filevar,nline_ini,nline_end, &
                            "_symmetry_equiv_pos_as_xyz",string)
       !nline_ini=np1
       ! TR  feb. 2015 .(re-reading the same item with another name)
       !if(len_trim(string) == 0) then
       if(nline_ini == 1) then   ! TR june 2016
        nline_ini=np1
        call Read_Key_StrVal(filevar,nline_ini,nline_end, "_space_group_symop_operation_xyz",string)
       end if

       if (len_trim(string) /=0) then
          string=adjustl(string)

          if (string(1:1) /="#" .and. string(1:1) /= "?") then      ! Comentario
             np1=index(string,"'")
             np2=index(string,"'",back=.true.)
             if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
                n_oper=n_oper+1
                oper_symm(n_oper)=string(np1+1:np2-1)
             else
                np1=index(string,'"')
                np2=index(string,'"',back=.true.)
                if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
                   n_oper=n_oper+1
                   oper_symm(n_oper)=string(np1+1:np2-1)
                end if
             end if
          end if
       end if

       do i=nline_ini+1,nline_end
          string=adjustl(filevar(i))
          if (len_trim(string) /=0) then
             if (string(1:1) /="#" .and. string(1:1) /= "?") then      ! Comentario o Vacio
                np1=index(string,"'")
                np2=index(string,"'",back=.true.)
                if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
                   n_oper=n_oper+1
                   oper_symm(n_oper)=string(np1+1:np2-1)
                else
                   np1=index(string,'"')
                   np2=index(string,'"',back=.true.)
                   if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
                      n_oper=n_oper+1
                      oper_symm(n_oper)=string(np1+1:np2-1)
                   end if
                end if
             end if
          else
             nline_ini=i+1
             exit
          end if
       end do

       return
    End Subroutine Read_Cif_Symm

    !!----
    !!---- Subroutine Read_Cif_Title(Filevar,Nline_Ini,Nline_End,Title)
    !!----    character(len=*),  dimension(:), intent(in) :: filevar      !  In -> String vector
    !!----    integer,           intent(in out)           :: nline_ini    !  In -> Line to start the search
    !!----                                                                  Out -> Actual line on Filevar
    !!----    integer,           intent(in)               :: nline_end    !  In -> Line to finish the search
    !!----    character(len=*),  intent(out)              :: title        ! Out -> Title string
    !!----
    !!----    Obtaining Title from Cif file
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Cif_Title(Filevar,Nline_Ini,Nline_End,title)
       !---- Arguments ----!
       character(len=*),  dimension(:), intent(in) :: filevar
       integer,           intent(in out)           :: nline_ini
       integer,           intent(in)               :: nline_end
       character(len=*),  intent(out)              :: title

       !---- Local variables ----!
       integer :: np, np1, np2

       title=" "
       call Read_Key_StrVal(filevar,nline_ini,nline_end, &
                            "_publ_section_title",title)

       if (len_trim(title) ==0 ) title=adjustl(filevar(nline_ini+1))
       if (title =="; ?" .or. title=="#") then
          title=" "
       else
          np=len_trim(title)
          if (np <= 3) title=adjustl(filevar(nline_ini+2))
          np1=index(title,"'")
          np2=index(title,"'",back=.true.)
          if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
             title=title(np1+1:np2-1)
          else
             np1=index(title,'"')
             np2=index(title,'"',back=.true.)
             if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
                title=title(np1+1:np2-1)
             end if
          end if
       end if

       return
    End Subroutine Read_Cif_Title

    !!----
    !!---- Subroutine Read_Cif_Z(Filevar,Nline_Ini,Nline_End,Z)
    !!----    character(len=*), dimension(:), intent(in) :: filevar     !  In -> String vector
    !!----    integer,           intent(in out)          :: nline_ini   !  In -> Line to start the search
    !!----                                                                Out -> Actual line on Filevar
    !!----    integer,           intent(in)              :: nline_end   !  In -> Line to finish the search
    !!----    integer,           intent(out)             :: Z           ! Out -> Z value
    !!----
    !!----    Unit formula from Cif file
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Cif_Z(filevar,nline_ini,nline_end,z)
       !---- Arguments ----!
       character(len=*),  dimension(:), intent(in) :: filevar
       integer,           intent(in out)           :: nline_ini
       integer,           intent(in)               :: nline_end
       integer,           intent(out)              :: z

       !---- Local Variables ----!
       integer                     :: iv
       integer,dimension(1)        :: ivet
       real(kind=cp), dimension(1) :: vet

       z=0
       call read_key_value(filevar,nline_ini,nline_end, &
                           "_cell_formula_units_Z",vet,ivet,iv)
       if (iv == 1) then
          z=ivet(1)
       end if

       return
    End Subroutine Read_Cif_Z

    !!----
    !!---- Subroutine Read_File_Atom(Filevar,Nline_Ini,Nline_End,Atomos)
    !!----    character(len=*),dimension(:), intent(in)       :: filevar     !  In -> String vector
    !!----    integer,                       intent(in)       :: nline_ini   !  In -> Line to start the search
    !!----                                                                     Out -> Actual line on Filevar
    !!----    integer,                       intent(in)       :: nline_end   !  In -> Line to finish the search
    !!----    type (atom_list_type),        intent(out)      :: Atomos      ! Out -> Atom list
    !!----           or
    !!----    type (Point_list_Type),        intent(out)      :: Atomos      ! Out -> point list
    !!----
    !!----     Subroutine to read an atom (or point) list from a file. Atomos should be previously allocated.
    !!----     Control of error is present.
    !!----
    !!---- Update: June - 2005
    !!

    !!--++
    !!--++ Subroutine Read_File_Atomlist(Filevar,Nline_Ini,Nline_End,Atomos)
    !!--++    character(len=*),dimension(:), intent(in)       :: filevar     !  In -> String vector
    !!--++    integer,                       intent(in)       :: nline_ini   !  In -> Line to start the search
    !!--++                                                                     Out -> Actual line on Filevar
    !!--++    integer,                       intent(in)       :: nline_end   !  In -> Line to finish the search
    !!--++    type (atom_list_type),        intent(out)       :: Atomos      ! Out -> Atom list
    !!--++
    !!--++     Subroutine to read an atom list from a file. Atomos should be previously allocated.
    !!--++     Control of error is present
    !!--++
    !!--++ Update: June - 2005
    !!
    Subroutine Read_File_Atomlist(filevar,nline_ini,nline_end,Atomos)
       !---- Arguments ----!
       character(len=*), dimension(:),   intent(in)      :: filevar
       integer,                          intent(in out)  :: nline_ini
       integer,                          intent(in)      :: nline_end
       type (atom_list_type),            intent(in out)  :: Atomos

       !---- Local variables -----!
       character(len=len(filevar(1))) :: line
       character(len=4)               :: dire
       integer                        :: i,na
       type (Atom_Type)               :: Atomo

       !---- Initial Values ----!
       na=0
       do i=nline_ini,nline_end
          dire=adjustl(u_case(filevar(i)(1:4)))
          if (dire /= "ATOM") cycle
          line=adjustl(filevar(i))
          call read_atom(line,atomo)
          if (err_form) cycle

          !---- Trial to read anisotropic thermal parameters ----!
          if( i < size(filevar) ) then
           line=adjustl(filevar(i+1))
           select case (u_case(line(1:4)))
             case ("U_IJ")
                call read_uvals(line,atomo, "u_ij")
             case ("B_IJ")
                call read_uvals(line,atomo, "b_ij")
             case ("BETA")
                call read_uvals(line,atomo, "beta")
           end select
           if (err_form) cycle
          end if
          na=na+1
          Atomos%atom(na)=atomo
       end do

       Atomos%natoms=na

       return
    End Subroutine Read_File_Atomlist

    !!----
    !!---- Subroutine Read_File_PointList(Filevar,Nline_Ini,Nline_End,Atomos)
    !!----    character(len=*),dimension(:), intent(in)       :: filevar     !  In -> String vector
    !!----    integer,                       intent(in)       :: nline_ini   !  In -> Line to start the search
    !!----                                                                     Out -> Actual line on Filevar
    !!----    integer,                       intent(in)       :: nline_end   !  In -> Line to finish the search
    !!----    type (Point_List_Type),        intent(out)      :: Atomos      ! Out -> point list
    !!----
    !!----     Subroutine to read an point list from a file. Atomos should be previously allocated.
    !!----     Control of error is present
    !!----
    !!---- Update: June - 2005
    !!
    Subroutine Read_File_PointList(filevar,nline_ini,nline_end,Atomos)
       !---- Arguments ----!
       character(len=*), dimension(:),   intent(in)      :: filevar
       integer,                          intent(in out)  :: nline_ini
       integer,                          intent(in)      :: nline_end
       type (Point_List_Type),           intent(in out)  :: Atomos

       !---- Local variables -----!
       character(len=len(filevar(1))) :: line
       character(len=4)               :: dire
       integer                        :: i,na
       type (Atom_Type)               :: Atomo

       !---- Initial Values ----!
       na=0

       do i=nline_ini,nline_end
          dire=adjustl(u_case(filevar(i)(1:4)))
          if (dire /= "ATOM") cycle
          line=adjustl(filevar(i))
          call read_atom(line,atomo)
          if (err_form) cycle
          na=na+1
          Atomos%x(:,na) =atomo%x(:)
          Atomos%p(na)   = 0
          Atomos%nam(na) = atomo%lab
       end do

       Atomos%np=na

       return
    End Subroutine Read_File_PointList

    !!----
    !!---- Subroutine Read_File_Cell(Filevar,Nline_Ini,Nline_End,Celda)
    !!----    character(len=*), dimension(:), intent(in) :: filevar      !  In -> String Vector
    !!----    integer,           intent(in out)          :: nline_ini    !  In -> Line to start the search
    !!----                                                                 Out -> Atual line on Filevar
    !!----    integer,           intent(in)              :: nline_end    !  In -> line to finish the search
    !!----
    !!----    real(kind=cp),dimension(6), intent (out)   :: Celda        ! Out -> Cell variable
    !!----                            or
    !!----    type (Crystal_Cell_Type), intent (out)     :: Celda        ! Out -> Cell variable
    !!----
    !!----    Read Cell Parameters from file. Control error is present
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Read_File_Cellc(Filevar,Nline_Ini,Nline_End,Celda)
    !!--++    character(len=*), dimension(:), intent(in) :: filevar      !  In -> String Vector
    !!--++    integer,           intent(in out)          :: nline_ini    !  In -> Line to start the search
    !!--++                                                                 Out -> Atual line on Filevar
    !!--++    integer,           intent(in)              :: nline_end    !  In -> line to finish the search
    !!--++    real(kind=cp),dimension(6), intent (out)   :: Celda        ! Out -> Cell variable
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Read Cell Parameters from file. Control error is present
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Read_File_Cellc(filevar,nline_ini,nline_end,Celda)
       !---- Arguments ----!
       character(len=*),  dimension(:), intent(in)     :: filevar
       integer,                         intent(in)     :: nline_ini
       integer,                         intent(in)     :: nline_end
       real(kind=cp),dimension(6),      intent(out)    :: Celda

       !---- Local Variables ----!
       integer                     :: iv, i,j
       integer, dimension(6)       :: ivet
       real(kind=cp), dimension(6) :: vet

       !---- Valores iniciales ----!
       call init_err_form()

       i=nline_ini
       j=nline_end

       !---- Celda ----!
       call read_key_value(filevar,i,j,"cell",vet,ivet,iv)
       if (iv /=6) then
          err_form=.true.
          ERR_Form_Mess=" Bad Cell Parameters..."
          return
       else
          celda=vet(:)
       end if

       return
    End Subroutine Read_File_Cellc

    !!--++
    !!--++ Subroutine Read_File_Cellt(Filevar,Nline_Ini,Nline_End,Celda,CFrame)
    !!--++    character(len=*),  dimension(:), intent(in)     :: filevar     !  In -> String Vector
    !!--++    integer,                         intent(in)     :: nline_ini   !  In -> Line to start the search
    !!--++    integer,                         intent(in)     :: nline_end   !  In -> line to finish the search
    !!--++    type (Crystal_Cell_Type),        intent(out)    :: Celda       ! Out -> Cell structure
    !!--++    character(len=*),  optional,     intent(in)     :: CFrame      !  Cartesian Frame "A" or "C" (if absent -> "A")
    !!--++          ! Out -> Cell variable
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Read Cell Parameters from file. Control error is present
    !!--++    The object Celda is constructed just after reading the cell parameters.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Read_File_Cellt(filevar,nline_ini,nline_end,Celda,CFrame)
       !---- Arguments ----!
       character(len=*),  dimension(:), intent(in)     :: filevar
       integer,                         intent(in)     :: nline_ini
       integer,                         intent(in)     :: nline_end
       type (Crystal_Cell_Type),        intent(out)    :: Celda
       character(len=*),  optional,     intent(in)     :: CFrame

       !---- Local Variables ----!
       integer                     :: iv, i,j
       real(kind=cp), dimension(6) :: vet1,vet2

       !---- Valores iniciales ----!
       call init_err_form()

       i=nline_ini
       j=nline_end

       !---- Celda ----!

       call read_key_valueSTD(filevar,i,j,"cell",vet1,vet2,iv)
       if (iv /=6) then
          err_form=.true.
          ERR_Form_Mess=" Bad Cell Parameters..."
          return
       end if
       if(present(CFrame)) then
         call Set_Crystal_Cell(vet1(1:3),vet1(4:6),Celda,CFrame)
       else
         call Set_Crystal_Cell(vet1(1:3),vet1(4:6),Celda,"A")
       end if
       celda%cell_std=vet2(1:3)
       celda%ang_std=vet2(4:6)

       return
    End Subroutine Read_File_Cellt

    !!----
    !!---- Subroutine Read_File_lambda(Filevar,Nline_Ini,Nline_End,v1,v2,v3)
    !!----    character(len=*), dimension(:), intent(in)     :: filevar   !  In -> String Vector
    !!----    integer,                        intent(in out) :: nline_ini !  In -> Line to start the search
    !!----                                                                  Out -> Atual line on Filevar
    !!----    integer,                        intent(in)     :: nline_end !  In -> line to finish the search
    !!----    real(kind=cp),                  intent(   out) :: v1,v2,v3  ! Out -> Lambda1,lambda2,ratio
    !!----
    !!----    Read wavelengths and ratio.
    !!----    If no value is read, Lambda1=Lambda2=1.54056 Angstroms, ratio=0.0
    !!----    If only one value is read Lambda1=Lambda2=v1, ratio=0
    !!----    If only two values iare read Lambda1=v1, Lambda2=v2, ratio=0.5
    !!----    In other cases Lambda1=v1, Lambda2=v2, ratio=v3
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_File_Lambda(Filevar,Nline_Ini,Nline_End,v1,v2,v3)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in)     :: filevar
       integer,                        intent(in out) :: nline_ini
       integer,                        intent(in)     :: nline_end
       real(kind=cp),                  intent(   out) :: v1,v2,v3

       !---- Local Variables ----!
       integer                    :: iv, i,j
       integer, dimension(3)      :: ivet
       real(kind=cp), dimension(3):: vet

       !---- Valores iniciales ----!
       call init_err_form()

       i=nline_ini
       j=nline_end

       v3=0.0
       v1=1.54056
       !---- Read Lambda ----!
       call read_key_value(filevar,i,j,"wave",vet,ivet,iv)
       if      (iv == 0) then
         v2=1.54056
       else if (iv == 1) then
         v1=vet(1)
         v2=vet(1)
       else if (iv == 2) then
         v1=vet(1)
         v2=vet(2)
         v3=0.5
       else if (iv == 3) then
         v1=vet(1)
         v2=vet(2)
         v3=vet(3)
       end if

       return
    End Subroutine Read_File_Lambda

    !!----
    !!---- Subroutine Read_File_RngSintL(Filevar,Nline_Ini,Nline_End,v1,v2)
    !!----    character(len=*), dimension(:), intent(in)     :: filevar   !  In -> String Vector
    !!----    integer,                        intent(in out) :: nline_ini !  In -> Line to start the search
    !!----                                                                  Out -> Atual line on Filevar
    !!----    integer,                        intent(in)     :: nline_end !  In -> line to finish the search
    !!----    real(kind=cp),                  intent(   out) :: v1,v2     ! Out -> Interval [v1,v2] in sinT/Lambda
    !!----
    !!----    Read range for sintheta/lambda.
    !!----    If only one value is read v1=0 and v2= read value
    !!----    If the keyword RNGSL is not given in the file, the default
    !!----    values are v1=0.0, v2=1.0
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_File_RngSintL(Filevar,Nline_Ini,Nline_End,v1,v2)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in)     :: filevar
       integer,                        intent(in out) :: nline_ini
       integer,                        intent(in)     :: nline_end
       real(kind=cp),                  intent(   out) :: v1,v2

       !---- Local Variables ----!
       integer                     :: iv, i,j
       integer,       dimension(2) :: ivet
       real(kind=cp), dimension(2) :: vet

       !---- Valores iniciales ----!
       call init_err_form()

       i=nline_ini
       j=nline_end

       !---- Range in sinTheta/Lambda ----!
       call read_key_value(filevar,i,j,"rngsl",vet,ivet,iv)
       if      (iv == 0) then
         v1=0.0
         v2=1.0
       else if (iv == 1) then
         v1=0.0
         v2=vet(1)
       else if (iv == 2) then
         v1=vet(1)
         v2=vet(2)
       end if

       return
    End Subroutine Read_File_RngSintL

    !!----
    !!---- Subroutine Read_File_Spg (Filevar,Nline_Ini,Nline_End,Spg,Sub)
    !!----    character(len=*),  dimension(:), intent(in) :: filevar       !  In -> String vector
    !!----    integer,           intent(in)               :: nline_ini     !  In -> Line to start the search
    !!----                                                                   Out -> Actual line on Filevar
    !!----    integer,           intent(in)               :: nline_end     !  In -> Line to Finish the search
    !!----    character(len=*),  intent(out)              :: spg           ! Out -> Space Group symbol
    !!----    character(len=*),  intent(in ),optional     :: sub           ! in  -> The space sroup symbol is a subgroup
    !!----                                                                 !        of an already given space group
    !!----    Reads the cards "SPGR", "SPACEG" or "SUBG" in filvar. Control of error is present
    !!----
    !!---- Update: February - 2011
    !!
    Subroutine Read_File_Spg(filevar,nline_ini,nline_end,Spg,sub)
       !---- Arguments ----!
       character(len=*),  dimension(:), intent(in) :: filevar   ! Variable
       integer,           intent(in)               :: nline_ini
       integer,           intent(in)               :: nline_end
       character(len=*),  intent(out)              :: spg
       character(len=*),  intent(in),  optional    :: sub

       !--Local variables--!
       integer  :: i

       call init_err_form()
       i=nline_ini
       if(present(sub)) then
         call Read_Key_StrVal(filevar,i,nline_end, "subg",spg)
       else
         call Read_Key_StrVal(filevar,i,nline_end, "spgr",spg)
       end if
       if (len_trim(spg) == 0 ) then
         call Read_Key_StrVal(filevar,i,nline_end, "spaceg",spg)
         if (len_trim(spg) == 0 ) then
            call Read_Key_StrVal(filevar,i,nline_end, "shub",spg)
            if (len_trim(spg) == 0 ) then
              err_form=.true.
              ERR_Form_Mess=" Problems reading the Space Group symbol/number"
              return
            end if
         end if
       end if
       return
    End Subroutine Read_File_Spg

    !!----
    !!---- Read_File_Transf(Filevar,Nline_Ini,Nline_End,Transf,Orig)
    !!----    character(len=*), dimension(:), intent(in)     :: filevar      !  In -> String Vector
    !!----    integer,                        intent(in out) :: nline_ini    !  In -> Line to start the search
    !!----                                                                     Out -> Atual line on Filevar
    !!----    integer,                        intent(in)     :: nline_end    !  In -> line to finish the search
    !!----    real(kind=cp),dimension(3,3),   intent(out)    :: transf       ! Out -> Cell variable
    !!----    real(kind=cp),dimension(3  ),   intent(out)    :: orig
    !!----
    !!----    Read transformation matrix for changing the space group or cell setting.
    !!----    First the matrix M is read row by row and then the origin in the old setting
    !!----    is finally read. A single line with 12 real numbers should be given.
    !!--<<
    !!----    e.g.: TRANS  m11 m12 m13  m21 m22 m33  m31 m32 m33   o1 o2 o3
    !!----
    !!----    That means       a'=m11 a + m12 b + m13 c
    !!----                     b'=m21 a + m22 b + m23 c
    !!----                     c'=m31 a + m32 b + m33 c
    !!----
    !!----                     X' = inv(Mt) (X-O)
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_File_transf(filevar,nline_ini,nline_end,trans,orig)
       !---- Arguments ----!
       character(len=*),  dimension(:), intent(in)     :: filevar
       integer,                         intent(in)     :: nline_ini
       integer,                         intent(in)     :: nline_end
       real(kind=cp),dimension(3,3),    intent(out)    :: trans
       real(kind=cp),dimension(3  ),    intent(out)    :: orig

       !---- Local Variables ----!
       integer                      :: iv, i,j
       integer,       dimension(12) :: ivet
       real(kind=cp), dimension(12) :: vet
       character(len=80)            :: transf_key

       !---- Initial values ----!
       call init_err_form()

       i=nline_ini
       j=nline_end

       !---- transformation matrix ----!
       call read_key_value(filevar,i,j,"trans",vet,ivet,iv,"#",transf_key)
       if (iv /= 12 .or. err_string) then
          !Try to read the transformation from transf_key
          if(len_trim(transf_key) /= 0) then
            call Get_Transf(transf_key,trans,orig)
            if(err_string) then
               err_form=.true.
               ERR_Form_Mess=" Bad matrix/origin setting in string: "//trim(transf_key)//" -> "//trim(Err_String_Mess)
               return
            end if
          else
               err_form=.true.
               ERR_Form_Mess=" Bad matrix/origin setting..."
               return
          end if

       else
          trans(1,1:3)=vet(1:3)
          trans(2,1:3)=vet(4:6)
          trans(3,1:3)=vet(7:9)
          orig(1:3) = vet(10:12)
       end if

       return
    End Subroutine Read_File_transf

    !!----
    !!---- Subroutine Read_Shx_Atom(Filevar,Nline_Ini,Nline_End,N_Fvar,Fvar,Elem_Type,Celda,Atm_List)
    !!----    character(len=*), dimension(:), intent(in)      :: filevar        !  In -> String vector
    !!----    integer,                        intent(in out)  :: nline_ini      !  In -> Line to start the search
    !!----                                                                         Out -> Actual line on Filevar
    !!----    integer,                        intent(in)      :: nline_end      !  In -> Line to finish the search
    !!----    integer,                        intent(in)      :: n_fvar         !  In -> Number of parameters on FVAR
    !!----    real(kind=cp), dimension(:),    intent(in)      :: fvar           !  In -> Values for FVAR
    !!----    character(len=*), dimension(:), intent(in)      :: elem_type      !  In -> type of elements
    !!----    type (Crystal_Cell_Type),       intent(in)      :: Celda          !  In -> Cell type variable
    !!----    type (Atom_list_type),          intent(out)     :: Atm_List       ! Out -> number of atoms
    !!----         ! Out -> Atom List
    !!----
    !!----    Obtaining Atoms parameters from Shelx file (.ins or .res)
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Shx_Atom(filevar,nline_ini,nline_end,n_fvar,fvar,elem_type,celda,Atm_List)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in)      :: filevar
       integer,                        intent(in out)  :: nline_ini
       integer,                        intent(in)      :: nline_end
       integer,                        intent(in)      :: n_fvar
       real(kind=cp), dimension(:),    intent(in)      :: fvar
       character(len=*), dimension(:), intent(in)      :: elem_type
       type (Crystal_Cell_Type),       intent(in)      :: Celda
       type (Atom_list_type),          intent(out)     :: Atm_List

       !---- Local Variables ----!
       character(len=80)               :: string
       character(len=30),dimension(15) :: label
       character(len=2)                :: el
       integer                         :: i, nc, iv
       integer                         :: j, n_atom
       integer, dimension(15)          :: ivet
       real(kind=cp)                   :: x, p, u
       real(kind=cp), dimension(15)    :: vet
       type(atom_list_type)            :: Atm

       call allocate_atom_list(nline_end-nline_ini+1,Atm)
       n_atom=0

       do i=nline_ini,nline_end
          string=filevar(i)
          if (len_trim(string) == 0) cycle
          call getword(string,label,nc)
          select case (nc)
             case (5) ! Atomname Sfac X Y Z
                call getnum(label(2),vet,ivet,iv)   ! Is Sfac integer?
                if (iv /= 1) cycle
                call getnum(label(3),vet,ivet,iv)   ! Is X real?
                if (iv /= 1) cycle
                call getnum(label(4),vet,ivet,iv)   ! Is Y real?
                if (iv /= 1) cycle
                call getnum(label(5),vet,ivet,iv)   ! Is Z real?
                if (iv /= 1) cycle

                n_atom=n_atom+1
                atm%atom(n_atom)%lab=label(1)(1:4)
                call getnum(label(2),vet,ivet,iv)
                el=elem_type(ivet(1))
                atm%atom(n_atom)%chemSymb=U_case(el(1:1))//L_case(el(2:2))
                call getnum(label(3),vet,ivet,iv)
                atm%atom(n_atom)%x(1)=vet(1)
                call getnum(label(4),vet,ivet,iv)
                atm%atom(n_atom)%x(2)=vet(1)
                call getnum(label(5),vet,ivet,iv)
                atm%atom(n_atom)%x(3)=vet(1)
                atm%atom(n_atom)%utype="u_ij"

             case (6) ! Atomname Sfac X Y Z Occ
                call getnum(label(2),vet,ivet,iv)   ! Is Sfac integer?
                if (iv /= 1) cycle
                call getnum(label(3),vet,ivet,iv)   ! Is X real?
                if (iv /= 1) cycle
                call getnum(label(4),vet,ivet,iv)   ! Is Y real?
                if (iv /= 1) cycle
                call getnum(label(5),vet,ivet,iv)   ! Is Z real?
                if (iv /= 1) cycle
                call getnum(label(6),vet,ivet,iv)   ! Is Occ real?
                if (iv /= 1) cycle

                n_atom=n_atom+1
                atm%atom(n_atom)%lab=trim(label(1))
                call getnum(label(2),vet,ivet,iv)
                el=elem_type(ivet(1))
                atm%atom(n_atom)%chemSymb=U_case(el(1:1))//L_case(el(2:2))
                call getnum(label(3),vet,ivet,iv)
                atm%atom(n_atom)%x(1)=vet(1)
                call getnum(label(4),vet,ivet,iv)
                atm%atom(n_atom)%x(2)=vet(1)
                call getnum(label(5),vet,ivet,iv)
                atm%atom(n_atom)%x(3)=vet(1)
                call getnum(label(6),vet,ivet,iv)
                atm%atom(n_atom)%occ=vet(1)
                atm%atom(n_atom)%utype="u_ij"

             case (7,8) ! Atomname Sfac X Y Z Occ Uiso   (TR: item 8 can be electronic density created by SHELXS)
                call getnum(label(2),vet,ivet,iv)   ! Is Sfac integer?
                if (iv /= 1) cycle
                call getnum(label(3),vet,ivet,iv)   ! Is X real?
                if (iv /= 1) cycle
                call getnum(label(4),vet,ivet,iv)   ! Is Y real?
                if (iv /= 1) cycle
                call getnum(label(5),vet,ivet,iv)   ! Is Z real?
                if (iv /= 1) cycle
                call getnum(label(6),vet,ivet,iv)   ! Is Occ real?
                if (iv /= 1) cycle
                call getnum(label(7),vet,ivet,iv)   ! Is Uiso real?
                if (iv /= 1) cycle

                n_atom=n_atom+1
                atm%atom(n_atom)%lab=trim(label(1))
                call getnum(label(2),vet,ivet,iv)
                el=elem_type(ivet(1))
                atm%atom(n_atom)%chemSymb=U_case(el(1:1))//L_case(el(2:2))
                call getnum(label(3),vet,ivet,iv)
                atm%atom(n_atom)%x(1)=vet(1)
                call getnum(label(4),vet,ivet,iv)
                atm%atom(n_atom)%x(2)=vet(1)
                call getnum(label(5),vet,ivet,iv)
                atm%atom(n_atom)%x(3)=vet(1)
                call getnum(label(6),vet,ivet,iv)
                atm%atom(n_atom)%occ=vet(1)
                call getnum(label(7),vet,ivet,iv)
                atm%atom(n_atom)%ueq=vet(1)
                atm%atom(n_atom)%utype="u_ij"
                atm%atom(n_atom)%thtype="isotr"

          case (9) ! Atomname Sfac X Y Z Occ U11 U22 = U33 U23 U13 U12
                call getnum(label(2),vet,ivet,iv)   ! Is Sfac integer?
                if (iv /= 1) cycle
                call getnum(label(3),vet,ivet,iv)   ! Is X real?
                if (iv /= 1) cycle
                call getnum(label(4),vet,ivet,iv)   ! Is Y real?
                if (iv /= 1) cycle
                call getnum(label(5),vet,ivet,iv)   ! Is Z real?
                if (iv /= 1) cycle
                call getnum(label(6),vet,ivet,iv)   ! Is Occ real?
                if (iv /= 1) cycle
                call getnum(label(7),vet,ivet,iv)   ! Is U11 real?
                if (iv /= 1) cycle
                call getnum(label(8),vet,ivet,iv)   ! Is U22 real?
                if (iv /= 1) cycle
                call getnum(filevar(i+1),vet,ivet,iv) ! Are U33 U23 U13 U12?
                if (iv /= 4) cycle

                n_atom=n_atom+1
                atm%atom(n_atom)%lab=trim(label(1))
                call getnum(label(2),vet,ivet,iv)
                el=elem_type(ivet(1))
                atm%atom(n_atom)%chemSymb=U_case(el(1:1))//L_case(el(2:2))
                call getnum(label(3),vet,ivet,iv)
                atm%atom(n_atom)%x(1)=vet(1)
                call getnum(label(4),vet,ivet,iv)
                atm%atom(n_atom)%x(2)=vet(1)
                call getnum(label(5),vet,ivet,iv)
                atm%atom(n_atom)%x(3)=vet(1)
                call getnum(label(6),vet,ivet,iv)
                atm%atom(n_atom)%occ=vet(1)
                !---- U11 U22 U33 U12 U13 U23 Order ----!
                call getnum(label(7),vet,ivet,iv)
                atm%atom(n_atom)%u(1)=vet(1)
                call getnum(label(8),vet,ivet,iv)
                atm%atom(n_atom)%u(2)=vet(1)
                call getnum(filevar(i+1),vet,ivet,iv)
                atm%atom(n_atom)%u(3)=vet(1)
                atm%atom(n_atom)%u(4)=vet(4)
                atm%atom(n_atom)%u(5)=vet(3)
                atm%atom(n_atom)%u(6)=vet(2)
                atm%atom(n_atom)%utype="u_ij"
                atm%atom(n_atom)%thtype="aniso"
             case default
                cycle
          end select
       end do

       !---- Adjusting ... ----!
       call allocate_atom_list(n_atom,Atm_list)
       do i=1,n_atom
          atm_list%atom(i)=atm%atom(i)
       end do
       call Deallocate_atom_list(atm)

       !---- Tratamiento de Datos del Shelx ----!
       do i=1,n_atom
          !---- coordinates ----!
          if (atm_list%atom(i)%x(1) >= 10.0) atm_list%atom(i)%x(1)=atm_list%atom(i)%x(1)-10.0
          if (atm_list%atom(i)%x(2) >= 10.0) atm_list%atom(i)%x(2)=atm_list%atom(i)%x(2)-10.0
          if (atm_list%atom(i)%x(3) >= 10.0) atm_list%atom(i)%x(3)=atm_list%atom(i)%x(3)-10.0

          !---- ocupancy ----!
          if (abs(atm_list%atom(i)%occ)  > 10.0) then
             x=atm_list%atom(i)%occ
             if (x > 10.0) then
                atm_list%atom(i)%occ=x-10.0
             else
                x=abs(atm_list%atom(i)%occ)
                do j=2,n_fvar
                   if (x > 10.0*real(j) .and. x < 10.0*real(j+1)) then
                      p=x-10.0*real(j)
                      if (atm_list%atom(i)%occ > 0.0) then
                         atm_list%atom(i)%occ=p*fvar(j)
                      else
                         atm_list%atom(i)%occ=p*(fvar(j)-1.0)
                      end if
                   end if
                end do
             end if
          end if

          !---- Thermal factors ----!
          if (atm_list%atom(i)%thtype == "aniso") then
             atm_list%atom(i)%ueq=U_Equiv(celda,atm_list%atom(i)%u(1:6))  ! Uequi
             atm_list%atom(i)%biso= atm_list%atom(i)%ueq*78.95683521
          else
             if (atm_list%atom(i)%ueq < 0.0) then
                u=-atm_list%atom(i)%ueq
                if (u <= 5.0 .and. u >= 0.5) then
                   do j=i-1,1,-1
                      if (atm_list%atom(j)%ChemSymb == "H " .or. atm_list%atom(j)%ChemSymb == "h " ) cycle
                      atm_list%atom(i)%ueq=u*U_Equiv(celda,atm_list%atom(j)%u(1:6))  ! Uequi
                      atm_list%atom(i)%biso= atm_list%atom(i)%ueq*78.95683521
                   end do
                end if
             end if
          end if

       end do

       return
    End Subroutine Read_Shx_Atom

    !!----
    !!---- Subroutine Read_Shx_Cell(Filevar,Nline_Ini,Nline_End,Celda,Stdcelda,Lambda,Z)
    !!----    character(len=*), dimension(:), intent(in)     :: filevar       !  In -> String vector
    !!----    integer,                        intent(in out) :: nline_ini     !  In -> Line to start the search
    !!----                                                                      Out -> Actual line on Filevar
    !!----    integer,                        intent(in)     :: nline_end     !  In -> Line to finish the search
    !!----    real(kind=cp),dimension(6),     intent(out)    :: celda         ! Out -> Cell Parameters
    !!----    real(kind=cp),dimension(6),     intent(out)    :: Stdcelda      ! Out -> Std Cell Parameters
    !!----    real(kind=cp),                  intent(out)    :: lambda        ! Out -> Lambda
    !!----    integer,                        intent(out)    :: Z             ! Out -> Z
    !!----
    !!----    Obtaining Cell Parameter from Shelx file
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Shx_Cell(filevar,nline_ini,nline_end,Celda,StdCelda,lambda,z)
       !---- Arguments ----!
       character(len=*), dimension(:),     intent(in)     :: filevar
       integer,                            intent(in out) :: nline_ini
       integer,                            intent(in)     :: nline_end
       real(kind=cp),dimension(6),         intent(out)    :: Celda
       real(kind=cp),dimension(6),optional,intent(out)    :: StdCelda
       real(kind=cp),             optional,intent(out)    :: lambda
       integer,          optional,         intent(out)    :: z

       !---- Local Variables ----!
       integer                      :: iv,z_shx
       integer, dimension(10)       :: ivet
       real(kind=cp), dimension(10) :: vet
       real(kind=cp)                :: lambda_shx
       real(kind=cp),dimension(6)   :: std

       !---- Valores iniciales ----!
       celda=0.0
       if (present(stdcelda)) stdcelda=0.0
       if (present(Lambda))   lambda=0.0
       if (present(z))        z=0

       !---- CELL ----!
       call read_key_value(filevar,nline_ini,nline_end,"CELL",vet,ivet,iv)
       if (iv == 7) then
          lambda_shx = vet(1)
          celda      = vet(2:7)
       end if

       !---- Z, STD ----!
       call read_key_value(filevar,nline_ini,nline_end,"ZERR",vet,ivet,iv)
       if (iv == 7) then
          z_shx= ivet(1)
          std  = vet(2:7)
       end if

       if (present(stdcelda)) stdcelda=std
       if (present(lambda)) lambda=lambda_shx
       if (present(z)) z=z_shx

       return
    End Subroutine Read_Shx_Cell

    !!----
    !!---- Subroutine Read_Shx_Cont(Filevar,Nline_Ini,Nline_End,N_Elem_Type,Elem_Type,N_Elem)
    !!----    character(len=*),  dimension(:),    intent(in)    :: filevar       !  In -> String Vector
    !!----    integer,                            intent(in out):: nline_ini     !  In -> Line to start the search
    !!----                                                                         Out -> Actual Line on Filevar
    !!----    integer,                            intent(in)    :: nline_end     !  In -> Line to finish the search
    !!----    integer,                            intent(out)   :: n_elem_type   ! Out -> N. of different species
    !!----    character(len=*), dimension(:),     intent(out)   :: elem_type     ! Out -> Character to identify the specie
    !!----    real(kind=cp),dimension(:),optional,intent(out)   :: n_elem        ! Out -> Number of elements into the same species
    !!----
    !!----    Obtaining Chemical contents from Shelx file (.ins or .res)
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Shx_Cont(filevar,nline_ini,nline_end,n_elem_type,elem_type,n_elem)
       !---- Arguments ----!
       character(len=*), dimension(:),     intent(in)      :: filevar
       integer,                            intent(in out)  :: nline_ini
       integer,                            intent(in)      :: nline_end
       integer,                            intent(out)     :: n_elem_type
       character(len=*), dimension(:),     intent(out)     :: elem_type
       real(kind=cp),dimension(:),optional,intent(out)     :: n_elem

       !---- Local  variables ----!
       character(len=len(filevar(1)))      :: string
       integer                     :: iv
       integer,      dimension(15) :: ivet
       real(kind=cp),dimension(15) :: vet

       n_elem_type = 0
       elem_type   = " "
       if (present(n_elem)) n_elem = 0.0

       call Read_Key_StrVal(filevar,nline_ini,nline_end,"SFAC",string)
       if (len_trim(string) /=0) then
          call getword(string,elem_type,n_elem_type)
       end if

       if (present(n_elem)) then
          call read_key_value(filevar,nline_ini,nline_end,"UNIT",vet,ivet,iv)
          if (iv /= 0) n_elem=vet
       end if

       return
    End Subroutine Read_Shx_Cont

    !!----
    !!---- Subroutine Read_Shx_Fvar(Filevar,Nline_Ini,Nline_End,N_Fvar,Fvar)
    !!----    character(len=*), dimension(:), intent(in)    :: filevar       !  In -> String vector
    !!----    integer,                        intent(in out):: nline_ini     !  In -> Line to start the search
    !!----                                                                   ! Out -> Actual line on Filevar
    !!----    integer,                        intent(in)    :: nline_end     !  In -> Line to finish the search
    !!----    integer,                        intent(out)   :: n_fvar        ! Out -> N. of parameters on FVAR
    !!----    real(kind=cp), dimension(:),    intent(out)   :: fvar          ! Out -> values of FVAR
    !!----
    !!----    Obtaining Fvar parameters from Shelx file (.ins or .res)
    !!----
    !!---- Update: February - 2003
    !!
    Subroutine Read_Shx_Fvar(filevar,nline_ini,nline_end,n_fvar,fvar)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in)    :: filevar
       integer,                        intent(in out):: nline_ini
       integer,                        intent(in)    :: nline_end
       integer,                        intent(out)   :: n_fvar
       real(kind=cp), dimension(:),    intent(out)   :: fvar

       !---- Local  variables ----!
       integer                      :: iv
       integer,       dimension(15) :: ivet
       real(kind=cp), dimension(15) :: vet

       n_fvar = 1
       fvar   = 1.0

       call read_key_value(filevar,nline_ini,nline_end,"FVAR",vet,ivet,iv)
       if (iv /= 0) then
          n_fvar=iv
          fvar=vet
       end if

       return
    End Subroutine Read_Shx_Fvar

    !!----
    !!---- Subroutine Read_Shx_Latt(Filevar,Nline_Ini,Nline_End,Latt)
    !!----    character(len=*), dimension(:), intent(in) :: filevar     !  In -> String Vector
    !!----    integer,           intent(in out)          :: nline_ini   !  In -> Line to start the search
    !!----                                                                Out -> Actual line on Filevar
    !!----    integer,           intent(in)              :: nline_end   !  In -> Line to finish the search
    !!----    integer,           intent(out)             :: latt        ! Out -> Lattice number
    !!----
    !!----    Obtaining lattice from Shelx file (.ins or .res)
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Shx_Latt(filevar,nline_ini,nline_end,latt)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in) :: filevar
       integer,           intent(in out)          :: nline_ini
       integer,           intent(in)              :: nline_end
       integer,           intent(out)             :: latt

       !---- Local Variables ----!
       integer                     :: iv
       integer,       dimension(2) :: ivet
       real(kind=cp), dimension(2) :: vet

       latt=1
       call read_key_value(filevar,nline_ini,nline_end,"LATT",vet,ivet,iv)
       if (iv == 1) latt = ivet(1)

       return
    End Subroutine Read_Shx_Latt

    !!----
    !!---- Subroutine Read_Shx_Symm(Filevar,Nline_Ini,Nline_End,N_Oper,Oper_Symm)
    !!----    character(len=*), dimension(:), intent(in) :: filevar       !  In -> String Vector
    !!----    integer,           intent(in out)          :: nline_ini     !  In -> Line to start the search
    !!----                                                                  Out -> Actual Line on Filevar
    !!----    integer,           intent(in)              :: nline_end     !  In -> Line to finish the search
    !!----    integer,           intent(out)             :: n_oper        ! Out -> Number of Operators
    !!----    character(len=*), dimension(:),intent(out) :: oper_symm     ! Out -> String for Symmetry Operators
    !!----
    !!----    Obtaining Symmetry Operators from Shelx file (.ins or .res)
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Shx_Symm(filevar,nline_ini,nline_end,n_oper,oper_symm)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in) :: filevar
       integer,          intent(in out)           :: nline_ini
       integer,          intent(in)               :: nline_end
       integer,          intent(out)              :: n_oper
       character(len=*), dimension(:),intent(out) :: oper_symm

       !---- Local variables ----!
       character(len=80) :: string
       integer           :: nline

       n_oper=0
       oper_symm=" "

       do
          call Read_Key_StrVal(filevar,nline_ini,nline_end,"SYMM",string)
          if (len_trim(string) /=0) then
             n_oper=n_oper+1
             oper_symm(n_oper)=string
             nline_ini=nline_ini+1
             nline=nline_ini
          else
             exit
          end if
       end do
       nline_ini=nline

       return
    End Subroutine Read_Shx_Symm

    !!----
    !!---- Subroutine Read_Shx_Titl(Filevar,Nline_Ini,Nline_End,Title)
    !!----    character(len=*),dimension(:), intent(in)     :: filevar      !  In -> String Vector
    !!----    integer,                       intent(in out) :: nline_ini    !  In -> Line to start the search
    !!----                                                                    Out -> Actual Line on Filevar
    !!----    integer,                       intent(in)     :: nline_end    !  In -> Line to finish the search
    !!----    character(len=*),              intent(out)    :: title        ! Out -> Title
    !!----
    !!----    Obtaining Title from Shelx file
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Shx_Titl(filevar,nline_ini,nline_end,Title)
       !---- Arguments ----!
       character(len=*),dimension(:), intent(in)     :: filevar
       integer,                       intent(in out) :: nline_ini
       integer,                       intent(in)     :: nline_end
       character(len=*),              intent(out)    :: title

       call Read_Key_StrVal(filevar,nline_ini,nline_end,"TITL",title)

       return
    End Subroutine Read_Shx_Titl

    !!----
    !!---- Subroutine Read_Uvals(Line,Atomo,Ulabel)
    !!----    character(len=*),  intent(in out)  :: line      !  In -> String
    !!----    Type (Atom_Type),  intent(in out)  :: Atomo     !  In -> Atomo variable
    !!----                                                      Out ->
    !!----    character(len=4),  intent(in)      :: ulabel    !  In -> u_ij, b_ij, beta
    !!----
    !!----    Subroutine to read the anisotropic thermal parameters from a given Line
    !!----    it complets the object Atomo of type Atom.
    !!----    Assumes the string Line has been read from a file and
    !!----    starts with one of the words (u_ij, b_ij or beta), that is removed before reading
    !!----    the values of the parameters.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Uvals(Line,Atomo,Ulabel)
       !---- Arguments ----!
       character(len=*),  intent(in )     :: line
       Type (Atom_Type),  intent(in out)  :: Atomo
       character(len=4),  intent(in)      :: ulabel

       !---- Local variables -----!
       character(len=len(line)),dimension(1):: line2
       real(kind=cp), dimension (6)         :: vet1,vet2
       integer                              :: iv,n

       call init_err_form()

       atomo%utype    = ulabel
       line2(1)=line
       n=1
       call cutst(line2(1))
       line2(1)="Uval "//line2(1)(1:len(line2(1))-5)  !this form of writing is to avoid gfortran warning -Wstring-overflow
       call Read_Key_ValueSTD(line2,n,n,"Uval",vet1,vet2,iv)

        if (iv /= 6) then
          err_form=.true.
          ERR_Form_Mess="  Error reading the anisotropic thermal parameters of atom:"//atomo%lab
          return
       end if
       atomo%U(1:6)=vet1(1:6)
       atomo%U_std(1:6)=vet2(1:6)
       atomo%thtype="aniso"

       return
    End Subroutine Read_Uvals

    !!---- Subroutine Readn_Set_Magnetic_Space_Group(file_line,n_ini,n_end,MGp,mode,uvw)
    !!----    character(len=*),dimension(:),  intent (in)  :: file_line
    !!----    integer,                        intent (in)  :: n_ini,n_end
    !!----    type(Magnetic_Space_Group_Type),intent (out) :: MGp
    !!----    character(len=*),               intent (in)  :: mode
    !!----    character(len=*), optional,     intent (in)  :: uvw
    !!----
    !!----  This subroutine reads the lattice centring and anti-centring vectors
    !!----  as well as the symmetry operators of a magnetic space group in an
    !!----  arbitrary BNS setting. It construct the relevant magnetic space group
    !!----  components necessary for magnetic structure factor calculations.
    !!----  It may be used for reading from a CFL or a PCR file.
    !!----
    !!----  Created: March 2016.
    !!----
    !!----
    Subroutine Readn_Set_Magnetic_Space_Group(file_line,n_ini,n_end,MGp,mode,uvw)
       character(len=*),dimension(:),  intent (in)  :: file_line
       integer,                        intent (in)  :: n_ini,n_end
       type(Magnetic_Space_Group_Type),intent (out) :: MGp
       character(len=*),               intent (in)  :: mode
       character(len=*), optional,     intent (in)  :: uvw
       !
       ! --- Local variables ---!
       character(len=8)                 :: typ
       character(len=40)                :: Symbol
       integer                          :: i,j,ind,Nsym, Cen, N_Clat, N_Ant,ini, &
                                           num_sym,m,nop,k,n,L,ier,icount
       integer, dimension(3,3)          :: isim,msim
       real(kind=cp)                    :: p_mag
       real(kind=cp), dimension(3)      :: tr,v
       character(len=132)               :: line,ShOp_symb,setting,Parent
       character(len=40),dimension(10)  :: words
       logical                          :: u_type,m_type,inv_type, ttst,nonmag

       typ=l_case(adjustl(mode))

       call Init_Magnetic_Space_Group_Type(MGp)

       !Check if the database has to be read.
       nonmag=.false.; ttst=.false.
       if(typ /= "database") then
          do i=n_ini,n_end
           line=l_case(adjustl(file_line(i)))
           ind=index(line,"transform to standard:")
           if(ind /= 0) ttst=.true.
           ind=index(line,"<--nonmagnetic")
           if(ind /= 0) nonmag=.true.
           if(nonmag .and. ttst) then
             typ="database"
             exit
           end if
          end do
       end if

       Select Case(trim(typ))

          Case("pcr")
             line=adjustl(file_line(n_ini))
             ind=index(line,"Magnetic Space")
             if(ind == 0) then
               Err_Form=.true.
               Err_Form_Mess=" The Magnetic Space Group symbol is not provided in the PCR file! "
               return
             else
               j=index(line,"number:")
               MGp%BNS_symbol=trim(line(1:j-1))
               MGp%BNS_number=trim(line(j+7:ind-4))
             end if
             ini=n_ini+1
             Nsym=0; Cen=0; N_Clat=0;  N_Ant=0
             do i=ini,N_end
               line=adjustl(file_line(i))
               ind=index(line,"Transform to standard:")
               if(ind /= 0) then
                 MGp%trn_to_standard=adjustl(line(ind+22:))
               end if
               ind=index(line,"Parent Space Group:")
               if(ind /= 0) then
                 j=index(line,"IT_number:")
                 if( j /= 0) then
                   MGp%Parent_spg=adjustl(line(ind+19:j-1))
                   read(unit=line(j+10:),fmt=*,iostat=ier) MGp%Parent_num
                   if(ier /= 0) MGp%Parent_num=0
                 else
                   MGp%Parent_spg=adjustl(line(ind+19:))
                 end if
               end if
               ind=index(line,"Transform from Parent:")
               if(ind /= 0) then
                 MGp%trn_from_parent=adjustl(line(ind+22:))
               end if
               ind=index(line,"N_Clat")
               if(ind == 0) cycle
               read(unit=file_line(i+1),fmt=*) Nsym, Cen, N_Clat, N_Ant
               ini=i+2
               exit
             end do
             if(Nsym == 0) then
               Err_Form=.true.
               Err_Form_Mess=" The number of symmetry operators is not provided in the PCR file! "
               return
             end if
             !Allocate components of the magnetic space group
             MGp%Num_aLat=0
             allocate(MGp%Latt_trans(3,N_Clat+1))
             MGp%Latt_trans=0.0
             MGp%Num_Lat=N_Clat+1
             if(N_Ant > 0) then
               allocate(MGp%aLatt_trans(3,N_Ant))
               MGp%aLatt_trans=0.0
               MGp%Num_aLat=N_Ant
               MGp%MagType=4
             end if
             MGp%Numops = Nsym
             MGp%Centred= max(1,Cen)
             MGp%Multip = MGp%Numops * MGp%Centred * (MGp%Num_Lat + MGp%Num_aLat)
             num_sym=MGp%Multip
             allocate(Mgp%SymopSymb(num_sym))
             allocate(Mgp%Symop(num_sym))
             allocate(Mgp%MSymopSymb(num_sym))
             allocate(Mgp%MSymop(num_sym))
             if(N_Clat > 0) then
               do i=ini,N_end
                 line=adjustl(file_line(i))
                 ind=index(line,"Centring vectors")
                 if(ind == 0) cycle
                 ini=i+1
                 exit
               end do
               if(ind == 0) then
                 Err_Form=.true.
                 Err_Form_Mess=" 'Centring vectors' line is not provided in the PCR file! "
                 return
               end if
               m=1
               do i=ini,ini+N_Clat-1
                 m=m+1
                 read(unit=file_line(i),fmt=*) MGp%Latt_trans(:,m)
               end do
               ini=ini+N_Clat
             end if
             if(N_Ant > 0) then
               do i=ini,N_end
                 line=adjustl(file_line(i))
                 ind=index(line,"Anti-Centring vectors")
                 if(ind == 0) cycle
                 ini=i+1
                 exit
               end do
               if(ind == 0) then
                 Err_Form=.true.
                 Err_Form_Mess=" 'Anti-Centring vectors' line is not provided in the PCR file! "
                 return
               end if
               m=0
               do i=ini,ini+N_Ant-1
                 m=m+1
                 read(unit=file_line(i),fmt=*) MGp%aLatt_trans(:,m)
               end do
               ini=ini+N_Ant
             end if
             !Check the type of symmetry operators given
             do i=ini,N_end
                line=adjustl(file_line(i))
                if(line(1:1) == "!") cycle
                j=index(line,"!")
                if( j > 1) line=line(1:j-1)  !remove comments
                call Getword(line, words, icount)
                ! Icount=2 => SHSYM  x,-y,z+1/2,-1    <= This type
                ! Icount=3 => SHSYM  x,-y,z+1/2  -1   <= This type or these types => SHSYM x,-y,z+1/2  -u,v,-w  or SHSYM x,-y,z+1/2  -mx,my,-mz
                ! Icount=4 => SHSYM x,-y,z+1/2  -u,v,-w -1    <= This type or this type =>  SHSYM  x,-y,z+1/2  -mx,my,-mz  -1
                if( icount < 2 .or. icount > 4) then
                 Err_Form=.true.
                 Err_Form_Mess=" Error in Shubnikov operator: "//trim(line)
                 return
                end if
                u_type=(index(line,"u") /= 0) .or. (index(line,"U") /= 0)
                m_type=(index(line,"mx") /= 0) .or. (index(line,"MX") /= 0)
                if(.not. (u_type .or. m_type)) inv_type=.true.
                exit
             end do

             !Reading reduced set of symmetry operators
             m=0
             do i=ini,N_end
               line=adjustl(file_line(i))
               if(line(1:1) == "!") cycle
               j=index(line,"!")
               if( j > 1) line=line(1:j-1)  !remove comments
               j=index(line," ")
               line=adjustl(line(j:))
               m=m+1
               if(m > Nsym) exit
               call Getword(line, words, j)
               Select Case (icount)
                 Case(2)
                    j=index(line,",",back=.true.)
                    MGp%SymopSymb(m)=line(1:j-1)
                    read(unit=line(j+1:),fmt=*,iostat=ier) n
                    if(ier /= 0) then
                       Err_Form=.true.
                       Err_Form_Mess=" Error reading the time inversion in line: "//trim(file_line(i))
                       return
                    else
                       MGp%MSymOp(m)%phas=real(n)
                    end if
                    !write(*,"(a,i3)") trim(MGp%SymopSymb(m)),n
                 Case(3)
                    MGp%SymopSymb(m)=words(1)
                    MGp%MSymopSymb(m)=words(2)  !u,v,w or mx,my,mz or +/-1

                 Case(4)
                    MGp%SymopSymb(m)=words(1)
                    MGp%MSymopSymb(m)=words(2)  !u,v,w or mx,my,mz
                    read(unit=words(3),fmt=*,iostat=ier) n
                    if(ier /= 0) then
                       Err_Form=.true.
                       Err_Form_Mess=" Error reading the time inversion in line: "//trim(file_line(i))
                       return
                    else
                       MGp%MSymOp(m)%phas=real(n)
                    end if

               End Select
               call Read_Xsym(MGp%SymopSymb(m),1,isim,tr)
               MGp%Symop(m)%Rot=isim
               MGp%Symop(m)%tr=tr
               if(inv_type) then
                 j=determ_a(isim)
                 msim=nint(MGp%MSymOp(m)%phas)*j*isim
               else if (u_type) then
                 line=trim(MGp%MSymopSymb(m))//",0.0"
                 CALL read_msymm(line,msim,p_mag)
               else !should be mx,my,mz
                 line=trim(MGp%MSymopSymb(m))
                 do j=1,len_trim(line)
                    if(line(j:j) == "m" .or. line(j:j) == "M") line(j:j)=" "
                    if(line(j:j) == "x" .or. line(j:j) == "X") line(j:j)="u"
                    if(line(j:j) == "y" .or. line(j:j) == "Y") line(j:j)="v"
                    if(line(j:j) == "z" .or. line(j:j) == "Z") line(j:j)="w"
                 end do
                 line=pack_string(line)//",0.0"
                 CALL read_msymm(line,msim,p_mag)
               end if
               MGp%MSymop(m)%Rot=msim
               if(m_type .and. .not. present(uvw)) then
                 call Get_Shubnikov_Operator_Symbol(isim,msim,tr,ShOp_symb,.true.,invt=j)
                 MGp%mcif=.true.
               else
                 call Get_Shubnikov_Operator_Symbol(isim,msim,tr,ShOp_symb,invt=j)
                 MGp%mcif=.false.
               end if
               !write(*,"(a,i3)") trim(ShOp_symb),j
               MGp%MSymOp(m)%phas=j
               if(m_type .and. .not. present(uvw)) then
                 call Getword(ShOp_symb, words, j)
                 MGp%MSymopSymb(m)=words(2)
               else
                 j=index(ShOp_symb,";")
                 k=index(ShOp_symb,")")
                 MGp%MSymopSymb(m)=ShOp_symb(j+1:k-1)
               end if
             end do

          Case("cfl") !Only standard symbol plus an eventual setting change is allowed

               do i=ini,N_end
                 line=adjustl(file_line(i))
                 ind=index(l_case(line),"shubnikov")
                 if(ind == 0) cycle
                 ini=i+1
                 exit
               end do
               if(ind == 0) then
                  Err_Form=.true.
                  Err_Form_Mess=" Error reading the keyword SHUBNIKOV: the keyword is not found! "
                  return
               end if
               i=index(l_case(line),"setting:")
               symbol=trim(line(ind+9:i-1))
               j=index(line,"!")
               if(j /= 0) then
                 setting=line(i+8:j-1)
               else
                 setting="a,b,c;0,0,0"
               end if
               Parent=" "
               if(present(uvw)) then
                 if(index(uvw,"mx") /= 0) then
                    call Set_Magnetic_Space_Group(symbol,setting,MGp,mcif=.true.,trn_to=.true.)
                 else
                    call Set_Magnetic_Space_Group(symbol,setting,MGp,trn_to=.true.)
                 end if
               else
                 call Set_Magnetic_Space_Group(symbol,setting,MGp,trn_to=.true.)
               end if
               return

          Case("database")
             line=adjustl(file_line(n_ini))
             ind=index(line,"Magnetic Space")
             if(ind == 0) then
               Err_Form=.true.
               Err_Form_Mess=" The Magnetic Space Group symbol is not provided in the PCR/CFL file! "
               return
             else
               j=index(line," ")
               symbol=trim(line(1:j-1))
             end if
             ini=n_ini+1
             line=adjustl(file_line(ini))
                       !     123456789012345678901234567890
             ind=index(line,"Transform to standard:")
             if(ind == 0) then
               Err_Form=.true.
               Err_Form_Mess=" The transformation to standard is needed even if it is: a,b,c;0,0,0 "
               return
             else
               ind=index(line,"<--")
               if( ind == 0) then
                 line=adjustl(line(23:))
                 ind=index(line," ")
                 setting=line(23:ind)
               else
                 setting=line(23:ind-1)
               end if
             end if
             ini=ini+1
             line=adjustl(file_line(ini))
             Parent=" "      !12345678901234567890
             ind= index(line,"Parent space group:")
             j  = index(line,"IT_number:")
             if(ind /= 0 .and. j /= 0) then
               Parent= adjustl(line(20:j-1))
               ind=index(line,"<--")
               Parent=trim(Parent)//" "//line(j+10:ind-1)
             end if
             ini=ini+1
             line=adjustl(file_line(ini))
             ind= index(line,"Transform from Parent:")
             if(ind /= 0) then
               j=index(line,"<--")
               Parent=trim(Parent)//"  "//line(23:j-1)
             end if
!C_ac  number: "9.41"                           <--Magnetic Space Group (BNS symbol and number)
!Transform to standard:  c,-b,a;0,0,0           <--Basis transformation from current setting to standard BNS
!Parent Space Group: Pna2_1  IT_number:   33    <--Non-magnetic Parent Group
!123456789012345678901234567890
!Transform from Parent:   a,2b,2c;0,0,0         <--Basis transformation from parent to current setting
             !write(*,"(a)") trim(symbol)//" "//trim(setting)//" "//trim(parent)
             ! trn_to=.true. always because magCIF considers thre transformation from the current
             ! setting to the standard setting
             if(len_trim(Parent) /= 0) then
               call Set_Magnetic_Space_Group(symbol,setting,MGp,parent,trn_to=.true.)
             else
               call Set_Magnetic_Space_Group(symbol,setting,MGp,trn_to=.true.)
             end if
             return       !The clean-up of operators is not needed
       End Select

       !Expand symmetry operators if Cen=2 (centre of symmetry at the origin)
       m=MGp%Numops
       if(Cen == 2) then
          do i=1,MGp%Numops
            m=m+1
            MGp%SymOp(m)%Rot  = -MGp%SymOp(i)%Rot
            MGp%SymOp(m)%tr   =  modulo_lat(-MGp%SymOp(i)%tr)
            MGp%MSymOp(m)%phas= MGp%MSymOp(i)%phas
            MGp%MSymOp(m)%Rot = MGp%MSymOp(i)%Rot
            call Get_Symsymb(MGp%SymOp(m)%Rot,MGp%SymOp(m)%tr,MGp%SymopSymb(m))
            call Get_Symsymb(MGp%MSymOp(m)%Rot,(/0.0,0.0,0.0/),line)
            !Expand the operator "line" to convert it to mx,my,mz like
            MGp%MSymopSymb(m)=Get_MagMatSymb(line,MGp%mcif)
          end do
       end if
       nop=m
       !Expand symmetry operators for lattice centrings
       do L=1,N_clat
         tr=MGp%Latt_trans(:,L+1)
         do j=1,nop
           m=m+1
           v=MGp%SymOp(j)%tr(:) + tr
           MGp%SymOp(m)%Rot  = MGp%SymOp(j)%Rot
           MGp%SymOp(m)%tr   = modulo_lat(v)
           MGp%MSymOp(m)%Rot = MGp%MSymOp(j)%Rot
           MGp%MSymOp(m)%phas= MGp%MSymOp(j)%phas
           call Get_Symsymb(MGp%SymOp(m)%Rot,MGp%SymOp(m)%tr,MGp%SymopSymb(m))
           call Get_Symsymb(MGp%MSymOp(m)%Rot,(/0.0,0.0,0.0/),line)
           !Expand the operator "line" to convert it to mx,my,mz like
           MGp%MSymopSymb(m)=Get_MagMatSymb(line,MGp%mcif)
         end do
       end do
       !Expand symmetry operators for lattice anti-centrings
       do L=1,N_Ant
         tr=MGp%aLatt_trans(:,L)
         do j=1,nop
           m=m+1
           v=MGp%SymOp(j)%tr(:) + tr
           MGp%SymOp(m)%Rot  = MGp%SymOp(j)%Rot
           MGp%SymOp(m)%tr   = modulo_lat(v)
           MGp%MSymOp(m)%Rot = -MGp%MSymOp(j)%Rot
           MGp%MSymOp(m)%phas= -MGp%MSymOp(j)%phas
           call Get_Symsymb(MGp%SymOp(m)%Rot,MGp%SymOp(m)%tr,MGp%SymopSymb(m))
           call Get_Symsymb(MGp%MSymOp(m)%Rot,(/0.0,0.0,0.0/),line)
           !Expand the operator "line" to convert it to mx,my,mz like
           MGp%MSymopSymb(m)=Get_MagMatSymb(line,MGp%mcif)
         end do
       end do
       ! Symmetry operators treatment done!

    End Subroutine Readn_Set_Magnetic_Space_Group

    !!----
    !!---- Subroutine Readn_Set_Magnetic_Structure_MCIF(file_mcif,mCell,MGp,Am)
    !!----    character(len=*),               intent (in)  :: file_mcif
    !!----    type(Crystal_Cell_type),        intent (out) :: mCell
    !!----    type(Magnetic_Space_Group_Type),intent (out) :: MGp
    !!----    type(Atom_List_Type),           intent (out) :: Am
    !!----
    !!----    Subroutine for reading and construct a magnetic structure.
    !!----    The atom list and the unit cell reading an mCIF file.
    !!----
    !!----  Created: January-2014 (JRC)
    !!----  Updated: August-2014 (JRC), January 2020
    !!
    Subroutine Readn_Set_Magnetic_Structure_MCIF(file_mcif,mCell,MGp,Am)
       character(len=*),               intent (in)  :: file_mcif
       type(Crystal_Cell_type),        intent (out) :: mCell
       type(Magnetic_Space_Group_Type),intent (out) :: MGp
       type(Atom_List_Type),           intent (out) :: Am

       !---- Local Variables ----!
       integer :: i,num_sym, num_constr, num_kvs,num_matom, num_mom, num_magscat, ier, j, m, n, k, L,   &
                  ncar,mult,nitems,iv, num_irreps, nitems_irreps, num_rsym, num_centering,det,kfin
       integer,          dimension(10)     :: lugar
       integer,          dimension(7)      :: irrep_pos
       integer,          dimension(5)      :: pos
       integer,          dimension(3,3)    :: Rot
       real(kind=cp),    dimension(3)      :: cel,ang,cel_std,ang_std,tr,v
       real(kind=cp),    dimension(6)      :: values,std
       real(kind=cp),    dimension(3,3)    :: matr
       real(kind=cp),    dimension(3,384)  :: orb
       character(len=132)                  :: lowline,keyword,line, mxmymz_op,linat
       character(len=132),dimension(384)   :: sym_strings, cent_strings
       character(len=132),dimension(384)   :: atm_strings
       character(len=132),dimension(384)   :: mom_strings
       character(len=132),dimension(30)    :: constr_strings, mag_scatt_string
       character(len=132),dimension(30)    :: irreps_strings
       character(len=132),dimension(30)    :: kv_strings
       character(len=20), dimension(15)    :: lab_items
       character(len=50)                   :: shubk
       character(len=2)                    :: chars
       character(len=10)                   :: label
       character(len=4)                    :: symbcar
       logical                             :: ktag,no_symop_mxmymz,no_cent_mxmymz,mom_symmform

       !type(Magnetic_Group_Type)  :: SG
       type(file_list_type)       :: mcif

       call init_err_Form()
       call File_To_FileList(file_mcif,mcif)
       !Remove all possible tabs and non-ASCII characters in the CIF
       do i=1,mcif%nlines
         do j=1,len_trim(mcif%line(i))
           if(mcif%line(i)(j:j) == char(9)) mcif%line(i)(j:j)=" "
         end do
       end do
       num_constr=0; num_kvs=0; num_matom=0; num_mom=0; num_sym=0; num_magscat=0; num_rsym=0; num_centering=0
       cel=0.0; ang=0.0; num_irreps=0; nitems_irreps=0
       i=0
       call Init_Magnetic_Space_Group_Type(MGp)
       ktag=.false.
       no_symop_mxmymz=.false.
       no_cent_mxmymz=.false.
       mom_symmform=.false.

       do
          i=i+1
          if(i > mcif%nlines) exit
          if (index(mcif%line(i)(1:1),"!")/=0 .or. index(mcif%line(i)(1:1),"#")/=0 .or. len_trim(mcif%line(i)) == 0) cycle
          line=adjustl(mcif%line(i))
          lowline=l_case(line)
          j=index(lowline," ")
          keyword=lowline(1:j-1)
          !write(*,"(a)") " Keyword: "//trim(keyword)

          Select Case (trim(keyword))

             Case("_magnetic_space_group_standard_setting","_magnetic_space_group.standard_setting")
                chars=adjustl(line(j+1:))
                if(chars(2:2) == "y" .or. chars(2:2) == "Y") MGp%standard_setting=.true.
                !write(unit=*,fmt="(a)") "  Treating item: _magnetic_space_group_standard_setting -> "//trim(chars)

             Case("_parent_space_group.name_h-m", "_parent_space_group_name_h-m","_parent_space_group.name_h-m_alt")
                shubk=adjustl(line(j+1:))
                m=len_trim(shubk)
                MGp%Parent_spg=shubk(2:m-1)
                !write(unit=*,fmt="(a)") "  Treating item: _parent_space_group_name_h-m -> "// MGp%Parent_spg

             Case("_parent_space_group.it_number","_parent_space_group_it_number")
                read(unit=lowline(j:),fmt=*,iostat=ier) m
                if(ier /= 0) then
                  Err_Form=.true.
                  Err_Form_Mess=" Error reading the number of the parent space group"
                  return
                end if
                MGp%Parent_num=m
                !write(unit=*,fmt="(a,i4)") "  Treating item: _parent_space_group_it_number -> ", MGp%Parent_num

             Case("_magnetic_space_group_bns_number","_space_group.magn_number_bns","_space_group_magn.number_bns")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=adjustl(shubk(2:k-1))
                MGp%BNS_number=shubk
                !write(unit=*,fmt="(a)") "  Treating item: _space_group.magn_number_bns -> "//trim(MGp%BNS_number)

             Case("_magnetic_space_group_bns_name","_space_group_magn.name_bns","_space_group.magn_name_bns")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=adjustl(shubk(2:k-1))
                MGp%BNS_symbol=pack_string(shubk)
                !write(unit=*,fmt="(a)") "  Treating item: _space_group.magn_name_bns -> "//trim(MGp%BNS_symbol)

             Case("_magnetic_space_group_og_number","_space_group_magn.number_og","_space_group.magn_number_og")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=adjustl(shubk(2:k-1))
                MGp%OG_number=shubk
                !write(unit=*,fmt="(a)") "  Treating item: _space_group.magn_number_og -> "//trim(MGp%OG_number)

             Case("_magnetic_space_group_point_group","_space_group_magn.point_group","_space_group.magn_point_group")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=adjustl(shubk(2:k-1))
                MGp%PG_symbol=pack_string(shubk)
                !write(unit=*,fmt="(a)") "  Treating item: _space_group_magn.point_group -> "//trim(MGp%PG_symbol)

             Case("_magnetic_space_group_og_name","_space_group_magn.name_og","_space_group.magn_name_og")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=adjustl(shubk(2:k-1))
                MGp%OG_symbol=pack_string(shubk)
                !write(unit=*,fmt="(a)") "  Treating item: _space_group.magn_name_og -> "//trim(MGp%OG_symbol)

             Case("_magnetic_space_group.transform_from_parent_pp_abc","_magnetic_space_group_transform_from_parent_pp_abc", &
                   "_parent_space_group.child_transform_pp_abc")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=adjustl(shubk(2:k-1))
                MGp%trn_from_parent=pack_string(shubk)
                !write(unit=*,fmt="(a)") "  Treating item: _magnetic_space_group_transform_from_parent_pp_abc -> "//trim(MGp%trn_from_parent)

             Case("_parent_space_group.transform_pp_abc")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=adjustl(shubk(2:k-1))
                MGp%trn_to_parent=pack_string(shubk)
                !write(unit=*,fmt="(a)") "  Treating item: _magnetic_space_group_transform_from_parent_pp_abc -> "//trim(MGp%trn_from_parent)

             Case("_magnetic_space_group.transform_to_standard_pp_abc","_magnetic_space_group_transform_to_standard_pp_abc")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=adjustl(shubk(2:k-1))
                MGp%trn_to_standard=pack_string(shubk)

             Case("_space_group_magn.transform_bns_pp_abc")
                shubk=adjustl(line(j+1:))
                k=len_trim(shubk)
                if(shubk(1:1) == '"' .or. shubk(1:1) == "'") shubk=adjustl(shubk(2:k-1))
                MGp%trn_to_standard=pack_string(shubk)

             Case("_magnetic_cell_length_a","_cell_length_a")
                call getnum_std(lowline(j:),values,std,iv)
                if(err_string) then
                  Err_Form=.true.
                  Err_Form_Mess=" Error reading the magnetic unit cell parameter 'a' -> "//trim(err_string_mess)
                  return
                end if
                cel(1)=values(1)
                cel_std(1)=std(1)
                MGp%m_cell=.true.
                !write(unit=*,fmt="(a)") "  Treating item: _cell_length_a"

             Case("_magnetic_cell_length_b","_cell_length_b")
                call getnum_std(lowline(j:),values,std,iv)
                if(err_string) then
                  Err_Form=.true.
                  Err_Form_Mess=" Error reading the magnetic unit cell parameter 'b' -> "//trim(err_string_mess)
                  return
                end if
                cel(2)=values(1)
                cel_std(2)=std(1)
                !write(unit=*,fmt="(a)") "  Treating item: _cell_length_b"

             Case("_magnetic_cell_length_c","_cell_length_c")
                call getnum_std(lowline(j:),values,std,iv)
                if(err_string) then
                  Err_Form=.true.
                  Err_Form_Mess=" Error reading the magnetic unit cell parameter 'c' -> "//trim(err_string_mess)
                  return
                end if
                cel(3)=values(1)
                cel_std(3)=std(1)
                !write(unit=*,fmt="(a)") "  Treating item: _cell_length_c"

             Case("_magnetic_cell_angle_alpha","_cell_angle_alpha")
                call getnum_std(lowline(j:),values,std,iv)
                if(err_string) then
                  Err_Form=.true.
                  Err_Form_Mess=" Error reading the magnetic unit cell parameter 'alpha' -> "//trim(err_string_mess)
                  return
                end if
                ang(1)=values(1)
                ang_std(1)=std(1)
                !write(unit=*,fmt="(a)") "  Treating item: _cell_angle_alpha"

             Case("_magnetic_cell_angle_beta","_cell_angle_beta")
                call getnum_std(lowline(j:),values,std,iv)
                if(err_string) then
                  Err_Form=.true.
                  Err_Form_Mess=" Error reading the magnetic unit cell parameter 'beta' -> "//trim(err_string_mess)
                  return
                end if
                ang(2)=values(1)
                ang_std(2)=std(1)
                !write(unit=*,fmt="(a)") "  Treating item: _cell_angle_beta"

             Case("_magnetic_cell_angle_gamma","_cell_angle_gamma")
                call getnum_std(lowline(j:),values,std,iv)
                if(err_string) then
                  Err_Form=.true.
                  Err_Form_Mess=" Error reading the magnetic unit cell parameter 'gamma' -> "//trim(err_string_mess)
                  return
                end if
                ang(3)=values(1)
                ang_std(3)=std(1)
                !write(unit=*,fmt="(a)") "  Treating item: _cell_angle_gamma"

             Case("loop_")
                 i=i+1
                 line=adjustl(mcif%line(i))
                 lowline=l_case(line)
                 j=index(lowline," ")
                 keyword=lowline(1:j-1)
                 !write(*,"(a)") "         Loop_Keyword: "//trim(keyword)
                 Select Case(trim(keyword))

                   Case("_space_group_magn_transforms.id")
                      !write(*,"(a)") "         Loop_Keyword: "//trim(keyword)

                      do k=1,2
                        i=i+1
                        if(index(mcif%line(i),"_space_group_magn_transforms") == 0) then
                          Err_Form=.true.
                          Err_Form_Mess=" Error reading _space_group_magn_transforms in loop"
                          return
                        end if
                      end do
                      i=i+1
                      call getword(mcif%line(i),lab_items,iv)
                      !write(unit=*,fmt="(3a)")  (lab_items(k),k=1,3)
                      if(lab_items(3)(1:3) == "BNS") then
                        MGp%trn_to_standard=lab_items(2)
                      end if
                      i=i+1
                      call getword(mcif%line(i),lab_items,iv)
                      !write(unit=*,fmt="(3a)")  (lab_items(k),k=1,3)
                      if(lab_items(3)(1:2) == "OG") then
                        !nothing to do
                      end if

                   Case("_irrep_id")
                      irrep_pos=0
                      irrep_pos(1)=1
                      j=1
                      do k=1,6
                         i=i+1
                         if(index(mcif%line(i),"_irrep_dimension") /= 0) then
                            j=j+1
                            irrep_pos(2)=j
                            cycle
                         end if
                         if(index(mcif%line(i),"_small_irrep_dimension") /= 0 .or.  &
                            index(mcif%line(i),"_irrep_small_dimension") /= 0) then
                            j=j+1
                            irrep_pos(3)=j
                            cycle
                         end if
                         if(index(mcif%line(i),"_irrep_direction_type") /= 0) then
                            j=j+1
                            irrep_pos(4)=j
                            cycle
                         end if
                         if(index(mcif%line(i),"_irrep_action") /= 0) then
                            j=j+1
                            irrep_pos(5)=j
                            cycle
                         end if
                         if(index(mcif%line(i),"_irrep_modes_number") /= 0) then
                            j=j+1
                            irrep_pos(6)=j
                            cycle
                         end if
                         if(index(mcif%line(i),"_irrep_presence") /= 0) then
                            j=j+1
                            irrep_pos(7)=j
                            cycle
                         end if
                         exit
                      end do

                      i=i-1
                      nitems_irreps=count(irrep_pos > 0)

                      k=0
                      do
                        i=i+1
                        if(i > mcif%nlines) exit
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        irreps_strings(k)=mcif%line(i)
                      end do
                      num_irreps=k
                      !Treat later the list of irreps

                   Case("_magnetic_propagation_vector_seq_id")
                      do k=1,3
                        i=i+1
                        if(index(mcif%line(i),"_magnetic_propagation_vector") == 0) then
                          Err_Form=.true.
                          Err_Form_Mess=" Error reading the propagation vector loop"
                          return
                        end if
                        if(index(mcif%line(i),"_magnetic_propagation_vector_kxkykz") /= 0) then
                          ktag=.true.  !new format for k-vector klabel '0,1/2,0'
                          exit
                        end if
                      end do
                      k=0
                      do
                        i=i+1
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        kv_strings(k)=mcif%line(i)
                      end do
                      num_kvs=k
                      MGp%n_kv=k
                      if(allocated(Mgp%kv)) deallocate(Mgp%kv)
                      allocate(Mgp%kv(3,k))
                      if(allocated(Mgp%kv_label)) deallocate(Mgp%kv_label)
                      allocate(Mgp%kv_label(k))
                      !Treat later the propagation vectors

                   Case("_atom_type_symbol")
                      !write(unit=*,fmt="(a)") "  Treating item: _atom_type_symbol"
                      do k=1,3
                        i=i+1
                        if(index(mcif%line(i),"_atom_type_symbol") == 0) then
                          Err_Form=.true.
                          Err_Form_Mess=" Error reading the _atom_type_symbol in loop"
                          return
                        end if
                      end do
                      k=0
                      do
                        i=i+1
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        mag_scatt_string(k)=mcif%line(i)
                      end do
                      num_magscat=k
                      !Treat later the scattering factor

                   Case("_magnetic_atom_site_moment_symmetry_constraints_label")
                      !write(unit=*,fmt="(a)") "  Treating item: _magnetic_atom_site_moment_symmetry_constraints_label"
                      i=i+1
                      if(index(mcif%line(i),"_atom_site_magnetic_moment_symmetry_constraints_mxmymz") == 0) then
                        Err_Form=.true.
                        Err_Form_Mess=" Error reading the magnetic_atom_site_moment_symmetry_constraints loop"
                        return
                      end if
                      k=0
                      do
                        i=i+1
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        constr_strings(k)=mcif%line(i)
                      end do
                      num_constr=k
                      MGp%m_constr=.true.
                      !Treat later the constraints

                   Case("_magnetic_space_group_symop_id")
                      !write(unit=*,fmt="(a)") "  Treating item: _space_group_symop_magn_operation.id"
                      do k=1,3
                        i=i+1
                        j=index(mcif%line(i),"_magnetic_space_group_symop_operation")
                        if(j == 0 ) then
                          Err_Form=.true.
                          Err_Form_Mess=" Error reading the _magnetic_space_group_symop_operation loop"
                          return
                        end if
                      end do
                      k=0
                      do
                        i=i+1
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        sym_strings(k)=mcif%line(i)
                      end do
                      !now allocate the list of symmetry operators
                      num_sym=k
                      MGp%Multip=k

                   Case("_space_group_symop_magn_operation.id","_space_group_symop_magn.id") !The second item is added to be compatible with BCS error
                      !write(unit=*,fmt="(a)") "  Treating item: _space_group_symop_magn_operation.id"

                      i=i+1
                      j=index(mcif%line(i),"_space_group_symop_magn_operation.xyz")
                      if(j == 0 ) then
                        Err_Form=.true.
                        Err_Form_Mess=" Error reading the _space_group_symop_magn_operation loop"
                        return
                      end if

                      k=0
                      do
                        i=i+1
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        sym_strings(k)=mcif%line(i)
                      end do
                      !now allocate the list of symmetry operators
                      num_sym=k
                      MGp%Multip=k

                   Case("_space_group_symop.magn_id")
                      !write(unit=*,fmt="(a)") "  Treating item: _space_group_symop_magn_id"
                      do k=1,2
                        i=i+1
                        if(index(mcif%line(i),"_space_group_symop.magn_operation") == 0 .and. &
                           index(mcif%line(i),"_space_group_symop_magn_operation") == 0) then
                          Err_Form=.true.
                          Err_Form_Mess=" Error reading the _space_group_symop_magn_operation loop"
                          return
                        end if
                      end do
                      if(index(mcif%line(i),"_space_group_symop.magn_operation_mxmymz") == 0 .and. &
                         index(mcif%line(i),"_space_group_symop_magn_operation_mxmymz") == 0) then
                         i=i-1
                         no_symop_mxmymz=.true.
                      end if
                      k=0
                      do
                        i=i+1
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        sym_strings(k)=mcif%line(i)
                      end do

                      num_rsym=k

                   Case("_space_group_symop_magn_id")   !here the symmetry operators are separated from the translations
                      !write(unit=*,fmt="(a)") "  Treating item: _space_group_symop_magn_id"
                      i=i+1
                      if(index(mcif%line(i),"_space_group_symop.magn_operation") == 0 .and. &
                         index(mcif%line(i),"_space_group_symop_magn_operation") == 0) then
                        Err_Form=.true.
                        Err_Form_Mess=" Error reading the _space_group_symop.magn_operation loop"
                        return
                      end if
                      if(index(mcif%line(i),"_space_group_symop.magn_operation_mxmymz") == 0 .and. &
                         index(mcif%line(i),"_space_group_symop_magn_operation_mxmymz") == 0) then
                         no_symop_mxmymz=.true.
                      end if
                      k=0
                      do
                        i=i+1
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        sym_strings(k)=mcif%line(i)
                      end do

                      num_rsym=k

                   Case("_space_group_symop.magn_centering_id")   !here we read the translations and anti-translations
                      !write(unit=*,fmt="(a)") "  Treating item: _space_group_symop_magn_centering_id"
                      do k=1,2
                        i=i+1
                        if(index(mcif%line(i),"_space_group_symop.magn_centering") == 0 .and. &
                           index(mcif%line(i),"_space_group_symop_magn_centering") == 0 ) then
                          Err_Form=.true.
                          Err_Form_Mess=" Error reading the _space_group_symop_magn_centering loop"
                          return
                        end if
                      end do
                      if(index(mcif%line(i),"_space_group_symop.magn_centering_mxmymz") == 0 .and. &
                         index(mcif%line(i),"_space_group_symop_magn_centering_mxmymz") == 0 ) then
                         i=i-1
                         no_cent_mxmymz=.true.
                      end if
                      k=0
                      do
                        i=i+1
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        cent_strings(k)=mcif%line(i)
                      end do
                      num_centering=k

                   Case("_space_group_symop_magn_centering.id")   !here we read the translations and anti-translations
                      !write(unit=*,fmt="(a)") "  Treating item: _space_group_symop_magn_centering_id"
                      i=i+1
                      if(index(mcif%line(i),"_space_group_symop_magn_centering.xyz") == 0) then
                        Err_Form=.true.
                        Err_Form_Mess=" Error reading the _space_group_symop_magn_centering.xyz loop"
                        return
                      end if
                      k=0
                      do
                        i=i+1
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        cent_strings(k)=mcif%line(i)
                      end do
                      num_centering=k

                   Case("_magnetic_atom_site_label","_atom_site_label")
                      !write(unit=*,fmt="(a)") "  Treating item: _atom_site_label"
                      !Count the number of keywords following the _loop
                      do k=1,10
                        linat=adjustl(mcif%line(i+k))
                        if(linat(1:1) /=  "_") then
                          kfin=k+1
                          iv=i+k
                          exit
                        end if
                      end do
                      lugar=0
                      lugar(1)=1
                      j=1
                      do k=1,kfin
                         i=i+1
                         if(index(mcif%line(i),"_atom_site_type_symbol") /= 0) then
                            j=j+1
                            lugar(2)=j
                            cycle
                         end if
                         if(index(mcif%line(i),"_atom_site_fract_x") /= 0) then
                            j=j+1
                            lugar(3)=j
                            cycle
                         end if
                         if(index(mcif%line(i),"_atom_site_fract_y") /= 0) then
                            j=j+1
                            lugar(4)=j
                            cycle
                         end if
                         if(index(mcif%line(i),"_atom_site_fract_z") /= 0) then
                            j=j+1
                            lugar(5)=j
                            cycle
                         end if
                         if (index(mcif%line(i),"_atom_site_U_iso_or_equiv") /= 0) then
                            j=j+1
                            lugar(6)=j
                            cycle
                         end if
                         if (index(mcif%line(i),"_atom_site_B_iso_or_equiv") /= 0) then
                            j=j+1
                            lugar(10)=j
                            cycle
                         end if
                         if (index(mcif%line(i),"_atom_site_occupancy") /= 0) then
                            j=j+1
                            lugar(7)=j
                            cycle
                         end if
                         if (index(mcif%line(i),"_atom_site_symmetry_multiplicity") /= 0) then
                            j=j+1
                            lugar(8)=j
                            cycle
                         end if
                         if (index(mcif%line(i),"_atom_site_Wyckoff_label") /= 0) then
                            j=j+1
                            lugar(9)=j
                            cycle
                         end if
                         exit
                      end do

                      if (any(lugar(3:5) == 0)) then
                          Err_Form=.true.
                          Err_Form_Mess=" Error reading the asymmetric unit of magnetic atoms"
                          return
                      end if

                      i=iv-1
                      nitems=count(lugar > 0)

                      k=0
                      do
                        i=i+1
                        if(i > mcif%nlines) exit
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        atm_strings(k)=adjustl(mcif%line(i))
                      end do
                      num_matom=k
                      !Treat late the list atoms

                   Case("_magnetic_atom_site_moment_label","_atom_site_moment_label","_atom_site_moment.label")
                      !write(unit=*,fmt="(a)") "  Treating item: _atom_site_moment_label"
                      do k=1,3
                        i=i+1
                        if(index(mcif%line(i),"_atom_site_moment_crystalaxis") == 0 .and. &
                           index(mcif%line(i),"_atom_site_moment.crystalaxis") == 0) then
                          Err_Form=.true.
                          Err_Form_Mess=" Error reading the magnetic_atom_site_moment loop"
                          return
                        end if
                      end do
                      i=i+1
                      if(index(mcif%line(i),"_atom_site_moment.symmform") /= 0) then
                        !write(*,*) " _atom_site_moment.symmform FOUND"
                        mom_symmform=.true.
                      else
                        i=i-1
                      end if
                      k=0
                      do
                        i=i+1
                        if(i > mcif%nlines) exit
                        if(len_trim(mcif%line(i)) == 0) exit
                        k=k+1
                        mom_strings(k)=mcif%line(i)
                      end do
                      num_mom=k
                      !Treat later the magnetic moment of the atoms
                 End Select
          End Select
       end do

       if(MGp%m_cell) then
         call Set_Crystal_Cell(cel,ang,mCell)
         mCell%cell_std=cel_std
         mCell%ang_std=ang_std
       end if

       !Treat symmetry operators
       !write(unit=*,fmt="(a,2i4)") " num_sym, num_rsym :",num_sym,num_rsym
       if(num_sym == 0 .and. num_rsym == 0) then
          Err_Form=.true.
          Err_Form_Mess=" No symmetry operators have been provided in the MCIF file "//trim(file_mcif)
          return
       else
          if(no_cent_mxmymz) then  !Full number of symmetry operators is not separated from the centering

            if(allocated(Mgp%SymopSymb)) deallocate(Mgp%SymopSymb)
            allocate(Mgp%SymopSymb(num_sym))
            if(allocated(Mgp%Symop)) deallocate(Mgp%Symop)
            allocate(Mgp%Symop(num_sym))
            if(allocated(Mgp%MSymopSymb)) deallocate(Mgp%MSymopSymb)
            allocate(Mgp%MSymopSymb(num_sym))
            if(allocated(Mgp%MSymop)) deallocate(Mgp%MSymop)
            allocate(Mgp%MSymop(num_sym))
            !write(unit=*,fmt="(a)") "  Decoding symmetry operators 1"

            ! Decode the symmetry operators
            do i=1,num_sym
              line=adjustl(sym_strings(i))
              j=index(line," ")
              line=adjustl(line(j+1:))
              j=index(line," ")
              MGp%SymopSymb(i)=line(1:j-1)
              line=adjustl(line(j+1:))
              j=index(line," ")
              MGp%MSymopSymb(i)=line(1:j-1)
              read(unit=line(j:),fmt=*,iostat=ier) n
              if(ier /= 0) then
                 Err_Form=.true.
                 Err_Form_Mess=" Error reading the time inversion in line: "//trim(sym_strings(i))
                 return
              else
                 MGp%MSymOp(i)%phas=real(n)
              end if
              call Read_Xsym(MGp%SymopSymb(i),1,MGp%Symop(i)%Rot,MGp%Symop(i)%tr)
              line=MGp%MSymopSymb(i)
              do k=1,len_trim(line)
                if(line(k:k) == "m") line(k:k)=" "
              end do
              line=Pack_String(line)
              call Read_Xsym(line,1,MGp%MSymop(i)%Rot)
            end do

          else

            if( num_rsym == 0) num_rsym=num_sym
            ! First allocate the full number of symmetry operators after decoding if centering lattice
            ! have been provided and if the group is centred or not
            if(num_centering == 0) then
               MGp%Multip=num_rsym
            else
               MGp%Multip=num_rsym*num_centering
            end if

            num_sym=MGp%Multip
            if(allocated(Mgp%SymopSymb)) deallocate(Mgp%SymopSymb)
            allocate(Mgp%SymopSymb(num_sym))
            if(allocated(Mgp%Symop)) deallocate(Mgp%Symop)
            allocate(Mgp%Symop(num_sym))
            if(allocated(Mgp%MSymopSymb)) deallocate(Mgp%MSymopSymb)
            allocate(Mgp%MSymopSymb(num_sym))
            if(allocated(Mgp%MSymop)) deallocate(Mgp%MSymop)
            allocate(Mgp%MSymop(num_sym))
            ! Decode the symmetry operators
            !write(unit=*,fmt="(a)") "  Decoding symmetry operators 2"
            do i=1,num_rsym
              line=adjustl(sym_strings(i))
              j=index(line," ")
              line=adjustl(line(j+1:))
              j=index(line," ")
              MGp%SymopSymb(i)=line(1:j-1)
              k=index(MGp%SymopSymb(i),",",back=.true.)
              read(unit=MGp%SymopSymb(i)(k+1:),fmt=*,iostat=ier) n
              if(ier /= 0) then
                 Err_Form=.true.
                 Err_Form_Mess=" Error reading the time inversion in line: "//trim(sym_strings(i))
                 return
              else
                 MGp%MSymOp(i)%phas=real(n)
              end if
              MGp%SymopSymb(i)=MGp%SymopSymb(i)(1:k-1)
              call Read_Xsym(MGp%SymopSymb(i),1,MGp%Symop(i)%Rot,MGp%Symop(i)%tr)

              !Now construc the magnetic rotation symbols
              line=adjustl(line(j+1:))
              if(len_trim(line) /= 0) then
                j=index(line," ")
                MGp%MSymopSymb(i)=line(1:j-1)
                line=MGp%MSymopSymb(i)
                do k=1,len_trim(line)
                  if(line(k:k) == "m") line(k:k)=" "
                end do
                line=Pack_String(line)
                call Read_Xsym(line,1,MGp%MSymop(i)%Rot)
              else
                det=determ_a(MGp%Symop(i)%Rot)
                MGp%MSymop(i)%Rot=MGp%Symop(i)%Rot*det*nint(MGp%MSymOp(i)%phas)
                call Get_Symsymb(MGp%MSymOp(i)%Rot,(/0.0,0.0,0.0/),line)
                !Expand the operator "line" to convert it to mx,my,mz like
                mxmymz_op=" "
                do j=1,len_trim(line)
                  Select Case(line(j:j))
                    case("x")
                       mxmymz_op=trim(mxmymz_op)//"mx"
                    case("y")
                       mxmymz_op=trim(mxmymz_op)//"my"
                    case("z")
                       mxmymz_op=trim(mxmymz_op)//"mz"
                    case default
                       mxmymz_op=trim(mxmymz_op)//line(j:j)
                  End Select
                end do
                MGp%MSymopSymb(i)=trim(mxmymz_op)
              end if


            end do
            !Decode lattice translations and anti-translations

            !write(unit=*,fmt="(a)") "  Decoding lattice translations and anti-translations"
            m=num_rsym
            do L=2,num_centering
              line=adjustl(cent_strings(L))
              j=index(line," ")
              line=adjustl(line(j+1:))
              j=index(line," ")
              line=line(1:j-1)
              k=index(line,",",back=.true.)
              read(unit=line(k+1:),fmt=*,iostat=ier) n
              if(ier /= 0) then
                 Err_Form=.true.
                 Err_Form_Mess=" Error reading the time inversion in line: "//trim(cent_strings(i))
                 return
              end if
              line=line(1:k-1)
              call Read_Xsym(line,1,Rot,tr)

              do j=1,num_rsym
                m=m+1
                v=MGp%SymOp(j)%tr(:) + tr
                MGp%SymOp(m)%Rot  = MGp%SymOp(j)%Rot
                MGp%SymOp(m)%tr   = modulo_lat(v)
                MGp%MSymOp(m)%Rot = n*MGp%MSymOp(j)%Rot
                MGp%MSymOp(m)%phas= n*MGp%MSymOp(j)%phas
                call Get_Symsymb(MGp%SymOp(m)%Rot,MGp%SymOp(m)%tr,MGp%SymopSymb(m))
                call Get_Symsymb(MGp%MSymOp(m)%Rot,(/0.0,0.0,0.0/),line)
                !Expand the operator "line" to convert it to mx,my,mz like
                mxmymz_op=" "
                do i=1,len_trim(line)
                  Select Case(line(i:i))
                    case("x")
                       mxmymz_op=trim(mxmymz_op)//"mx"
                    case("y")
                       mxmymz_op=trim(mxmymz_op)//"my"
                    case("z")
                       mxmymz_op=trim(mxmymz_op)//"mz"
                    case default
                       mxmymz_op=trim(mxmymz_op)//line(i:i)
                  End Select
                end do
                MGp%MSymopSymb(m)=trim(mxmymz_op)
              end do
            end do
          end if
       end if
       ! Symmetry operators treatment done
       Call cleanup_symmetry_operators(MGp)
       if(Err_Form) then
          return
          !write(unit=*,fmt="(a)") " => "//trim(Err_Form)
       end if

       !Treating irreps

       if(num_irreps == 0) then

          MGp%n_irreps=0

       else
          !write(*,"(a,i3)") " Treating irreps: ",num_irreps
          MGp%n_irreps=num_irreps
          if(allocated(MGp%irrep_dim))          deallocate(MGp%irrep_dim)
          if(allocated(MGp%small_irrep_dim))    deallocate(MGp%small_irrep_dim)
          if(allocated(MGp%irrep_id))           deallocate(MGp%irrep_id)
          if(allocated(MGp%irrep_direction))    deallocate(MGp%irrep_direction)
          if(allocated(MGp%irrep_action))       deallocate(MGp%irrep_action)
          if(allocated(MGp%irrep_modes_number)) deallocate(MGp%irrep_modes_number)
          allocate(MGp%irrep_dim(num_irreps),MGp%small_irrep_dim(num_irreps),MGp%irrep_id(num_irreps), &
                   MGp%irrep_direction(num_irreps),MGp%irrep_action(num_irreps),MGp%irrep_modes_number(num_irreps))

          MGp%irrep_dim=0; MGp%small_irrep_dim=0; MGp%irrep_id=" "; MGp%irrep_direction=" "; MGp%irrep_action=" "
          MGp%irrep_modes_number=0

          do i=1,MGp%n_irreps

            call getword(irreps_strings(i),lab_items,iv)

            !if(iv /= nitems_irreps) write(*,"(2(a,i2))") " => Warning irreps_nitems=",nitems_irreps," /= items read=",iv

            MGp%irrep_id(i)=lab_items(irrep_pos(1))
            if(MGp%irrep_id(i) == "?") then
               MGp%n_irreps=0
               exit
            end if

            if (irrep_pos(2) /= 0) then
               read(unit=lab_items(irrep_pos(2)),fmt=*,iostat=ier) MGp%irrep_dim(i)
               if(ier /= 0) MGp%irrep_dim(i)=0
            end if

            if (irrep_pos(3) /= 0) then
               read(unit=lab_items(irrep_pos(3)),fmt=*,iostat=ier) MGp%small_irrep_dim(i)
               if(ier /= 0) MGp%small_irrep_dim(i)=0
            end if

            if (irrep_pos(4) /= 0) then
               MGp%irrep_direction(i)=lab_items(irrep_pos(4))
            end if

            if (irrep_pos(5) /= 0) then
               MGp%irrep_action(i)=lab_items(irrep_pos(5))
            end if

            if (irrep_pos(6) /= 0) then
               read(unit=lab_items(irrep_pos(6)),fmt=*,iostat=ier) MGp%irrep_modes_number(i)
               if(ier /= 0) MGp%irrep_modes_number(i)=0
            end if

          end do
       end if
       ! End treatment of irreps

       ! Treating propagation vectors
       if(num_kvs == 0) then
         MGp%n_kv=0
       else
         !write(*,"(a,i3)") " Treating propagation vectors: ",num_kvs
         do i=1,MGp%n_kv
            line=adjustl(kv_strings(i))
            j=index(line," ")
            MGp%kv_label(i)=line(1:j-1)
            line=adjustl(line(j+1:))
            n=len_trim(line)
            if(ktag) then
              line=adjustl(line(2:n-1))
              n=n-2
              Call Get_Separator_Pos(line,",",pos,ncar)
            else
              Call Get_Separator_Pos(line," ",pos,ncar)
            end if
            keyword=line(1:pos(1)-1)//"a,"//line(pos(1)+1:pos(2)-1)//"b,"//trim(line(pos(2)+1:))//"c"
            keyword=Pack_String(keyword)
            call Get_Mat_From_Symb(keyword,Matr, (/"a","b","c"/) )
            do k=1,3
               MGp%kv(k,i)=Matr(k,k)
            end do
         end do
       end if
       ! Propagation vectors treatment done!

       !Treating magnetic atoms
       if(num_matom == 0) then
          Am%natoms = 0
          return
       else
          !write(*,"(a,i4)") " Treating magnetic atoms:  ",num_matom
          Call Allocate_Atom_list(num_matom,Am)

          do i=1,Am%natoms

            call getword(atm_strings(i),lab_items,iv)
            !if(iv /= nitems) write(*,"(2(a,i2))") " => Warning nitems=",nitems," /= items read=",iv
            Am%atom(i)%lab=lab_items(lugar(1))
            if (lugar(2) /= 0) then
               Am%atom(i)%SfacSymb=lab_items(lugar(2))(1:4)
               if(index("1234567890+-",lab_items(lugar(2))(2:2)) /= 0 ) then
                  Am%atom(i)%chemSymb=U_case(lab_items(lugar(2))(1:1))
               else
                  Am%atom(i)%chemSymb=U_case(lab_items(lugar(2))(1:1))//L_case(lab_items(lugar(2))(2:2))
               end if
            else
               if(index("1234567890+-",lab_items(lugar(1))(2:2)) /= 0 ) then
                  Am%atom(i)%chemSymb=U_case(lab_items(lugar(1))(1:1))
               else
                  Am%atom(i)%chemSymb=U_case(lab_items(lugar(1))(1:1))//L_case(lab_items(lugar(1))(2:2))
               end if
               Am%atom(i)%SfacSymb=Am%atom(i)%chemSymb
            end if
            call getnum_std(lab_items(lugar(3)),values,std,iv)    ! _atom_site_fract_x
            Am%atom(i)%x(1)=values(1)
            Am%atom(i)%x_std(1)=std(1)
            call getnum_std(lab_items(lugar(4)),values,std,iv)    ! _atom_site_fract_y
            Am%atom(i)%x(2)=values(1)
            Am%atom(i)%x_std(2)=std(1)
            call getnum_std(lab_items(lugar(5)),values,std,iv)    ! _atom_site_fract_z
            Am%atom(i)%x(3)=values(1)
            Am%atom(i)%x_std(3)=std(1)

            if (lugar(6) /= 0) then  ! _atom_site_U_iso_or_equiv
               call getnum_std(lab_items(lugar(6)),values,std,iv)
               Am%atom(i)%ueq=values(1)
               Am%atom(i)%Biso=values(1)*78.95683521     !If anisotropic they
               Am%atom(i)%Biso_std=std(1)*78.95683521    !will be put to zero
            else if (lugar(10) /= 0) then    ! _atom_site_B_iso_or_equiv
               call getnum_std(lab_items(lugar(10)),values,std,iv)
               Am%atom(i)%ueq=values(1)/78.95683521
               Am%atom(i)%Biso=values(1)     !If anisotropic they
               Am%atom(i)%Biso_std=std(1)    !will be put to zero
            else
               Am%atom(i)%ueq=0.0
               Am%atom(i)%Biso=0.0
               Am%atom(i)%Biso_std=0.0
            end if
            Am%atom(i)%utype="u_ij"

            if (lugar(7) /= 0) then ! _atom_site_occupancy
               call getnum_std(lab_items(lugar(7)),values,std,iv)
            else
               values=1.0
               std=0.0
            end if
            Am%atom(i)%occ=values(1)
            Am%atom(i)%occ_std=std(1)

            if(lugar(8) /= 0) then
              read(unit=lab_items(lugar(8)),fmt=*) Mult
              Am%atom(i)%mult=Mult
            else
              Call Get_mOrbit(Am%atom(i)%x,MGp,Mult,orb)
              Am%atom(i)%mult=Mult
            end if
            !Conversion from occupancy to occupation factor
            Am%atom(i)%occ=Am%atom(i)%occ*real(Mult)/real(MGp%Multip)

            if(lugar(9) /= 0) then
               Am%atom(i)%wyck=adjustl(trim(lab_items(lugar(9))))
            end if

          end do
       end if

       !Treating moments of magnetic atoms
       if(num_mom /= 0) then
          !write(*,"(a,i4)") " Treating magnetic moments:  ",num_mom
          m=4
          if(mom_symmform) m=5
          do i=1,num_mom
            call getword(mom_strings(i),lab_items,iv)
            !write(*,"(2i6,tr4,5(a,tr3))") k,iv,lab_items(1:iv)
            if(iv /= m) then
               Err_Form=.true.
               write(unit=Err_Form_Mess,fmt="(a,i4)")" Error reading magnetic moment #",i
               Err_Form_Mess=trim(Err_Form_Mess)//" -> 4-5 items expected in this line: 'Label mx my mz', read: "// &
                                                      trim(mom_strings(i))
               return
            end if
            label=Lab_items(1)
            do j=1,Am%natoms
               if(label == Am%Atom(j)%lab) then
                 do k=1,3
                     call getnum_std(lab_items(1+k),values,std,iv)
                     Am%Atom(j)%M_xyz(k)=values(1)
                     Am%Atom(j)%sM_xyz(k)=std(1)
                 end do
                 Am%Atom(j)%moment=99.0  !used for indicating that this atom is susceptible to bring a magnetic moment
               end if
            end do
          end do
       end if

       if(num_constr /= 0) then

         !write(*,"(a,i4)") " Treating constraints:  ",num_constr
         do i=1,num_constr
           line=adjustl(constr_strings(i))
           j=index(line," ")
           label=line(1:j-1)
           keyword=adjustl(line(j+1:))
           Call Get_Separator_Pos(keyword,",",pos,ncar)
           if(ncar == 0) then !There are no ","
             j=index(keyword," ")
             shubk=keyword(1:j-1)//","
             keyword=adjustl(keyword(j+1:))
             j=index(keyword," ")
             shubk=trim(shubk)//keyword(1:j-1)//","
             keyword=trim(shubk)//trim(adjustl(keyword(j+1:)))
           end if
           do j=1,len_trim(keyword)
             if(keyword(j:j) == "m") keyword(j:j) = " "
           end do
           keyword=Pack_String(keyword)
           !write(*,"(a)") "  constr_string: "//trim(line)
           !write(*,"(a)") "        keyword: "//trim(keyword)
           call Get_Mat_From_Symb(keyword,Matr, (/"x","y","z"/) )
           !write(*,"(9f10.3)") Matr
           do j=1,Am%natoms
             if(label == Am%Atom(j)%lab) then
                Am%Atom(j)%M_xyz=matmul(Matr,Am%Atom(j)%M_xyz)
                Am%Atom(j)%AtmInfo=constr_strings(i)
                Am%Atom(j)%moment=99.0  !used for indicating that this atom is susceptible to bring a magnetic moment
                exit
             end if
           end do
           !The treatment of the codes will be done in the future
         end do
       end if

       if(num_magscat > 0) then !Reading the valence for determining the magnetic form factor
         do i=1,num_magscat
           call getword(mag_scatt_string(i),lab_items,iv)
           do j=1,Am%natoms
             if(Am%atom(j)%chemSymb == lab_items(1)) then
               Am%atom(j)%SfacSymb=lab_items(2)
               if(lab_items(2) /= ".") then !magnetic atoms
                  Am%Atom(j)%moment=99.0  !used for indicating that this atom is susceptible to bring a magnetic moment
               end if
             end if
           end do
         end do
       end if

       !Get pointers to the magnetic form factors
       !Stored for each atom in the component ind(1)
       call Set_Magnetic_Form()

       !---- Find Species in Magnetic_Form ----!
       do i=1,Am%natoms
          symbcar=get_magnetic_form_factor(Am%atom(i)%SfacSymb)
          do j=1,num_mag_form
             if (symbcar /= Magnetic_Form(j)%Symb) cycle
             Am%atom(i)%ind(2)=j
             Am%atom(i)%SfacSymb=symbcar
             exit
          end do
       end do

       return
    End Subroutine Readn_Set_Magnetic_Structure_MCIF

    !!--++
    !!--++ Subroutine Readn_Set_XTal_CFL(file_dat,nlines,Cell,SpG,A,CFrame,NPhase,Job_Info)
    !!--++    character(len=*),dimension(:),intent(in)   :: file_dat
    !!--++    integer,                      intent(in)   :: nlines
    !!--++    Type (Crystal_Cell_Type),     intent(out)  :: Cell
    !!--++    Type (Space_Group_Type),      intent(out)  :: SpG
    !!--++    Type (atom_list_type),        intent(out)  :: A
    !!--++    character(len=*),    optional,intent(in)   :: CFrame
    !!--++    Integer,             optional,intent( in)  :: Nphase
    !!--++    Type(Job_Info_type), optional,intent(out)  :: Job_Info
    !!--++
    !!--++ (Private)
    !!--++ Read and Set Crystal Information in a CFL File
    !!--++
    !!--++ Update: April - 2005
    !!
    Subroutine Readn_Set_XTal_CFL(file_dat,nlines,Cell,SpG,A,CFrame,NPhase,Job_Info)
       !---- Arguments ----!
       character(len=*),dimension(:),intent(in)   :: file_dat
       integer,                      intent(in)   :: nlines
       Type (Crystal_Cell_Type),     intent(out)  :: Cell
       Type (Space_Group_Type),      intent(out)  :: SpG
       Type (atom_list_type),        intent(out)  :: A
       character(len=*),    optional,intent(in)   :: CFrame
       Integer,             optional,intent( in)  :: Nphase
       Type(Job_Info_type), optional,intent(out)  :: Job_Info

       !---- Local variables ----!
       character(len=132)               :: line
       character(len= 20)               :: Spp
       character(len= 40),dimension(192):: gen
       integer                          :: i, nauas, ndata, iph, n_ini,n_end,ngen,k,nsym
       integer, parameter               :: maxph=21  !Maximum number of phases "maxph-1"
       integer, dimension(maxph)        :: ip

       real(kind=cp),dimension(3):: vet

       !---- Standard CrysFML file *.CFL ----!
       nauas=0
       ndata=0
       ip=nlines
       ip(1)=1

       !---- Calculating number of Phases ----!
       do i=1,nlines
          line=adjustl(file_dat(i))
          if (l_case(line(1:6)) == "phase_")  then
             ndata=ndata+1
             ip(ndata)=i
          end if
       end do

       !---- Reading Phase Information ----!
       iph=1
       if (present(nphase)) iph=nphase
       if (present(Job_Info)) then
          n_ini=ip(iph)           !Updated values to handle non-conventional order
          n_end=ip(iph+1)
          call Get_Job_Info(file_dat,n_ini,n_end,Job_info)
       end if

       !---- Reading Cell Parameters ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)
       if(present(CFrame)) then
         call read_File_Cell(file_dat,n_ini,n_end,Cell,CFrame) !Read and construct Cell
       else
         call read_File_Cell(file_dat,n_ini,n_end,Cell) !Read and construct Cell
       end if
       if (err_form) return

       !---- Reading Space Group Information ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)
       call read_File_Spg (file_dat,n_ini,n_end,Spp)
       if (err_form) then !Try to read symmetry operators or generators
         ngen=0
         nsym=0
         do i=n_ini, n_end
           line=l_case(adjustl(file_dat(i)))
           if(line(1:4) == "symm") nsym=nsym+1
           if(line(1:3) == "gen")  ngen=ngen+1
         end do
         if(ngen > 0) then
           k=0
           do i=n_ini, n_end
             line=l_case(adjustl(file_dat(i)))
             if(line(1:3) == "gen")  then
              k=k+1
              gen(k)=adjustl(line(5:))
             end if
           end do
           call Set_SpaceGroup(" ",SpG,gen,ngen,"gen")   !Construct the space group from generators
         else if (nsym > 0) then
           k=0
           do i=n_ini, n_end
             line=l_case(adjustl(file_dat(i)))
             if(line(1:4) == "symm")  then
              k=k+1
              gen(k)=adjustl(line(6:))
             end if
           end do
           call Set_SpaceGroup(" ",SpG,gen,nsym,"fix")  !Construct the space group from fixed symmetry elements
         else
           return
         end if
       else
          call Set_SpaceGroup(Spp,SpG) !Construct the space group
       end if
       !---- Read Atoms Information ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)

       !---- Calculating number of Atoms in the Phase ----!
       do i=n_ini,n_end
          line=adjustl(file_dat(i))
          if (l_case(line(1:4)) == "atom")  nauas=nauas+1
       end do

       if (nauas > 0) then
          call Allocate_atom_list(nauas,A)  !allocation space for Atom list
          call read_File_Atom(file_dat,n_ini,n_end,A)
          if (err_form) return

          do i=1,A%natoms
             vet=A%atom(i)%x
             A%atom(i)%Mult=Get_Multip_Pos(vet,SpG)
             if(A%atom(i)%occ < epsv) A%atom(i)%occ=real(A%atom(i)%Mult)/real(SpG%Multip)
             if (A%atom(i)%thtype == "aniso") then
                select case (A%atom(i)%Utype)
                   case ("u_ij")
                      A%atom(i)%u(1:6) =  Convert_U_Betas(A%atom(i)%u(1:6),Cell)
                   case ("b_ij")
                      A%atom(i)%u(1:6) =  Convert_B_Betas(A%atom(i)%u(1:6),Cell)
                end select
                A%atom(i)%Utype="beta"
             end if
          end do
       end if

       return
    End Subroutine Readn_Set_XTal_CFL

    !!--++
    !!--++ Subroutine Readn_Set_XTal_CFL_Molec(file_dat, nlines, Molcrys)
    !!--++    character(len=*),dimension(:),  intent(in)     :: file_dat
    !!--++    integer,                        intent(in)     :: nlines
    !!--++    Type (Molecular_Crystal_Type),  intent(in out) :: Molcrys
    !!--++
    !!--++ (Private)
    !!--++ Read Molecule Information in a CFL
    !!--++
    !!--++ Update: April - 2005
    !!
    Subroutine Readn_Set_XTal_CFL_Molec(file_dat, nlines, Molcrys, Nphase)
       !---- Arguments ----!
       character(len=*),dimension(:),  intent(in)     :: file_dat
       integer,                        intent(in)     :: nlines
       type (Molecular_Crystal_Type),  intent(in out) :: Molcrys
       Integer, optional,              intent(in)     :: Nphase

       !---- Local variables ----!
       character(len=132)            :: line
       integer                       :: i,n,nmol,npos,n_ini,n_end,ierr,nauas, iph, ndata
       integer, parameter               :: maxph=21  !Maximum number of phases "maxph-1"
       integer, dimension(maxph)        :: ip
       real(kind=cp)                 :: theta,phi,chi
       real(kind=cp), dimension(3)   :: x1f,x2f,x3f
       real(kind=cp), dimension(3,3) :: EuM

       !---- Standard CrysFML file *.CFL ----!
       nauas=0
       ndata=0
       ip=nlines
       ip(1)=1

       !---- Calculating number of Phases ----!
       do i=1,nlines
          line=adjustl(file_dat(i))
          if (l_case(line(1:6)) == "phase_")  then
             ndata=ndata+1
             ip(ndata)=i
          end if
       end do

       !---- Reading Phase Information ----!

       if (present(nphase)) then
           iph=nphase
       else
           iph=1
       end if

       n_ini=ip(iph)
       n_end=ip(iph+1)

       !---- Detecting the Molecules defined in the file ----!
       nmol=0
       do i=n_ini,n_end
          line=u_case(adjustl(file_dat(i)))
          if (line(1:1) == " ") cycle
          if (line(1:1) == "!") cycle
          npos=index(line,"MOLE")
          if (npos /= 0) nmol=nmol+1
       end do
       if (nmol==0) return

       !---- Allocating Memory for all molecules ----!
       if (allocated(molcrys%mol)) deallocate(molcrys%mol)
       molcrys%n_mol=nmol
       allocate(molcrys%mol(nmol))

       !---- Reading Molecules ----!

       do n=1,nmol
          !---- Read ----!
          do i=n_ini,n_end
             line=u_case(adjustl(file_dat(i)))
             if (line(1:1) == " ") cycle
             if (line(1:1) == "!") cycle
             npos=index(line,"MOLE")
             if (npos == 0) cycle
             call read_molecule(file_dat,n_ini,n_end,molcrys%mol(n))
             err_form=err_molec
             ERR_Form_Mess=err_molec_mess
             if (err_form) then
                molcrys%n_mol=n-1
                return
             end if
             exit
          end do

          !---- Search for three points (fractional coordinates) ----!
          !---- defining a Cartesian frame                       ----!
          do
             if (n_ini > n_end) exit
             line=adjustl(file_dat(n_ini))
             if (u_case(line(1:9)) == "XYZ_FRAME") then
                read(unit=line(10:),fmt=*,iostat=ierr) x1f,x2f,x3f
                if (ierr == 0) then
                   call get_euler_from_fract(x1f,x2f,x3f,molcrys%Cell%Cr_Orth_cel,phi,theta,chi,EuM, Code="D")
                   molcrys%mol(n)%orient(1)= phi
                   molcrys%mol(n)%orient(2)= theta
                   molcrys%mol(n)%orient(3)= chi
                   molcrys%mol(n)%xcentre= x3f
                   call Set_euler_matrix(molcrys%mol(n)%rot_type, phi,theta,chi,EuM)
                   molcrys%mol(n)%Euler=EuM
                   molcrys%mol(n)%is_EulerMat=.true.
                   molcrys%mol(n)%in_Xtal=.true.
                end if
                n_ini=n_ini+1
                exit
             else
                if (u_case(line(1:4)) =="MOLE") exit
                n_ini=n_ini+1
             end if
          end do

       end do

       return
    End Subroutine Readn_Set_XTal_CFL_Molec

    !!--++
    !!--++ Subroutine Readn_Set_XTal_CFL_Shub(file_dat,nlines,Cell,SpG,A,CFrame,NPhase,Job_Info)
    !!--++    character(len=*),dimension(:),intent(in)   :: file_dat
    !!--++    integer,                      intent(in)   :: nlines
    !!--++    Type (Crystal_Cell_Type),     intent(out)  :: Cell
    !!--++    Type (Magnetic_Space_Group_Type), intent(out)  :: SpG
    !!--++    Type (atom_list_type),        intent(out)  :: A
    !!--++    character(len=*),    optional,intent(in)   :: CFrame
    !!--++    Integer,             optional,intent( in)  :: Nphase
    !!--++    Type(Job_Info_type), optional,intent(out)  :: Job_Info
    !!--++
    !!--++ (Private)
    !!--++ Read and Set Crystal Information in a CFL File
    !!--++
    !!--++ Update: April - 2005
    !!
    Subroutine Readn_Set_XTal_CFL_Shub(file_dat,nlines,Cell,SpG,A,CFrame,NPhase,Job_Info)
       !---- Arguments ----!
       character(len=*),dimension(:),    intent(in)   :: file_dat
       integer,                          intent(in)   :: nlines
       Type (Crystal_Cell_Type),         intent(out)  :: Cell
       Type (Magnetic_Space_Group_Type), intent(out)  :: SpG
       Type (atom_list_type),            intent(out)  :: A
       character(len=*),        optional,intent(in)   :: CFrame
       Integer,                 optional,intent( in)  :: Nphase
       Type(Job_Info_type),     optional,intent(out)  :: Job_Info

       !---- Local variables ----!
       character(len=132)               :: line
       character(len= 50)               :: Spp,setting
       !character(len= 40),dimension(192):: gen
       integer                          :: i, nauas, ndata, iph, n_ini,n_end,k !,ngen
       integer, parameter               :: maxph=21  !Maximum number of phases "maxph-1"
       integer, dimension(maxph)        :: ip

       real(kind=cp),dimension(3):: vet

       !---- Standard CrysFML file *.CFL ----!
       nauas=0
       ndata=0
       ip=nlines
       ip(1)=1

       !---- Calculating number of Phases ----!
       do i=1,nlines
          line=adjustl(file_dat(i))
          if (l_case(line(1:6)) == "phase_")  then
             ndata=ndata+1
             ip(ndata)=i
          end if
       end do

       !---- Reading Phase Information ----!
       iph=1
       if (present(nphase)) iph=nphase
       if (present(Job_Info)) then
          n_ini=ip(iph)           !Updated values to handle non-conventional order
          n_end=ip(iph+1)
          call Get_Job_Info(file_dat,n_ini,n_end,Job_info)
       end if

       !---- Reading Cell Parameters ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)
       if(present(CFrame)) then
         call read_File_Cell(file_dat,n_ini,n_end,Cell,CFrame) !Read and construct Cell
       else
         call read_File_Cell(file_dat,n_ini,n_end,Cell) !Read and construct Cell
       end if
       if (err_form) return

       !---- Reading Space Group Information ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)
       call read_File_Spg (file_dat,n_ini,n_end,Spp)
       i=index(Spp,"{")
       k=len_trim(Spp)
       setting=" "
       if(i /= 0) then
         setting=Spp(i+1:k-1)
         Spp=Spp(1:i-1)
       end if
       call Set_Magnetic_Space_Group(Spp,setting,Spg) !Construct the magnetic space group

       !---- Read Atoms Information ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)

       !---- Calculating number of Atoms in the Phase ----!
       do i=n_ini,n_end
          line=adjustl(file_dat(i))
          if (l_case(line(1:4)) == "atom")  nauas=nauas+1
       end do

       if (nauas > 0) then
          call Allocate_atom_list(nauas,A)  !allocation space for Atom list
          call read_File_Atom(file_dat,n_ini,n_end,A)
          if (err_form) return

          do i=1,A%natoms
             vet=A%atom(i)%x
             A%atom(i)%Mult=Get_Multip_Pos(vet,SpG)
             if(A%atom(i)%occ < epsv) A%atom(i)%occ=real(A%atom(i)%Mult)/real(SpG%Multip)
             if (A%atom(i)%thtype == "aniso") then
                select case (A%atom(i)%Utype)
                   case ("u_ij")
                      A%atom(i)%u(1:6) =  Convert_U_Betas(A%atom(i)%u(1:6),Cell)
                   case ("b_ij")
                      A%atom(i)%u(1:6) =  Convert_B_Betas(A%atom(i)%u(1:6),Cell)
                end select
                A%atom(i)%Utype="beta"
             end if
          end do
       end if

       return
    End Subroutine Readn_Set_XTal_CFL_Shub
    !!--++
    !!--++ Subroutine Readn_Set_XTal_CIF(file_dat, nlines, Cell, Spg, A, CFrame, NPhase)
    !!--++    character(len=*),dimension(:),intent(in)   :: file_dat
    !!--++    integer,                      intent(in)   :: nlines
    !!--++    Type (Crystal_Cell_Type),     intent(out)  :: Cell
    !!--++    Type (Space_Group_Type),      intent(out)  :: SpG
    !!--++    Type (atom_list_type),        intent(out)  :: A
    !!--++    Character(len=*),    optional,intent( in)  :: CFrame
    !!--++    Integer,             optional,intent( in)  :: Nphase
    !!--++
    !!--++ (Private)
    !!--++ Read and Set Crystal Information in a CIF File
    !!--++
    !!--++ Update: April - 2005
    !!
    Subroutine Readn_Set_XTal_CIF(file_dat, nlines, Cell, Spg, A, CFrame, NPhase)
       !---- Arguments ----!
       character(len=*),dimension(:),intent(in)   :: file_dat
       integer,                      intent(in)   :: nlines
       Type (Crystal_Cell_Type),     intent(out)  :: Cell
       Type (Space_Group_Type),      intent(out)  :: SpG
       Type (atom_list_type),        intent(out)  :: A
       Character(len=*),    optional,intent( in)  :: CFrame
       Integer,             optional,intent( in)  :: Nphase

       !---- Local Variables ----!
       character(len=132)                :: line
       character(len= 20)                :: Spp
       character(len=60), dimension(192) :: symm_car

       integer                   :: i, nauas, ndata, iph, n_ini,n_end,noper
       integer, parameter        :: maxph=250  !Maximum number of phases "maxph-1"
       integer, dimension(maxph) :: ip

       real(kind=cp),dimension(6):: vet,vet2

       ip=nlines
       ip(1)=1

       !---- First determine if there is more than one structure ----!
       do i=1,nlines
          line=adjustl(file_dat(i))
          if (l_case(line(1:5)) == "data_" .and. l_case(line(1:11)) /= "data_global" )  then
             n_ini=i
             ip(1)=i
             exit
          end if
       end do

       ndata=0
       do i=n_ini,nlines
          line=adjustl(file_dat(i))
          if (l_case(line(1:5)) == "data_")  then
             ndata=ndata+1
             if (ndata > maxph-1) then
                err_form=.true.
                ERR_Form_Mess=" => Too many phases in this file "
                return
             end if
             ip(ndata)=i   !Pointer to the number of the line starting a single phase
          end if
       end do

       iph=1
       if (present(nphase)) iph=nphase

       !---- Read Cell Parameters ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)
       call Read_Cif_Cell(file_dat,n_ini,n_end,vet,vet2)
       if (err_form) return
       if(present(CFrame)) then
         call Set_Crystal_Cell(vet(1:3),vet(4:6),Cell,CFrame,vet2(1:3),vet2(4:6))
       else
         call Set_Crystal_Cell(vet(1:3),vet(4:6),Cell,"A",vet2(1:3),vet2(4:6))
       end if
       !---- Read Atoms Information ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)
       call Read_Cif_Atom(file_dat,n_ini,n_end,nauas,A)
       if (err_form) return

       !---- SpaceGroup Information ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)
       call Read_Cif_Hm(file_dat,n_ini,n_end,Spp)

       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)
       if (len_trim(Spp) == 0) call Read_Cif_Hall(file_dat,n_ini,n_end,Spp)

       if (len_trim(Spp) == 0) then
          n_ini=ip(iph)           !Updated values to handle non-conventional order
          n_end=ip(iph+1)
          call Read_Cif_Symm(file_dat,n_ini,n_end,noper,symm_car)

          if (noper ==0) then
             err_form=.true.
             ERR_Form_Mess=" => No Space Group/No Symmetry information in this file "
             return
          else
             call Set_SpaceGroup("  ",SpG,symm_car,noper,"GEN")
          end if
       else
          call Set_SpaceGroup(Spp,SpG) !Construct the space group
       end if

       !---- Modify occupation factors and set multiplicity of atoms
       !---- in order to be in agreement with the definitions of Sfac in CrysFML
       !---- Convert Us to Betas and Uiso to Biso
       do i=1,A%natoms
          vet(1:3)=A%atom(i)%x
          A%atom(i)%Mult=Get_Multip_Pos(vet(1:3),SpG)
          A%atom(i)%Occ=A%atom(i)%Occ*real(A%atom(i)%Mult)/max(1.0,real(SpG%Multip))
          if(A%atom(i)%occ < epsv) A%atom(i)%occ=real(A%atom(i)%Mult)/max(1.0,real(SpG%Multip))

          select case (A%atom(i)%thtype)
             case ("isotr")
                A%atom(i)%biso= A%atom(i)%ueq*78.95683521

             case ("aniso")
                select case (A%atom(i)%Utype)
                   case ("u_ij")
                      A%atom(i)%u(1:6) =  Convert_U_Betas(A%atom(i)%u(1:6),Cell)
                   case ("b_ij")
                      A%atom(i)%u(1:6) = Convert_B_Betas(A%atom(i)%u(1:6),Cell)
                end select
                A%atom(i)%Utype="beta"

             case default
                A%atom(i)%biso = A%atom(i)%ueq*78.95683521
                A%atom(i)%thtype = "isotr"
          end select
       end do

       return
    End Subroutine Readn_Set_XTal_CIF

    !!--++
    !!--++ Subroutine Readn_Set_XTal_PCR(file_dat, nlines, Cell, Spg, A, CFrame, NPhase)
    !!--++    character(len=*),dimension(:),intent(in)   :: file_dat
    !!--++    integer,                      intent(in)   :: nlines
    !!--++    Type (Crystal_Cell_Type),     intent(out)  :: Cell
    !!--++    Type (Space_Group_Type),      intent(out)  :: SpG
    !!--++    Type (atom_list_type),        intent(out)  :: A
    !!--++    character(len=*),    optional,intent(in)   :: CFrame
    !!--++    Integer,             optional,intent( in)  :: Nphase
    !!--++
    !!--++ (Private)
    !!--++ Read and Set Crystal Information in a PCR File
    !!--++
    !!--++ Update: 17/05/2010
    !!
    Subroutine Readn_Set_XTal_PCR(file_dat, nlines, Cell, Spg, A, CFrame, NPhase)
       !---- Arguments ----!
       character(len=*),dimension(:),intent(in)   :: file_dat
       integer,                      intent(in)   :: nlines
       Type (Crystal_Cell_Type),     intent(out)  :: Cell
       Type (Space_Group_Type),      intent(out)  :: SpG
       Type (atom_list_type),        intent(out)  :: A
       character(len=*),    optional,intent(in)   :: CFrame
       Integer,             optional,intent(in)   :: Nphase

       !---- Local Variables ----!
       logical                           :: multi,ask_phase,is_codewords
       character(len=132)                :: line
       character(len= 20)                :: Spp, label
       integer                           :: i,j, k,iv, nauas, ndata, iph, n_ini,n_end, nlong1
       integer, parameter                :: maxph=21  !Maximum number of phases "maxph-1"
       integer, dimension(maxph)         :: ip
       integer, dimension(30)            :: ivet

       real(kind=cp),dimension(30)       :: vet

       ip=nlines
       ip(1)=1

       !> Simple / Multi format
       multi=.false.
       do i=1,nlines
          line=adjustl(file_dat(i))
          if (line(1:1) =='!' .or. line(1:1)==' ') cycle
          if (index(line,'NPATT ') <=0) cycle
          multi=.true.
       end do

       !> Number of Phases
       if (.not. multi) then
          do i=2,nlines
             line=adjustl(file_dat(i))
             if (line(1:1) =='!' .or. line(1:1)==' ') cycle
             call getnum(line,vet,ivet,iv)
             if (iv > 3) then
                iph=ivet(3)
                exit
             end if
          end do

       else
          do i=1,nlines
             line=adjustl(file_dat(i))
             if (line(1:4) /='!Nph') cycle

             line=adjustl(file_dat(i+1))
             call getnum(line,vet,ivet,iv)
             if (iv > 1) then
                iph=ivet(1)
                exit
             end if
          end do
       end if
       if (iph == 0) then
          err_form=.true.
          ERR_Form_Mess=" No Phase information was found in this PCR file. Please, check it! "
          return
       end if

       !> Locate where begin each Phase
       k=0
       ask_phase=.true.

       do i=1,nlines
          line=adjustl(file_dat(i))
          if (ask_phase) then
             if (index(line,'Data for PHASE') <= 0) cycle
          else
             if (line(1:1) /='!') then
                k=k+1
                ip(k)=i
                if (k == iph) exit

                ask_phase=.true.
             end if
             cycle
          end if
          ask_phase=.false.
       end do
       if (iph /= k) then
          err_form=.true.
          ERR_Form_Mess=" Locating Phases failed in this PCR. Please, check it!"
          return
       end if

       !> Select the Phase
       iph=1
       if (present(nphase)) iph=nphase
       n_ini=ip(iph)
       n_end=ip(iph+1)

       !---- Read Cell Parameters ----!
       do i=n_ini,n_end
          if (index(file_dat(i),'alpha') /=0 .and. index(file_dat(i),'gamma') /=0) then
             do j=i+1,n_end
                line=adjustl(file_dat(j))
                if (line(1:1) == '!' .or. line(1:1) == ' ') cycle
                iv=index(line,'#')
                if (iv > 1) line=line(1:iv-1)

                call getnum(line, vet, ivet,iv)
                if (iv /= 6) then
                   err_form=.true.
                   ERR_Form_Mess=" => Problems reading Cell Parameters on PCR file "
                   return
                end if
                if(present(CFrame)) then
                  call Set_Crystal_Cell(vet(1:3),vet(4:6),Cell,CFrame)
                else
                  call Set_Crystal_Cell(vet(1:3),vet(4:6),Cell)
                end if
                exit
             end do
             exit
          end if
       end do

       !---- SpaceGroup Information ----!
       Spp=' '
       do i=n_ini,n_end
          line=adjustl(file_dat(i))
          if (line(1:1) == '!' .or. line(1:1)==' ') cycle
          if (index(file_dat(i),'<--Space') /=0) then
             j=index(file_dat(i),'<--Space')
             Spp=adjustl(file_dat(i)(1:j-1))
             if (len_trim(Spp) <= 0) then
                err_form=.true.
                ERR_Form_Mess=" => Problems reading Space group on PCR file "
                return
             end if
             call Set_SpaceGroup(Spp,SpG) !Construct the space group
             exit
          end if
       end do

       !---- Read Atoms Information ----!
       do i=n_ini,n_end
          line=adjustl(file_dat(i))
          if (line(1:4) /= '!Nat') cycle
          do j=i+1,n_end
             line=adjustl(file_dat(j))
             if (line(1:1) == '!' .or. line(1:1)==' ') cycle
             call getnum(line(1:5),vet,ivet,iv)
             ndata=ivet(1)
             exit
          end do
          exit
       end do

       if (ndata > 0) then
          call allocate_atom_list(ndata,A)

          is_codewords=.false.
          nauas=0

          do i=n_ini,n_end
             line=adjustl(file_dat(i))
             if (index(line,'!Atom') == 0 .or. index(line,'Typ') == 0) cycle

             do j=i+1,n_end
                line=adjustl(file_dat(j))
                if (line(1:1) == '!' .or. line(1:1)==' ') cycle
                if (is_codewords) then
                   is_codewords=.false.
                   cycle
                end if

                iv=index(line,'#')
                if (iv > 1) line=line(1:iv-1)

                nauas=nauas+1
                ! Atom Label
                call cutst(line,nlong1,label)
                A%atom(nauas)%lab=trim(label)

                ! Atom Type
                call cutst(line,nlong1,label)
                A%Atom(nauas)%chemsymb=U_case(label(1:1))//L_case(label(2:2))

                ! Atom Coordinates,Biso and Occ
                call getnum(line,vet,ivet,iv)
                if (iv < 5) then    !Line reading for the second time anisotropic temperature factors
                   nauas = nauas -1 !see below
                   is_codewords=.true.
                   cycle
                end if

                A%atom(nauas)%x=vet(1:3)
                A%atom(nauas)%Mult=Get_Multip_Pos(vet(1:3),SpG)
                A%atom(nauas)%biso=vet(4)
                A%atom(nauas)%occ=vet(5)
                A%atom(nauas)%thtype='isotr'
                A%atom(nauas)%Utype="beta"
                if (ivet(8) == 2) then    ! Anisotropic reading
                   A%atom(nauas)%thtype='aniso'
                   call getnum(file_dat(j+2),vet,ivet,iv)
                   A%atom(nauas)%u(1:6)=vet(1:6)
                end if
                is_codewords=.true.
                if (nauas == ndata) exit
             end do
             exit
          end do
       end if

       return
    End Subroutine Readn_Set_XTal_PCR

    !!--++
    !!--++ Subroutine Readn_Set_XTal_SHX(file_dat,nlines,Cell,SpG,A,CFrame)
    !!--++    character(len=*),dimension(:),intent(in)   :: file_dat
    !!--++    integer,                      intent(in)   :: nlines
    !!--++    Type (Crystal_Cell_Type),     intent(out)  :: Cell
    !!--++    Type (Space_Group_Type),      intent(out)  :: SpG
    !!--++    Type (Atom_list_type),        intent(out)  :: A
    !!--++    Character(len=*), optional,   intent(in)   :: CFrame
    !!--++
    !!--++ (Private)
    !!--++ Read and Set Crystal Information in a Shelx File
    !!--++
    !!--++ Update: April - 2005
    !!
    Subroutine Readn_Set_XTal_SHX(file_dat,nlines,Cell,SpG,A,CFrame)
       !---- Arguments ----!
       character(len=*),dimension(:),intent(in)   :: file_dat
       integer,                      intent(in)   :: nlines
       Type (Crystal_Cell_Type),     intent(out)  :: Cell
       Type (Space_Group_Type),      intent(out)  :: SpG
       Type (Atom_list_type),        intent(out)  :: A
       Character(len=*), optional,   intent(in)   :: CFrame

       !---- Local Variables ----!
       character(len=60), dimension(192) :: symm_car
       character(len=2),  dimension(15)  :: elem_atm
       integer                           :: i,n_ini, n_end, nl, noper
       integer                           :: n_elem_atm, n_fvar
       real(kind=cp), dimension(6)       :: vet,vet2
       real(kind=cp), dimension(10)      :: fvar

       n_ini=1
       n_end=nlines

       !---- CELL / ZERR ----!
       call Read_Shx_Cell(file_dat,n_ini,n_end,vet,vet2)
       if(present(CFrame)) then
         call Set_Crystal_Cell(vet(1:3),vet(4:6),Cell,CFrame,vet2(1:3),vet2(4:6))
       else
         call Set_Crystal_Cell(vet(1:3),vet(4:6),Cell,"A",vet2(1:3),vet2(4:6))
       end if

       !---- OBTAIN SPACE GROUP (LATT / SYMM) ----!
       call Read_Shx_Latt(file_dat,n_ini,n_end,nl)
       call Read_Shx_Symm(file_dat,n_ini,n_end,noper,symm_car)
       if (nl > 0) then
          noper=noper+1
          symm_car(noper)="-X,-Y,-Z"
       end if
       select case (abs(nl))
          case (2) ! I
             noper=noper+1
             symm_car(noper)="X+1/2,Y+1/2,Z+1/2"
          case (3) ! Rom, Hex
             noper=noper+1
             symm_car(noper)="X+2/3,Y+1/3,Z+1/3"
             noper=noper+1
             symm_car(noper)="X+1/3,Y+2/3,Z+2/3"
          case (4) ! F
             noper=noper+1
             symm_car(noper)="X,Y+1/2,Z+1/2"
          case (5) ! A
             noper=noper+1
             symm_car(noper)="X,Y+1/2,Z+1/2"
             noper=noper+1
             symm_car(noper)="X+1/2,Y,Z+1/2"
             noper=noper+1
             symm_car(noper)="X+1/2,Y+1/2,Z"
          case (6) ! B
             noper=noper+1
             symm_car(noper)="X+1/2,Y,Z+1/2"
          case (7) ! C
             noper=noper+1
             symm_car(noper)="X+1/2,Y+1/2,Z"
       end select ! nl
       call set_spacegroup(" ",SPG,symm_car,noper,"gen")

       !---- ATOMS ----!
       call Read_Shx_Cont(file_dat,n_ini,n_end,n_elem_atm,elem_atm)
       call Read_Shx_Fvar(file_dat,n_ini,n_end,n_fvar,fvar)
       call Read_Shx_Atom(file_dat,n_ini,n_end,n_fvar,fvar,elem_atm,cell,A)
       if (err_form) return

       !---- Convert Us to Betas and Uiso to Biso
       do i=1,A%natoms
          vet(1:3)=A%atom(i)%x
          A%atom(i)%Mult=Get_Multip_Pos(vet(1:3),SpG)

          select case (A%atom(i)%thtype)
             case ("isotr")
                A%atom(i)%biso= A%atom(i)%ueq*78.95683521

             case ("aniso")
                A%atom(i)%ueq=U_Equiv(cell,a%atom(i)%u(1:6))  ! Uequi
                A%atom(i)%biso= A%atom(i)%ueq*78.95683521
                select case (A%atom(i)%Utype)
                   case ("u_ij")
                      A%atom(i)%u(1:6) =  Convert_U_Betas(A%atom(i)%u(1:6),Cell)
                   case ("b_ij")
                      A%atom(i)%u(1:6) = Convert_B_Betas(A%atom(i)%u(1:6),Cell)
                end select
                A%atom(i)%Utype="beta"

             case default
                A%atom(i)%ueq=0.05
                A%atom(i)%biso = A%atom(i)%ueq*78.95683521
                A%atom(i)%thtype = "isotr"
          end select
       end do

       return
    End Subroutine Readn_Set_XTal_SHX

    !!--++
    !!--++ Subroutine Readn_Set_Xtal_Structure_Molcr(filenam,Molcrys,Mode,Iphase, Job_Info, file_list,CFrame)
    !!--++    character(len=*),              intent( in)     :: filenam  ! In -> Name of the file
    !!--++    Type (Molecular_Crystal_Type), intent(out)     :: Molcrys  ! Molecular crytal
    !!--++    Character(len=*),    optional, intent( in)     :: Mode     ! In -> if Mode="CIF" filenam
    !!--++                                                                       is of CIF type format
    !!--++    Integer,             optional, intent( in)     :: Iphase   ! Number of the phase.
    !!--++    Type(Job_Info_type), optional, intent(out)     :: Job_Info ! Diffaction conditions
    !!--++    Type(file_list_type),optional, intent(in out)  :: file_list! Complete file to be used by
    !!--++                                                              the calling program or other procedures
    !!--++    Character(len=*),    optional, intent(in)      :: CFrame
    !!--++    Overloaded
    !!--++    Subroutine to read and input file and construct the crystal structure
    !!--++    in terms of the ofjects Cell, SpG and A. The optional argument Iphase is an integer
    !!--++    telling to the program to read the phase number Iphase in the case of the presence
    !!--++    of more than one phase. If absent only the first phase is read.
    !!--++
    !!--++ Update: April - 2005
    !!
    Subroutine Readn_Set_Xtal_Structure_Molcr(filenam,Molcrys,Mode,Iphase,Job_Info,file_list,CFrame)
       !---- Arguments ----!
       character(len=*),              intent( in)     :: filenam
       Type (Molecular_Crystal_Type), intent(out)     :: Molcrys
       Character(len=*),     optional,intent( in)     :: Mode
       Integer,              optional,intent( in)     :: Iphase
       Type(Job_Info_type),  optional,intent(out)     :: Job_Info
       Type(file_list_type), optional,intent(in out)  :: file_list
       Character(len=*),     optional,intent(in)      :: CFrame
       !---- Local variables -----!
       Type (Atom_list_type)                         :: A
       character(len=132), allocatable, dimension(:) :: file_dat
       character(len=3)                              :: modec
       integer                                       :: i,nlines


       call init_err_form()

       nlines=0
       if (present(file_list)) nlines=file_list%nlines

       !---- Number of Lines in the input file ----!
       if(nlines == 0) then
           call Number_Lines(trim(filenam), nlines)
           if (nlines==0) then
              err_form=.true.
              ERR_Form_Mess="The file "//trim(filenam)//" contains nothing"
              return
           else
              if (allocated(file_dat)) deallocate( file_dat)
              allocate( file_dat(nlines))
              call reading_Lines(trim(filenam),nlines,file_dat)
           end if
           if (present(file_list)) then
              file_list%nlines=nlines
              if (allocated(file_list%line)) deallocate(file_list%line)
              allocate(file_list%line(nlines))
              file_list%line=file_dat
           end if
       else
           if (allocated(file_dat)) deallocate( file_dat)
           allocate( file_dat(nlines))
           file_dat=file_list%line
       end if


       !---- Define the type of file: CIF, CFL, RES,... ----!
       modec=" "
       if (present(mode)) modec=l_case(mode(1:3))

       select case(modec)
           case("cif")
              if (present(iphase)) then
                 if(present(CFrame)) then
                   call readn_set_xtal_cif(file_dat,nlines,molcrys%Cell,molcrys%Spg, A,CFrame,NPhase=IPhase)
                 else
                   call readn_set_xtal_cif(file_dat,nlines,molcrys%Cell,molcrys%Spg, A,NPhase=IPhase)
                 end if
              else
                 if(present(CFrame)) then
                   call readn_set_xtal_cif(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,CFrame)
                 else
                   call readn_set_xtal_cif(file_dat,nlines,molcrys%Cell,molcrys%Spg,A)
                 end if
              end if

           case("pcr")
              if (present(iphase)) then
                 if(present(CFrame)) then
                   call readn_set_xtal_pcr(file_dat,nlines,molcrys%Cell,molcrys%Spg, A,CFrame,NPhase=IPhase)
                 else
                   call readn_set_xtal_pcr(file_dat,nlines,molcrys%Cell,molcrys%Spg, A,NPhase=IPhase)
                 end if
              else
                 if(present(CFrame)) then
                   call readn_set_xtal_pcr(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,CFrame)
                 else
                   call readn_set_xtal_pcr(file_dat,nlines,molcrys%Cell,molcrys%Spg,A)
                 end if
              end if

           case("shx")
              if(present(CFrame)) then
                call readn_set_xtal_shx(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,CFrame)
              else
                call readn_set_xtal_shx(file_dat,nlines,molcrys%Cell,molcrys%Spg,A)
              end if
           case default
              !---- CFL Format ----!
              if (present(Job_Info)) then
                 if (present(iphase)) then
                    if(present(CFrame)) then
                      call readn_set_xtal_cfl(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,CFrame,NPhase=IPhase,Job_Info=Job_Info)
                    else
                      call readn_set_xtal_cfl(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,NPhase=IPhase,Job_Info=Job_Info)
                    end if
                 else
                    if(present(CFrame)) then
                      call readn_set_xtal_cfl(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,CFrame,Job_Info=Job_Info)
                    else
                      call readn_set_xtal_cfl(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,Job_Info=Job_Info)
                    end if
                 end if
              else
                 if (present(iphase)) then
                    if(present(CFrame)) then
                      call readn_set_xtal_cfl(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,CFrame,NPhase=IPhase)
                    else
                      call readn_set_xtal_cfl(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,NPhase=IPhase)
                    end if
                 else
                    if(present(CFrame)) then
                      call readn_set_xtal_cfl(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,CFrame)
                    else
                      call readn_set_xtal_cfl(file_dat,nlines,molcrys%Cell,molcrys%Spg,A)
                    end if
                 end if
              end if
              !---- Reading molecules ----!
              if (present(iphase)) then
                call readn_set_xtal_cfl_molec(file_dat,nlines,molcrys,NPhase=IPhase)
              else
                call readn_set_xtal_cfl_molec(file_dat,nlines,molcrys)
              end if

       end select
       if (err_form) return

       !---- Passing from Atom_List_Type -> Molcrys ----!
       molcrys%n_free=A%natoms
       if (A%natoms > 0) then
          if (allocated(molcrys%Atm)) deallocate(molcrys%Atm)
          allocate(molcrys%Atm(A%natoms))
          molcrys%Atm=A%Atom
       end if

       call deallocate_atom_list(A)

       !---- Testing if Xtal was defined ----!
       if (all(molcrys%cell%cell > 0.0)) then
          do i=1,molcrys%n_mol
             if (.not. molcrys%mol(i)%in_xtal) then
                 molcrys%mol(i)%in_xtal=.true.
             end if
          end do
       end if

       return
    End Subroutine Readn_Set_Xtal_Structure_Molcr

    !!--++
    !!--++ Subroutine Readn_Set_Xtal_Structure_Split(filenam,Cell,SpG,A,Mode,Iphase,Job_Type,File_List,CFrame)
    !!--++    character(len=*),              intent( in)     :: filenam  ! In -> Name of the file
    !!--++    Type (Crystal_Cell_Type),      intent(out)     :: Cell     ! Out -> Cell object
    !!--++    Type (Space_Group_Type),       intent(out)     :: SpG      ! Out -> Space Group object
    !!--++    Type (atom_list_type),         intent(out)     :: A        ! Out -> Atom_List object
    !!--++    Character(len=*),    optional, intent( in)     :: Mode     ! In -> if Mode="CIF" filenam
    !!--++                                                                       is of CIF type format
    !!--++    Integer,             optional, intent( in)     :: Iphase   ! Number of the phase.
    !!--++    Type(Job_Info_type), optional, intent(out)     :: Job_Info ! Diffaction conditions
    !!--++    Type(file_list_type),optional, intent(in out)  :: file_list! Complete file to be used by
    !!--++                                                                 the calling program or other procedures
    !!--++    Character(len=*),    optional, intent( in)     :: CFrame   !Cartesian Frame
    !!--++
    !!--++    Overloaded
    !!--++    Subroutine to read and input file and construct the crystal structure
    !!--++    in terms of the ofjects Cell, SpG and A. The optional argument Iphase is an integer
    !!--++    telling to the program to read the phase number Iphase in the case of the presence
    !!--++    of more than one phase. If absent only the first phase is read.
    !!--++
    !!--++ Update: April - 2005
    !!
    Subroutine Readn_Set_Xtal_Structure_Split(filenam,Cell,SpG,A,Mode,Iphase,Job_Info,file_list,CFrame)
       !---- Arguments ----!
       character(len=*),             intent( in)     :: filenam
       Type (Crystal_Cell_Type),     intent(out)     :: Cell
       Type (Space_Group_Type),      intent(out)     :: SpG
       Type (atom_list_type),        intent(out)     :: A
       Character(len=*),    optional,intent( in)     :: Mode
       Integer,             optional,intent( in)     :: Iphase
       Type(Job_Info_type), optional,intent(out)     :: Job_Info
       Type(file_list_type),optional,intent(in out)  :: file_list
       Character(len=*),    optional,intent( in)     :: CFrame

       !---- Local variables -----!
       character(len=132), allocatable, dimension(:) :: file_dat
       character(len=3)                              :: modec
       integer                                       :: nlines

       call init_err_form()

       nlines=0
       if (present(file_list)) nlines=file_list%nlines

       !---- Number of Lines in the input file ----!
       if(nlines == 0) then
           call Number_Lines(trim(filenam), nlines)
           if (nlines==0) then
              err_form=.true.
              ERR_Form_Mess="The file "//trim(filenam)//" contains nothing"
              return
           else
              if (allocated(file_dat)) deallocate( file_dat)
              allocate( file_dat(nlines))
              call reading_Lines(trim(filenam),nlines,file_dat)
           end if
           if (present(file_list)) then
              file_list%nlines=nlines
              if (allocated(file_list%line)) deallocate(file_list%line)
              allocate(file_list%line(nlines))
              file_list%line=file_dat
           end if
       else
           if (allocated(file_dat)) deallocate( file_dat)
           allocate( file_dat(nlines))
           file_dat=file_list%line
       end if

       !---- Define the type of file: CIF, CFL, RES,... ----!
       modec=" "
       if (present(mode)) modec=l_case(mode(1:3))

       select case(modec)
           case("cif")
              if (present(iphase)) then
                 if(present(CFrame)) then
                   call readn_set_xtal_cif(file_dat,nlines,Cell,Spg, A,CFrame,NPhase=IPhase)
                 else
                   call readn_set_xtal_cif(file_dat,nlines,Cell,Spg, A,NPhase=IPhase)
                 end if
              else
                 if(present(CFrame)) then
                   call readn_set_xtal_cif(file_dat,nlines,Cell,Spg,A,CFrame)
                 else
                   call readn_set_xtal_cif(file_dat,nlines,Cell,Spg,A)
                 end if
              end if

           case("pcr")
              if (present(iphase)) then
                 if(present(CFrame)) then
                   call readn_set_xtal_pcr(file_dat,nlines,Cell,Spg, A,CFrame,NPhase=IPhase)
                 else
                   call readn_set_xtal_pcr(file_dat,nlines,Cell,Spg, A,NPhase=IPhase)
                 end if
              else
                 if(present(CFrame)) then
                   call readn_set_xtal_pcr(file_dat,nlines,Cell,Spg,A,CFrame)
                 else
                   call readn_set_xtal_pcr(file_dat,nlines,Cell,Spg,A)
                 end if
              end if

           case("shx")
              if(present(CFrame)) then
                call readn_set_xtal_shx(file_dat,nlines,Cell,Spg,A,CFrame)
              else
                call readn_set_xtal_shx(file_dat,nlines,Cell,Spg,A)
              end if

           case default
              !---- CFL Format ----!
              if (present(Job_Info)) then
                 if (present(iphase)) then
                    if(present(CFrame)) then
                      call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,CFrame,NPhase=IPhase,Job_Info=Job_Info)
                    else
                      call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,NPhase=IPhase,Job_Info=Job_Info)
                    end if
                 else
                    if(present(CFrame)) then
                      call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,CFrame,Job_Info=Job_Info)
                    else
                      call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,Job_Info=Job_Info)
                    end if
                 end if
              else
                 if (present(iphase)) then
                    if(present(CFrame)) then
                      call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,CFrame,NPhase=IPhase)
                    else
                      call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,NPhase=IPhase)
                    end if
                 else
                    if(present(CFrame)) then
                      call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,CFrame)
                    else
                      call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A)
                    end if
                 end if
              end if

       end select

    End Subroutine Readn_Set_Xtal_Structure_Split

    Subroutine Readn_Set_Xtal_Structure_Magn(filenam,Cell,SpG,A,Mode,Iphase,Job_Info,file_list,CFrame)
       !---- Arguments ----!
       character(len=*),                 intent( in)     :: filenam
       Type (Crystal_Cell_Type),         intent(out)     :: Cell
       Type (Magnetic_Space_Group_Type), intent(out)     :: SpG
       Type (atom_list_type),            intent(out)     :: A
       Character(len=*),    optional,    intent( in)     :: Mode
       Integer,             optional,    intent( in)     :: Iphase
       Type(Job_Info_type), optional,    intent(out)     :: Job_Info
       Type(file_list_type),optional,    intent(in out)  :: file_list
       Character(len=*),    optional,    intent( in)     :: CFrame
       !
       character(len=132), allocatable, dimension(:) :: file_dat
       character(len=3)                              :: modec
       integer                                       :: nlines

       call init_err_form()

       nlines=0
       if (present(file_list)) nlines=file_list%nlines

       !---- Number of Lines in the input file ----!
       if(nlines == 0) then
           call Number_Lines(trim(filenam), nlines)
           if (nlines==0) then
              err_form=.true.
              ERR_Form_Mess="The file "//trim(filenam)//" contains nothing"
              return
           else
              if (allocated(file_dat)) deallocate( file_dat)
              allocate( file_dat(nlines))
              call reading_Lines(trim(filenam),nlines,file_dat)
           end if
           if (present(file_list)) then
              file_list%nlines=nlines
              if (allocated(file_list%line)) deallocate(file_list%line)
              allocate(file_list%line(nlines))
              file_list%line=file_dat
           end if
       else
           if (allocated(file_dat)) deallocate( file_dat)
           allocate( file_dat(nlines))
           file_dat=file_list%line
       end if

       !---- Define the type of file: CIF, CFL, RES,... ----!
       modec=" "
       if (present(mode)) modec=l_case(mode(1:3))

       select case(modec)
           case("cif")

              call Readn_Set_Magnetic_Structure_MCIF(filenam,Cell,Spg,A)

           case default
              !---- CFL Format ----!
              if (present(Job_Info)) then
                 if (present(iphase)) then
                    if(present(CFrame)) then
                      call Readn_Set_XTal_CFL_Shub(file_dat,nlines,Cell,SpG,A,CFrame,NPhase=IPhase,Job_Info=Job_Info)
                    else
                      call Readn_Set_XTal_CFL_Shub(file_dat,nlines,Cell,Spg,A,NPhase=IPhase,Job_Info=Job_Info)
                    end if
                 else
                    if(present(CFrame)) then
                      call Readn_Set_XTal_CFL_Shub(file_dat,nlines,Cell,Spg,A,CFrame,Job_Info=Job_Info)
                    else
                      call Readn_Set_XTal_CFL_Shub(file_dat,nlines,Cell,Spg,A,Job_Info=Job_Info)
                    end if
                 end if
              else
                 if (present(iphase)) then
                    if(present(CFrame)) then
                      call Readn_Set_XTal_CFL_Shub(file_dat,nlines,Cell,Spg,A,CFrame,NPhase=IPhase)
                    else
                      call Readn_Set_XTal_CFL_Shub(file_dat,nlines,Cell,Spg,A,NPhase=IPhase)
                    end if
                 else
                    if(present(CFrame)) then
                      call Readn_Set_XTal_CFL_Shub(file_dat,nlines,Cell,Spg,A,CFrame)
                    else
                      call Readn_Set_XTal_CFL_Shub(file_dat,nlines,Cell,Spg,A)
                    end if
                 end if
              end if

       end select
    End Subroutine Readn_Set_Xtal_Structure_Magn
    !!----
    !!---- Subroutine Set_Magnetic_Space_Group(symb,setting,MSpg,parent,mcif,keepd,trn_to)
    !!----    character (len=*),                intent(in) :: symb        !  In -> String with the BNS symbol of the Shubnikov Group
    !!----    character (len=*),                intent(in ):: setting     !  In -> setting in the form -a,c,2b;1/2,0,0 (if empty no transformation is performed)
    !!----    Type (Magnetic_Space_Group_Type), intent(out):: MGp         ! Out -> Magnetic Space Group object
    !!----    character (len=*), optional,      intent(in ):: Parent      !  In -> Parent crystallographic group
    !!----    logical,  optional,               intent(in ):: mcif        !  In -> True if one wants to store the symbols as mx,my,mz
    !!----    logical,  optional,               intent(in ):: keepd       !  In -> True if one wants to keep the database allocated
    !!----    logical,  optional,               intent(in ):: trn_to      !  In -> True if the setting is from current TO standard setting
    !!----
    !!----    Subroutine constructing the object MGp from the BNS symbol by
    !!----    reading the database compiled by Harold T. Stokes and Branton J. Campbell
    !!----
    !!---- Created: November - 2016 (JRC)
    !!
    Subroutine Set_Magnetic_Space_Group(symb,setting,MSpg,parent,mcif,keepd,trn_to)
      character(len=*),               intent (in)  :: symb,setting
      type(Magnetic_Space_Group_Type),intent (out) :: MSpg
      character(len=*),optional,      intent (in)  :: parent
      logical,         optional,      intent (in)  :: mcif
      logical,         optional,      intent (in)  :: keepd
      logical,         optional,      intent (in)  :: trn_to
      !--- Local variables ---!
      integer                          :: i,j,m,k,n,L,ier,num,idem !,inv_time
      real(kind=cp)                    :: det
      !real(kind=cp), dimension(3)      :: orig
      real(kind=cp), dimension(3,3)    :: e !,S,Sinv
      integer, dimension(3,3)          :: identity
      character(len=256)               :: line,ShOp_symb
      logical                          :: change_setting,centring
      type(Magnetic_Space_Group_Type)  :: MGp

      call Init_Err_Form()
      call Allocate_DataBase()
      call read_magnetic_data()
      identity=0
      do i=1,3
        identity(i,i)=1
      end do
      e=identity
      !write(*,"(a)") trim(symb)//"  "//trim(setting)
      !if(present(parent)) write(*,"(a)") trim(Parent)
      !Check if the number of the magnetic group has been given
      !instead of the symbol
      read(unit=symb,fmt=*,iostat=ier) num
      if(ier /= 0) then
        num=0 !It is supposed that a symbol has been provided
        do i=1,magcount
          !write(*,"(i5,tr5,a)") i, spacegroup_label_bns(i)
          if(trim(symb) == trim(spacegroup_label_bns(i)) .or. &
             trim(symb) == trim(spacegroup_label_og(i))) then
            num=i
            exit
          end if
        end do
        if(num == 0) then
           write(unit=Err_Form_Mess,fmt="(a)") " => The BNS symbol: "//trim(symb)//" is illegal! "
           Err_Form=.true.
           if(.not. present(keepd)) call deAllocate_DataBase()
           return
        end if
      else
        if(num < 1 .or. num > magcount) then !magcount=1651
           write(unit=Err_Form_Mess,fmt="(a,i4,a)") " => The number of the Shubnikov group: ",num," is illegal!"
           Err_Form=.true.
           if(.not. present(keepd)) call deAllocate_DataBase()
           return
        end if
      end if
      if(len_trim(setting) == 0 .or. setting =='a,b,c;0,0,0') then
        change_setting=.false.
      else
        change_setting=.true.
      end if

      MGp%Sh_number=num
      MGp%BNS_number=nlabel_bns(num)
      MGp%OG_number= nlabel_og(num)
      MGp%BNS_symbol=spacegroup_label_bns(num)
      MGp%OG_symbol=spacegroup_label_og(num)
      MGp%MagType=magtype(num)
      !Setting the magnetic point group symbol from the BNS label
      m=0
      Select Case (MGp%MagType)
         Case(1,2,3)
           MGp%PG_Symbol=MGp%BNS_symbol(2:) !Remove the type of lattice
           do i=2,len_trim(MGp%BNS_symbol)
             m=m+1
             if(MGp%BNS_symbol(i:i) == "a" .or. MGp%BNS_symbol(i:i) == "b"  &
           .or. MGp%BNS_symbol(i:i) == "c" .or. MGp%BNS_symbol(i:i) == "d"  &
           .or. MGp%BNS_symbol(i:i) == "e" .or. MGp%BNS_symbol(i:i) == "g"  &
           .or. MGp%BNS_symbol(i:i) == "n") MGp%PG_Symbol(m:m)="m"
             if(MGp%BNS_symbol(i:i) == "_") MGp%PG_Symbol(m:m+1)=" "
           end do
           MGp%PG_Symbol=pack_string(MGp%PG_Symbol)

         Case(4)
           MGp%PG_Symbol=MGp%BNS_symbol(4:) !Remove the type of lattice
           do i=4,len_trim(MGp%BNS_symbol)
             m=m+1
             if(MGp%BNS_symbol(i:i) == "a" .or. MGp%BNS_symbol(i:i) == "b"  &
           .or. MGp%BNS_symbol(i:i) == "c" .or. MGp%BNS_symbol(i:i) == "d"  &
           .or. MGp%BNS_symbol(i:i) == "e" .or. MGp%BNS_symbol(i:i) == "g"  &
           .or. MGp%BNS_symbol(i:i) == "n") MGp%PG_Symbol(m:m)="m"
             if(MGp%BNS_symbol(i:i) == "_") MGp%PG_Symbol(m:m+1)=" "
           end do
           MGp%PG_Symbol=pack_string(MGp%PG_Symbol//"1'")
      End Select

      if(len_trim(setting) == 0 .or. setting =='a,b,c;0,0,0') then
        MGp%standard_setting=.true.
      else
        MGp%standard_setting=.false.
      end if
      MGp%mcif=.false.     !true if mx,my,mz notation is used , false is u,v,w notation is used
      if(present(mcif)) MGp%mcif=mcif
      MGp%m_cell=.true.    !true if magnetic cell is used for symmetry operators
      MGp%m_constr=.false. !true if constraints have been provided
      MGp%trn_from_parent=" "
      MGp%trn_to_standard="a,b,c;0,0,0"
      MGp%trn_from_standard="a,b,c;0,0,0"
      !Info about Parent Crystallographic Space Group
      if(present(parent)) then
        !Parent should be of the form  Xnnn  num  trn_from_parent
        line=adjustl(parent)
        i=index(line," ")
        MGp%Parent_spg=parent(1:i-1)
        line=adjustl(line(i:))
        i=index(line," ")
        read(unit=line(1:i),fmt=*,iostat=ier) MGp%Parent_num
        if(ier /= 0) then
           MGp%Parent_num=0
           MGp%trn_from_parent=line(1:i)
        else
           line=adjustl(line(i:))
           i=index(line," ")
           MGp%trn_from_parent=line(1:i-1)
        end if
      else
        !Try to deduce the parent space group from the BNS/OG numbers
        line=MGp%BNS_number
        i=index(line,".")
        line=line(1:i-1)
        read(unit=line,fmt=*) MGp%Parent_num
        if(MGp%MagType < 4) then
          MGp%Parent_spg=MGp%BNS_symbol
          if(MGp%MagType == 2) MGp%Parent_spg=MGp%Parent_spg(1:len_trim(MGp%Parent_spg)-2)
          do i=1,len_trim(MGp%Parent_spg)
            if(MGp%Parent_spg(i:i) == "'") MGp%Parent_spg(i:i) = " "
          end do
          MGp%Parent_spg=Pack_String(MGp%Parent_spg)
        else
          line=MGp%OG_number
          i=index(line,".")
          line=line(1:i-1)
          read(unit=line,fmt=*) MGp%Parent_num
          MGp%Parent_spg=MGp%OG_symbol
          if(MGp%Parent_spg(3:3) == "2") then
             MGp%Parent_spg(2:4)=" "
          else
             MGp%Parent_spg(2:3)=" "
          end if
          do i=1,len_trim(MGp%Parent_spg)
            if(MGp%Parent_spg(i:i) == "'") MGp%Parent_spg(i:i) = " "
          end do
          MGp%Parent_spg=Pack_String(MGp%Parent_spg)
        end if
      end if
      MGp%standard_setting = .true.
      ! Crystal system
      Select Case (num)
        case(1:7)
          MGp%CrystalSys="Triclinic"
        case(8:98)
          MGp%CrystalSys="Monoclinic"
        case(99:660)
          MGp%CrystalSys="Orthorhombic"
        case(661:1230)
          MGp%CrystalSys="Tetragonal"
        case(1231:1338)
          MGp%CrystalSys="Trigonal"
        case(1339:1502)
          MGp%CrystalSys="Hexagonal"
        case(1503:1651)
          MGp%CrystalSys="Cubic"
        case default
          MGp%CrystalSys="Unknown"
      End Select
      if(MGp%MagType == 4) then
        MGp%SPG_lat=spacegroup_label_bns(num)(1:3)
      else
        MGp%SPG_lat=spacegroup_label_bns(num)(1:1)
      end if
      MGp%SPG_latsy=MGp%SPG_lat !provisional before knowing the crystal system

      MGp%Num_Lat=lattice_bns_vectors_count(num)-2         ! Number of lattice points in a cell
      if(allocated(MGp%Latt_trans)) deallocate(MGp%Latt_trans)
      allocate(MGp%Latt_trans(3,MGp%Num_Lat))
      MGp%Latt_trans=0.0
      centring=.false.
      if(MGp%Num_Lat > 1) centring=.true.
      m=1
      do j=4,lattice_bns_vectors_count(num)
         m=m+1
         MGp%Latt_trans(:,m)= real(lattice_bns_vectors(:,j,num))/real(lattice_bns_vectors_denom(j,num))
      end do

      j=1
      MGp%Multip=wyckoff_mult(j,num)
      if(allocated(MGp%SymopSymb)) deallocate(MGp%SymopSymb)  ! Alphanumeric Symbols for SYMM
      if(allocated(MGp%SymOp))     deallocate(MGp%SymOp)      ! Crystallographic symmetry operators
      if(allocated(MGp%MSymopSymb))deallocate(MGp%MSymopSymb) ! Alphanumeric Symbols for MSYMM
      if(allocated(MGp%MSymOp))    deallocate(MGp%MSymOp)     ! Magnetic symmetry operators
      allocate(MGp%SymOp(MGp%Multip))
      allocate(MGp%SymopSymb(MGp%Multip))
      allocate(MGp%MSymOp(MGp%Multip))
      allocate(MGp%MSymopSymb(MGp%Multip))

      m=0
      !write(*,"(3(a,i5))") "Shubnikov number: ",num,"Wyckoff position count: ",wyckoff_pos_count(j,num)," Multiplicity: ",MGp%Multip
      Do k=1,wyckoff_pos_count(j,num)
        idem=wyckoff_bns_fract_denom(k,j,num)
        MGp%SymOp(k)%tr=real(wyckoff_bns_fract(:,k,j,num))/real(idem)
        MGp%SymOp(k)%Rot = wyckoff_bns_xyz(:,:,k,j,num)
        MGp%MSymOp(k)%Rot = wyckoff_bns_mag(:,:,k,j,num)
        !inv_time=ops_bns_timeinv(k,num)  !Errors in the Database ... to be explored
        !MGp%MSymOp(k)%Phas=inv_time
        det=determ_a(MGp%SymOp(k)%Rot)
        if(det > 0.0) then
           if(equal_matrix(MGp%MSymOp(k)%Rot,MGp%SymOp(k)%Rot,3)) then
              MGp%MSymOp(k)%Phas=1.0
           else
              MGp%MSymOp(k)%Phas=-1.0
           end if
        else
           if(equal_matrix(MGp%MSymOp(k)%Rot,-MGp%SymOp(k)%Rot,3)) then
              MGp%MSymOp(k)%Phas=1.0
           else
              MGp%MSymOp(k)%Phas=-1.0
           end if
        end if
        if(MGp%mcif) then
           Call Get_Shubnikov_Operator_Symbol(MGp%SymOp(k)%Rot,MGp%MSymOp(k)%Rot,MGp%SymOp(k)%tr,ShOp_symb,MGp%mcif)
        else
           Call Get_Shubnikov_Operator_Symbol(MGp%SymOp(k)%Rot,MGp%MSymOp(k)%Rot,MGp%SymOp(k)%tr,ShOp_symb)
        end if
        !write(*,"(a)") trim(ShOp_symb)
        i=index(ShOp_symb,";")
        MGp%SymopSymb(k)=ShOp_symb(2:i-1)
        MGp%MSymopSymb(k)=ShOp_symb(i+1:len_trim(ShOp_symb)-1)
        if(MGp%MagType == 2) cycle
        if(equal_matrix(MGp%SymOp(k)%Rot,identity,3) .and. MGp%MSymOp(k)%Phas < 0.0) m=m+1 !counting anti-translations
      End Do

      if(centring) then
        n=wyckoff_pos_count(j,num)
        m=m*MGp%Num_Lat
        do L=2,MGp%Num_Lat
         do k=1,wyckoff_pos_count(j,num)
           MGp%SymOp(k+n)%Rot=MGp%SymOp(k)%Rot
           MGp%SymOp(k+n)%tr=Modulo_Lat(MGp%SymOp(k)%tr+MGp%Latt_trans(:,L))
           MGp%MSymOp(k+n)%Rot=MGp%MSymOp(k)%Rot
           MGp%MSymOp(k+n)%Phas=MGp%MSymOp(k)%Phas
           MGp%MSymopSymb(k+n)=MGp%MSymopSymb(k)
           Call Get_Shubnikov_Operator_Symbol(MGp%SymOp(k+n)%Rot,MGp%MSymOp(k+n)%Rot,MGp%SymOp(k+n)%tr,ShOp_symb)
           i=index(ShOp_symb,";")
           MGp%SymopSymb(k+n)=ShOp_symb(2:i-1)
         end do
         n=n+wyckoff_pos_count(j,num)
        end do
      end if

      MGp%Num_aLat=m       ! Number of anti-lattice points in a cell
      if(allocated(MGp%aLatt_trans)) deallocate(MGp%aLatt_trans)
      allocate(MGp%aLatt_trans(3,m))     ! Lattice anti-translations

      m=0
      if(MGp%MagType /= 2) then
        do k=1,MGp%multip
          if(equal_matrix(MGp%SymOp(k)%Rot,identity,3) .and. MGp%MSymOp(k)%Phas < 0) then
            m=m+1
            MGp%aLatt_trans(:,m) = MGp%SymOp(k)%tr
          end if
        end do
      end if
      MGp%Centred=0        ! Centric or Acentric [ =0 Centric(-1 no at origin),=1 Acentric,=2 Centric(-1 at origin)]
      MGp%Centre_coord=0.0 ! Fractional coordinates of the inversion centre
      do k=1,wyckoff_pos_count(j,num) !j=1 multiplicity of the general position
        if(equal_matrix(MGp%SymOp(k)%Rot,-identity,3) .and. MGp%MSymOp(k)%Phas > 0) then
          m=k
          MGp%Centred=max(MGp%Centred,1)
          if(sum(abs(MGp%SymOp(k)%tr)) < 0.001) then
            MGp%Centred=2
            exit
          end if
        end if
      end do
      MGp%NumOps=wyckoff_pos_count(j,num)
      MGp%Centre="Non-Centrosymmetric"    ! Alphanumeric information about the center of symmetry
      if(MGp%Centred == 1) then
        MGp%Centre="Centrosymmetric, -1 not @the origin "       ! Alphanumeric information about the center of symmetry
        MGp%Centre_coord=0.5*MGp%SymOp(m)%tr
      else if(MGp%Centred == 2) then
        MGp%Centre="Centrosymmetric, -1@the origin "       ! Alphanumeric information about the center of symmetry
        MGp%NumOps=MGp%NumOps/2
      end if
      !write(*,"(a)")    "  "//trim(MGp%Centre)
      !write(*,"(a,i4)") " Number of minimal S.O. (Numops): ",MGp%NumOps
      if(change_setting) then
        if(present(trn_to)) then
          call Setting_Change(setting,MGp,MSpg,trn_to)
        else
          call Setting_Change(setting,MGp,MSpg)
        end if
        if(Err_Form) then
          if(.not. present(keepd)) call deAllocate_DataBase()
          return
        end if
      else
        MSpg=MGp !everything is allocated in the assignement (Fortran 2003)
      end if
      if(.not. present(keepd)) call deAllocate_DataBase()
    End Subroutine Set_Magnetic_Space_Group

    !!----
    !!---- Subroutine Write_Cif_Powder_Profile(Filename,Pat,r_facts)
    !!----    character(len=*),                    intent(in) :: filename     !  In -> Name of File
    !!----    type(Diffraction_Pattern_Type),      intent(in) :: Pat
    !!----    real(kind=cp), dimension(4),optional,intent(in) :: r_facts      !R_patt,R_wpatt,R_exp, Chi2
    !!----
    !!----    Write a Cif Powder Profile file (converted from FullProf)
    !!----
    !!---- Update: January - 2020
    !!
    Subroutine Write_Cif_Powder_Profile(filename,Pat,r_facts)
       !---- Arguments ----!
       character(len=*),                    intent(in) :: filename
       type(Diffraction_Pattern_Type),      intent(in) :: Pat
       real(kind=cp), dimension(4),optional,intent(in) :: r_facts

       !---- Local Variables ----!
       logical             :: info
       character(len=132)  :: line
       character(len=30)   :: comm,date_time
       !character(len=1)    :: statut
       integer,save        :: iunit
       integer             :: i,j,n !,mirf,mult,iph,ivk,ix,icz,irc
       real(kind=cp)       :: an, R_patt,R_wpatt,R_exp, Chi2 !, phas, dspac, fobs2, fcal2
       !integer,      dimension(3):: hi
       !real(kind=cp),dimension(3):: hr

       !---- Inicialization of variables ----!
       info=.false.
       if(present(r_facts)) then
          R_patt = r_facts(1)
          R_wpatt= r_facts(2)
          R_exp  = r_facts(3)
          Chi2   = r_facts(4)
       end if

      !---- Is this file opened? ----!
       inquire(file=filename,opened=info)
       if (info) then
          inquire(file=filename,number=iunit)
          close(unit=iunit)
          open(iunit,file=filename,status="unknown",action="write",position="append")
       else
          iunit=61
          open(iunit,file=filename,status="replace",action="write",position="rewind")
       end if

       !---- Writing ----!
       !---- Head ----!
       write(unit=iunit,fmt='(a)')    " "
       write(unit=iunit,fmt='(a)')    "#==========================="
       write(unit=iunit,fmt='(a,i3)') "# Powder diffraction pattern "
       write(unit=iunit,fmt='(a)')    "#==========================="
       call Write_Date_Time(dtim=date_time)
       write(unit=iunit,fmt='(a)')    "#  "//trim(date_time)
       write(unit=iunit,fmt='(a)')    " "

       write(iunit,'(a)') "data_profile"
       j=index(date_time, "Time:")-1
       i=index(date_time, "Date:")+5
       write(unit=iunit,fmt="(a)")"_audit_creation_date "//date_time(i:j)
       write(unit=iunit,fmt="(a)")'_audit_creation_method  "CrysFML"'
       write(unit=iunit,fmt='(a)') " "
       write(unit=iunit,fmt='(a)') "_pd_block_id      ?"

       write(iunit,'(a)')"#==============================================================================             "
       write(unit=iunit,fmt='(a)')"# 9. INSTRUMENT CHARACTERIZATION                                                            "
       write(unit=iunit,fmt='(a)')"                                                                                            "
       write(unit=iunit,fmt='(a)')"_exptl_special_details                                                                      "
       write(unit=iunit,fmt='(a)')"; ?                                                                                         "
       write(unit=iunit,fmt='(a)')";                                                                                           "
       write(unit=iunit,fmt='(a)')"                                                                                            "
       write(unit=iunit,fmt='(a)')"# if regions of the data are excluded, the reason(s) are supplied here:                     "
       write(unit=iunit,fmt='(a)')"_pd_proc_info_excluded_regions                                                              "
       write(unit=iunit,fmt='(a)')"; ?                                                                                         "
       write(unit=iunit,fmt='(a)')";                                                                                           "
       write(unit=iunit,fmt='(a)')"                                                                                            "
       write(unit=iunit,fmt='(a)')"# The following item is used to identify the equipment used to record                       "
       write(unit=iunit,fmt='(a)')"# the powder pattern when the diffractogram was measured at a laboratory                    "
       write(unit=iunit,fmt='(a)')"# other than the authors' home institution, e.g. when neutron or synchrotron                "
       write(unit=iunit,fmt='(a)')"# radiation is used.                                                                        "
       write(unit=iunit,fmt='(a)')"                                                                                            "
       write(unit=iunit,fmt='(a)')"_pd_instr_location                                                                          "
       write(unit=iunit,fmt='(a)')"; ?                                                                                         "
       write(unit=iunit,fmt='(a)')";                                                                                           "
       write(unit=iunit,fmt='(a)')"_pd_calibration_special_details           # description of the method used                  "
       write(unit=iunit,fmt='(a)')"                                          # to calibrate the instrument                     "
       write(unit=iunit,fmt='(a)')"; ?                                                                                         "
       write(unit=iunit,fmt='(a)')";                                                                                           "
       write(unit=iunit,fmt='(a)')"                                                                                            "
       write(unit=iunit,fmt='(a)')"_diffrn_ambient_temperature    ?                                                            "
       write(unit=iunit,fmt='(a)')"_diffrn_source                 ?                                                            "
       write(unit=iunit,fmt='(a)')"_diffrn_source_target          ?                                                            "
       write(unit=iunit,fmt='(a)')"_diffrn_source_type            ?                                                            "
       write(unit=iunit,fmt='(a)')"_diffrn_measurement_device_type?                                                            "
       write(unit=iunit,fmt='(a)')"_diffrn_detector               ?                                                            "
       write(unit=iunit,fmt='(a)')"_diffrn_detector_type          ?  # make or model of detector                               "
       write(unit=iunit,fmt='(a)')"                                                                                            "
       write(unit=iunit,fmt='(a)')"_pd_meas_scan_method           ?  # options are 'step', 'cont',                             "
       write(unit=iunit,fmt='(a)')"                                  # 'tof', 'fixed' or                                       "
       write(unit=iunit,fmt='(a)')"                                  # 'disp' (= dispersive)                                   "
       write(unit=iunit,fmt='(a)')"_pd_meas_special_details                                                                    "
       write(unit=iunit,fmt='(a)')";  ?                                                                                        "
       write(unit=iunit,fmt='(a)')";                                                                                           "
       write(unit=iunit,fmt='(a)')"                                                                                            "
       write(unit=iunit,fmt='(a)')"# The following two items identify the program(s) used (if appropriate).                    "
       write(unit=iunit,fmt='(a)')"_computing_data_collection        ?                                                         "
       write(unit=iunit,fmt='(a)')"_computing_data_reduction         ?                                                         "
       write(unit=iunit,fmt='(a)')"                                                                                            "
       write(unit=iunit,fmt='(a)')"# Describe any processing performed on the data, prior to refinement.                       "
       write(unit=iunit,fmt='(a)')"# For example: a manual Lp correction or a precomputed absorption correction                "
       write(unit=iunit,fmt='(a)')"_pd_proc_info_data_reduction      ?                                                         "
       write(unit=iunit,fmt='(a)')"                                                                                            "
       write(unit=iunit,fmt='(a)')"# The following item is used for angular dispersive measurements only.                      "
       write(unit=iunit,fmt='(a)')"                                                                                            "
       write(unit=iunit,fmt='(a)')"_diffrn_radiation_monochromator   ?                                                         "
       write(unit=iunit,fmt='(a)')"                                                                                            "
       write(unit=iunit,fmt='(a)')"# The following items are used to define the size of the instrument.                        "
       write(unit=iunit,fmt='(a)')"# Not all distances are appropriate for all instrument types.                               "
       write(unit=iunit,fmt='(a)')"                                                                                            "
       write(unit=iunit,fmt='(a)')"_pd_instr_dist_src/mono           ?                                                         "
       write(unit=iunit,fmt='(a)')"_pd_instr_dist_mono/spec          ?                                                         "
       write(unit=iunit,fmt='(a)')"_pd_instr_dist_src/spec           ?                                                         "
       write(unit=iunit,fmt='(a)')"_pd_instr_dist_spec/anal          ?                                                         "
       write(unit=iunit,fmt='(a)')"_pd_instr_dist_anal/detc          ?                                                         "
       write(unit=iunit,fmt='(a)')"_pd_instr_dist_spec/detc          ?                                                         "
       write(unit=iunit,fmt='(a)')"                                                                                            "
       write(unit=iunit,fmt='(a)')"# 10. Specimen size and mounting information                                                "
       write(unit=iunit,fmt='(a)')"                                                                                            "
       write(unit=iunit,fmt='(a)')"# The next three fields give the specimen dimensions in mm.  The equatorial                 "
       write(unit=iunit,fmt='(a)')"# plane contains the incident and diffracted beam.                                          "
       write(unit=iunit,fmt='(a)')"                                                                                            "
       write(unit=iunit,fmt='(a)')"_pd_spec_size_axial               ?       # perpendicular to                                "
       write(unit=iunit,fmt='(a)')"                                          # equatorial plane                                "
       write(unit=iunit,fmt='(a)')"                                                                                            "
       write(unit=iunit,fmt='(a)')"_pd_spec_size_equat               ?       # parallel to                                     "
       write(unit=iunit,fmt='(a)')"                                          # scattering vector                               "
       write(unit=iunit,fmt='(a)')"                                          # in transmission                                 "
       write(unit=iunit,fmt='(a)')"                                                                                            "
       write(unit=iunit,fmt='(a)')"_pd_spec_size_thick               ?       # parallel to                                     "
       write(unit=iunit,fmt='(a)')"                                          # scattering vector                               "
       write(unit=iunit,fmt='(a)')"                                          # in reflection                                   "
       write(unit=iunit,fmt='(a)')"                                                                                            "
       write(unit=iunit,fmt='(a)')"_pd_spec_mounting                         # This field should be                            "
       write(unit=iunit,fmt='(a)')"                                          # used to give details of the                     "
       write(unit=iunit,fmt='(a)')"                                          # container.                                      "
       write(unit=iunit,fmt='(a)')"; ?                                                                                         "
       write(unit=iunit,fmt='(a)')";                                                                                           "
       write(unit=iunit,fmt='(a)')"                                                                                            "
       write(unit=iunit,fmt='(a)')"_pd_spec_mount_mode               ?       # options are 'reflection'                        "
       write(unit=iunit,fmt='(a)')"                                          # or 'transmission'                               "
       write(unit=iunit,fmt='(a)')"                                                                                            "
       write(unit=iunit,fmt='(a)')"_pd_spec_shape                    ?       # options are 'cylinder'                          "
       write(unit=iunit,fmt='(a)')"                                          # 'flat_sheet' or 'irregular'                     "
       write(unit=iunit,fmt='(a)')"                                                                                            "
       write(unit=iunit,fmt='(a)')"     "
       write(unit=iunit,fmt='(a)')"_diffrn_radiation_probe   "//trim(pat%diff_kind)
       if(Pat%scat_var=="2theta") then
          write(unit=iunit,fmt='(a,f12.6)') "_diffrn_radiation_wavelength ",Pat%conv(1)
       end if
       if(present(r_facts)) then
          write(unit=iunit,fmt='(a)')"     "
          write(unit=iunit,fmt='(a)') "#  The following profile R-factors are NOT CORRECTED for background"
          write(unit=iunit,fmt='(a)') "#  The sum is extended to all non-excluded points."
          write(unit=iunit,fmt='(a)') "#  These are the current CIF standard"
          write(unit=iunit,fmt='(a)') " "
          write(unit=iunit,fmt='(a,f12.4)') "_pd_proc_ls_prof_R_factor          ",R_patt
          write(unit=iunit,fmt='(a,f12.4)') "_pd_proc_ls_prof_wR_factor         ",R_wpatt
          write(unit=iunit,fmt='(a,f12.4)') "_pd_proc_ls_prof_wR_expected       ",R_exp
          write(unit=iunit,fmt='(a,f12.4)') "_pd_proc_ls_prof_chi2              ",chi2
       end if
       write(unit=iunit,fmt='(a)')"  "
       write(unit=iunit,fmt='(a)')"_pd_proc_ls_background_function   "
       write(unit=iunit,fmt='(a)')";   Background function description  "
      !write(unit=iunit,fmt='(a)')" Shifted Chebyshev function of 1st kind                                       "
      !write(unit=iunit,fmt='(a)')"      1:    61.5838     2:    15.7371     3:    40.4163     4:    9.42001     "
      !write(unit=iunit,fmt='(a)')"      5:    20.9238     6:    6.19575     7:    10.2828     8:    3.66734     "
      !write(unit=iunit,fmt='(a)')"      9:    2.70015    10:  -0.618124                                         "
       write(unit=iunit,fmt='(a)')";                                                                              "
       write(unit=iunit,fmt='(a)')"                                                                               "
       write(unit=iunit,fmt='(a)')"_exptl_absorpt_process_details                                                 "
       write(unit=iunit,fmt='(a)')";   Absorption/surface roughness correction description    "
       write(unit=iunit,fmt='(a)')" No correction is applied ?.                                                   "
       write(unit=iunit,fmt='(a)')";                                                                              "
      !write(unit=iunit,fmt='(a)')"_exptl_absorpt_correction_T_min        1.00000   "
      !write(unit=iunit,fmt='(a)')"_exptl_absorpt_correction_T_max        1.00000   "
       write(unit=iunit,fmt='(a)')"                                                                               "
       write(unit=iunit,fmt='(a)')"_pd_proc_ls_profile_function                                                   "
       write(unit=iunit,fmt='(a)')";   Profile function description                                              "
      !write(unit=iunit,fmt='(a)')" CW Profile function number 1 with   6 terms                                  "
      !write(unit=iunit,fmt='(a)')" Profile coefficients for Simpson's rule integration of Gaussian function     "
      !write(unit=iunit,fmt='(a)')" C.J. Howard (1982). J. Appl. Cryst.,15,615-620.                              "
      !write(unit=iunit,fmt='(a)')" Cooper & Sayer, J. Appl. Cryst., 8, 615-618 (1975).                          "
      !write(unit=iunit,fmt='(a)')" Thomas, J. Appl. Cryst., 10, 12-13(1977).                                    "
      !write(unit=iunit,fmt='(a)')" #1(U)    =  218.631 #2(V)    = -251.292 #3(W)    =  159.521                  "
      !write(unit=iunit,fmt='(a)')" #4(asym) =   8.7128 #5(F1)   =    0.000 #6(F2)   =    0.000                  "
      !write(unit=iunit,fmt='(a)')" Peak tails are ignored  where the intensity is below 0.0050 times the peak   "
      !write(unit=iunit,fmt='(a)')"   Aniso. broadening axis   0.0   0.0   1.0                                   "
       write(unit=iunit,fmt='(a)')";                                                                              "
       write(unit=iunit,fmt='(a)')"_pd_proc_ls_peak_cutoff 0.00500                                                "
     ! write(unit=iunit,fmt='(a)')"_pd_proc_info_datetime                 2002-12-21T19:04:06                    "
       write(unit=iunit,fmt='(a,a)')'_pd_calc_method  "   Rietveld Refinement" '
       write(unit=iunit,fmt='(a)')"                                                                               "
       write(unit=iunit,fmt='(a)')"#---- raw/calc data loop -----   "
       select case (trim(l_case(Pat%scat_var)))
          case ("2theta")   ! 2_Theta
               write(unit=iunit,fmt='(a,f14.6)')"_pd_meas_2theta_range_min " , Pat%xmin
               write(unit=iunit,fmt='(a,f14.6)')"_pd_meas_2theta_range_max " , Pat%xmax
               write(unit=iunit,fmt='(a,f14.6)')"_pd_meas_2theta_range_inc " , Pat%step
           case ("tof")   ! T.O.F.
             write(unit=iunit,fmt='(a)') "_pd_proc_d_spacing "
       end select


       !---- Profile ----!
       write(unit=iunit,fmt='(a)') " "

       write(unit=iunit,fmt='(a)') "loop_"
       write(unit=iunit,fmt='(a)') "_pd_proc_point_id"
       select case (trim(l_case(Pat%scat_var)))
          case ("2theta")   !
             write(unit=iunit,fmt='(a)') "_pd_proc_2theta_corrected   "
          case ("tof")   ! T.O.F.
             write(unit=iunit,fmt='(a)') "_pd_proc_d_spacing "
          case ("energy")   ! Energy
             write(unit=iunit,fmt='(a)') "_pd_proc_energy_incident  "
       end select
       write(iunit,'(a)') "_pd_proc_intensity_total"
       write(iunit,'(a)') "_pd_calc_intensity_total"
       write(iunit,'(a)') "_pd_proc_intensity_bkg_calc"
       write(iunit,'(a)') " "
       n=0
       do_poi: do i=1,Pat%npts
          if(Pat%istat(i) == 0) cycle
          an=Pat%x(i)
          line=" "
          write(line(1:6),'(i6)') i
          n=n+1
          select case (trim(l_case(Pat%scat_var)))
             case ("2theta")   ! 2_Theta
                write(line(10:),'(f8.4,5x,a)') an-Pat%zerop,'.    .'
             case ("tof")   ! T.O.F.
                write(line(10:),'(f8.4,5x,a)') (an-Pat%zerop)/pat%conv(1),'.    .'  !dtt1
             case ("energy")   ! Energy
                write(line(10:),'(a,f15.4,2x,a)') '. ',1000.0*(an-Pat%zerop),'.'
          end select
          call setnum_std(Pat%y(i),sqrt(Pat%sigma(i)),comm)
          write(line(21:),'(a,2f18.4)') trim(comm)//" ", Pat%ycalc(i),Pat%bgr(i)
          write(iunit,'(a)') line
       end do do_poi

       write(iunit,'(a,i7)')  "_pd_proc_number_of_points",n
       write(iunit,'(a)') " "

      !Writing the reflections of the current pattern
      !write(iunit,'(a)') "loop_                          "
      !write(iunit,'(a)') "      _refln_index_h           "
      !write(iunit,'(a)') "      _refln_index_k           "
      !write(iunit,'(a)') "      _refln_index_l           "
      !write(iunit,'(a)') "      _refln_mult              "
      !write(iunit,'(a)') "      _pd_refln_phase_id       "
      !write(iunit,'(a)') "      _refln_observed_status   "
      !write(iunit,'(a)') "      _refln_F_squared_meas    "
      !write(iunit,'(a)') "      _refln_F_squared_calc    "
      !write(iunit,'(a)') "      _refln_phase_calc        "
      !write(iunit,'(a)') "      _refln_d_spacing         "
      !
      !icz=sum(Num_refl(:,n_pat))
      !statut="o"
      !do_ref: do ix=1,icz
      !   iph=irefs(ix,n_pat)/256/256/256/8
      !   if(iph < 0 ) iph = 8 - iph
      !   mirf=abs(irefs(ix,n_pat))
      !   irc=MOD(mirf/256/256/256,8)
      !   if(irc /= 1) cycle
      !   hi(3)=MOD(mirf,256)-128
      !   hi(2)=MOD(mirf/256,256)-128
      !   hi(1)=MOD(mirf/256/256,256)-128
      !   dspac=refs(ix,4,n_pat)
      !   an=refs(ix,2,n_pat)
      !   DO j=1,nexcrg(n_pat)
      !      IF((an >= alow(j,n_pat) .AND. an <= ahigh(j,n_pat)) &
      !         .or. an < thmin(n_pat)-glb(1,n_pat) .or.  an > thmax(n_pat)-glb(1,n_pat) ) cycle do_ref
      !   END DO
      !   fcal2=ff(ix,n_pat)
      !   fobs2=fobs(ix,n_pat)
      !   phas=phasen(ix,n_pat)*57.29577951308
      !   mult=mlt(ix,n_pat)
      !   ivk=ihkl(ix,n_pat)
      !   IF(nvk(iph) /= 0 .AND. ivk /= 0) THEN
      !     hr(:)=REAL(hi(:))+pvk(iph,ivk,:)
      !     write(unit=iunit,fmt="(3f7.2,2i3,tr2,a,tr2,2f14.4,f8.2,f10.5)") &
      !                           hr,mult,iph,statut,fobs2,fcal2,phas,dspac
      !   else
      !     write(unit=iunit,fmt="(3i4,2i3,tr2,a,tr2,2f14.4,f8.2,f10.5)") &
      !                           hi,mult,iph,statut,fobs2,fcal2,phas,dspac
      !   END IF
      !
      !end do do_ref


       write(iunit,'(a)') " "
       write(iunit,'(a)') "# The following lines are used to test the character set of files sent by     "
       write(iunit,'(a)') "# network email or other means. They are not part of the CIF data set.        "
       write(iunit,'(a)') "# abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789              "
       write(iunit,'(a)') "# !@#$%^&*()_+{}:""~<>?|\-=[];'`,./ "


       return
    End Subroutine Write_Cif_Powder_Profile

    !!----
    !!---- Subroutine Write_Cif_Template(filename,type_data,code,cell,SpG,A)
    !!----    character(len=*),        intent(in) :: filename   !  In -> Filename
    !!----    integer,                 intent(in) :: type_data  !  In -> 0: Single Crystal, 1: Powder Data, 2:Only structural data
    !!----    character(len=*),        intent(in) :: code       !  In -> Code name for the data set
    !!----    Type (Crystal_Cell_Type),intent(in) :: Cell       ! Cell type to be output
    !!----    Type (Space_Group_Type), intent(in) :: SpG        ! Space group type to be output
    !!----    Type (Atom_List_Type),   intent(in) :: A          ! Atom list type to be output
    !!----
    !!----    Write a Cif File
    !!----
    !!---- Updated: February - 2005, January 2015
    !!
    Subroutine Write_Cif_Template(filename,type_data,code,cell,SpG,A)
       !---- Arguments ----!
       character(len=*),        intent(in) :: filename
       integer,                 intent(in) :: type_data
       character(len=*),        intent(in) :: code
       Type (Crystal_Cell_Type),intent(in) :: Cell
       Type (Space_Group_Type), intent(in) :: SpG
       Type (Atom_List_Type),   intent(in) :: A

       !---- Local Variables ----!
       logical                           :: info,aniso
       character(len=132)                :: line
       character(len=1), parameter       :: qmark='?'
       character(len=30)                 :: comm,adptyp
       character(len=30),dimension(6)    :: text
       real(kind=cp)                     :: u,su, ocf
       real(kind=cp), dimension(6)       :: Ua,sua,aux
       real(kind=cp), dimension(A%natoms):: occup,soccup
       integer                           :: iunit,i, j

       !---- Initialization of variables ----!
       info=.false.
       iunit=0

       !---- Is this file opened? ----!
       inquire(file=filename,opened=info)
       if (info) then
          inquire(file=filename,number=iunit)
          close(unit=iunit)
       end if

       !---- Writing ----!
       if (iunit==0) iunit=61
       open(unit=iunit,file=filename,status="unknown",action="write")
       rewind(unit=iunit)

       !---- Head Information ----!
       if(type_data == 0) then
           write(unit=iunit,fmt="(a)") "##############################################################################"
           write(unit=iunit,fmt="(a)") "###    CIF submission form for molecular structure report (Acta Cryst. C)  ###"
           write(unit=iunit,fmt="(a)") "##############################################################################"
           write(unit=iunit,fmt="(a)") " "
           write(unit=iunit,fmt="(a)") "#============================================================================="
           write(unit=iunit,fmt="(a)") "data_global"
           write(unit=iunit,fmt="(a)") "#============================================================================="
           write(unit=iunit,fmt="(a)") " "
          else if(type_data > 1) then
           write(unit=iunit,fmt="(a)") "##################################################################"
           write(unit=iunit,fmt="(a)") "###    CIF file from CrysFML, contains only structural data    ###"
           write(unit=iunit,fmt="(a)") "##################################################################"

       end if

       !---- Processing Summary ----!
       if(type_data < 2) then
         write(unit=iunit,fmt="(a)") "# PROCESSING SUMMARY (IUCr Office Use Only)"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_journal_data_validation_number      ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_journal_date_recd_electronic        ?"
         write(unit=iunit,fmt="(a)") "_journal_date_to_coeditor            ?"
         write(unit=iunit,fmt="(a)") "_journal_date_from_coeditor          ?"
         write(unit=iunit,fmt="(a)") "_journal_date_accepted               ?"
         write(unit=iunit,fmt="(a)") "_journal_date_printers_first         ?"
         write(unit=iunit,fmt="(a)") "_journal_date_printers_final         ?"
         write(unit=iunit,fmt="(a)") "_journal_date_proofs_out             ?"
         write(unit=iunit,fmt="(a)") "_journal_date_proofs_in              ?"
         write(unit=iunit,fmt="(a)") "_journal_coeditor_name               ?"
         write(unit=iunit,fmt="(a)") "_journal_coeditor_code               ?"
         write(unit=iunit,fmt="(a)") "_journal_coeditor_notes"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_journal_techeditor_code             ?"
         write(unit=iunit,fmt="(a)") "_journal_techeditor_notes"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_journal_coden_ASTM                  ?"
         write(unit=iunit,fmt="(a)") "_journal_name_full                   ?"
         write(unit=iunit,fmt="(a)") "_journal_year                        ?"
         write(unit=iunit,fmt="(a)") "_journal_volume                      ?"
         write(unit=iunit,fmt="(a)") "_journal_issue                       ?"
         write(unit=iunit,fmt="(a)") "_journal_page_first                  ?"
         write(unit=iunit,fmt="(a)") "_journal_page_last                   ?"
         write(unit=iunit,fmt="(a)") "_journal_paper_category              ?"
         write(unit=iunit,fmt="(a)") "_journal_suppl_publ_number           ?"
         write(unit=iunit,fmt="(a)") "_journal_suppl_publ_pages            ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !---- Submission details ----!
         write(unit=iunit,fmt="(a)") "# 1. SUBMISSION DETAILS"
         write(unit=iunit,fmt="(a)") " "

         write(unit=iunit,fmt="(a)") "_publ_contact_author_name            ?   # Name of author for correspondence"
         write(unit=iunit,fmt="(a)") "_publ_contact_author_address             # Address of author for correspondence"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_contact_author_email           ?"
         write(unit=iunit,fmt="(a)") "_publ_contact_author_fax             ?"
         write(unit=iunit,fmt="(a)") "_publ_contact_author_phone           ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_publ_contact_letter"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_publ_requested_journal              ?"
         write(unit=iunit,fmt="(a)") "_publ_requested_coeditor_name        ?"
         write(unit=iunit,fmt="(a)") "_publ_requested_category             ?   # Acta C: one of CI/CM/CO/FI/FM/FO"

         write(unit=iunit,fmt="(a)") "#=============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !---- Title  and Author List ----!
         write(unit=iunit,fmt="(a)") "# 3. TITLE AND AUTHOR LIST"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_publ_section_title"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_title_footnote"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") ";"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "# The loop structure below should contain the names and addresses of all "
         write(unit=iunit,fmt="(a)") "# authors, in the required order of publication. Repeat as necessary."

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "loop_"
         write(unit=iunit,fmt="(a)") "    _publ_author_name"
         write(unit=iunit,fmt="(a)") "    _publ_author_footnote"
         write(unit=iunit,fmt="(a)") "    _publ_author_address"
         write(unit=iunit,fmt="(a)") "?                                   #<--'Last name, first name' "
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !---- Text ----!
         write(unit=iunit,fmt="(a)") "# 4. TEXT"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_publ_section_synopsis"
         write(unit=iunit,fmt="(a)") ";  ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_abstract"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";          "
         write(unit=iunit,fmt="(a)") "_publ_section_comment"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_exptl_prep      # Details of the preparation of the sample(s)"
         write(unit=iunit,fmt="(a)") "                              # should be given here. "
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_exptl_refinement"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_references"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_figure_captions"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_acknowledgements"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !---- Identifier ----!
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") "# If more than one structure is reported, the remaining sections should be "
         write(unit=iunit,fmt="(a)") "# completed per structure. For each data set, replace the '?' in the"
         write(unit=iunit,fmt="(a)") "# data_? line below by a unique identifier."
       end if !type_data < 2

       write(unit=iunit,fmt="(a)") " "
       if (len_trim(code) == 0) then
          write(unit=iunit,fmt="(a)") "data_?"
       else
          write(unit=iunit,fmt="(a)") "data_"//code(1:len_trim(code))
       end if
       write(unit=iunit,fmt="(a)") " "
       if(type_data < 2) then
         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !---- Chemical Data ----!
         write(unit=iunit,fmt="(a)") "# 5. CHEMICAL DATA"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_chemical_name_systematic"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_chemical_name_common             ?"
         write(unit=iunit,fmt="(a)") "_chemical_formula_moiety          ?"
         write(unit=iunit,fmt="(a)") "_chemical_formula_structural      ?"
         write(unit=iunit,fmt="(a)") "_chemical_formula_analytical      ?"
         write(unit=iunit,fmt="(a)") "_chemical_formula_iupac           ?"
         write(unit=iunit,fmt="(a)") "_chemical_formula_sum             ?"
         write(unit=iunit,fmt="(a)") "_chemical_formula_weight          ?"
         write(unit=iunit,fmt="(a)") "_chemical_melting_point           ?"
         write(unit=iunit,fmt="(a)") "_chemical_compound_source         ?       # for minerals and "
         write(unit=iunit,fmt="(a)") "                                          # natural products"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "loop_"
         write(unit=iunit,fmt="(a)") "    _atom_type_symbol               "
         write(unit=iunit,fmt="(a)") "    _atom_type_description          "
         write(unit=iunit,fmt="(a)") "    _atom_type_scat_dispersion_real "
         write(unit=iunit,fmt="(a)") "    _atom_type_scat_dispersion_imag "
         write(unit=iunit,fmt="(a)") "    _atom_type_scat_source          "
         write(unit=iunit,fmt="(a)") "    _atom_type_scat_length_neutron       # include if applicable"
         write(unit=iunit,fmt="(a)") "    ?    ?    ?    ?    ?      ?    "

       end if !type_data < 2
       write(unit=iunit,fmt="(a)") " "
       write(unit=iunit,fmt="(a)") "#============================================================================="
       write(unit=iunit,fmt="(a)") " "
       !---- Crystal Data ----!
       select case (type_data)
          case (0,2) ! Single Crystal or structural data only
             write(unit=iunit,fmt="(a)") "# 6. CRYSTAL DATA"
          case (1) ! Powder Data + Crystal Data
             write(unit=iunit,fmt="(a)") "# 6. POWDER SPECIMEN AND CRYSTAL DATA"
       end select

       write(unit=iunit,fmt="(a)") " "
       write(unit=iunit,fmt="(a)") "_symmetry_cell_setting               ?"
       line=SpG%SPG_Symb
       write(unit=iunit,fmt="(a)") "_symmetry_space_group_name_H-M       '"//trim(line)//"'"
       line=SpG%Hall
       write(unit=iunit,fmt="(a)") "_symmetry_space_group_name_Hall      '"//trim(line)//"'"

       write(unit=iunit,fmt="(a)") " "
       write(unit=iunit,fmt="(a)") "loop_"
       write(unit=iunit,fmt="(a)") "    _symmetry_equiv_pos_as_xyz"
       do i=1,SpG%multip
          line="'"//trim(SpG%SymopSymb(i))//"'"
          write(iunit,'(a)') trim(line)
       end do

       write(unit=iunit,fmt="(a)") " "
       do i=1,3
          call setnum_std(Cell%cell(i),Cell%cell_std(i),text(i))
       end do
       do i=1,3
          call setnum_std(Cell%ang(i),Cell%ang_std(i),text(i+3))
       end do
       write(iunit,'(a)') "_cell_length_a                       "//trim(adjustl(text(1)))
       write(iunit,'(a)') "_cell_length_b                       "//trim(adjustl(text(2)))
       write(iunit,'(a)') "_cell_length_c                       "//trim(adjustl(text(3)))
       write(iunit,'(a)') "_cell_angle_alpha                    "//trim(adjustl(text(4)))
       write(iunit,'(a)') "_cell_angle_beta                     "//trim(adjustl(text(5)))
       write(iunit,'(a)') "_cell_angle_gamma                    "//trim(adjustl(text(6)))

       write(unit=iunit,fmt="(a,f14.4)") "_cell_volume                   ",Cell%CellVol
       if(type_data < 2) then
         write(unit=iunit,fmt="(a)") "_cell_formula_units_Z                ?"
         write(unit=iunit,fmt="(a)") "_cell_measurement_temperature        ?"
         write(unit=iunit,fmt="(a)") "_cell_special_details"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
       end if

       select case (type_data)
          case (0) ! Single Crystal
             write(unit=iunit,fmt="(a)") "_cell_measurement_reflns_used        ?"
             write(unit=iunit,fmt="(a)") "_cell_measurement_theta_min          ?"
             write(unit=iunit,fmt="(a)") "_cell_measurement_theta_max          ?"

             write(unit=iunit,fmt="(a)") " "
             write(unit=iunit,fmt="(a)") "_exptl_crystal_description           ?"
             write(unit=iunit,fmt="(a)") "_exptl_crystal_colour                ?"
             write(unit=iunit,fmt="(a)") "_exptl_crystal_size_max              ?"
             write(unit=iunit,fmt="(a)") "_exptl_crystal_size_mid              ?"
             write(unit=iunit,fmt="(a)") "_exptl_crystal_size_min              ?"
             write(unit=iunit,fmt="(a)") "_exptl_crystal_size_rad              ?"
             write(unit=iunit,fmt="(a)") "_exptl_crystal_density_diffrn        ?"
             write(unit=iunit,fmt="(a)") "_exptl_crystal_density_meas          ?"
             write(unit=iunit,fmt="(a)") "_exptl_crystal_density_method        ?"
             write(unit=iunit,fmt="(a)") "_exptl_crystal_F_000                 ?"

          case (1) ! Powder Data
             write(unit=iunit,fmt="(a)") "# The next three fields give the specimen dimensions in mm.  The equatorial"
             write(unit=iunit,fmt="(a)") "# plane contains the incident and diffracted beam."

             write(unit=iunit,fmt="(a)") " "
             write(unit=iunit,fmt="(a)") "_pd_spec_size_axial               ?       # perpendicular to "
             write(unit=iunit,fmt="(a)") "                                          # equatorial plane"

             write(unit=iunit,fmt="(a)") "_pd_spec_size_equat               ?       # parallel to "
             write(unit=iunit,fmt="(a)") "                                          # scattering vector"
             write(unit=iunit,fmt="(a)") "                                          # in transmission"
             write(unit=iunit,fmt="(a)") "_pd_spec_size_thick               ?       # parallel to "
             write(unit=iunit,fmt="(a)") "                                          # scattering vector"
             write(unit=iunit,fmt="(a)") "                                          # in reflection"

             write(unit=iunit,fmt="(a)") " "
             write(unit=iunit,fmt="(a)") "# The next five fields are character fields that describe the specimen."

             write(unit=iunit,fmt="(a)") " "
             write(unit=iunit,fmt="(a)") "_pd_spec_mounting                         # This field should be"
             write(unit=iunit,fmt="(a)") "                                          # used to give details of the "
             write(unit=iunit,fmt="(a)") "                                          # container."
             write(unit=iunit,fmt="(a)") "; ?"
             write(unit=iunit,fmt="(a)") ";"
             write(unit=iunit,fmt="(a)") "_pd_spec_mount_mode               ?       # options are 'reflection'"
             write(unit=iunit,fmt="(a)") "                                          # or 'transmission'"
             write(unit=iunit,fmt="(a)") "_pd_spec_shape                    ?       # options are 'cylinder' "
             write(unit=iunit,fmt="(a)") "                                          # 'flat_sheet' or 'irregular'"
             write(unit=iunit,fmt="(a)") "_pd_char_particle_morphology      ?"
             write(unit=iunit,fmt="(a)") "_pd_char_colour                   ?       # use ICDD colour descriptions"

             write(unit=iunit,fmt="(a)") " "
             write(unit=iunit,fmt="(a)") "# The following three fields describe the preparation of the specimen."
             write(unit=iunit,fmt="(a)") "# The cooling rate is in K/min.  The pressure at which the sample was "
             write(unit=iunit,fmt="(a)") "# prepared is in kPa.  The temperature of preparation is in K.        "

             write(unit=iunit,fmt="(a)") " "
             write(unit=iunit,fmt="(a)") "_pd_prep_cool_rate                ?"
             write(unit=iunit,fmt="(a)") "_pd_prep_pressure                 ?"
             write(unit=iunit,fmt="(a)") "_pd_prep_temperature              ?"
       end select
       if(type_data < 2) then
         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "# The next four fields are normally only needed for transmission experiments."
         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_exptl_absorpt_coefficient_mu        ?"
         write(unit=iunit,fmt="(a)") "_exptl_absorpt_correction_type       ?"
         write(unit=iunit,fmt="(a)") "_exptl_absorpt_process_details       ?"
         write(unit=iunit,fmt="(a)") "_exptl_absorpt_correction_T_min      ?"
         write(unit=iunit,fmt="(a)") "_exptl_absorpt_correction_T_max      ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !---- Experimental Data ----!
         write(unit=iunit,fmt="(a)") "# 7. EXPERIMENTAL DATA"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_exptl_special_details"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"

         if (type_data == 1) then
            write(unit=iunit,fmt="(a)") " "
            write(unit=iunit,fmt="(a)") "# The following item is used to identify the equipment used to record "
            write(unit=iunit,fmt="(a)") "# the powder pattern when the diffractogram was measured at a laboratory "
            write(unit=iunit,fmt="(a)") "# other than the authors' home institution, e.g. when neutron or synchrotron"
            write(unit=iunit,fmt="(a)") "# radiation is used."

            write(unit=iunit,fmt="(a)") " "
            write(unit=iunit,fmt="(a)") "_pd_instr_location"
            write(unit=iunit,fmt="(a)") "; ?"
            write(unit=iunit,fmt="(a)") ";"
            write(unit=iunit,fmt="(a)") "_pd_calibration_special_details           # description of the method used"
            write(unit=iunit,fmt="(a)") "                                          # to calibrate the instrument"
            write(unit=iunit,fmt="(a)") "; ?"
            write(unit=iunit,fmt="(a)") ";"
         end if

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_diffrn_ambient_temperature          ?"
         write(unit=iunit,fmt="(a)") "_diffrn_radiation_type               ?"
         write(unit=iunit,fmt="(a)") "_diffrn_radiation_wavelength         ?"
         write(unit=iunit,fmt="(a)") "_diffrn_radiation_source             ?"
         write(unit=iunit,fmt="(a)") "_diffrn_source                       ?"
         write(unit=iunit,fmt="(a)") "_diffrn_source_target                ?"
         write(unit=iunit,fmt="(a)") "_diffrn_source_type                  ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_diffrn_radiation_monochromator      ?"
         write(unit=iunit,fmt="(a)") "_diffrn_measurement_device_type      ?"
         write(unit=iunit,fmt="(a)") "_diffrn_measurement_method           ?"
         write(unit=iunit,fmt="(a)") "_diffrn_detector_area_resol_mean     ?   # Not in version 2.0.1"
         write(unit=iunit,fmt="(a)") "_diffrn_detector                     ?"
         write(unit=iunit,fmt="(a)") "_diffrn_detector_type                ?   # make or model of detector"
         if (type_data == 1) then
            write(unit=iunit,fmt="(a)") "_pd_meas_scan_method                 ?   # options are 'step', 'cont',"
            write(unit=iunit,fmt="(a)") "                                         # 'tof', 'fixed' or"
            write(unit=iunit,fmt="(a)") "                                         # 'disp' (= dispersive)"
            write(unit=iunit,fmt="(a)") "_pd_meas_special_details"
            write(unit=iunit,fmt="(a)") ";  ?"
            write(unit=iunit,fmt="(a)") ";"
         end if

         select case (type_data)
            case (0)
               write(unit=iunit,fmt="(a)") " "
               write(unit=iunit,fmt="(a)") "_diffrn_reflns_number                ?"
               write(unit=iunit,fmt="(a)") "_diffrn_reflns_av_R_equivalents      ?"
               write(unit=iunit,fmt="(a)") "_diffrn_reflns_av_sigmaI/netI        ?"
               write(unit=iunit,fmt="(a)") "_diffrn_reflns_theta_min             ?"
               write(unit=iunit,fmt="(a)") "_diffrn_reflns_theta_max             ?"
               write(unit=iunit,fmt="(a)") "_diffrn_reflns_theta_full            ?"
               write(unit=iunit,fmt="(a)") "_diffrn_measured_fraction_theta_max  ?"
               write(unit=iunit,fmt="(a)") "_diffrn_measured_fraction_theta_full ?"
               write(unit=iunit,fmt="(a)") "_diffrn_reflns_limit_h_min           ?"
               write(unit=iunit,fmt="(a)") "_diffrn_reflns_limit_h_max           ?"
               write(unit=iunit,fmt="(a)") "_diffrn_reflns_limit_k_min           ?"
               write(unit=iunit,fmt="(a)") "_diffrn_reflns_limit_k_max           ?"
               write(unit=iunit,fmt="(a)") "_diffrn_reflns_limit_l_min           ?"
               write(unit=iunit,fmt="(a)") "_diffrn_reflns_limit_l_max           ?"
               write(unit=iunit,fmt="(a)") "_diffrn_reflns_reduction_process     ?"

               write(unit=iunit,fmt="(a)") " "
               write(unit=iunit,fmt="(a)") "_diffrn_standards_number             ?"
               write(unit=iunit,fmt="(a)") "_diffrn_standards_interval_count     ?"
               write(unit=iunit,fmt="(a)") "_diffrn_standards_interval_time      ?"
               write(unit=iunit,fmt="(a)") "_diffrn_standards_decay_%            ?"
               write(unit=iunit,fmt="(a)") "loop_"
               write(unit=iunit,fmt="(a)") "    _diffrn_standard_refln_index_h"
               write(unit=iunit,fmt="(a)") "    _diffrn_standard_refln_index_k"
               write(unit=iunit,fmt="(a)") "    _diffrn_standard_refln_index_l"
               write(unit=iunit,fmt="(a)") "?   ?   ?"

            case (1)
               write(unit=iunit,fmt="(a)") " "
               write(unit=iunit,fmt="(a)") "#  The following four items give details of the measured (not processed)"
               write(unit=iunit,fmt="(a)") "#  powder pattern.  Angles are in degrees."

               write(unit=iunit,fmt="(a)") " "
               write(unit=iunit,fmt="(a)") "_pd_meas_number_of_points         ?"
               write(unit=iunit,fmt="(a)") "_pd_meas_2theta_range_min         ?"
               write(unit=iunit,fmt="(a)") "_pd_meas_2theta_range_max         ?"
               write(unit=iunit,fmt="(a)") "_pd_meas_2theta_range_inc         ?"

               write(unit=iunit,fmt="(a)") " "
               write(unit=iunit,fmt="(a)") "# The following three items are used for time-of-flight measurements only."

               write(unit=iunit,fmt="(a)") " "
               write(unit=iunit,fmt="(a)") "_pd_instr_dist_src/spec           ?"
               write(unit=iunit,fmt="(a)") "_pd_instr_dist_spec/detc          ?"
               write(unit=iunit,fmt="(a)") "_pd_meas_2theta_fixed             ?"

         end select

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !---- Refinement Data ----!
         write(unit=iunit,fmt="(a)") "# 8. REFINEMENT DATA"

         write(unit=iunit,fmt="(a)") " "

         write(unit=iunit,fmt="(a)") "_refine_special_details"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"

         if (type_data == 1) then
            write(unit=iunit,fmt="(a)") " "
            write(unit=iunit,fmt="(a)") "# Use the next field to give any special details about the fitting of the"
            write(unit=iunit,fmt="(a)") "# powder pattern."

            write(unit=iunit,fmt="(a)") " "
            write(unit=iunit,fmt="(a)") "_pd_proc_ls_special_details"
            write(unit=iunit,fmt="(a)") "; ?"
            write(unit=iunit,fmt="(a)") ";"

            write(unit=iunit,fmt="(a)") " "
            write(unit=iunit,fmt="(a)") "# The next three items are given as text."
            write(unit=iunit,fmt="(a)") " "

            write(unit=iunit,fmt="(a)") "_pd_proc_ls_profile_function      ?"
            write(unit=iunit,fmt="(a)") "_pd_proc_ls_background_function   ?"
            write(unit=iunit,fmt="(a)") "_pd_proc_ls_pref_orient_corr"
            write(unit=iunit,fmt="(a)") "; ?"
            write(unit=iunit,fmt="(a)") ";"
         end if

         select case (type_data)
            case (0)
               write(unit=iunit,fmt="(a)") " "
               write(unit=iunit,fmt="(a)") "_reflns_number_total                 ?"
               write(unit=iunit,fmt="(a)") "_reflns_number_gt                    ?"
               write(unit=iunit,fmt="(a)") "_reflns_threshold_expression         ?"

            case (1)
               write(unit=iunit,fmt="(a)") " "
               write(unit=iunit,fmt="(a)") "_pd_proc_ls_prof_R_factor         ?"
               write(unit=iunit,fmt="(a)") "_pd_proc_ls_prof_wR_factor        ?"
               write(unit=iunit,fmt="(a)") "_pd_proc_ls_prof_wR_expected      ?"

              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "# The following four items apply to angular dispersive measurements."
              write(unit=iunit,fmt="(a)") "# 2theta minimum, maximum and increment (in degrees) are for the "
              write(unit=iunit,fmt="(a)") "# intensities used in the refinement."

              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "_pd_proc_2theta_range_min         ?"
              write(unit=iunit,fmt="(a)") "_pd_proc_2theta_range_max         ?"
              write(unit=iunit,fmt="(a)") "_pd_proc_2theta_range_inc         ?"
              write(unit=iunit,fmt="(a)") "_pd_proc_wavelength               ?"

              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "_pd_block_diffractogram_id        ?  # The id used for the block containing"
              write(unit=iunit,fmt="(a)") "                                     # the powder pattern profile (section 11)."

              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "# Give appropriate details in the next two text fields."
              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "_pd_proc_info_excluded_regions    ?"
              write(unit=iunit,fmt="(a)") "_pd_proc_info_data_reduction      ?"

         end select

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_refine_ls_structure_factor_coef     ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_matrix_type               ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_R_I_factor                ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_R_Fsqd_factor             ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_R_factor_all              ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_R_factor_gt               ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_wR_factor_all             ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_wR_factor_ref             ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_goodness_of_fit_all       ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_goodness_of_fit_ref       ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_restrained_S_all          ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_restrained_S_obs          ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_number_reflns             ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_number_parameters         ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_number_restraints         ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_number_constraints        ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_hydrogen_treatment        ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_weighting_scheme          ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_weighting_details         ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_shift/su_max              ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_shift/su_mean             ?"
         write(unit=iunit,fmt="(a)") "_refine_diff_density_max             ?"
         write(unit=iunit,fmt="(a)") "_refine_diff_density_min             ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_extinction_method         ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_extinction_coef           ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_abs_structure_details     ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_abs_structure_Flack       ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_abs_structure_Rogers      ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "# The following items are used to identify the programs used."
         write(unit=iunit,fmt="(a)") " "

         write(unit=iunit,fmt="(a)") "_computing_data_collection           ?"
         write(unit=iunit,fmt="(a)") "_computing_cell_refinement           ?"
         write(unit=iunit,fmt="(a)") "_computing_data_reduction            ?"
         write(unit=iunit,fmt="(a)") "_computing_structure_solution        ?"
         write(unit=iunit,fmt="(a)") "_computing_structure_refinement      ?"
         write(unit=iunit,fmt="(a)") "_computing_molecular_graphics        ?"
         write(unit=iunit,fmt="(a)") "_computing_publication_material      ?"

       End if  !(type_data < 2) then
       write(unit=iunit,fmt="(a)") " "
       write(unit=iunit,fmt="(a)") "#============================================================================="
       write(unit=iunit,fmt="(a)") " "
       !---- Atomic Coordinates and Displacement Parameters ----!
       write(unit=iunit,fmt="(a)") "# 9. ATOMIC COORDINATES AND DISPLACEMENT PARAMETERS"

       write(unit=iunit,fmt="(a)") " "

       write(unit=iunit,fmt="(a)") "loop_"
       write(unit=iunit,fmt='(a)') "    _atom_site_label"
       write(unit=iunit,fmt='(a)') "    _atom_site_type_symbol"
       write(unit=iunit,fmt='(a)') "    _atom_site_fract_x"
       write(unit=iunit,fmt='(a)') "    _atom_site_fract_y"
       write(unit=iunit,fmt='(a)') "    _atom_site_fract_z"
       write(unit=iunit,fmt='(a)') "    _atom_site_U_iso_or_equiv"
       write(unit=iunit,fmt='(a)') "    _atom_site_occupancy"
       write(unit=iunit,fmt='(a)') "    _atom_site_adp_type"
       write(unit=iunit,fmt='(a)') "    _atom_site_type_symbol"

       !Calculation of the factor corresponding to the occupation factor provided in A
       do i=1,A%natoms
         occup(i)=A%Atom(i)%occ/(real(A%Atom(i)%mult)/real(SpG%multip))
         soccup(i)=A%Atom(i)%occ_std/(real(A%Atom(i)%mult)/real(SpG%multip))
       end do
       ocf=sum(abs(A%atom(1)%x-A%atom(2)%x))
       if( ocf < 0.001) then
         ocf=occup(1)+occup(2)
       else
         ocf=occup(1)
       end if
       occup=occup/ocf; soccup=soccup/ocf
       aniso=.false.
       do i=1,A%natoms
          line(1:132)=" "
          line(2:)= A%Atom(i)%Lab//"  "//A%Atom(i)%SfacSymb
           ! _atom_site_fract_x, _atom_site_fract_y, _atom_site_fract_z
          do j=1,3
            comm=" "
            call setnum_std(A%Atom(i)%x(j),A%Atom(i)%x_std(j),comm)
            line=trim(line)//" "//trim(comm)
          end do
          ! _atom_site_U_iso_or_equiv
            comm=" "
          if(A%Atom(i)%thtype == "isotr") then
            adptyp='Uiso'
            u=A%Atom(i)%Biso/(8.0*pi*pi)
            su=A%Atom(i)%Biso_std/(8.0*pi*pi)
            call setnum_std(u,su,comm)

          else if(A%Atom(i)%thtype == "aniso") then
            aniso=.true.
            adptyp='Uani'
            if(A%atom(i)%Utype == "beta") then
               aux=A%atom(i)%u
               ua=convert_betas_u(aux,cell)
               aux=A%atom(i)%u_std
               sua=convert_betas_u(aux,cell)
            else if(A%atom(i)%Utype == "u_ij") then
               ua=A%atom(i)%u
               sua=A%atom(i)%u_std
            end if
            u=(ua(1)+ua(2)+ua(3))/3.0
            su=(ua(1)+ua(2)+ua(3))/3.0
            call setnum_std(u,su,comm)
          else
            adptyp='.'
          end if
            line=trim(line)//" "//trim(comm)

           !_atom_site_occupancy
            comm=" "
            call setnum_std(occup(i),soccup(i),comm)
            line=trim(line)//" "//trim(comm)

          WRITE(iunit,"(a)") trim(line)//" "//trim(adptyp)//" "//A%atom(i)%SfacSymb
       end do
       if(aniso) then
          write(unit=iunit,fmt="(a)") " "
          write(unit=iunit,fmt="(a)") "loop_"
          write(unit=iunit,fmt="(a)") "    _atom_site_aniso_label "
          write(unit=iunit,fmt="(a)") "    _atom_site_aniso_U_11  "
          write(unit=iunit,fmt="(a)") "    _atom_site_aniso_U_22  "
          write(unit=iunit,fmt="(a)") "    _atom_site_aniso_U_33  "
          write(unit=iunit,fmt="(a)") "    _atom_site_aniso_U_12  "
          write(unit=iunit,fmt="(a)") "    _atom_site_aniso_U_13  "
          write(unit=iunit,fmt="(a)") "    _atom_site_aniso_U_23  "
          write(unit=iunit,fmt="(a)") "    _atom_site_aniso_type_symbol"
          do i=1,A%natoms
             if(A%Atom(i)%thtype /= "aniso") cycle
             line(1:132)=" "
             line(2:)= A%Atom(i)%Lab
             if(A%atom(i)%Utype == "beta") then
                aux=A%atom(i)%u
                ua=convert_betas_u(aux,cell)
                aux=A%atom(i)%u_std
                sua=convert_betas_u(aux,cell)
             else if(A%atom(i)%Utype == "u_ij") then
                ua=A%atom(i)%u
                sua=A%atom(i)%u_std
             end if
             do j=1,6
               comm=" "
               call setnum_std(ua(j),sua(j),comm)
               line=trim(line)//" "//trim(comm)
             end do
              WRITE(iunit,"(a)") trim(line)//"  "//A%atom(i)%SfacSymb
          end do
       end if

       if(type_data < 2) then
          write(unit=iunit,fmt="(a)") " "
          write(unit=iunit,fmt="(a)") "# Note: if the displacement parameters were refined anisotropically"
          write(unit=iunit,fmt="(a)") "# the U matrices should be given as for single-crystal studies."

          write(unit=iunit,fmt="(a)") " "
          write(unit=iunit,fmt="(a)") "#============================================================================="
          write(unit=iunit,fmt="(a)") " "

          !---- Molecular Geometry ----!
          write(unit=iunit,fmt="(a)") "# 10. MOLECULAR GEOMETRY"

          write(unit=iunit,fmt="(a)") " "


          write(unit=iunit,fmt="(a)") "_geom_special_details                ?"

          write(unit=iunit,fmt="(a)") " "
          write(unit=iunit,fmt="(a)") "loop_"
          write(unit=iunit,fmt="(a)") "    _geom_bond_atom_site_label_1  "
          write(unit=iunit,fmt="(a)") "    _geom_bond_atom_site_label_2  "
          write(unit=iunit,fmt="(a)") "    _geom_bond_site_symmetry_1    "
          write(unit=iunit,fmt="(a)") "    _geom_bond_site_symmetry_2    "
          write(unit=iunit,fmt="(a)") "    _geom_bond_distance           "
          write(unit=iunit,fmt="(a)") "    _geom_bond_publ_flag          "
          write(unit=iunit,fmt="(a)") "    ?   ?   ?   ?   ?   ?"

          write(unit=iunit,fmt="(a)") " "
          write(unit=iunit,fmt="(a)") "loop_"
          write(unit=iunit,fmt="(a)") "    _geom_contact_atom_site_label_1 "
          write(unit=iunit,fmt="(a)") "    _geom_contact_atom_site_label_2 "
          write(unit=iunit,fmt="(a)") "    _geom_contact_distance          "
          write(unit=iunit,fmt="(a)") "    _geom_contact_site_symmetry_1   "
          write(unit=iunit,fmt="(a)") "    _geom_contact_site_symmetry_2   "
          write(unit=iunit,fmt="(a)") "    _geom_contact_publ_flag         "
          write(unit=iunit,fmt="(a)") "    ?   ?   ?   ?   ?   ?"

          write(unit=iunit,fmt="(a)") " "
          write(unit=iunit,fmt="(a)") "loop_"
          write(unit=iunit,fmt="(a)") "_geom_angle_atom_site_label_1 "
          write(unit=iunit,fmt="(a)") "_geom_angle_atom_site_label_2 "
          write(unit=iunit,fmt="(a)") "_geom_angle_atom_site_label_3 "
          write(unit=iunit,fmt="(a)") "_geom_angle_site_symmetry_1   "
          write(unit=iunit,fmt="(a)") "_geom_angle_site_symmetry_2   "
          write(unit=iunit,fmt="(a)") "_geom_angle_site_symmetry_3   "
          write(unit=iunit,fmt="(a)") "_geom_angle                   "
          write(unit=iunit,fmt="(a)") "_geom_angle_publ_flag         "
          write(unit=iunit,fmt="(a)") "?   ?   ?   ?   ?   ?   ?   ?"

          write(unit=iunit,fmt="(a)") " "
          write(unit=iunit,fmt="(a)") "loop_"
          write(unit=iunit,fmt="(a)") "_geom_torsion_atom_site_label_1 "
          write(unit=iunit,fmt="(a)") "_geom_torsion_atom_site_label_2 "
          write(unit=iunit,fmt="(a)") "_geom_torsion_atom_site_label_3 "
          write(unit=iunit,fmt="(a)") "_geom_torsion_atom_site_label_4 "
          write(unit=iunit,fmt="(a)") "_geom_torsion_site_symmetry_1   "
          write(unit=iunit,fmt="(a)") "_geom_torsion_site_symmetry_2   "
          write(unit=iunit,fmt="(a)") "_geom_torsion_site_symmetry_3   "
          write(unit=iunit,fmt="(a)") "_geom_torsion_site_symmetry_4   "
          write(unit=iunit,fmt="(a)") "_geom_torsion                   "
          write(unit=iunit,fmt="(a)") "_geom_torsion_publ_flag         "
          write(unit=iunit,fmt="(a)") "?   ?   ?   ?   ?   ?   ?   ?   ?   ?"

          write(unit=iunit,fmt="(a)") " "
          write(unit=iunit,fmt="(a)") "loop_"
          write(unit=iunit,fmt="(a)") "_geom_hbond_atom_site_label_D "
          write(unit=iunit,fmt="(a)") "_geom_hbond_atom_site_label_H "
          write(unit=iunit,fmt="(a)") "_geom_hbond_atom_site_label_A "
          write(unit=iunit,fmt="(a)") "_geom_hbond_site_symmetry_D   "
          write(unit=iunit,fmt="(a)") "_geom_hbond_site_symmetry_H   "
          write(unit=iunit,fmt="(a)") "_geom_hbond_site_symmetry_A   "
          write(unit=iunit,fmt="(a)") "_geom_hbond_distance_DH       "
          write(unit=iunit,fmt="(a)") "_geom_hbond_distance_HA       "
          write(unit=iunit,fmt="(a)") "_geom_hbond_distance_DA       "
          write(unit=iunit,fmt="(a)") "_geom_hbond_angle_DHA         "
          write(unit=iunit,fmt="(a)") "_geom_hbond_publ_flag         "
          write(unit=iunit,fmt="(a)") "?   ?   ?   ?   ?   ?   ?   ?   ?   ?   ?"

          write(unit=iunit,fmt="(a)") " "
          write(unit=iunit,fmt="(a)") "#============================================================================="
          write(unit=iunit,fmt="(a)") " "


          !---- Final Informations ----!
          write(unit=iunit,fmt="(a)") "#============================================================================="
          write(unit=iunit,fmt="(a)") "# Additional structures (last six sections and associated data_? identifiers) "
          write(unit=iunit,fmt="(a)") "# may be added at this point.                                                 "
          write(unit=iunit,fmt="(a)") "#============================================================================="

          write(unit=iunit,fmt="(a)") " "
          write(unit=iunit,fmt="(a)") "# The following lines are used to test the character set of files sent by     "
          write(unit=iunit,fmt="(a)") "# network email or other means. They are not part of the CIF data set.        "
          write(unit=iunit,fmt="(a)") "# abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789              "
          write(unit=iunit,fmt="(a)") "# !@#$%^&*()_+{}:"//""""//"~<>?|\-=[];'`,./ "
       end if

       close(unit=iunit)

       return
    End Subroutine Write_Cif_Template

    !!----
    !!---- Subroutine Write_Shx_Template(Filename,Code,Title,Lambda,Z,Celda,Space,Atomos)
    !!----    character(len=*),        intent(in) :: filename  !  In -> Filename
    !!----    integer,                 intent(in) :: code      !  In -> 0 Shelxs-Patterson
    !!----                                                              1 Shelxs-Direct Methods
    !!----                                                              2 Shelxl-Refinement
    !!----    character(len=*),        intent(in) :: title     !  In -> Title
    !!----    real(kind=cp),           intent(in) :: lambda    !  In -> Lambda
    !!----    integer,                 intent(in) :: z         !  In -> Z
    !!----    type(Crystal_cell_Type), intent(in) :: celda     !  In -> Cell variable
    !!----    type(Space_Group_Type),  intent(in) :: Space     !  In -> SpaceGroup variable
    !!----    type(atom_list_type),    intent(in) :: atomos    !  In -> Atom List
    !!----
    !!----    Write a Shelx File
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Write_Shx_Template(filename,code,title,lambda,z,celda,space,atomos)
       !---- Arguments ----!
       character(len=*),        intent(in) :: filename
       integer,                 intent(in) :: code
       character(len=*),        intent(in) :: title
       real(kind=cp),           intent(in) :: lambda
       integer,                 intent(in) :: z
       type(Crystal_cell_Type), intent(in) :: celda
       type(Space_Group_Type),  intent(in) :: Space
       type(atom_list_type),    intent(in) :: atomos

       !---- Local Variables ----!
       logical                :: info

       integer                :: i,j,k,nc,iunit !,nlong
       integer                :: nlat
       integer, dimension(15) :: z_cont

       !---- Inicializacion de variables ----!
       info=.false.
       iunit=0
       z_cont=0
       nc=0  !this depends on scattering factor?

       !---- Esta abierto este Fichero? ----!
       inquire(file=filename,opened=info)
       if (info) then
          inquire(file=filename,number=iunit)
          close(unit=iunit)
       end if

       !---- Escritura ----!
       if (iunit == 0) iunit=61
       open(unit=iunit,file=filename,status="unknown",action="write")
       rewind(unit=iunit)

       !---- Title ----!
       write(unit=iunit,fmt="(a)") "TITL "//title(1:len_trim(title))

       !---- Lambda, Cell ----!
       write(unit=iunit,fmt="(a,f8.5,3f9.4,3f7.3)") "CELL ",lambda,celda%cell,celda%ang

       !---- Z, Std ----!
       write(unit=iunit,fmt="(a,i3,a,3f8.4,3f7.3)") "ZERR ",z,"     ",celda%cell_std,celda%ang_std

       !---- Latt ----!
       nlat=1
       select case (space%centred)
          case (0) ! Centric

          case (1) ! Acentric
             nlat=-1

          case (2) ! Not used in Shelx
             write(unit=iunit,fmt="(a)") " ERROR: Origin not at -1 "
             close(unit=iunit)
             return

       end select
       select case (space%spg_lat)
          case ("P")

          case ("I")
             nlat=2*nlat

          case ("R")
             nlat=3*nlat

          case ("F")
             nlat=4*nlat

          case ("A")
             nlat=5*nlat

          case ("B")
             nlat=6*nlat

          case ("C")
             nlat=7*nlat

       end select
       write(unit=iunit,fmt="(a,i2)") "LATT ",nlat

       !---- Symm ----!
       do i=2,space%numops
          write(unit=iunit,fmt="(a)") "SYMM "//u_case(space%symopsymb(i))
       end do

       !---- Sfac ----!
       j=0
       do i=1,atomos%natoms
          if (j == 0) then
             j=1
             z_cont(j)=atomos%atom(i)%z
          else
             do k=1,j
                if (z_cont(k) == atomos%atom(i)%z) exit
             end do
             if (z_cont(k) /= atomos%atom(i)%z) then
                j=j+1
                z_cont(j)=atomos%atom(i)%z
             end if
          end if
       end do


       write(unit=iunit,fmt="(a)") "SFAC "

       !---- Unit ----!
       write(unit=iunit,fmt="(a)") "UNIT "

       select case (code)
          case (0) ! Shelxs - Patterson
             write(unit=iunit,fmt="(a)") "PATT "

          case (1) ! Shelxs - Direct Methods
             write(unit=iunit,fmt="(a)") "TREF "

          case (2) ! Shelxl - Refinement
             !---- L.S. ----!
             write(unit=iunit,fmt="(a)") "L.S. 10"

             !---- Fvar ----!
             write(unit=iunit,fmt="(a)") "FVAR 1.0"

             !---- Weight ----!
             write(unit=iunit,fmt="(a)") "WGHT 0.2"

             !---- Fmap ----!
             write(unit=iunit,fmt="(a)") "FMAP 2"

             !---- Atoms ----!
             do i=1,atomos%natoms
                write(unit=iunit,fmt="(a4,i3,4f11.5)") &
                     atomos%atom(i)%lab, nc, atomos%atom(i)%x, atomos%atom(i)%occ+10.0
             end do
       end select

       !---- Format ----!
       write(unit=iunit,fmt="(a)") "HKLF 4"

       !---- End ----!
       write(unit=iunit,fmt="(a)") "END "

       return
    End Subroutine Write_Shx_Template

    !!----
    !!---- Subroutine Get_Phases_File(filecode, Nphas, PhasesName,ILines)
    !!----    character(len=*),                intent(in)   :: filecode
    !!----    Integer,                         intent(out)  :: Nphas
    !!----    Character(len=80), dimension(:), intent(out)  :: PhasesName
    !!----    Integer,dimension(2,:),          intent(out)  :: ILines
    !!----
    !!---- Determine how many phases there are in a CIF or PCR file and
    !!---- give the lines to locate
    !!----
    !!---- Update: 01/05/2013
    !!
    Subroutine Get_Phases_File(filecode, NPhas, PhasesName,ILines)
       !---- Arguments ----!
       character(len=*),             intent(in)   :: filecode
       integer,                      intent(out)  :: Nphas
       character(len=*),dimension(:),intent(out)  :: PhasesName
       integer,dimension(:,:),       intent(out)  :: ILines

       !---- Local Variables ----!
       character(len=3) :: ext
       integer          :: npos

       !> Error
       call init_err_form()

       !> Init
       Nphas=0
       PhasesName=' '
       Ilines=0

       !> PCR or CIF file
       npos=index(filecode,'.',back=.true.)
       if (npos <=0) then
          err_form=.true.
          err_form_mess='No extension was found in the name of the file!'
          return
       end if

       ext=filecode(npos+1:)
       ext=u_case(ext)
       select case (ext)
          case ('CIF')
             call get_nphases_ciffile(filecode, NPhas, PhasesName,ILines)
          case ('PCR')
             call get_nphases_pcrfile(filecode, NPhas, PhasesName,ILines)
          case default
             err_form=.true.
             err_form_mess='Extension for this file not valid!'
       end select

       return
    End Subroutine Get_Phases_File

    !!--++
    !!--++ Subroutine Get_NPhases_CIFFile(Filecode,NPhas,PhasesName,ILines)
    !!--++    character(len=*),                 intent(in)  :: Filecode    ! Filename
    !!--++    integer,                          intent(out) :: NPhas       ! Number of Phases in the file
    !!--++    character(len=*), dimension(:),   intent(out) :: PhasesName     ! Name of Phases in the file
    !!--++    integer,          dimension(:,:), intent(out) :: ILines        ! Index for lines for each Phase
    !!--++
    !!--++ Determine the number of phases are included into the file
    !!--++
    !!--++ Date: 01/05/2013
    !!
    Subroutine Get_NPhases_CIFFile(Filecode,NPhas,PhasesName,ILines)
       !---- Arguments ----!
       character(len=*),                 intent(in)  :: Filecode    ! Filename
       integer,                          intent(out) :: NPhas       ! Number of Phases in the file
       character(len=*), dimension(:),   intent(out) :: PhasesName     ! Name of Phases in the file
       integer,          dimension(:,:), intent(out) :: ILines        ! Index for lines for each Phase

       !---- Local Variables ----!
       character(len=150), dimension(:), allocatable :: filen
       character(len=150)                            :: line
       integer                                       :: i,j,nl

       !> Error
       call init_err_form()

       !> Initialize
       NPhas=0
       PhasesName=' '
       ILines=0

       !> Reading file
       nl=0
       call number_lines(trim(filecode),nl)
       if (nl <=0) then
          err_form=.true.
          err_form_mess='No lines were read for '//trim(filecode)//' !!'
          return
       end if
       allocate(filen(nl))
       call reading_lines(trim(filecode),nl,filen)

       !> Number of Phases
       do i=1,nl
          line=adjustl(filen(i))

          !> empty line
          if (len_trim(line) <= 0) cycle

          !> comment line
          if (line(1:1) =='#') cycle

          !> No data_global
          j=index(line,'data_global')
          if (j > 0) cycle

          !> Just only lines beginning with data...
          j=index(line,'data_')
          if (j /= 1) cycle

          nphas=nphas+1
          ILines(1,Nphas)=i
          PhasesName(nphas)=trim(line(j+5:))
          if (nphas > 1) ILines(2,nphas-1)=i-1
       end do
       if (nphas > 0 .and. ILines(2,nphas)==0) ILines(2,nphas)=nl

       if (allocated(filen)) deallocate(filen)

       return
    End Subroutine  Get_NPhases_CIFFile

    !!--++
    !!--++ Subroutine Get_NPhases_PCRFile(filecode, Nphas,PhasesName,ILines)
    !!--++    character(len=*),                intent(in)   :: filecode
    !!--++    Integer,                         intent(out)  :: Nphas
    !!--++    Character(len=80), dimension(:), intent(out)  :: PhasesName
    !!--++    Integer,dimension(2,:),          intent(out)  :: ILines
    !!--++
    !!--++ Determine how many phases and where there in a PCR file
    !!--++
    !!--++ Update: 01/05/2013
    !!
    Subroutine Get_NPhases_PCRFile(filecode, NPhas, PhasesName,ILines)
       !---- Arguments ----!
       character(len=*),             intent(in)   :: filecode
       integer,                      intent(out)  :: Nphas
       character(len=*),dimension(:),intent(out)  :: PhasesName
       integer,dimension(:,:),       intent(out)  :: ILines

       !---- Local Variables ----!
       logical                                      :: multi, ask_phase
       character(len=80), dimension(:), allocatable :: file_dat
       character(len=80)                            :: line
       integer                                      :: i,k,iv,nlines
       integer, dimension(30)                       :: ivet
       real(kind=cp), dimension(30)                 :: vet

       !> Err
       call init_err_form()

       !> Init
       NPhas=0
       PhasesName=' '
       ILines=0

       !> Reading file
       nlines=0
       call number_lines(trim(filecode),nlines)
       if (nlines <=0) then
          err_form=.true.
          err_form_mess='No lines were read for '//trim(filecode)//' !!'
          return
       end if
       allocate(file_dat(nlines))
       call reading_lines(trim(filecode),nlines,file_dat)

       ILines(1,:)=1
       ILines(2,:)=nlines

       !> Simple / Multi format
       multi=.false.
       do i=1,nlines
          line=adjustl(file_dat(i))
          if (line(1:1) =='!' .or. line(1:1)==' ') cycle
          if (index(line,'NPATT ') <=0) cycle
          multi=.true.
       end do

       !> Number of Phases
       if (.not. multi) then
          do i=2,nlines
             line=adjustl(file_dat(i))
             if (line(1:1) =='!' .or. line(1:1)==' ') cycle
             call getnum(line,vet,ivet,iv)
             if (iv > 3) then
                NPhas=ivet(3)
                exit
             end if
          end do

       else
          do i=1,nlines
             line=adjustl(file_dat(i))
             if (line(1:4) /='!Nph') cycle

             line=adjustl(file_dat(i+1))
             call getnum(line,vet,ivet,iv)
             if (iv > 1) then
                NPhas=ivet(1)
                exit
             end if
          end do
       end if

       if (NPhas == 0) then
          err_form=.true.
          err_form_mess=" No Phase information was found in this PCR file. Please, check it! "
          return
       end if

       !> Locate where begin each Phase
       k=0
       ask_phase=.true.

       do i=1,nlines
          line=adjustl(file_dat(i))
          if (ask_phase) then
             if (index(line,'Data for PHASE') <= 0) cycle
          else
             if (line(1:1) /='!') then
                k=k+1
                ILines(1,k)=i
                PhasesName(k)=trim(adjustl(line))
                if (k == NPhas) exit

                ask_phase=.true.
             end if
             cycle
          end if
          ask_phase=.false.
       end do

       if (NPhas /= k) then
          err_form=.true.
          err_form_mess=" Locating Phases failed in this PCR. Please, check it!"
          return
       end if

       do i=1,Nphas
          if (nphas > 1) then
             ilines(2,i)=ilines(1,i+1)-1
          end if
       end do

       return
    End Subroutine Get_NPhases_PCRFile
    !!----
    !!---- Subroutine Write_CFL(lun,Cel,SpG,Atm,comment)
    !!----    integer,                  intent(in)    :: lun
    !!----    type (Space_Group_Type),  intent(in)    :: SpG
    !!----    type (Crystal_Cell_Type), intent(in)    :: Cel
    !!----    type (atom_list_type),    intent(in)    :: Atm
    !!----    character(len=*),optional,intent(in)    :: comment
    !!----
    !!----    (OVERLOADED)
    !!----
    !!----    Write a CFL-file with atom_list_type
    !!----
    !!---- Update: July - 2014
    !!
    Subroutine Write_CFL_Atom_List_Type(lun,Cel,SpG,Atm,comment)
       !---- Arguments ----!
       integer,                  intent(in)    :: lun
       type (Space_Group_Type),  intent(in)    :: SpG
       type (Crystal_Cell_Type), intent(in)    :: Cel
       type (atom_list_type),    intent(in)    :: Atm
       character(len=*),optional,intent(in)    :: comment

       !----- Local variables -----!
       integer                         :: j !,loc
       real(kind=cp), dimension(6)     :: a,sa
       character(len=30), dimension(6) :: text

       if(present(comment)) write(unit=lun,fmt="(a)") "TITLE "//trim(comment)
       write(unit=lun,fmt="(a)") "!  Automatically generated CFL file (Write_CFL)"

       a(1:3)=Cel%Cell
       a(4:6)=Cel%ang
       sa(1:3)=Cel%Cell_std
       sa(4:6)=Cel%ang_std
       do j=1,6
          call SetNum_Std(a(j), sa(j), text(j))
       end do
       write(unit=lun,fmt="(a)") "!         a               b               c            alpha           beta            gamma"
       write(unit=lun,fmt="(a,6a16)") "Cell ",text
       write(unit=lun,fmt="(a,i3)")"!     Space Group # ",SpG%NumSpg
       write(unit=lun,fmt="(a,a)") "Spgr  ",SpG%SPG_Symb
       call Write_Atoms_CFL(Atm,Lun,cel)

       return
    End Subroutine Write_CFL_Atom_List_Type
    !!----
    !!---- Subroutine Write_CFL(lun,Molx,comment)
    !!----    integer,                       intent(in) :: lun
    !!----    type (Molecular_Crystal_Type), intent(in) :: Molx
    !!----    character(len=*),optional,     intent(in) :: comment
    !!----
    !!----    (OVERLOADED)
    !!----
    !!----    Write a CFL-file with molecular_crystal_type
    !!----
    !!---- Update: July - 2014
    !!
    Subroutine Write_CFL_Molcrys(lun,Molx,comment)
       !---- Arguments ----!
       integer,                       intent(in) :: lun
       type (Molecular_Crystal_Type), intent(in) :: Molx
       character(len=*),optional,     intent(in) :: comment

       !----- Local variables -----!
       integer                         :: j !,loc
       real(kind=cp), dimension(6)     :: a,sa
       character(len=30), dimension(6) :: text

       if(present(comment)) write(unit=lun,fmt="(a)") "TITLE "//trim(comment)
       write(unit=lun,fmt="(a)") "!  Automatically generated CFL file (Write_CFL)"

       a(1:3)=molx%cell%Cell
       a(4:6)=molx%cell%ang
       sa(1:3)=molx%cell%Cell_std
       sa(4:6)=molx%cell%ang_std
       do j=1,6
          call SetNum_Std(a(j), sa(j), text(j))
       end do
       write(unit=lun,fmt="(a)") "!         a               b               c            alpha           beta            gamma"
       write(unit=lun,fmt="(a,6a16)") "Cell ",text
       write(unit=lun,fmt="(a,i3)")"!     Space Group # ",molx%spg%NumSpg
       write(unit=lun,fmt="(a,a)") "Spgr  ",molx%spg%SPG_Symb
       call Write_Atoms_CFL(Molx,Lun)

       return
    End Subroutine Write_CFL_Molcrys
    !!----
    !!---- Subroutine Write_Atoms_CFL(Ats,Lun,Cell)
    !!----    Type (atom_list_type),dimension(:),  intent(in) :: Ats     !  In -> Atom List
    !!----    integer, optional,                   intent(in) :: lun     !  In -> Unit to write
    !!----    Type(Crystal_Cell_Type), optional,   intent(in) :: Cell    !  In -> Transform to thermal parameters
    !!----
    !!----    Write the atoms in the asymmetric unit for a CFL file
    !!----
    !!---- Update: February - 2003
    !!
    Subroutine Write_Atoms_CFL_ATM(Ats,Lun,cell)
       !---- Arguments ----!
       type (atom_list_type),            intent(in) :: Ats
       integer, optional,                intent(in) :: Lun
       Type(Crystal_Cell_Type), optional,intent(in) :: Cell

       !---- Local Variables ----!
       character(len=30),dimension(6) :: text
       character(len=36)              :: forma,fom
       integer                        :: i, j, iunit, leng, maxl,ish
       real(kind=cp), dimension(6)    :: u,bet,sb

       iunit=6
       if (present(lun)) iunit=lun

       if(ats%natoms == 0) then
         write (unit=iunit,fmt="(a)") "!  No atoms ..."
         return
       end if
       !Determine the maximum length of the atom labels
       maxl=0
       do i=1,ats%natoms
         leng=len_trim(ats%atom(i)%lab)
         if(leng > maxl) maxl=leng
       end do
       maxl=max(maxl,4)+1
       ish=maxl-4
       fom   ="(a,tr  ,a)"
       Select Case(ish)
          Case(:9)
            write(unit=fom(6:6),fmt="(i1)") ish
          Case(10:)
            write(unit=fom(6:7),fmt="(i2)") ish
       End Select
       forma="(a,a  ,tr2,a,tr3,5a14,2f8.2,tr3,a)"
       Select Case(maxl)
         Case(:9)
             write(unit=forma(5:5),fmt="(i1)") maxl
         Case(10:)
             write(unit=forma(5:6),fmt="(i2)") maxl
       End Select
       write (unit=iunit,fmt=fom) "!     ", &
             "Atom  Type     x/a           y/b           z/c           Biso          Occ           Spin    Charge    Info"
       do i=1,ats%natoms

          do j=1,3
             call SetNum_Std(ats%atom(i)%x(j), ats%atom(i)%x_std(j), text(j))
          end do
          call SetNum_Std(ats%atom(i)%Biso, ats%atom(i)%Biso_std, text(4))
          call SetNum_Std(ats%atom(i)%Occ, ats%atom(i)%Occ_std, text(5))

          write (unit=iunit,fmt=forma) &
                "Atom   ",trim(ats%atom(i)%lab),ats%atom(i)%chemsymb, (text(j),j=1,5), &
                 ats%atom(i)%moment,ats%atom(i)%charge,"# "//ats%atom(i)%AtmInfo

          if (ats%atom(i)%thtype == "aniso") then

             if (ats%atom(i)%utype == "beta") then
                bet=ats%atom(i)%u(1:6)
                sb=ats%atom(i)%u_std(1:6)
                do j=1,6
                   call SetNum_Std(bet(j), sb(j), text(j))
                end do
                write (unit=iunit,fmt="(a,tr1,6a14)") "Beta  ", text
                if (present(Cell)) then
                   u=convert_betas_u(bet,cell)
                   sb=convert_betas_u(ats%atom(i)%u_std,cell)
                   do j=1,6
                      call SetNum_Std(u(j), sb(j), text(j))
                   end do
                   write(unit=iunit,fmt="(a,6a14)") "!U_ij  ", text
                end if

             else if(ats%atom(i)%thtype == "u_ij") then
                u=ats%atom(i)%u(1:6)
                sb=ats%atom(i)%u_std(1:6)
                do j=1,6
                   call SetNum_Std(u(j), sb(j), text(j))
                end do
                write(unit=iunit,fmt="(a,6a14)") "U_ij  ", text
                if (present(Cell)) then
                   bet=convert_u_betas(u,cell)
                   sb=convert_u_betas(ats%atom(i)%u_std,cell)
                   do j=1,6
                      call SetNum_Std(bet(j), sb(j), text(j))
                   end do
                   write(unit=iunit,fmt="(a,6a14)") "!Beta  ", text
                end if
             end if

          end if
       end do

       return
    End Subroutine Write_Atoms_CFL_ATM

    !!----
    !!---- Subroutine Write_Atoms_CFL(Ats,Lun,Cell)
    !!----    Type (atom_list_type),dimension(:),  intent(in) :: Ats     !  In -> Atom List
    !!----    integer, optional,                   intent(in) :: lun     !  In -> Unit to write
    !!----    Type(Crystal_Cell_Type), optional,   intent(in) :: Cell    !  In -> Transform to thermal parameters
    !!----
    !!----    Write the atoms in the asymmetric unit for a CFL file
    !!----
    !!---- Update: February - 2003
    !!
    Subroutine Write_Atoms_CFL_MOLX(Molx,Lun)
        !---- Arguments ----!
        type (Molecular_Crystal_Type), intent(in) :: Molx
        integer, optional,             intent(in) :: Lun

        !---- Local Variables ----!
        character(len=30),dimension(6) :: text
        character(len=36)              :: forma,fom
        integer                        :: i, j, iunit, leng, maxl,ish
        real(kind=cp), dimension(6)    :: u,bet,sb

        iunit=6
        if (present(lun)) iunit=lun

        if(molx%n_free > 0) then
            !Determine the maximum length of the atom labels
            maxl=0
            do i=1,molx%n_free
                leng=len_trim(molx%atm(i)%lab)
                if(leng > maxl) maxl=leng
            end do
            maxl=max(maxl,4)+1
            ish=maxl-4
            fom   ="(a,tr  ,a)"
            Select Case(ish)
                Case(:9)
                    write(unit=fom(6:6),fmt="(i1)") ish
                Case(10:)
                    write(unit=fom(6:7),fmt="(i2)") ish
            End Select
            forma="(a,a  ,tr2,a,tr3,5a14,2f8.2,tr3,a)"
            Select Case(maxl)
                Case(:9)
                    write(unit=forma(5:5),fmt="(i1)") maxl
                Case(10:)
                    write(unit=forma(5:6),fmt="(i2)") maxl
            End Select
            write (unit=iunit,fmt=fom) "!     ", &
                  "Atom  Type     x/a           y/b           z/c           Biso          Occ           Spin    Charge    Info"
            do i=1,molx%n_free

                do j=1,3
                   call SetNum_Std(molx%atm(i)%x(j), molx%atm(i)%x_std(j), text(j))
                end do
                call SetNum_Std(molx%atm(i)%Biso, molx%atm(i)%Biso_std, text(4))
                call SetNum_Std(molx%atm(i)%Occ, molx%atm(i)%Occ_std, text(5))

                write (unit=iunit,fmt=forma) &
                      "Atom   ",trim(molx%atm(i)%lab),molx%atm(i)%chemsymb, (text(j),j=1,5), &
                       molx%atm(i)%moment,molx%atm(i)%charge,"# "//molx%atm(i)%AtmInfo

                if (molx%atm(i)%thtype == "aniso") then

                    if (molx%atm(i)%utype == "beta") then
                        bet=molx%atm(i)%u(1:6)
                        sb=molx%atm(i)%u_std(1:6)
                        do j=1,6
                            call SetNum_Std(bet(j), sb(j), text(j))
                        end do
                        write (unit=iunit,fmt="(a,tr1,6a14)") "Beta  ", text
                        u=convert_betas_u(bet,molx%cell)
                        sb=convert_betas_u(molx%atm(i)%u_std,molx%cell)
                        do j=1,6
                            call SetNum_Std(u(j), sb(j), text(j))
                        end do
                        write(unit=iunit,fmt="(a,6a14)") "!U_ij  ", text
                    else if(molx%atm(i)%thtype == "u_ij") then
                        u=molx%atm(i)%u(1:6)
                        sb=molx%atm(i)%u_std(1:6)
                        do j=1,6
                            call SetNum_Std(u(j), sb(j), text(j))
                        end do
                        write(unit=iunit,fmt="(a,6a14)") "U_ij  ", text
                        bet=convert_u_betas(u,molx%cell)
                        sb=convert_u_betas(molx%atm(i)%u_std,molx%cell)
                        do j=1,6
                            call SetNum_Std(bet(j), sb(j), text(j))
                        end do
                        write(unit=iunit,fmt="(a,6a14)") "!Beta  ", text
                    end if
                end if
            end do ! i=1,molx%n_free
        end if ! molx%n_free > 0

        if (molx%n_mol > 0) then
            do i=1,molx%n_mol
                write(unit=iunit,fmt="(/,a,tr2,i3,tr2,a,tr2,a)") &
                     "MOLEX",molx%mol(i)%natoms,trim(molx%mol(i)%Name_mol),molx%mol(i)%coor_type
                write(unit=iunit,fmt="(a)") &
                     "!    Xc         Yc          Zc        Phi        Theta      Chi     TypeAngles TypeThermal"
                write(unit=iunit,fmt="(6f11.5,tr6,a,tr10,a)") &
                     molx%mol(i)%xcentre,molx%mol(i)%orient,molx%mol(i)%rot_type,molx%mol(i)%therm_type
                write(unit=iunit,fmt="(t1,6i10,tr2,a)") &
                     molx%mol(i)%lxcentre,molx%mol(i)%lorient," ! Refinemencodes"

                select case (molx%mol(i)%coor_type)
                    case ("C","c")
                        write(unit=iunit,fmt="(a)") &
                        "!Atom   Type        XC          YC          ZC    N1  N2  N3      Biso        Occ "
                    case ("F","f")
                        write(unit=iunit,fmt="(a)") &
                        "!Atom   Type        X           Y           Z     N1  N2  N3      Biso        Occ "
                    case ("S","s")
                        write(unit=iunit,fmt="(a)") &
                        "!Atom   Type    distance      Theta       Phi     N1  N2  N3      Biso        Occ "
                    case ("Z","z")
                        write(unit=iunit,fmt="(a)") &
                        "!Atom   Type    distance  Bond-Angle Torsion-Ang  N1  N2  N3      Biso        Occ "
                    case default
                        write(unit=iunit,fmt="(a)") &
                        "!Atom   Type      Coor1       Coor2       Coor3   N1  N2  N3      Biso        Occ "
                end select ! molx%mol(i)%coor_type

                do j=1,molx%mol(i)%natoms
                    write(unit=iunit,fmt="(a,tr2,a,3f12.5,3i4,2f12.5)")  &
                          molx%mol(i)%AtName(j), molx%mol(i)%AtSymb(j),molx%mol(i)%I_Coor(:,j),  &
                          molx%mol(i)%Conn(:,j), molx%mol(i)%Biso(j),  molx%mol(i)%Occ(j)
                end do ! j = molx%mol(i)%natoms
            end do ! i = 1,molx%n_mol
        end if ! molx%n_mol > 0
        return
    End Subroutine Write_Atoms_CFL_MOLX
    !!----
    !!---- Subroutine Write_Atoms_CFL(Ats,Lun,Cell)
    !!----    Type (atom_list_type),dimension(:),  intent(in) :: Ats     !  In -> Atom List
    !!----    integer, optional,                   intent(in) :: lun     !  In -> Unit to write
    !!----    Type(Crystal_Cell_Type), optional,   intent(in) :: Cell    !  In -> Transform to thermal parameters
    !!----
    !!----    Write the atoms in the asymmetric unit for a CFL file
    !!----
    !!---- Update: February - 2003
    !!
    Subroutine Write_Atoms_CFL_MOLX_orig(Molx,Lun)
       !---- Arguments ----!
       type (Molecular_Crystal_Type), intent(in) :: Molx
       integer, optional,             intent(in) :: Lun

       !---- Local Variables ----!
       character(len=30),dimension(6) :: text
       character(len=36)              :: forma,fom
       integer                        :: i, j, iunit, leng, maxl,ish
       real(kind=cp), dimension(6)    :: u,bet,sb

       iunit=6
       if (present(lun)) iunit=lun

       if(molx%n_free == 0) then
         write (unit=iunit,fmt="(a)") "!  No atoms ..."
         return
       end if
       !Determine the maximum length of the atom labels
       maxl=0
       do i=1,molx%n_free
         leng=len_trim(molx%atm(i)%lab)
         if(leng > maxl) maxl=leng
       end do
       maxl=max(maxl,4)+1
       ish=maxl-4
       fom   ="(a,tr  ,a)"
       Select Case(ish)
          Case(:9)
            write(unit=fom(6:6),fmt="(i1)") ish
          Case(10:)
            write(unit=fom(6:7),fmt="(i2)") ish
       End Select
       forma="(a,a  ,tr2,a,tr3,5a14,2f8.2,tr3,a)"
       Select Case(maxl)
         Case(:9)
             write(unit=forma(5:5),fmt="(i1)") maxl
         Case(10:)
             write(unit=forma(5:6),fmt="(i2)") maxl
       End Select
       write (unit=iunit,fmt=fom) "!     ", &
             "Atom  Type     x/a           y/b           z/c           Biso          Occ           Spin    Charge    Info"
       do i=1,molx%n_free

          do j=1,3
             call SetNum_Std(molx%atm(i)%x(j), molx%atm(i)%x_std(j), text(j))
          end do
          call SetNum_Std(molx%atm(i)%Biso, molx%atm(i)%Biso_std, text(4))
          call SetNum_Std(molx%atm(i)%Occ, molx%atm(i)%Occ_std, text(5))

          write (unit=iunit,fmt=forma) &
                "Atom   ",trim(molx%atm(i)%lab),molx%atm(i)%chemsymb, (text(j),j=1,5), &
                 molx%atm(i)%moment,molx%atm(i)%charge,"# "//molx%atm(i)%AtmInfo

          if (molx%atm(i)%thtype == "aniso") then

             if (molx%atm(i)%utype == "beta") then
                bet=molx%atm(i)%u(1:6)
                sb=molx%atm(i)%u_std(1:6)
                do j=1,6
                   call SetNum_Std(bet(j), sb(j), text(j))
                end do
                write (unit=iunit,fmt="(a,tr1,6a14)") "Beta  ", text
                u=convert_betas_u(bet,molx%cell)
                sb=convert_betas_u(molx%atm(i)%u_std,molx%cell)
                do j=1,6
                    call SetNum_Std(u(j), sb(j), text(j))
                end do
                write(unit=iunit,fmt="(a,6a14)") "!U_ij  ", text
             else if(molx%atm(i)%thtype == "u_ij") then
                u=molx%atm(i)%u(1:6)
                sb=molx%atm(i)%u_std(1:6)
                do j=1,6
                   call SetNum_Std(u(j), sb(j), text(j))
                end do
                write(unit=iunit,fmt="(a,6a14)") "U_ij  ", text
                bet=convert_u_betas(u,molx%cell)
                sb=convert_u_betas(molx%atm(i)%u_std,molx%cell)
                do j=1,6
                    call SetNum_Std(bet(j), sb(j), text(j))
                end do
                write(unit=iunit,fmt="(a,6a14)") "!Beta  ", text
             end if
          end if
       end do

       return
    End Subroutine Write_Atoms_CFL_MOLX_orig

    Subroutine Write_MCIF(Ipr,mCell,MSGp,Am,Cell)
       Integer,                         intent(in)           :: Ipr
       type(Magnetic_Space_Group_Type), intent(in)           :: MSGp
       type(Crystal_Cell_Type),         intent(in)           :: mCell
       type(Atom_List_Type),            intent(in)           :: Am
       type(Crystal_Cell_Type),optional,intent(in)           :: Cell
       !
       Character(len=132)             :: line
       character(len=80),dimension(6) :: text
       character(len=2)               :: invc
       real(kind=cp)                  :: occ,occ_std,uiso,uiso_std
       integer :: i,j

       write(unit=Ipr,fmt="(a)") "#  --------------------------------------"
       write(unit=Ipr,fmt="(a)") "#  Magnetic CIF file generated by CrysFML"
       write(unit=Ipr,fmt="(a)") "#  --------------------------------------"
       write(unit=Ipr,fmt="(a)") "# https://forge.epn-campus.eu/projects/crysfml/repository"
       call Write_Date_Time(dtim=line)
       write(unit=Ipr,fmt="(a)") trim(line)
       write(unit=Ipr,fmt="(a)") " "

       write(unit=Ipr,fmt="(a)") "data_"
       write(unit=Ipr,fmt="(a)") "_citation_journal_abbrev ?"
       write(unit=Ipr,fmt="(a)") "_citation_journal_volume ?"
       write(unit=Ipr,fmt="(a)") "_citation_page_first     ?"
       write(unit=Ipr,fmt="(a)") "_citation_page_last      ?"
       write(unit=Ipr,fmt="(a)") "_citation_article_id     ?"
       write(unit=Ipr,fmt="(a)") "_citation_year           ?"
       write(unit=Ipr,fmt="(a)") "_loop "
       write(unit=Ipr,fmt="(a)") "_citation_author_name"
       write(unit=Ipr,fmt="(a)") "?"
       write(unit=Ipr,fmt="(a)")
       write(unit=Ipr,fmt="(a)") "_atomic_positions_source_database_code_ICSD  ?"
       write(unit=Ipr,fmt="(a)") "_atomic_positions_source_other    .  "
       write(unit=Ipr,fmt="(a)")
       write(unit=Ipr,fmt="(a)") "_Neel_temperature  ?"
       write(unit=Ipr,fmt="(a)") "_magn_diffrn_temperature  ?"
       write(unit=Ipr,fmt="(a)") "_exptl_crystal_magnetic_properties_details"
       write(unit=Ipr,fmt="(a)") ";"
       write(unit=Ipr,fmt="(a)") ";"
       write(unit=Ipr,fmt="(a)") "_active_magnetic_irreps_details"
       write(unit=Ipr,fmt="(a)") ";"
       write(unit=Ipr,fmt="(a)") ";"
       write(unit=Ipr,fmt="(a)") " "
       if(MSGp%standard_setting) then
          write(unit=Ipr,fmt="(a)") "_magnetic_space_group_standard_setting  'yes'"
       else
          write(unit=Ipr,fmt="(a)") "_magnetic_space_group_standard_setting  'no'"
       end if
       write(unit=Ipr,fmt="(a)")    '_parent_space_group.name_H-M  "'//trim(MSGp%Parent_spg)//'"'
       write(unit=Ipr,fmt="(a,i3)") "_parent_space_group.IT_number  ",MSGp%Parent_num
       write(unit=Ipr,fmt="(a)")    "_magnetic_space_group.transform_from_parent_Pp_abc  '"//trim(MSGp%trn_from_parent)//"'"
       write(unit=Ipr,fmt="(a)")    "_magnetic_space_group.transform_to_standard_Pp_abc  '"//trim(MSGp%trn_to_standard)//"'"
       write(unit=Ipr,fmt="(a)")
       if(len_trim(MSGp%BNS_number) /= 0) &
       write(unit=Ipr,fmt="(a)") "_space_group.magn_number_BNS  "//trim(MSGp%BNS_number)
       if(len_trim(MSGp%BNS_symbol) /= 0) &
       write(unit=Ipr,fmt="(a)") '_space_group.magn_name_BNS  "'//trim(MSGp%BNS_symbol)//'"'
       if(len_trim(MSGp%OG_number) /= 0) &
       write(unit=Ipr,fmt="(a)") '_space_group.magn_number_OG '//trim(MSGp%OG_number)
       if(len_trim(MSGp%OG_symbol) /= 0) &
       write(unit=Ipr,fmt="(a)") '_space_group.magn_name_OG  "'//trim(MSGp%OG_symbol)//'"'
       write(unit=Ipr,fmt="(a)")

       if(MSGp%n_irreps /= 0) then
          write(unit=Ipr,fmt="(a)") "loop_"
          write(unit=Ipr,fmt="(a)") "_irrep_id"
          write(unit=Ipr,fmt="(a)") "_irrep_dimension"
          if( any(MSGp%small_irrep_dim > 0) ) write(unit=Ipr,fmt="(a)") "_small_irrep_dimension"
          write(unit=Ipr,fmt="(a)") "_irrep_direction_type"
          write(unit=Ipr,fmt="(a)") "_irrep_action"
          if( any(MSGp%irrep_modes_number > 0) ) write(unit=Ipr,fmt="(a)") "_irrep_modes_number"
          do i=1,MSGp%n_irreps
            if(MSGp%small_irrep_dim(i) > 0) then
               write(unit=line,fmt=("(2i4)"))  MSGp%irrep_dim(i), MSGp%small_irrep_dim(i)
            else
               write(unit=line,fmt=("(i4)"))  MSGp%irrep_dim(i)
            end if
            line= trim(MSGp%irrep_id(i))//"  "//trim(line)//"   "// &
                                      trim(MSGp%irrep_direction(i))//"  "//trim(MSGp%irrep_action(i))
            if( MSGp%irrep_modes_number(i) > 0) then
               j=len_trim(line)
              write(unit=line(j+1:),fmt="(i4)") MSGp%irrep_modes_number(i)
            end if
            write(unit=Ipr,fmt="(a)") trim(line)
          end do
          write(unit=Ipr,fmt="(a)")
       else
          write(unit=Ipr,fmt="(a)") "loop_"
          write(unit=Ipr,fmt="(a)") "_irrep_id"
          write(unit=Ipr,fmt="(a)") "_irrep_dimension"
          write(unit=Ipr,fmt="(a)") "_small_irrep_dimension"
          write(unit=Ipr,fmt="(a)") "_irrep_direction_type"
          write(unit=Ipr,fmt="(a)") "_irrep_action"
          write(unit=Ipr,fmt="(a)") "_irrep_modes_number"
          write(unit=Ipr,fmt="(a)") " ?  ?  ?  ?  ?  ?"
          write(unit=Ipr,fmt="(a)")
       end if

       if(MSGp%m_cell) then
          do i=1,3
            call setnum_std(mCell%Cell(i),mCell%cell_std(i),text(i))
            call setnum_std(mCell%ang(i),mCell%ang_std(i),text(i+3))
          end do
          write(unit=Ipr,fmt="(a)") "_cell_length_a    "//trim(text(1))
          write(unit=Ipr,fmt="(a)") "_cell_length_b    "//trim(text(2))
          write(unit=Ipr,fmt="(a)") "_cell_length_c    "//trim(text(3))
          write(unit=Ipr,fmt="(a)") "_cell_angle_alpha "//trim(text(4))
          write(unit=Ipr,fmt="(a)") "_cell_angle_beta  "//trim(text(5))
          write(unit=Ipr,fmt="(a)") "_cell_angle_gamma "//trim(text(6))
          write(unit=Ipr,fmt="(a)")
       else
          if(present(Cell)) then
             do i=1,3
               call setnum_std(Cell%Cell(i),Cell%cell_std(i),text(i))
               call setnum_std(Cell%ang(i),Cell%ang_std(i),text(i+3))
             end do
             write(unit=Ipr,fmt="(a)") "_cell_length_a    "//trim(text(1))
             write(unit=Ipr,fmt="(a)") "_cell_length_b    "//trim(text(2))
             write(unit=Ipr,fmt="(a)") "_cell_length_c    "//trim(text(3))
             write(unit=Ipr,fmt="(a)") "_cell_angle_alpha "//trim(text(4))
             write(unit=Ipr,fmt="(a)") "_cell_angle_beta  "//trim(text(5))
             write(unit=Ipr,fmt="(a)") "_cell_angle_gamma "//trim(text(6))
             write(unit=Ipr,fmt="(a)")
          end if
       end if
       if(MSGp%n_kv > 0) then
          write(unit=Ipr,fmt="(a)") "loop_"
          write(unit=Ipr,fmt="(a)") "_magnetic_propagation_vector_seq_id"
          write(unit=Ipr,fmt="(a)") "_magnetic_propagation_vector_kxkykz"
          do i=1,MSGp%n_kv
            call Frac_Trans_2Dig(MSGp%kv(:,i),line)
            line=adjustl(line(2:len_trim(line)-1))
            write(unit=Ipr,fmt="(a)") trim(MSGp%kv_label(i))//"  '"//trim(line)//"'"
          end do
       end if
       if(MSGp%m_constr) then
          write(unit=Ipr,fmt="(a)")
          write(unit=Ipr,fmt="(a)") "loop_"
          write(unit=Ipr,fmt="(a)") "_magnetic_atom_site_moment_symmetry_constraints_label"
          write(unit=Ipr,fmt="(a)") "_atom_site_magnetic_moment_symmetry_constraints_mxmymz"
          do i=1,Am%natoms
            line=Am%Atom(i)%AtmInfo
            if(len_trim(line) < 8) cycle
            write(unit=Ipr,fmt="(a)")trim(line)
          end do
       end if
       write(unit=Ipr,fmt="(a)")
       write(unit=Ipr,fmt="(a)")  "loop_"
       write(unit=Ipr,fmt="(a)")  "_space_group_magn_symop_operation.id"
       write(unit=Ipr,fmt="(a)")  "_space_group_magn_symop_operation.xyz"
       write(unit=Ipr,fmt="(a)")  "_space_group_magn_symop_operation.mxmymz"
       do i=1,MSGp%Multip            !New mCIF format
          write(unit=invc,fmt="(i2)") nint(MSgp%MSymop(i)%Phas)
          if(invc(1:1) == " ") invc(1:1)="+"
          write(unit=Ipr,fmt="(i3,a)") i," "//trim(MSgp%SymopSymb(i))//","//invc//" "//trim(MSgp%MSymopSymb(i))
       end do
       write(unit=Ipr,fmt="(a)")
       write(unit=Ipr,fmt="(a)") "loop_"
       write(unit=Ipr,fmt="(a)") "_atom_site_label"
       write(unit=Ipr,fmt="(a)") "_atom_site_type_symbol"
       write(unit=Ipr,fmt="(a)") "_atom_site_fract_x"
       write(unit=Ipr,fmt="(a)") "_atom_site_fract_y"
       write(unit=Ipr,fmt="(a)") "_atom_site_fract_z"
       write(unit=Ipr,fmt="(a)") "_atom_site_U_iso_or_equiv"
       write(unit=Ipr,fmt="(a)") "_atom_site_occupancy"
       write(unit=Ipr,fmt="(a)") "_atom_site_symmetry_multiplicity"
       write(unit=Ipr,fmt="(a)") "_atom_site_Wyckoff_label"
       line=" "
       do i=1,Am%natoms
          do j=1,3
            call setnum_std(Am%atom(i)%x(j),Am%atom(i)%x_std(j),text(j))
          end do
          occ=real(MSgp%Multip)/real(Am%atom(i)%Mult)*Am%atom(i)%occ
          occ_std=real(MSgp%Multip)/real(Am%atom(i)%Mult)*Am%atom(i)%occ_std
          call setnum_std(occ,occ_std,text(5))
          uiso=Am%atom(i)%biso/78.95683521
          uiso_std=Am%atom(i)%biso_std/78.95683521
          call setnum_std(uiso,uiso_std,text(4))
          write(unit=Ipr,fmt="(a6,a6,3a13,2a11,i4,a)") Am%Atom(i)%lab, Am%atom(i)%SfacSymb,(text(j),j=1,5),&
                                                       Am%atom(i)%Mult," "//Am%atom(i)%wyck
       end do
       write(unit=Ipr,fmt="(a)")
       write(unit=Ipr,fmt="(a)") "loop_"
       write(unit=Ipr,fmt="(a)") "_atom_site_moment_label"
       write(unit=Ipr,fmt="(a)") "_atom_site_moment_crystalaxis_x"
       write(unit=Ipr,fmt="(a)") "_atom_site_moment_crystalaxis_y"
       write(unit=Ipr,fmt="(a)") "_atom_site_moment_crystalaxis_z"
       do i=1,Am%natoms
          !if(sum(abs(Am%Atom(i)%Skr(:,1))) < 0.0001) cycle
          if(Am%Atom(i)%moment < 0.01) cycle
          do j=1,3
            call setnum_std(Am%atom(i)%M_xyz(j),Am%atom(i)%sM_xyz(j),text(j))
          end do
          write(unit=Ipr,fmt="(a8,3a12)") Am%Atom(i)%lab,(text(j),j=1,3)
       end do
       write(unit=Ipr,fmt="(a)")
       return
    End Subroutine Write_MCIF

 End Module CFML_IO_Formats

