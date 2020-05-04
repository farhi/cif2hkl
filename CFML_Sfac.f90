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
!!---- MODULE: CFML_Structure_Factors
!!----   INFO: Main module for Structure Factors Calculations
!!----
!!---- HISTORY
!!----    Update: 06/03/2011
!!----
!!----
!!---- DEPENDENCIES
!!----
!!--++     Use CFML_Scattering_Chemical_Tables
!!--++     Use CFML_Crystallographic_Symmetry,   only: Space_Group_Type,Magnetic_Space_Group_type
!!--++     Use CFML_Reflections_Utilities,       only: Reflection_List_Type, HKL_R
!!--++     Use CFML_Atom_TypeDef,                only: atom_list_type
!!--++     Use CFML_GlobalDeps,                  only: sp, tpi
!!--++     Use CFML_Math_General,                only: atan2d
!!--++     Use CFML_String_Utilities,            only: L_Case,U_Case
!!----
!!---- VARIABLES
!!--++    AF0                             [Private]
!!--++    AFP                             [Private]
!!--++    AFPP                            [Private]
!!--++    AJH                             [Private]
!!--++    BJH                             [Private]
!!----    ERR_SFAC
!!----    ERR_SFAC_MESS
!!--++    HR_TYPE                         [Private]
!!--++    HR                              [Private]
!!--++    HT                              [Private]
!!----    SCATTERING_SPECIES_TYPE
!!--++    SF_INITIALIZED                  [Private]
!!----    STRF_TYPE
!!----    STRF_LIST_TYPE
!!--++    TH                              [Private]
!!----
!!---- PUBLIC PROCEDURES
!!----    Functions:
!!--++       FJ                           [Private]
!!----
!!----    Subroutines:
!!----       ALLOCATE_SCATTERING_SPECIES
!!----       CALC_GENERAL_STRFACTOR
!!--++       CALC_TABLE_AB                [Private]
!!--++       CALC_TABLE_TH                [Private]
!!----       CALC_HKL_STRFACTOR
!!----       CALC_STRFACTOR
!!--++       CREATE_TABLE_AF0_ELECTRONS   [Private]
!!--++       CREATE_TABLE_AF0_XRAY        [Private]
!!--++       CREATE_TABLE_AFP_NEUTNUC     [Private]
!!--++       CREATE_TABLE_FABC_XRAY       [Private]
!!--++       CREATE_TABLE_HR_HT           [Private]
!!----       INIT_HKL_STRUCTURE_FACTORS
!!----       INIT_STRUCTURE_FACTORS
!!----       MODIFY_SF
!!--++       SET_FIXED_TABLES             [Private]
!!----       SET_FORM_FACTORS
!!----       STRUCTURE_FACTORS
!!--++       SUM_AB                       [Private]
!!--++       SUM_AB_NEUTNUC               [Private]
!!----       WRITE_STRUCTURE_FACTORS
!!----
!!
 Module CFML_Structure_Factors

    !---- Use Modules ----!
    Use CFML_IO_Formats,                  only: file_list_type
    Use CFML_GlobalDeps,                  only: cp,dp,tpi
    Use CFML_Math_General,                only: atan2d,cosd
    Use CFML_String_Utilities,            only: L_Case,U_Case
    Use CFML_Scattering_Chemical_Tables
    Use CFML_Crystal_Metrics,             only: Crystal_Cell_type
    Use CFML_Crystallographic_Symmetry,   only: Space_Group_Type,Magnetic_Space_Group_type
    Use CFML_Reflections_Utilities,       only: Reflection_List_Type, HKL_R, Reflect_Type,Reflect_List_Type,&
                                                hkl_gen_shub
    Use CFML_Atom_TypeDef,                only: atom_list_type

    !---- Variables ----!
    implicit none

    private

    !---- List of public functions ----!

    !---- List of public subroutines ----!
    public :: Init_Structure_Factors,Init_Calc_hkl_StrFactors, Structure_Factors,  &
              Modify_SF, Write_Structure_Factors,Calc_StrFactor, Calc_hkl_StrFactor, &
              Init_Calc_StrFactors, Allocate_Scattering_Species, Additional_Scattering_Factors, &
              Calc_General_StrFactor, Set_Form_Factors, Calc_Mag_Structure_Factor,Magnetic_Structure_Factors

    !---- List of private functions ----!
    private :: Fj

    !---- List of private subroutines ----!
    private :: Calc_Table_AB, Create_Table_AF0_Xray, Create_Table_AFP_NeutNuc, &
               Create_Table_HR_HT, Set_Fixed_Tables, Calc_Table_TH, Sum_AB,    &
               Sum_AB_NeutNuc, Create_Table_Fabc_Xray, Create_Table_AF0_Electrons, &
               Write_Structure_Factors_Crys,Write_Structure_Factors_Mag


    !---- Definitions ----!

    !!--++
    !!--++ AF0
    !!--++     real(kind=cp), dimension(:,:), allocatable, private :: AF0
    !!--++
    !!--++     Array for Atomic Factor. The dimensions are
    !!--++           AF0(Natoms,NRef)
    !!--++
    !!--++ Update: December - 2003
    !!
    real(kind=cp), dimension(:,:), allocatable, private :: AF0

    !!--++
    !!--++ AFP
    !!--++     real(kind=cp), dimension(:), allocatable, private :: AFP
    !!--++
    !!--++     Array for real part of anomalous scattering form factor.
    !!--++     The dimension is: AFP(Natoms)
    !!--++
    !!--++ Update: December - 2003
    !!
    real(kind=cp), dimension(:), allocatable, private :: AFP

    !!--++
    !!--++ AFPP
    !!--++     real(kind=cp), dimension(:), allocatable, private :: AFPP
    !!--++
    !!--++     Array for imaginary part of anomalous scattering form factor.
    !!--++     The dimension is: AFPP(Natoms)
    !!--++
    !!--++ Update: December - 2003
    !!
    real(kind=cp), dimension(:), allocatable, private :: AFPP

    !!--++
    !!--++ AJH
    !!--++     real(kind=cp), dimension(:,:), allocatable, private :: Ajh
    !!--++
    !!--++     Array for Aj(h). The dimensions are
    !!--++           Ajh(Natoms,Nref)
    !!--++     where
    !!--++           F(h)=Sum_j[Fj(h){Aj(h)+i Bj(h)}]
    !!--++
    !!--++ Update: December - 2003
    !!
    real(kind=cp), dimension(:,:), allocatable, private :: AJH

    !!--++
    !!--++ BJH
    !!--++     real(kind=cp), dimension(:,:), allocatable, private :: Bjh
    !!--++
    !!--++     Array for Bj(h). The dimensions are
    !!--++           Bjh(Natoms,Nref)
    !!--++     where
    !!--++           F(h)=Sum_j[Fj(h){Aj(h)+i Bj(h)}]
    !!--++
    !!--++ Update: December - 2003
    !!
    real(kind=cp), dimension(:,:), allocatable, private :: BJH

    !!----
    !!---- ERR_SFAC
    !!----    logical, public ::
    !!----
    !!----    Logical Variable in
    !!----
    !!---- Update: February - 200
    !!
    logical, public :: ERR_SFac

    !!----
    !!---- ERR_SFac_Mess
    !!----    character(len=150), public :: ERR_SFac_Mess
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: February - 2005
    !!
    character(len=150), public :: ERR_SFac_Mess


    !!--++
    !!--++ FF_A, FF_B, FF_C, FF_Z
    !!--++     real(kind=cp), dimension(:,:), allocatable, private :: FF_a,FF_b
    !!--++     real(kind=cp), dimension(  :), allocatable, private :: FF_c
    !!--++
    !!--++     Arrays for coefficients of X-rays scattering form factors.
    !!--++     The dimensions are: AFP(Nspecies)
    !!--++      FF_A(4,Nspecies), FF_B(4,Nspecies), FF_C(Nspecies), FF_Z(Nspecies)
    !!--++     Constructed in Create_Table_fabc_Xray(Atm,lambda,lun)
    !!--++     FF_Z contains atomic number Z (useful for electron diffraction)
    !!--++
    !!--++ Update: April - 2009
    !!
    real(kind=cp), dimension(:,:), allocatable, private :: FF_a, FF_b
    real(kind=cp), dimension(  :), allocatable, private :: FF_c
    real(kind=cp), dimension(  :), allocatable, private :: FF_Z

    !!--++
    !!--++    Type :: HR_Type
    !!--++       integer,dimension(3) :: H
    !!--++    End Type HR_Type
    !!--++
    !!--++    (Private)
    !!--++    Define a H vector
    !!--++
    !!--++ Update: February - 2005
    !!
    Type, Private :: HR_Type
       integer, dimension(3) :: H
    End Type HR_Type

    !!--++
    !!--++ HR
    !!--++     type(HR_Type), dimension(:,:), allocatable, private :: Hr
    !!--++
    !!--++     Array for HR Calculations. The dimension are
    !!--++           HR(Natoms,NRef)
    !!--++
    !!--++ Update: February - 2005
    !!
    type(HR_Type), dimension(:,:), allocatable, private :: HR

    !!----
    !!---- Type, public :: Scattering_Species_Type
    !!----    integer                                            :: Num_Species
    !!----    integer                                            :: Num_magspc !number of magnetic species
    !!----    character(len=6),        dimension(:), allocatable :: Symb,Symb_mag
    !!----    real(kind=cp),           dimension(:), allocatable :: br,bi
    !!----    real(kind=cp),           dimension(:), allocatable :: delta_fp,delta_fpp
    !!----    type(Xray_Form_Type),    dimension(:), allocatable :: Xcoef
    !!----    type(Magnetic_Form_Type),dimension(:), allocatable :: Mcoef
    !!---- End Type Scattering_Species_Type
    !!----
    !!---- Type encapsulting scattering factors for neutrons and X-rays. Constructed
    !!---- by calling Set_Form_Factors
    !!----
    Type, public :: Scattering_Species_Type
       integer                                            :: Num_Species
       integer                                            :: Num_magspc
       character(len=6),        dimension(:), allocatable :: Symb,Symb_mag
       real(kind=cp),           dimension(:), allocatable :: br,bi
       real(kind=cp),           dimension(:), allocatable :: delta_fp,delta_fpp
       type(Xray_Form_Type),    dimension(:), allocatable :: Xcoef
       type(Magnetic_Form_Type),dimension(:), allocatable :: Mcoef
    End Type Scattering_Species_Type

    !!----
    !!---- TYPE :: Strf_Type
    !!--..
    !!----    Type, Public  :: Strf_Type
    !!----       real(kind=cp)                    :: sqNuc       !Square of the nuclear structure factor
    !!----       real(kind=cp)                    :: sqMiV       !Square of the Magnetic Interaction vector
    !!----       complex(kind=cp)                 :: NsF         !Nuclear structure factor
    !!----       complex(kind=cp), dimension(3)   :: MsF         !Magnetic structure factor w.r.t. unitary Crystal Frame
    !!----       complex(kind=cp), dimension(3)   :: MiV         !Magnetic interaction vector w.r.t. unitary Crystal Frame
    !!----       complex(kind=cp), dimension(3)   :: MiVC        !Magnetic interaction vector in Cartesian components w.r.t. Crystal Frame
    !!----    End Type  Strf_Type
    !!----
    !!----
    !!---- Updated: January -2020
    !!
    Type, Public  :: Strf_Type
       real(kind=cp)                    :: sqNuc=0.0      !Square of the nuclear structure factor
       real(kind=cp)                    :: sqMiV=0.0      !Square of the Magnetic Interaction vector
       complex(kind=cp)                 :: NsF  =0.0      !Nuclear structure factor
       complex(kind=cp), dimension(3)   :: MsF  =0.0      !Magnetic structure factor w.r.t. unitary Crystal Frame
       complex(kind=cp), dimension(3)   :: MiV  =0.0      !Magnetic interaction vector w.r.t. unitary Crystal Frame
       complex(kind=cp), dimension(3)   :: MiVC =0.0      !Magnetic interaction vector in Cartesian components w.r.t. Crystal Frame
    End Type  Strf_Type

    !!----
    !!---- TYPE :: Strf_List_Type
    !!--..
    !!----    Type, Public  :: Strf_List_Type
    !!----       integer                                   :: Nref  !Number of reflections
    !!----       Type(Strf_Type), dimension(:),allocatable :: Strf
    !!----    End Type  Strf_Type
    !!----
    !!----
    !!---- Updated: January -2020
    !!
    Type, Public  :: Strf_List_Type
       integer                                   :: Nref
       Type(Strf_Type), dimension(:),allocatable :: Strf
    End Type  Strf_List_Type

    !!--++
    !!--++ HT
    !!--++    real(kind=cp), dimension(:,:), allocatable, private :: Ht
    !!--++
    !!--++    Array for HT Calculations. The dimension are
    !!--++          HT(Natoms,Nref)
    !!--++
    !!--++ Update: February - 2005
    !!
    real(kind=cp), dimension(:,:), allocatable, private :: HT

    !!--++
    !!--++ Nspecies
    !!--++     integer, private :: Nspecies
    !!--++
    !!--++     Number of chemical species for X-rays scattering form factors.
    !!--++  Constructed in Create_Table_fabc_Xray(Atm,lambda,lun)
    !!--++
    !!--++ Update: April - 2009
    !!
    integer, private :: Nspecies

    !!--++
    !!--++ P_A
    !!--++     integer, dimension(:), allocatable, private :: P_A
    !!--++
    !!--++     Integer pointer from atoms to species: P_A(Natoms), contains the species
    !!--++     of atom Natoms. Constructed in Create_Table_fabc_Xray(Atm,lambda,lun)
    !!--++
    !!--++ Update: April - 2009
    !!
    integer, dimension(:), allocatable, private :: P_A

    !!----
    !!---- SF_Initialized
    !!----    logical, private :: SF_Initialized
    !!----
    !!----  Logical Variable indicating if the module has been initialized.
    !!----
    !!---- Update: February - 2005
    !!
    logical, private :: SF_Initialized=.false.

    !!--++
    !!--++ TH
    !!--++    real(kind=cp), dimension(:,:), allocatable, private :: Th
    !!--++
    !!--++    Array for TH Calculations. The dimension are
    !!--++          TH(Natoms,Nref)
    !!--++
    !!--++ Update: February - 2005
    !!
    real(kind=cp), dimension(:,:), allocatable, private :: TH

    Interface Write_Structure_Factors
       Module Procedure Write_Structure_Factors_Crys
       Module Procedure Write_Structure_Factors_Mag
    End Interface Write_Structure_Factors

 Contains

    !---- Functions ----!

    !!--++
    !!--++ Pure Function Fj(s,a,b,c)
    !!--++    real(kind=cp),             intent(in) :: s
    !!--++    real(kind=cp),dimension(4),intent(in) :: a
    !!--++    real(kind=cp),dimension(4),intent(in) :: b
    !!--++    real(kind=cp),             intent(in) :: c
    !!--++
    !!--++    (Private)
    !!--++    Atomic scattering factor calculation according to:
    !!--++       Fj(s)=Sum_i[Ai*exp(-Bi*s*s)] + C (i=1..4)
    !!--++
    !!--++ Update: February - 2005
    !!
    Pure Function Fj(s,a,b,c) Result(res)
       !---- Arguments ----!
       real(kind=cp),             intent(in) :: s
       real(kind=cp),dimension(4),intent(in) :: a
       real(kind=cp),dimension(4),intent(in) :: b
       real(kind=cp),             intent(in) :: c
       real(kind=cp)                         :: res

       !---- Local variables ----!
       integer :: i

       res=0.0
       do i=1,4
          res=res + a(i)*exp(-b(i)*s*s)
       end do
       res=res+c

       return
    End Function Fj

    !---- Subroutines ----!

    !!----
    !!----  Subroutine Additional_Scattering_Factors(fil,add_Scatt,ok,mess)
    !!----    Type(File_List_Type),          intent(in)  :: fil
    !!----    Type(Scattering_Species_Type), intent(out) :: add_Scatt
    !!----    logical,                       intent(out) :: ok
    !!----    character(len=*),              intent(out) :: mess
    !!----
    !!----  Subroutine constructing add_Scatt by reading fil of File_List_Type
    !!----  These values replace those read from database CFML_Scattering_Chemical_Tables
    !!----
    !!----  Created (JRC): October-2015
    !!----
   Subroutine Additional_Scattering_Factors(fil,add_Scatt,ok,mess)
      Type(File_List_Type),          intent(in)  :: fil
      Type(Scattering_Species_Type), intent(out) :: add_Scatt
      logical,                       intent(out) :: ok
      character(len=*),              intent(out) :: mess
      !Local variables
      integer, parameter :: N_add = 20
      character(len=132) :: line
      character(len=4), dimension(N_add)   :: names
      real(kind=cp),    dimension(N_add)   :: b_real,b_imag
      real(kind=cp),    dimension(N_add)   :: d_fp,d_fpp,cc
      real(kind=cp),    dimension(4,N_add) :: ac,bc
      integer :: i,nsp,j,ier
      ok=.true.
      mess=" "
      b_real=0.0; b_imag=0.0; d_fp=0.0; d_fpp=0.0; cc=0.0
      ac=0.0; bc=0.0
      nsp=0
      do i=1,fil%nlines
        line=adjustl(fil%line(i))
        if(U_case(line(1:2)) == "B_") then
          nsp=nsp+1
          j=index(line," ")
          names(nsp)=line(3:j-1)
          read(unit=line(j:),fmt=*,iostat=ier) b_real(nsp),b_imag(nsp)
          if(ier /= 0) then
            ok=.false.
            mess="Error reading scattering length on line containing: "//trim(line)
            return
          end if
        end if
        if(U_case(line(1:5)) == "DELT_") then
          nsp=nsp+1
          j=index(line," ")
          names(nsp)=line(6:j-1)
          read(unit=line(j:),fmt=*,iostat=ier) d_fp(nsp),d_fpp(nsp)
          if(ier /= 0) then
            ok=.false.
            mess="Error reading anomalous scattering terms on line containing: "//trim(line)
            return
          end if
        end if
        if(U_case(line(1:7)) == "XCOEFF_") then
          nsp=nsp+1
          j=index(line," ")
          names(nsp)=line(8:j-1)
          read(unit=line(j:),fmt=*,iostat=ier) ac(:,nsp),bc(:,nsp),cc(nsp)
          if(ier /= 0) then
            ok=.false.
            mess="Error reading X-ray scattering coefficients on line containing: "//trim(line)
            return
          end if
        end if
      end do
      if(nsp > N_add) then
        ok=.false.
        write(unit=mess,fmt="(a,i3,a)") "The number of additional scattering factors is limited to ",N_add," !!!"
        return
      end if
      if(nsp > 0) then
        call Allocate_Scattering_Species(nsp,add_Scatt)
        do i=1,nsp
          add_Scatt%Symb(i)      = names(i)
          add_Scatt%br(i)        = b_real(i)
          add_Scatt%bi(i)        = b_imag(i)
          add_Scatt%delta_fp(i)  = d_fp(i)
          add_Scatt%delta_fpp(i) = d_fpp(i)
          add_Scatt%Xcoef(i)%Symb= names(i)
          add_Scatt%Xcoef(i)%a   = ac(:,i)
          add_Scatt%Xcoef(i)%b   = bc(:,i)
          add_Scatt%Xcoef(i)%c   = cc(i)
        end do
      else
        add_Scatt%Num_species=0
      end if

    End Subroutine Additional_Scattering_Factors


    !!----
    !!---- Subroutine Allocate_Scattering_Species(n,Scf)
    !!----   integer,                       intent(in)  :: n
    !!----   type(Scattering_Species_Type), intent(out) :: Scf
    !!----
    !!----  Allocates the components of Scf
    !!----
    !!----  Created (JRC): October-2015
    !!----
    Subroutine Allocate_Scattering_Species(n,Scf)
      integer,                       intent(in)  :: n
      type(Scattering_Species_Type), intent(out) :: Scf
      !--- Local variables ---!
      integer :: i

      Scf%Num_Species=n; Scf%Num_magspc=0
      allocate(Scf%br(n),Scf%bi(n),Scf%delta_fp(n),Scf%delta_fpp(n),Scf%symb(n))
      Scf%br=0.0; Scf%bi=0.0; Scf%delta_fp=0.0; Scf%delta_fpp=0.0; Scf%symb= " "
      allocate(Scf%Xcoef(n))

      do i=1,Scf%Num_Species
        Scf%Xcoef(i)%Z=0
        Scf%Xcoef(i)%a=0.0
        Scf%Xcoef(i)%b=0.0
        Scf%Xcoef(i)%c=0.0
      end do

    End Subroutine Allocate_Scattering_Species

    !!----
    !!---- Subroutine Calc_General_StrFactor(hn,sn,Atm,Grp,Scf,fn,fx,fe)
    !!----    !---- Arguments ----!
    !!----    real(kind=cp),dimension(3),         intent(in) :: hn
    !!----    real(kind=cp),                      intent(in) :: sn !(sinTheta/Lambda)**2
    !!----    type(atom_list_type),               intent(in) :: Atm
    !!----    type(space_group_type),             intent(in) :: Grp
    !!----    type(Scattering_Species_Type),      intent(in) :: Scf
    !!----    complex,                            intent(out):: fn,fx,fe
    !!----
    !!----   Calculates nuclear, x-rays and electrostatic structure factors from
    !!----   the list of atoms Atm, space group Grp and scattering species Scf.
    !!----
    !!----  Created (JRC): October-2015
    !!----
    Subroutine Calc_General_StrFactor(hn,sn,Atm,Grp,Scf,fn,fx,fe)
       !---- Arguments ----!
       real(kind=cp),dimension(3),         intent(in) :: hn
       real(kind=cp),                      intent(in) :: sn !(sinTheta/Lambda)**2
       type(atom_list_type),               intent(in) :: Atm
       type(space_group_type),             intent(in) :: Grp
       type(Scattering_Species_Type),      intent(in) :: Scf
       complex,                            intent(out):: fn,fx,fe

       !---- Local Variables ----!
       integer                               :: i,j,k !,m
       real(kind=cp)                         :: arg,anis,scosr,ssinr,b !,s
       real(kind=cp)                         :: a1,a3,b1,b3,av,bv,nffr,nffi        !fn
       real(kind=cp)                         :: xa1,xa3,xb1,xb3,xav,xbv,xffr,xffi  !fx
       real(kind=cp)                         :: ea1,ea3,eb1,eb3,effr !,eav,ebv     !fe
       real(kind=cp),dimension(3)            :: h
       real(kind=cp),dimension(6)            :: beta

       !--- Initialising local variables
       a1 =0.0;  a3 =0.0
       b1 =0.0;  b3 =0.0
       ea1=0.0;  ea3=0.0
       eb1=0.0;  eb3=0.0
       xa1=0.0;  xa3=0.0
       xb1=0.0;  xb3=0.0

       do i=1,Atm%natoms
          arg=0.0
          scosr=0.0
          ssinr=0.0
          do k=1,grp%NumOps
             h=matmul(hn,grp%Symop(k)%Rot)
             arg=tpi*(dot_product(h,Atm%atom(i)%x)+ dot_product(hn,grp%Symop(k)%tr))
             anis=1.0
             if(Atm%atom(i)%thtype == "aniso") then
               beta=Atm%atom(i)%u(1:6)
               anis=     h(1)*h(1)*beta(1)+     h(2)*h(2)*beta(2)+    h(3)*h(3)*beta(3) &
                    +2.0*h(1)*h(2)*beta(4)+ 2.0*h(1)*h(3)*beta(5)+2.0*h(2)*h(3)*beta(6)
               anis=exp(-anis)
             end if
             scosr=scosr+COS(arg)*anis  !Real part of geometrical structure factor for the current atom
             ssinr=ssinr+SIN(arg)*anis  !Imaginary part of geometrical structure factor for the current atom
          end do ! symmetry

          b= atm%atom(i)%occ * exp(-atm%atom(i)%biso*sn)
          !Calculation of scattering factors
          j=atm%atom(i)%ind(1)  !pointer to the form factor coefficients

          nffr = Scf%br(j)*b
          nffi = Scf%bi(j)*b

          xffr=Scf%Xcoef(j)%c
          do k=1,4
            xffr=xffr+Scf%Xcoef(j)%a(k)*exp(-Scf%Xcoef(j)%b(k)*sn)
          end do

          effr = (real(Scf%Xcoef(j)%Z)-xffr)*b  !<- Here delta_fp is not used ....
          !write(*,"(tr4,a,i4,3f10.4)")  Atm%atom(i)%chemsymb//": Z, xffr, b, effr", Scf%Xcoef(j)%Z, xffr, b, effr
          xffr = (xffr+Scf%delta_fp(j))*b       ! (f0+Deltaf')*OCC*Tiso
          xffi = Scf%delta_fpp(j)*b             !     Deltaf" *OCC*Tiso

          a1 = a1 + nffr*scosr  ! F=A+iB: components of A  and B (ai,bi)
          b1 = b1 + nffi*scosr  ! a2,b2,a4,b4 are components for anisotropic form factors
          a3 = a3 + nffi*ssinr  ! they are not used here
          b3 = b3 + nffr*ssinr  ! For general case: av = a1-a2-a3-a4, bv = b1-b2+b3+b4

          xa1 = xa1 + xffr*scosr  ! F=A+iB: components of A  and B (ai,bi)
          xb1 = xb1 + xffi*scosr  ! a2,b2,a4,b4 are components for anisotropic form factors
          xa3 = xa3 + xffi*ssinr  ! they are not used here
          xb3 = xb3 + xffr*ssinr  ! For general case: av = a1-a2-a3-a4, bv = b1-b2+b3+b4

          ea1 = ea1 + effr*scosr  ! No anomalous imaginary component is used here
          eb3 = eb3 + effr*ssinr  ! there is no anomalous scattering neutron + electrons
          !write(*,"(tr4,a,tr4,3f10.4)") Atm%atom(i)%chemsymb//": effr,   ea1, eb3",  effr, ea1, eb3

       end do ! Atoms

       av = a1-a3    !real part of the Nuclear structure factor
       bv = b1+b3    !imaginary part of the Nuclear structure factor
       fn=cmplx(av,bv) * Grp%Centred * Grp%NumLat

       xav = xa1-xa3    !real part of the X-rays structure factor
       xbv = xb1+xb3    !imaginary part of the X-rays structure factor
       fx=cmplx(xav,xbv) * Grp%Centred * Grp%NumLat

       fe=cmplx(ea1,eb3) * Grp%Centred * Grp%NumLat

    End Subroutine Calc_General_StrFactor

    !!----
    !!---- Subroutine Calc_Mag_Structure_Factor(mode,hm,Cell,Grp,Atm,Scf,Strf,magonly,mdom,tdom,twin)
    !!----   character(len=*),                   intent(in) :: mode !S-XTAL (S) or Powder (P)
    !!----   type(Reflect_Type),                 intent(in) :: hn     !Contains hkl,s,mult and imag within Reflect_List_Type
    !!----   type(Crystal_Cell_type),            intent(in) :: Cell
    !!----   type(Matom_list_type),              intent(in) :: Atm
    !!----   type(Magnetic_Space_Group_type),    intent(in) :: Grp
    !!----   Type(Scattering_Species_Type),      intent(in) :: Scf !species
    !!----   type(Strf_Type),                    intent(out):: Strf
    !!----   Integer,Dimension(3,3),    Optional,Intent(In) :: mdom !Matrix to be applied to all Fourier coefficients
    !!----   Real(Kind=Cp),Dimension(3),Optional,Intent(In) :: tdom !Translation to be applied to all atom positions together with mdom
    !!----   Character(Len=*),          Optional,Intent(In) :: twin !Representing a particular orientation domain
    !!----                                                          !Useful only for single crystals
    !!---
    !!----  Calculation of nuclear and magnetic structure factor, when
    !!----  the symmetry is given by a Magnetic Space Group, for reflection h
    !!----  Structure factor components are provided in Stf
    !!----
    !!----  Updated: January 2020
    !!----
    Subroutine Calc_Mag_Structure_Factor(mode,hm,Cell,Grp,Atm,Scf,Strf,magonly,mdom,tdom,twin)
      character(len=*),                   intent(in)  :: mode !S-XTAL (S) or Powder (P)
      type(Reflect_Type),                 intent(in)  :: hm
      type(Crystal_Cell_type),            intent(in)  :: Cell
      type(atom_list_type),               intent(in)  :: Atm
      type(Magnetic_Space_Group_type),    intent(in)  :: Grp
      Type(Scattering_Species_Type),      intent(in)  :: Scf !species
      type(Strf_Type),                    intent(out) :: Strf
      logical,                   optional,intent(in)  :: magonly
      Integer,Dimension(3,3),    optional,intent(In)  :: mdom !Matrix to be applied to all Fourier coefficients
      Real(Kind=Cp),Dimension(3),optional,intent(In)  :: tdom !Translation to be applied to all atom positions together with mdom
      Character(Len=*),          optional,intent(In)  :: twin !Representing a particular orientation domain
                                                             !Useful only for single crystals
      !-----------------------------------------------
      !   L o c a l   V a r i a b l e s
      !-----------------------------------------------
      integer :: i,ni, ii, ir
      real(kind=cp), dimension(Atm%natoms)     :: otr, oti, frc, frs
      real(kind=cp), dimension(3)              :: h, hnn, xi, cosa, side, aa, bb, t, ar, br,ed,ec
      real(kind=cp), dimension(6)              :: betas
      Real(Kind=Cp), Dimension(3,3)            :: sm,SMcos,SMsin
      complex(kind=cp),dimension(3)            :: Mc
      logical                                  :: mag,nuc,mag_only
      character(len=1)                         :: tw
      real(kind=cp) :: ffr, ffi, ffx, cosr, sinr, scosr, ssinr, temp,snexi !,x1, yy, z
      real(kind=cp) :: x, arg, arg2, exparg,ssnn
      real(kind=dp) :: a1, a3, b1, b3, av,bv
      real(kind=dp), parameter  :: pn=0.2695420113693928312
      tw="N"
      if(present(mdom)) then
         if(present(twin)) tw=twin(1:1)
      end if
      if(tw == "T") then
        hnn=matmul(hm%h,real(mdom))
      else
        hnn=hm%h
      end if

      mag_only=.false.
      a1=0.0; a3=0.0
      b1=0.0; b3=0.0
      aa=0.0; bb=0.0  !aa(1:3)
      cosa=cosd(Cell%Ang)
      side=Cell%Cell
      ssnn=hm%s*hm%s

      nuc= hm%imag == 0  .or. hm%imag == 2

      if(present(magonly)) mag_only=magonly
      mag= hm%imag /= 0
      if(mag_only) nuc=.false.

      !---------------------------
      DO i=1,Atm%natoms           !loop over atoms
      !---------------------------        ipiof=i+iof
        xi=Atm%atom(i)%x
        betas=Atm%atom(i)%U
        !Modify the first atom position according to the interpretation of domains with translations
        if(present(tdom)) xi(1:3) = matmul(real(mdom),xi(1:3))+tdom(1:3)
        temp=EXP(-Atm%atom(i)%Biso*ssnn)   !exp{-Bi (sintheta/Lambda)^2}

        !Nuclear form factor
        ffi=0.0; ffr=0.0
        if(nuc) then
           ni=Atm%atom(i)%ind(1)
           ffr=Scf%br(ni)
           ffi=0.0
        end if
        ffx=0.0

        IF(mag) Then
          ni=Atm%atom(i)%ind(2)
          ffx=Scf%Mcoef(ni)%SctM(7)
          DO ii=1,5,2
            ffx=ffx+Scf%Mcoef(ni)%SctM(ii)*EXP(-Scf%Mcoef(ni)%SctM(ii+1)*ssnn)  !Form factor value for Q=H+k
          END DO
          snexi=pn*ffx*temp*Atm%atom(i)%occ     ! 0.26954.f(Q).Temp(i).Occ
        END IF
        ! Loop over symmetry operators
        scosr=0.0
        ssinr=0.0
        ! Nuclear and Magnetic Structure Factor calculations
        !   Fm=(Ax,Ay,Az)+i(Bx,By,Bz) --->  aa(I),bb(I), I=1,2,3 for x,y,z
        SMcos=0.0; SMsin=0.0
       !+++++++++++++++++++++++++
        DO  ir=1,Grp%Multip   !Loop over symmetry operators
       !+++++++++++++++++++++++++
          sm(:,:)=Grp%Symop(ir)%Rot(:,:)
             t(:)=Grp%Symop(ir)%tr(:)
          x=dot_product(t,hnn)
          h=matmul(hnn,sm)
          arg=x+dot_product(h,xi(1:3))
          arg=tpi*arg
          arg2=    h(1)*h(1)*betas(1)+     h(2)*h(2)*betas(2)+    h(3)*h(3)*betas(3)+        &
               2.0*h(1)*h(2)*betas(4)+ 2.0*h(1)*h(3)*betas(5)+2.0*h(2)*h(3)*betas(6)
          exparg=EXP(-arg2)
          cosr=COS(arg)*exparg      !cos{2pi(hT Rs rj+ts)}*exp(-{hTRsBetaj RsTh})
          scosr=scosr+cosr          !FRC= SIG fr(j,s)cos{2pi(hT Rs rj+ts)}*Ta(s)
          sinr=SIN(arg)*exparg      !sin{2pi(hT Rs rj+ts)}*exp(-{hTRsBetaj RsTh})

          IF(Grp%centred == 1) then
           ssinr=ssinr+sinr          !FRS= SIG fr(j,s)sin{2pi(hT Rs rj+ts)}*Ta(s)
          END IF
          IF(mag) Then
            SMcos(:,:)=SMcos(:,:)+cosr*Grp%MSymop(ir)%Rot(:,:)
            SMsin(:,:)=SMsin(:,:)+sinr*Grp%MSymop(ir)%Rot(:,:)
          END IF
        !+++++++++++++++++++++++++
        END DO     !over symmetry operators !  END LOOP SYMM.OP.
        !+++++++++++++++++++++++++

        if(nuc) then
           frc(i)=scosr    !Components of geometrical struture factor of atom i
           frs(i)=ssinr
           otr(i)=ffr*Atm%atom(i)%occ*temp     ! (f0+Deltaf')*OCC*Tiso
           oti(i)=ffi*Atm%atom(i)%occ*temp     !     Deltaf" *OCC*Tiso

           !-----CALCULATE A AND B OF F
           a1 = a1 + otr(i)*frc(i)    ! components of A  and B
           b1 = b1 + oti(i)*frc(i)    ! A(h) = a1 - a3
                                      ! B(h) = b1 + b3
           IF(Grp%centred == 1) THEN
             a3 = a3 + oti(i)*frs(i)
             b3 = b3 + otr(i)*frs(i)
           END IF
        end if
        !Magnetic structure factor components
        IF(mag) Then
           ar = matmul(SMcos,Atm%atom(i)%M_xyz(:)/side)*side  !The introduction of Cell%cell
           br = matmul(SMsin,Atm%atom(i)%M_xyz(:)/side)*side  !of using non conventional settings for
           aa(:)= aa(:) + snexi*ar(:)
           bb(:)= bb(:) + snexi*br(:)
        END IF
      !---------------------------
      END DO  !over atoms

      if( nuc ) then
        !---------------------------
        !NUCLEAR STRUCTURE FACTOR
        !=========================
        av = a1-a3   !real part of the structure factor
        bv = b1+b3   !imaginary part of the structure factor
        Strf%NsF=cmplx(av,bv)
        ! For a powder h and -h cannot be measured independently, both kind
        ! of reflections contribute simultaneously to a peak, so the intensity
        ! is proportional to F^2(h)+ F^2(-h), the binary terms of the form:
        ! Fij(h)= 2.0 [a(i)*a(j) -b(i)*b(j)] = -Fij(-h)
        ! If the reflection -h is not generated, the calculation must be
        ! performed using only the diagonal terms.
        !     FNN = av*av + bv*bv !For a single crystal

        if(mode=="S") then
          Strf%sqNuc = av*av+bv*bv
        else
          Strf%sqNuc = a1*a1 + a3*a3 + b1*b1 + b3*b3
        end if
      else
        Strf%sqNuc = 0.0
        Strf%NsF=cmplx(0.0,0.0)
      end if

      if(mag) then
        !MAGNETIC STRUCTURE FACTOR
        !=========================
        !  Calcul of F2= Fm.Fm* - (e.Fm)*(e.Fm)
        !x=0.0
        !yy=0.0
        !DO i=1,3
        !  DO j=i,3
        !    k=6-i-j
        !    x1=1.0
        !    IF(i /= j) x1=2.0*cosa(k)
        !    x=x1*aa(i)*aa(j)+x   ! Fm.Fm with unit vectors along a,b,c metric tensor
        !    x1=2.0
        !    IF(i == j) x1=1.0
        !    yy=x1*hnn(i)*hnn(j)*aa(i)*aa(j)/(side(i)*side(j))+yy  !(e.Fm).(e.Fm)
        !    IF(icent == 1) THEN
        !      x1=1.0
        !      IF(i /= j) x1=2.0*cosa(k)
        !      x=x1*bb(i)*bb(j)+x
        !      x1=2.0
        !      IF(i == j) x1=1.0
        !      yy=x1*hnn(i)*hnn(j)*bb(i)*bb(j)/(side(i)*side(j))+yy
        !    END IF
        !  END DO
        !END DO
        !Strf%sqMiV=(x-0.25*yy/ssnn)!*nlatti(iph)      ! Halpern & Johnson F2= Fm.Fm* - (e.Fm)*(e.Fm)

        Strf%MsF=cmplx(aa,bb) !MsF in basis {e1,e2,e3}
        ed = matmul(cell%GR, 0.5*hm%h/hm%s)    ! unitary vector referred to the direct    "
        ec = matmul(Cell%Cr_Orth_cel,ed)       ! Cartesian
        Mc = Strf%MsF / Cell%cell              ! Magnetic structure factor in basis {a,b,c}
        Mc = matmul(Cell%Cr_Orth_cel,Mc)       !                            Cartesian
        Strf%MiVC = Mc - dot_product(ec,Mc) * ec      !Magnetic interaction vector in Cartesian components
        Strf%MiV = matmul(Cell%Orth_Cr_cel,Strf%MiVC)* Cell%cell       !Magnetic interaction vector in basis  {e1,e2,e3}
        Strf%sqMiV= dot_product(Strf%MiVC, Strf%MiVC)
      else
        Strf%sqMiV=0.0      ! Halpern & Johnson F2= Fm.Fm* - (e.Fm)*(e.Fm)
        Strf%MsF=cmplx(0.0,0.0)
        Strf%MiV=cmplx(0.0,0.0)
        Strf%MiVC=cmplx(0.0,0.0)
      end if

      Return
    End Subroutine Calc_Mag_Structure_Factor

    !!----
    !!---- Subroutine Calc_StrFactor(mode,rad,nn,sn,Atm,Grp,sf2,deriv,fc)
    !!----    character(len=*),                   intent(in) :: mode !S-XTAL (S) or Powder (P)
    !!----    character(len=*),                   intent(in) :: rad  !Radiation: X-rays, Neutrons
    !!----    integer,                            intent(in) :: nn
    !!----    real(kind=cp)                       intent(in) :: sn !(sinTheta/Lambda)**2
    !!----    type(atom_list_type),               intent(in) :: Atm
    !!----    type(space_group_type),             intent(in) :: Grp
    !!----    real(kind=cp)                       intent(out):: sf2
    !!----    real(kind=cp),dimension(:),optional,intent(out):: deriv
    !!----    complex, optional,                  intent(out):: fc
    !!----
    !!----    Calculate Structure Factor for reflection "nn" in the list
    !!----    and derivatives with respect to refined parameters
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Calc_StrFactor(mode,rad,nn,sn,Atm,Grp,sf2,deriv,fc)
       !---- Arguments ----!
       character(len=*),                   intent(in) :: mode
       character(len=*),                   intent(in) :: rad
       integer,                            intent(in) :: nn
       real(kind=cp),                      intent(in) :: sn !(sinTheta/Lambda)**2
       type(atom_list_type),               intent(in) :: Atm
       type(space_group_type),             intent(in) :: Grp
       real(kind=cp),                      intent(out):: sf2
       real(kind=cp),dimension(:),optional,intent(out):: deriv
       complex, optional,                  intent(out):: fc

       !---- Local Variables ----!
       character(len=1)                      :: modi
       integer                               :: i,j,k,m
       real(kind=cp)                         :: arg,anis,cosr,sinr,scosr,ssinr,fr,der !,fi
       real(kind=cp)                         :: a1,a2,a3,a4,b1,b2,b3,b4,av,bv,f
       real(kind=cp),dimension(3)            :: h
       real(kind=cp),dimension(6)            :: beta
       real(kind=cp),dimension(Atm%natoms)   :: frc,frs,otr,oti,afpxn
       real(kind=cp),dimension(9,Atm%natoms) :: drs,drc

       !--- Initialising local variables
       a1=0.0
       a2=0.0
       a3=0.0
       a4=0.0
       b1=0.0
       b2=0.0
       b3=0.0
       b4=0.0
       av=0.0
       bv=0.0
       fr=1.0
       !fi=0.0
       frc=0.0
       frs=0.0
       otr=0.0
       oti=0.0
       modi=u_case(mode(1:1))
       if(rad(1:1) == "N") then
         afpxn(:)=afp(:)
       else
         afpxn(:)=af0(:,nn)
       end if

       if(Grp%Centred == 2) then
            do i=1,Atm%natoms
               arg=0.0
               scosr=0.0
               ssinr=0.0
               drs(:,i)=0.0
               drc(:,i)=0.0
               do k=1,grp%NumOps
                  h=hr(k,nn)%h
                  arg=tpi*(dot_product(h,Atm%atom(i)%x)+ht(k,nn))
                  anis=1.0
                  if(Atm%atom(i)%thtype == "aniso") then
                    beta=Atm%atom(i)%u(1:6)
                    anis=     h(1)*h(1)*beta(1)+     h(2)*h(2)*beta(2)+    h(3)*h(3)*beta(3) &
                         +2.0*h(1)*h(2)*beta(4)+ 2.0*h(1)*h(3)*beta(5)+2.0*h(2)*h(3)*beta(6)
                    anis=exp(-anis)
                  end if
                  cosr=COS(arg)*anis*fr     !fr*cos{2pi(hT Rs rj+ts)}*exp(-{hTRsBetaj RsTh})
                  scosr=scosr+cosr          !FRC= SIG fr(j,s)cos{2pi(hT Rs rj+ts)}*Ta(s)

                  if(present(deriv)) then
                     sinr=SIN(arg)*anis*fr   !fr*sin{2pi(hT Rs rj+ts)}*exp(-{hTRsBetaj RsTh})
                     drc(1:3,i)=drc(1:3,i)+h(1:3)*sinr      ! -
                     drc(4,i)=drc(4,i)+h(1)*h(1)*cosr
                     drc(5,i)=drc(5,i)+h(2)*h(2)*cosr
                     drc(6,i)=drc(6,i)+h(3)*h(3)*cosr
                     drc(7,i)=drc(7,i)+h(1)*h(2)*cosr
                     drc(8,i)=drc(8,i)+h(1)*h(3)*cosr
                     drc(9,i)=drc(9,i)+h(2)*h(3)*cosr
                  end if

               end do ! symmetry

               frc(i) = scosr
               otr(i) = afpxn(i)*th(i,nn)
               oti(i) =  afpp(i)*th(i,nn)
               a1= a1 + otr(i)*frc(i)
               b1= b1 + oti(i)*frc(i)

            end do ! Atoms

            av = a1-a2-a3-a4    !real part of the structure factor
            bv = b1-b2+b3+b4    !imaginary part of the structure factor

       else

            do i=1,Atm%natoms
               arg=0.0
               scosr=0.0
               ssinr=0.0
               drs(:,i)=0.0
               drc(:,i)=0.0
               do k=1,grp%NumOps
                  h=hr(k,nn)%h
                  arg=tpi*(dot_product(h,Atm%atom(i)%x)+ht(k,nn))
                  anis=1.0
                  if(Atm%atom(i)%thtype == "aniso") then
                    beta=Atm%atom(i)%u(1:6)
                    anis=     h(1)*h(1)*beta(1)+     h(2)*h(2)*beta(2)+    h(3)*h(3)*beta(3) &
                         +2.0*h(1)*h(2)*beta(4)+ 2.0*h(1)*h(3)*beta(5)+2.0*h(2)*h(3)*beta(6)
                    anis=exp(-anis)
                  end if
                  cosr=COS(arg)*anis*fr     !fr*cos{2pi(hT Rs rj+ts)}*exp(-{hTRsBetaj RsTh})
                  sinr=SIN(arg)*anis*fr     !fr*sin{2pi(hT Rs rj+ts)}*exp(-{hTRsBetaj RsTh})
                  scosr=scosr+cosr          !FRC= SIG fr(j,s)cos{2pi(hT Rs rj+ts)}*Ta(s)
                  ssinr=ssinr+sinr          !FRS= SIG fr(j,s)sin{2pi(hT Rs rj+ts)}*Ta(s)

                  if(present(deriv)) then
                     drc(1:3,i)=drc(1:3,i)+h(1:3)*sinr      ! -
                     drs(1:3,i)=drs(1:3,i)+h(1:3)*cosr      ! +

                     drc(4,i)=drc(4,i)+h(1)*h(1)*cosr
                     drc(5,i)=drc(5,i)+h(2)*h(2)*cosr
                     drc(6,i)=drc(6,i)+h(3)*h(3)*cosr
                     drc(7,i)=drc(7,i)+h(1)*h(2)*cosr
                     drc(8,i)=drc(8,i)+h(1)*h(3)*cosr
                     drc(9,i)=drc(9,i)+h(2)*h(3)*cosr

                     drs(4,i)=drs(4,i)+h(1)*h(1)*sinr
                     drs(5,i)=drs(5,i)+h(2)*h(2)*sinr
                     drs(6,i)=drs(6,i)+h(3)*h(3)*sinr
                     drs(7,i)=drs(7,i)+h(1)*h(2)*sinr
                     drs(8,i)=drs(8,i)+h(1)*h(3)*sinr
                     drs(9,i)=drs(9,i)+h(2)*h(3)*sinr
                  end if

               end do ! symmetry

               frc(i) = scosr
               frs(i) = ssinr
               otr(i) = afpxn(i)*th(i,nn)
               oti(i) =  afpp(i)*th(i,nn)
               a1= a1 + otr(i)*frc(i)
               b1= b1 + oti(i)*frc(i)
               a3 = a3 + oti(i)*frs(i)
               b3 = b3 + otr(i)*frs(i)

            end do ! Atoms

            av = a1-a2-a3-a4    !real part of the structure factor
            bv = b1-b2+b3+b4    !imaginary part of the structure factor

       end if

       If(modi == "P") then
          sf2 = a1*a1 + a2*a2 + a3*a3 + a4*a4 + b1*b1 + b2*b2 + b3*b3 + b4*b4
          sf2 = sf2 + 2.0*(b1*b4 -  a1*a4 + a2*a3 - b2*b3)
       else
          sf2= av*av+bv*bv
       End if

       if(present(fc)) then
         fc=cmplx(av,bv)
       end if

       if(present(deriv)) then

         if(modi == "P") then

             do i=1,Atm%natoms
                !derivatives with respect to coordinates  POWDER
                do m=1,3
                   k= Atm%atom(i)%lx(m)
                   if(k /= 0) then
                     f=atm%atom(i)%mx(m)
                     der= otr(i)*(-a1*drc(m,i)+b3*drs(m,i))+oti(i)*(-b1*drc(m,i)+a3*drs(m,i))
                     der=2.0*der*tpi
                     deriv(k) = sign(1.0_cp,f)*der+deriv(k)
                   end if
                 end do

                 k=Atm%atom(i)%lbiso  !Derivatives w.r.t. Biso  POWDER
                 if(k /= 0) then
                   f=Atm%atom(i)%mbiso
                   der= otr(i)*(a1*frc(i) +b3*frs(i))+oti(i)*(b1*frc(i) +a3*frs(i))
                   der=-2.0*der*sn
                   deriv(k) = sign(1.0_cp,f)*der+deriv(k)
                 end if

                 k=Atm%atom(i)%locc    !Derivatives w.r.t. occupation factor   POWDER
                 if(k /= 0) then
                   f=Atm%atom(i)%mocc
                   der= otr(i)*(a1*frc(i)+b3*frs(i))+oti(i)*(b1*frc(i)+a3*frs(i))
                   der=2.0*der/atm%atom(i)%occ
                   deriv(k) = sign(1.0_cp,f)*der+deriv(k)
                 end if

                 do m=4,9      !Derivatives w.r.t. anisotropic temperature factors   POWDER
                    j=m-3
                    k=Atm%atom(i)%lu(j)
                    if(k /= 0) then
                      f=Atm%atom(i)%mu(j)
                      der=  otr(i)*(a1*drc(i,j)+b3*drs(m,i))+oti(i)*(b1*drc(m,i)+a3*drs(m,i))
                      der=-2.0*der
                      if(j > 3) der=2.0*der
                      deriv(k) = sign(1.0_cp,f)*der+deriv(k)
                    end if
                 end do

             end do

         else

             do i=1,Atm%natoms
                !derivatives with respect to coordinates  S-XTAL
                do m=1,3
                   k= Atm%atom(i)%lx(m)
                   if(k /= 0) then
                     f=atm%atom(i)%mx(m)
                     der=   -av*(otr(i)*drc(m,i) + oti(i)*drs(m,i))
                     der=der-bv*(oti(i)*drc(m,i) - otr(i)*drs(m,i))
                     der=2.0*der*tpi
                     deriv(k) = sign(1.0_cp,f)*der+deriv(k)
                   end if
                 end do

                 k=Atm%atom(i)%lbiso  !Derivatives w.r.t. Biso  S-XTAL
                 if(k /= 0) then
                   f=Atm%atom(i)%mbiso
                   der=   -av*( otr(i)*frc(i) - oti(i)*frs(i) )
                   der=der-bv*( oti(i)*frc(i) + otr(i)*frs(i) )
                   der=2.0*der*sn
                   deriv(k) = sign(1.0_cp,f)*der+deriv(k)
                 end if

                 k=Atm%atom(i)%locc    !Derivatives w.r.t. occupation factor  S-XTAL
                 if(k /= 0) then
                   f=Atm%atom(i)%mocc
                   der=    av*( otr(i)*frc(i) - oti(i)*frs(i) )
                   der=der+bv*( oti(i)*frc(i) + otr(i)*frs(i) )
                   der=2.0*der/atm%atom(i)%occ
                   deriv(k) = sign(1.0_cp,f)*der+deriv(k)
                 end if

                 do m=4,9        !Derivatives w.r.t. anisotropic temperature factors S-XTAL
                    j=m-3
                    k=Atm%atom(i)%lu(j)
                    if(k /= 0) then
                      f=Atm%atom(i)%mu(j)
                      der=   -av*(otr(i)*drc(m,i) - oti(i)*drs(m,i))
                      der=der-bv*(oti(i)*drc(m,i) + otr(i)*drs(m,i))
                      der=2.0*der
                      if(j > 3) der=2.0*der
                      deriv(k) = sign(1.0_cp,f)*der+deriv(k)
                    end if
                 end do

             end do
         end if !modi
       end if  !derivatives

       return
    End Subroutine Calc_StrFactor

    !!--++
    !!----
    !!---- Subroutine Calc_hkl_StrFactor(mode,rad,hn,sn,Atm,Grp,sf2,deriv,fc)
    !!----    character(len=*),                   intent(in) :: mode !S-XTAL (S) or Powder (P)
    !!----    character(len=*),                   intent(in) :: rad  !Radiation: X-rays, Neutrons
    !!----    integer, dimension(3)               intent(in) :: hn
    !!----    real(kind=cp)                       intent(in) :: sn !(sinTheta/Lambda)**2
    !!----    type(atom_list_type),               intent(in) :: Atm
    !!----    type(space_group_type),             intent(in) :: Grp
    !!----    real(kind=cp)                       intent(out):: sf2
    !!----    real(kind=cp),dimension(:),optional,intent(out):: deriv
    !!----    complex, optional,                  intent(out):: fc
    !!----
    !!----    Calculate Structure Factor for reflection "h=(hkl)" not related with
    !!----    previous lists and derivatives with respect to refined parameters.
    !!----    This subroutine calculates the form-factors internally without using
    !!----    global tables. The purpose of this procedure is to avoid the use of
    !!----    too much memory in tables.
    !!----
    !!---- Update: April - 2009
    !!
    Subroutine Calc_hkl_StrFactor(mode,rad,hn,sn,Atm,Grp,sf2,deriv,fc)
       !---- Arguments ----!
       character(len=*),                   intent(in) :: mode
       character(len=*),                   intent(in) :: rad
       integer,dimension(3),               intent(in) :: hn
       real(kind=cp),                      intent(in) :: sn !(sinTheta/Lambda)**2
       type(atom_list_type),               intent(in) :: Atm
       type(space_group_type),             intent(in) :: Grp
       real(kind=cp),                      intent(out):: sf2
       real(kind=cp),dimension(:),optional,intent(out):: deriv
       complex, optional,                  intent(out):: fc

       !---- Local Variables ----!
       character(len=1)                      :: modi
       integer                               :: i,j,k,m
       real(kind=cp)                         :: arg,anis,cosr,sinr,scosr,ssinr,fr,fi,der, hnt
       real(kind=cp)                         :: a1,a2,a3,a4,b1,b2,b3,b4,av,bv,f,occ,b, Tob
       real(kind=cp),dimension(3)            :: h
       real(kind=cp),dimension(6)            :: beta
       real(kind=cp),dimension(Atm%natoms)   :: frc,frs,otr,oti,afpxn,ff
       real(kind=cp),dimension(9,Atm%natoms) :: drs,drc


       !--- Initialising local variables
       a1=0.0
       a2=0.0
       a3=0.0
       a4=0.0
       b1=0.0
       b2=0.0
       b3=0.0
       b4=0.0
       av=0.0
       bv=0.0
       frc=0.0
       frs=0.0
       otr=0.0
       oti=0.0
       modi=u_case(mode(1:1))
       !Setting up the scattering form factors and multiply by group specific
       !coefficients for calculating structure factors per conventional cell
       !---- Modify the scattering factors to include the
       !---- multipliers factors concerning centre of symmetry and
       !---- centred translations
       fr=1.0; fi=1.0
       if (Grp%Centred == 2) fr=2.0
       if (Grp%NumLat  > 1)  fi=Grp%NumLat
       Select Case (rad(1:1))
          Case("N")
              afpxn(:)=fr*fi*afp(:)
          Case("X","E")
              do i=1,Nspecies
                ff(i)=FF_c(i)
                do j=1,4
                 ff(i)=ff(i)+FF_a(j,i)*exp(-sn*FF_b(j,i))
                end do
                 if (rad(1:1) == "E") ff(i)=0.023934*(FF_Z(i)-ff(i))/sn !Mott-Bethe formula fe=me^2/(8pi Eps0 h^2) (Z-fx(s))/s^2
              end do
              do i=1,Atm%natoms
                j=P_a(i)   !pointer has been set up in Initialization subroutine
                afpxn(i)= fr*fi*ff(j)
              end do
       End Select

       fr=1.0
       fi=0.0
       if(Grp%Centred == 2) then
            do i=1,Atm%natoms
               arg=0.0
               scosr=0.0
               ssinr=0.0
               drs(:,i)=0.0
               drc(:,i)=0.0
               do k=1,grp%NumOps
                  h=Hkl_R(hn,grp%symop(k))                   !Calculations in-lining
                  hnt=dot_product(real(hn),Grp%SymOp(k)%Tr)  !Calculations in-lining
                  arg=tpi*(dot_product(h,Atm%atom(i)%x)+hnt)
                  anis=1.0
                  if(Atm%atom(i)%thtype == "aniso") then
                    beta=Atm%atom(i)%u(1:6)
                    anis=     h(1)*h(1)*beta(1)+     h(2)*h(2)*beta(2)+    h(3)*h(3)*beta(3) &
                         +2.0*h(1)*h(2)*beta(4)+ 2.0*h(1)*h(3)*beta(5)+2.0*h(2)*h(3)*beta(6)
                    anis=exp(-anis)
                  end if
                  cosr=COS(arg)*anis*fr     !fr*cos{2pi(hT Rs rj+ts)}*exp(-{hTRsBetaj RsTh})
                  scosr=scosr+cosr          !FRC= SIG fr(j,s)cos{2pi(hT Rs rj+ts)}*Ta(s)

                  if(present(deriv)) then
                     sinr=SIN(arg)*anis*fr   !fr*sin{2pi(hT Rs rj+ts)}*exp(-{hTRsBetaj RsTh})
                     drc(1:3,i)=drc(1:3,i)+h(1:3)*sinr      ! -
                     drc(4,i)=drc(4,i)+h(1)*h(1)*cosr
                     drc(5,i)=drc(5,i)+h(2)*h(2)*cosr
                     drc(6,i)=drc(6,i)+h(3)*h(3)*cosr
                     drc(7,i)=drc(7,i)+h(1)*h(2)*cosr
                     drc(8,i)=drc(8,i)+h(1)*h(3)*cosr
                     drc(9,i)=drc(9,i)+h(2)*h(3)*cosr
                  end if

               end do ! symmetry
               occ= atm%atom(i)%occ
               b=atm%atom(i)%biso
               Tob= occ * exp(-b*sn)
               frc(i) = scosr
               otr(i) = afpxn(i)* Tob
               oti(i) =  afpp(i)* Tob
               a1= a1 + otr(i)*frc(i)
               b1= b1 + oti(i)*frc(i)

            end do ! Atoms

            av = a1-a2-a3-a4    !real part of the structure factor
            bv = b1-b2+b3+b4    !imaginary part of the structure factor

       else

            do i=1,Atm%natoms
               arg=0.0
               scosr=0.0
               ssinr=0.0
               drs(:,i)=0.0
               drc(:,i)=0.0
               do k=1,grp%NumOps
                  h=Hkl_R(hn,grp%symop(k))
                  hnt=dot_product(real(hn),Grp%SymOp(k)%Tr)
                  arg=tpi*(dot_product(h,Atm%atom(i)%x)+hnt)
                  anis=1.0
                  if(Atm%atom(i)%thtype == "aniso") then
                    beta=Atm%atom(i)%u(1:6)
                    anis=     h(1)*h(1)*beta(1)+     h(2)*h(2)*beta(2)+    h(3)*h(3)*beta(3) &
                         +2.0*h(1)*h(2)*beta(4)+ 2.0*h(1)*h(3)*beta(5)+2.0*h(2)*h(3)*beta(6)
                    anis=exp(-anis)
                  end if
                  cosr=COS(arg)*anis*fr     !fr*cos{2pi(hT Rs rj+ts)}*exp(-{hTRsBetaj RsTh})
                  sinr=SIN(arg)*anis*fr     !fr*sin{2pi(hT Rs rj+ts)}*exp(-{hTRsBetaj RsTh})
                  scosr=scosr+cosr          !FRC= SIG fr(j,s)cos{2pi(hT Rs rj+ts)}*Ta(s)
                  ssinr=ssinr+sinr          !FRS= SIG fr(j,s)sin{2pi(hT Rs rj+ts)}*Ta(s)

                  if(present(deriv)) then
                     drc(1:3,i)=drc(1:3,i)+h(1:3)*sinr      ! -
                     drs(1:3,i)=drs(1:3,i)+h(1:3)*cosr      ! +

                     drc(4,i)=drc(4,i)+h(1)*h(1)*cosr
                     drc(5,i)=drc(5,i)+h(2)*h(2)*cosr
                     drc(6,i)=drc(6,i)+h(3)*h(3)*cosr
                     drc(7,i)=drc(7,i)+h(1)*h(2)*cosr
                     drc(8,i)=drc(8,i)+h(1)*h(3)*cosr
                     drc(9,i)=drc(9,i)+h(2)*h(3)*cosr

                     drs(4,i)=drs(4,i)+h(1)*h(1)*sinr
                     drs(5,i)=drs(5,i)+h(2)*h(2)*sinr
                     drs(6,i)=drs(6,i)+h(3)*h(3)*sinr
                     drs(7,i)=drs(7,i)+h(1)*h(2)*sinr
                     drs(8,i)=drs(8,i)+h(1)*h(3)*sinr
                     drs(9,i)=drs(9,i)+h(2)*h(3)*sinr
                  end if

               end do ! symmetry
               occ= atm%atom(i)%occ
               b=atm%atom(i)%biso
               Tob= occ * exp(-b*sn)
               frc(i) = scosr
               frs(i) = ssinr
               otr(i) = afpxn(i)* Tob
               oti(i) =  afpp(i)* Tob
               a1= a1 + otr(i)*frc(i)
               b1= b1 + oti(i)*frc(i)
               a3 = a3 + oti(i)*frs(i)
               b3 = b3 + otr(i)*frs(i)

            end do ! Atoms

            av = a1-a2-a3-a4    !real part of the structure factor
            bv = b1-b2+b3+b4    !imaginary part of the structure factor

       end if

       If(modi == "P") then
          sf2 = a1*a1 + a2*a2 + a3*a3 + a4*a4 + b1*b1 + b2*b2 + b3*b3 + b4*b4
          sf2 = sf2 + 2.0*(b1*b4 -  a1*a4 + a2*a3 - b2*b3)
       else
          sf2= av*av+bv*bv
       End if

       if(present(fc)) then
         fc=cmplx(av,bv)
       end if

       if(present(deriv)) then

         if(modi == "P") then

             do i=1,Atm%natoms
                !derivatives with respect to coordinates  POWDER
                do m=1,3
                   k= Atm%atom(i)%lx(m)
                   if(k /= 0) then
                     f=atm%atom(i)%mx(m)
                     der= otr(i)*(-a1*drc(m,i)+b3*drs(m,i))+oti(i)*(-b1*drc(m,i)+a3*drs(m,i))
                     der=2.0*der*tpi
                     deriv(k) = sign(1.0_cp,f)*der+deriv(k)
                   end if
                 end do

                 k=Atm%atom(i)%lbiso  !Derivatives w.r.t. Biso  POWDER
                 if(k /= 0) then
                   f=Atm%atom(i)%mbiso
                   der= otr(i)*(a1*frc(i) +b3*frs(i))+oti(i)*(b1*frc(i) +a3*frs(i))
                   der=-2.0*der*sn
                   deriv(k) = sign(1.0_cp,f)*der+deriv(k)
                 end if

                 k=Atm%atom(i)%locc    !Derivatives w.r.t. occupation factor   POWDER
                 if(k /= 0) then
                   f=Atm%atom(i)%mocc
                   der= otr(i)*(a1*frc(i)+b3*frs(i))+oti(i)*(b1*frc(i)+a3*frs(i))
                   der=2.0*der/atm%atom(i)%occ
                   deriv(k) = sign(1.0_cp,f)*der+deriv(k)
                 end if

                 do m=4,9      !Derivatives w.r.t. anisotropic temperature factors   POWDER
                    j=m-3
                    k=Atm%atom(i)%lu(j)
                    if(k /= 0) then
                      f=Atm%atom(i)%mu(j)
                      der=  otr(i)*(a1*drc(i,j)+b3*drs(m,i))+oti(i)*(b1*drc(m,i)+a3*drs(m,i))
                      der=-2.0*der
                      if(j > 3) der=2.0*der
                      deriv(k) = sign(1.0_cp,f)*der+deriv(k)
                    end if
                 end do

             end do

         else

             do i=1,Atm%natoms
                !derivatives with respect to coordinates  S-XTAL
                do m=1,3
                   k= Atm%atom(i)%lx(m)
                   if(k /= 0) then
                     f=atm%atom(i)%mx(m)
                     der=   -av*(otr(i)*drc(m,i) + oti(i)*drs(m,i))
                     der=der-bv*(oti(i)*drc(m,i) - otr(i)*drs(m,i))
                     der=2.0*der*tpi
                     deriv(k) = sign(1.0_cp,f)*der+deriv(k)
                   end if
                 end do

                 k=Atm%atom(i)%lbiso  !Derivatives w.r.t. Biso  S-XTAL
                 if(k /= 0) then
                   f=Atm%atom(i)%mbiso
                   der=   -av*( otr(i)*frc(i) - oti(i)*frs(i) )
                   der=der-bv*( oti(i)*frc(i) + otr(i)*frs(i) )
                   der=2.0*der*sn
                   deriv(k) = sign(1.0_cp,f)*der+deriv(k)
                 end if

                 k=Atm%atom(i)%locc    !Derivatives w.r.t. occupation factor  S-XTAL
                 if(k /= 0) then
                   f=Atm%atom(i)%mocc
                   der=    av*( otr(i)*frc(i) - oti(i)*frs(i) )
                   der=der+bv*( oti(i)*frc(i) + otr(i)*frs(i) )
                   der=2.0*der/atm%atom(i)%occ
                   deriv(k) = sign(1.0_cp,f)*der+deriv(k)
                 end if

                 do m=4,9        !Derivatives w.r.t. anisotropic temperature factors S-XTAL
                    j=m-3
                    k=Atm%atom(i)%lu(j)
                    if(k /= 0) then
                      f=Atm%atom(i)%mu(j)
                      der=   -av*(otr(i)*drc(m,i) - oti(i)*drs(m,i))
                      der=der-bv*(oti(i)*drc(m,i) + otr(i)*drs(m,i))
                      der=2.0*der
                      if(j > 3) der=2.0*der
                      deriv(k) = sign(1.0_cp,f)*der+deriv(k)
                    end if
                 end do

             end do
         end if !modi
       end if  !derivatives

       return
    End Subroutine Calc_hkl_StrFactor

    !!--++
    !!--++ Subroutine Calc_Table_AB(Nref,Atm,Grp)
    !!--++    integer,                intent(in) :: Nref
    !!--++    type(atom_list_type),   intent(in) :: Atm
    !!--++    type(space_group_type), intent(in) :: Grp
    !!--++
    !!--++    (Private)
    !!--++    Calculate Table with Aj(h) and Bj(h) values
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Calc_Table_AB(Nref,Atm,Grp)
       !---- Arguments ----!
       integer,                intent(in) :: Nref
       type(atom_list_type),   intent(in) :: Atm
       type(space_group_type), intent(in) :: Grp

       !---- Local Variables ----!
       integer                       :: i,j,k
       real(kind=cp)                 :: arg,anis
       real(kind=cp),dimension(3)    :: h
       real(kind=cp),dimension(6)    :: beta

       Ajh=0.0
       Bjh=0.0
       if(Grp%Centred == 2) then
         do j=1,Nref
            do i=1,Atm%natoms
               arg=0.0
               do k=1,grp%NumOps
                  h=hr(k,j)%h
                  arg=tpi*(dot_product(h,Atm%atom(i)%x)+ht(k,j))
                  anis=1.0
                  if(Atm%atom(i)%thtype == "aniso") then
                    beta=Atm%atom(i)%u(1:6)
                    anis=     h(1)*h(1)*beta(1)+     h(2)*h(2)*beta(2)+    h(3)*h(3)*beta(3) &
                         +2.0*h(1)*h(2)*beta(4)+ 2.0*h(1)*h(3)*beta(5)+2.0*h(2)*h(3)*beta(6)
                    anis=exp(-anis)
                  end if
                  Ajh(i,j)=Ajh(i,j)+cos(arg)*anis
               end do ! symmetry
            end do ! Atoms
         end do ! Reflections
       else
         do j=1,Nref
            do i=1,Atm%natoms
               arg=0.0
               do k=1,grp%NumOps
                  h=hr(k,j)%h
                  arg=tpi*(dot_product(h,Atm%atom(i)%x)+ht(k,j))
                  anis=1.0
                  if(Atm%atom(i)%thtype == "aniso") then
                    beta=Atm%atom(i)%u(1:6)
                    anis=     h(1)*h(1)*beta(1)+     h(2)*h(2)*beta(2)+    h(3)*h(3)*beta(3) &
                         +2.0*h(1)*h(2)*beta(4)+ 2.0*h(1)*h(3)*beta(5)+2.0*h(2)*h(3)*beta(6)
                    anis=exp(-anis)
                  end if
                  Ajh(i,j)=Ajh(i,j)+cos(arg)*anis
                  Bjh(i,j)=Bjh(i,j)+sin(arg)*anis
               end do ! symmetry
            end do ! Atoms
         end do ! Reflections
       end if

       return
    End Subroutine Calc_Table_AB

    !!--++ Subroutine Calc_Table_TH(Reflex,Atm)
    !!--++    type(reflection_List_type),   intent(in) :: Reflex
    !!--++    type(atom_list_type),        intent(in) :: Atm
    !!--++
    !!--++    (Private)
    !!--++    Calculate a Table of Isotropinc Thermal contribution and occupation
    !!--..         TH(Natoms,Nref)
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Calc_Table_TH(Reflex,Atm)
       !---- Argument ----!
       type(reflection_list_type), intent(in) :: Reflex
       type(atom_list_type),      intent(in) :: Atm

       !---- Local variables ----!
       integer          :: i,j
       real(kind=cp)    :: b,s

       !---- Isotropic model ----!
       do j=1,reflex%nref
          s=reflex%ref(j)%s
          do i=1,atm%natoms
             b=atm%atom(i)%biso
             th(i,j)= atm%atom(i)%occ * exp(-b*s*s)
          end do
       end do

       return
    End Subroutine Calc_Table_TH

    !!--++
    !!--++ Subroutine Create_Table_fabc_Xray(Atm,lambda,lun)
    !!--++    type(atom_list_type),      intent(in) :: Atm
    !!--++    real(kind=cp), optiona      intent(in) :: lambda
    !!--++    integer, optional,          intent(in) :: lun
    !!--++
    !!--++    (Private)
    !!--++    Calculate a Table of Coefficients for Atomic Form Factors for X-Ray
    !!--++    ff_A(4,species),ff_B(4,Nspecies),ff_C(Nspecies), AFP(Nspecies), AFPP(Nspecies)
    !!--++    ff_z(Nspecies) contains the atomic number of the chemical species (useful for Electron diffraction)
    !!--++    p_a(Natoms) => pointer to the species of each atom
    !!--++
    !!--++ Update: April - 2009
    !!
    Subroutine Create_Table_fabc_Xray(Atm,lambda,elect,lun)
       !---- Arguments ----!
       type(atom_list_type),       intent(in) :: Atm
       real(kind=cp), optional,    intent(in) :: lambda
       integer, optional,          intent(in) :: elect
       integer, optional,          intent(in) :: lun

       !---- Local Variables ----!
       character(len=4)               :: symbcar
       integer                        :: i,j, k,n,L
       integer, dimension(atm%natoms) :: ix,jx,ia
       real(kind=cp)                  :: dmin,d

       !---- Init ----!
       err_sfac=.false.

       !---- Load form factor values for XRay ----!
       call Set_Xray_Form()

       !---- Found Species on Xray_Form ----!
       ix=0
       jx=0
       n=0
       if(allocated(P_a)) deallocate(P_a)
       allocate(P_a(atm%natoms))

       do i=1,atm%natoms
          symbcar=l_case(atm%atom(i)%SfacSymb)
          do j=1,Num_Xray_Form
             if (symbcar /= Xray_form(j)%Symb) cycle
             ix(i)=j
             if(any(jx == j) ) exit
             n=n+1
             jx(n)=j
             ia(n)=i
             exit
          end do
       end do

       if (any(ix==0)) then
          err_sfac=.true.
          ERR_SFac_Mess="The Species "//symbcar//" was not found"
          return
       end if

       do i=1,atm%natoms
         j=ix(i)
         do k=1,n
           if(jx(k) == j) then
             P_a(i)=k              !The atom i is of species k
             exit
           end if
         end do
       end do
       Nspecies=n !Global private variable (Total number of chemical species)
       if(allocated(FF_a)) deallocate (FF_a)
       if(allocated(FF_b)) deallocate (FF_b)
       if(allocated(FF_c)) deallocate (FF_c)
       if(allocated(FF_z)) deallocate (FF_z)
       allocate(FF_a(4,n),FF_b(4,n),FF_c(n),FF_z(n))
       do k=1,n
          j = jx(k)
          i = ia(k)
          FF_a(:,k)= xray_form(j)%a(:)
          FF_b(:,k)= xray_form(j)%b(:)
          FF_c(  k)= xray_form(j)%c
          FF_z(  k)= xray_form(j)%Z
       end do

       if (present(lun)) then
          if(present(elect)) then
            write(unit=lun,fmt="(/,a)") "  INFORMATION FROM TABULATED X-RAY SCATTERING FACTORS (For Electron Diffraction)"
            write(unit=lun,fmt="(a,/)") "  =============================================================================="
          else
            write(unit=lun,fmt="(/,a)") "  INFORMATION FROM TABULATED X-RAY SCATTERING FACTORS"
            write(unit=lun,fmt="(a,/)") "  ==================================================="
          end if
       End if
       if(.not. present(elect)) then
        if (present(lambda)) then
          !---- Load anomalous scattering form factor values for XRays ----!
          call Set_Delta_Fp_Fpp()

          !---- Select wavelength (by default is CuKalpha1: k=5 in the list) ----!
          dmin=1000.0
          do i=1,5
             d=abs(lambda-Xray_Wavelengths(i)%Kalfa(1))
             if (d < dmin) then
                dmin=d
                k=i        !Selection of the index for fp and fpp lists
             end if
          end do

          !---- Found Species on Anomalous_ScFac ----!
          do i=1,atm%natoms
             symbcar=l_case(atm%atom(i)%chemsymb)
             do j=1,Num_Delta_Fp
                if (symbcar /= Anomalous_ScFac(j)%Symb) cycle
                afp(i)=Anomalous_ScFac(j)%fp(k)
                afpp(i)=Anomalous_ScFac(j)%fpp(k)
                exit
             end do
          end do
          call Remove_Delta_Fp_Fpp()
        else
           if (present(lun)) then
             write(unit=lun,fmt="(a)")    "  Missed lambda, anomalous dipersion corrections not applied   "
             write(unit=lun,fmt="(a)")    "  The default wavelength is that of Cu-Kalpha1 spectral line  "
           end if
        end if
       end if !present(elect)


       !---- Printing Information ----!
       if (present(lun)) then
         if(present(elect)) then
          write(unit=lun,fmt="(/,a,/)")    "   ATOMIC SCATTERING FACTOR COEFFICIENTS: {A(i),B(i),I=1,4},C  and Atomic Number "
          write(unit=lun,fmt="(a,i3)")     "   Number of chemically different species: ",n
         write(unit=lun,fmt="(/,a)") &
               "   Atom     a1       b1       a2       b2       a3       b3       a4       b4        c       Z"
          do k=1,n
             j = jx(k)
             i = ia(k)
             write(unit=lun,fmt="(a,9F9.5,i7)")    &
                           "     "//atm%atom(i)%chemsymb, &
                           (xray_form(j)%a(L),xray_form(j)%b(L), L=1,4), xray_form(j)%c, &
                           xray_form(j)%z
          end do
        else
          write(unit=lun,fmt="(/,a,/)")    "   ATOMIC SCATTERING FACTOR COEFFICIENTS: {A(i),B(i),I=1,4},C  Dfp  Dfpp "
          write(unit=lun,fmt="(a,i3)")     "   Number of chemically different species: ",n
          write(unit=lun,fmt="(/,a)") &
               "   Atom     a1       b1       a2       b2       a3       b3       a4       b4        c      Dfp     Dfpp"
          do k=1,n
             j = jx(k)
             i = ia(k)
             write(unit=lun,fmt="(a,11F9.5)")    &
                           "     "//atm%atom(i)%chemsymb, &
                           (xray_form(j)%a(L),xray_form(j)%b(L), L=1,4), xray_form(j)%c, &
                           afp(i), afpp(i)
          end do
         end if
          write(unit=lun,fmt="(/,/)")
       end if

       call Remove_Xray_Form()

       return
    End Subroutine Create_Table_fabc_Xray

    !!--++
    !!--++ Subroutine Create_Table_AF0_Electrons(Reflex,Atm,lun)
    !!--++    type(reflection_List_type), intent(in) :: Reflex
    !!--++    type(atom_list_type),       intent(in) :: Atm
    !!--++    integer, optional,          intent(in) :: lun
    !!--++
    !!--++    (Private)
    !!--++    Calculate a Table of Atomic Factors for Electrons
    !!--++    applying the Mott-Bethe formula:
    !!--++    fe=me^2/(8pi Eps0 h^2) (Z-fx(s))/s^2
    !!--++
    !!--++ Update: April - 2009
    !!
    Subroutine Create_Table_AF0_Electrons(Reflex,Atm,lun)
       !---- Arguments ----!
       type(reflection_list_type), intent(in) :: Reflex
       type(atom_list_type),       intent(in) :: Atm
       integer, optional,          intent(in) :: lun

       !---- Local Variables ----!
       character(len=4)               :: symbcar
       integer                        :: i,j, k,n,L
       integer, dimension(atm%natoms) :: ix,jx,ia
       real(kind=cp)                  :: fx

       !---- Init ----!
       err_sfac=.false.

       !---- Load form factor values for XRay ----!
       call Set_Xray_Form()

       !---- Found Species on Xray_Form ----!
       ix=0
       jx=0
       n=0
       do i=1,atm%natoms
          symbcar=l_case(atm%atom(i)%SfacSymb)
          do j=1,Num_Xray_Form
             if (symbcar /= Xray_form(j)%Symb) cycle
             ix(i)=j
             if(any(jx == j) ) exit
             n=n+1
             jx(n)=j
             ia(n)=i
             exit
          end do
       end do

       if (present(lun)) then
         write(unit=lun,fmt="(/,a)") "  INFORMATION FROM TABULATED X-RAY SCATTERING FACTORS (For Electron Diffraction)"
          write(unit=lun,fmt="(a,/)") "  =============================================================================="
       End if

       if (any(ix==0)) then
          err_sfac=.true.
          ERR_SFac_Mess="The Species "//symbcar//" was not found"
       else
          !---- Fill AF Table ----!
          do j=1,reflex%nref
             do i=1,atm%natoms
                fx=fj(reflex%ref(j)%s,xray_form(ix(i))%a,xray_form(ix(i))%b,xray_form(ix(i))%c)+afp(i)
                !Mott-Bethe formula fe=me^2/(8pi Eps0 h^2) (Z-fx(s))/s^2
                af0(i,j)=0.023934*(xray_form(ix(i))%Z-fx)/(reflex%ref(j)%s*reflex%ref(j)%s)
             end do
          end do
       end if

       !---- Printing Information ----!
       if (present(lun)) then
          write(unit=lun,fmt="(/,a,/)")    "   ATOMIC SCATTERING FACTOR COEFFICIENTS: {A(i),B(i),I=1,4},C  and Atomic Number "
          write(unit=lun,fmt="(a,i3)")     "   Number of chemically different species: ",n
          write(unit=lun,fmt="(/,a)") &
               "   Atom     a1       b1       a2       b2       a3       b3       a4       b4        c       Z"
          do k=1,n
             j = jx(k)
             i = ia(k)
             write(unit=lun,fmt="(a,9F9.5,i7)")    &
                           "     "//atm%atom(i)%chemsymb, &
                           (xray_form(j)%a(L),xray_form(j)%b(L), L=1,4), xray_form(j)%c, &
                            xray_form(j)%Z
          end do
          write(unit=lun,fmt="(/,/)")
       end if

       call Remove_Xray_Form()

       return
    End Subroutine Create_Table_AF0_Electrons

    !!--++
    !!--++ Subroutine Create_Table_AF0_Xray(Reflex,Atm,lambda,lun)
    !!--++    type(reflection_List_type), intent(in) :: Reflex
    !!--++    type(atom_list_type),      intent(in) :: Atm
    !!--++    real(kind=cp), optiona      intent(in) :: lambda
    !!--++    integer, optional,          intent(in) :: lun
    !!--++
    !!--++    (Private)
    !!--++    Calculate a Table of Atomic Factors for X-Ray
    !!--..      AF0(Natoms,Nref), AFP(Natoms), AFPP(Natoms)
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Create_Table_AF0_Xray(Reflex,Atm,lambda,lun)
       !---- Arguments ----!
       type(reflection_list_type), intent(in) :: Reflex
       type(atom_list_type),      intent(in) :: Atm
       real(kind=cp), optional,    intent(in) :: lambda
       integer, optional,          intent(in) :: lun

       !---- Local Variables ----!
       character(len=4)               :: symbcar
       integer                        :: i,j, k,n,L
       integer, dimension(atm%natoms) :: ix,jx,ia
       real(kind=cp)                  :: dmin,d

       !---- Init ----!
       err_sfac=.false.

       !---- Load form factor values for XRay ----!
       call Set_Xray_Form()

       !---- Found Species on Xray_Form ----!
       ix=0
       jx=0
       n=0
       do i=1,atm%natoms
          symbcar=l_case(atm%atom(i)%SfacSymb)
          do j=1,Num_Xray_Form
             if (symbcar /= Xray_form(j)%Symb) cycle
             ix(i)=j
             if(any(jx == j) ) exit
             n=n+1
             jx(n)=j
             ia(n)=i
             exit
          end do
       end do

       if (present(lun)) then
          write(unit=lun,fmt="(/,a)") "  INFORMATION FROM TABULATED X-RAY SCATTERING FACTORS"
          write(unit=lun,fmt="(a,/)") "  ==================================================="
       End if
       if (present(lambda)) then
          !---- Load anomalous scattering form factor values for XRays ----!
          call Set_Delta_Fp_Fpp()

          !---- Select wavelength (by default is CuKalpha1: k=5 in the list) ----!
          dmin=1000.0
          do i=1,5
             d=abs(lambda-Xray_Wavelengths(i)%Kalfa(1))
             if (d < dmin) then
                dmin=d
                k=i        !Selection of the index for fp and fpp lists
             end if
          end do

          !---- Found Species on Anomalous_ScFac ----!
          do i=1,atm%natoms
             symbcar=l_case(atm%atom(i)%chemsymb)
             do j=1,Num_Delta_Fp
                if (symbcar /= Anomalous_ScFac(j)%Symb) cycle
                afp(i)=Anomalous_ScFac(j)%fp(k)
                afpp(i)=Anomalous_ScFac(j)%fpp(k)
                exit
             end do
          end do
          call Remove_Delta_Fp_Fpp()
       else
           if (present(lun)) then
             write(unit=lun,fmt="(a)")    "  Missed lambda, anomalous dipersion corrections not applied   "
             write(unit=lun,fmt="(a)")    "  The default wavelength is that of Cu-Kalpha1 spectral line  "
           end if
       end if

       if (any(ix==0)) then
          err_sfac=.true.
          ERR_SFac_Mess="The Species "//symbcar//" was not found"
       else
          !---- Fill AF Table ----!
          do j=1,reflex%nref
             do i=1,atm%natoms
                af0(i,j)=fj(reflex%ref(j)%s,xray_form(ix(i))%a,xray_form(ix(i))%b,xray_form(ix(i))%c)+afp(i)
             end do
          end do
       end if

       !---- Printing Information ----!
       if (present(lun)) then
          write(unit=lun,fmt="(/,a,/)")    "   ATOMIC SCATTERING FACTOR COEFFICIENTS: {A(i),B(i),I=1,4},C  Dfp  Dfpp "
          write(unit=lun,fmt="(a,i3)")     "   Number of chemically different species: ",n
          write(unit=lun,fmt="(/,a)") &
               "   Atom     a1       b1       a2       b2       a3       b3       a4       b4        c      Dfp     Dfpp"
          do k=1,n
             j = jx(k)
             i = ia(k)
             write(unit=lun,fmt="(a,11F9.5)")    &
                           "     "//atm%atom(i)%chemsymb, &
                           (xray_form(j)%a(L),xray_form(j)%b(L), L=1,4), xray_form(j)%c, &
                           afp(i), afpp(i)
          end do
          write(unit=lun,fmt="(/,/)")
       end if

       call Remove_Xray_Form()

       return
    End Subroutine Create_Table_AF0_Xray

    !!--++
    !!--++ Subroutine Create_Table_AFP_NeutNuc(Atm,lun)
    !!--++    type(atom_list_type),              intent(in) :: Atm
    !!--++    integer, optional,                  intent(in) :: lun
    !!--++
    !!--++    (Private)
    !!--++    Setting a Table of Fermi Lengths for Neutron Nuclear Scattering
    !!--..      AFP(Natoms)
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Create_Table_AFP_NeutNuc(Atm,lun)
       !---- Arguments ----!
       type(atom_list_type),    intent(in) :: Atm
       integer, optional,       intent(in) :: lun

       !---- Local Variables ----!
       character(len=4)                        :: symbcar
       integer                                 :: i,k,n
       character(len=4), dimension(atm%natoms) :: symb
       real(kind=cp),    dimension(atm%natoms) :: bs
       real(kind=cp)                           :: b

       !---- Init ----!
       err_sfac=.false.

       !---- Load chemical information ----!
       call set_chem_info()

       !---- Getting Fermi Lengths of atoms ----!
       symb="    "
       bs=0.0
       n=0
       do i=1,atm%natoms
          symbcar=u_case(atm%atom(i)%chemsymb)
          call Get_Fermi_Length(symbcar,b)
          if (abs(b) < 0.0001) then
             err_sfac=.true.
             ERR_SFac_Mess="The Fermi Length of Species "//symbcar//" was not found"
             return
          else
             afp(i) = b
             if(any(symb == symbcar)) cycle
             n=n+1
             symb(n)=symbcar
             bs(n) = b
          end if
       end do

       !---- Printing Information ----!
       if (present(lun)) then
          write(unit=lun,fmt="(/,a)")  "  INFORMATION FROM TABULATED NEUTRON SCATTERING FACTORS"
          write(unit=lun,fmt="(a,/)")  "  ==================================================="
          write(unit=lun,fmt="(a)")    "  FERMI LENGTHS "
          write(unit=lun,fmt="(a,i3)") "   Number of chemically different species: ",n
          write(unit=lun,fmt="(/,a)")  "   Atom     Fermi Length [10^(-12) cm]"
          do k=1,n
             write(unit=lun,fmt="(a,F15.6)")  "     "//symb(k), bs(k)
          end do
          write(unit=lun,fmt="(/,/)")
       end if

       call Remove_chem_info()

       return
    End Subroutine Create_Table_AFP_NeutNuc

    !!--++
    !!--++ Subroutine Create_Table_HR_HT(Reflex,Grp)
    !!--++    type(reflection_list_type), intent(in) :: Hkl
    !!--++    type(space_group_type),     intent(in) :: Grp
    !!--++
    !!--++    (Private)
    !!--++    Calculate a Table with HR and HT values
    !!--..       Hr(Grp%Numops,Reflex%Nref)
    !!--..       HT(Grp%Numops,Reflex%Nref)
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Create_Table_HR_HT(Reflex,Grp)
       !---- Arguments ----!
       type(reflection_list_type), intent(in) :: Reflex
       type(space_group_type),     intent(in) :: Grp

       !---- Local Variables ----!
       integer :: i,j

       do j=1,reflex%nref
          do i=1,grp%NumOps
             hr(i,j)%h=Hkl_R(reflex%ref(j)%h,grp%symop(i))
             ht(i,j)=dot_product(real(reflex%ref(j)%h),Grp%SymOp(i)%Tr)
          end do
       end do

       return
    End Subroutine Create_Table_HR_HT

    !!----
    !!---- Subroutine Init_Calc_hkl_StrFactors(Atm,Mode,lambda,lun)
    !!----    type(atom_list_type),                intent(in) :: Atm
    !!----    character(len=*),          optional, intent(in) :: Mode
    !!----    real(kind=cp),             optional, intent(in) :: lambda
    !!----    integer,                   optional, intent(in) :: lun  !Logical unit for writing scatt-factors
    !!----
    !!----    Allocates and initializes arrays for hkl - Structure Factors calculations.
    !!----    No calculation of fixed tables is performed. Should be called before using
    !!----    the subroutine Calc_hkl_StrFactor
    !!----
    !!---- Update: April - 2009
    !!
    Subroutine Init_Calc_hkl_StrFactors(Atm,Mode,lambda,lun)
       !---Arguments ---!
       type(atom_list_type),                intent(in) :: Atm
       character(len=*),          optional, intent(in) :: Mode
       real(kind=cp),             optional, intent(in) :: lambda
       integer,                   optional, intent(in) :: lun

       !--- Local variables ---!

       integer :: Natm
       integer :: ierr
       character(len=3) :: tipo

       tipo="XRA"
       if (present(mode)) tipo=adjustl(mode)
       tipo=U_Case(tipo)
       err_sfac=.false.
       Natm = Atm%natoms


       !---- Anomalous Scattering factor tables ----!
       if (allocated(AFP)) deallocate(AFP)
       allocate(AFP(Natm),stat=ierr)
       if (ierr /=0) then
          err_sfac=.true.
          ERR_SFac_Mess="Error on memory for AFP"
          return
       end if
       AFP=0.0

       if (allocated(AFPP)) deallocate(AFPP)
       allocate(AFPP(Natm),stat=ierr)
       if (ierr /=0) then
          err_sfac=.true.
          ERR_SFac_Mess="Error on memory for AFPP"
          return
       end if
       AFPP=0.0


       !---- Table Fabc ----!
       select case (tipo)

          case ("XRA")
             if (present(lambda)) then
                if (present(lun)) then
                   call Create_Table_Fabc_Xray(Atm,lambda,lun)
                else
                   call Create_Table_Fabc_Xray(Atm,lambda)
                end if
             else
                if (present(lun)) then
                   call Create_Table_Fabc_Xray(Atm,lun=lun)
                else
                   call Create_Table_Fabc_Xray(Atm)
                end if
             end if

          case ("ELE")

             if (present(lun)) then
                call Create_Table_Fabc_Xray(Atm,lun=lun)
             else
                call Create_Table_Fabc_Xray(Atm)
             end if


          case ("NUC")
             if (present(lun)) then
                call Create_Table_AFP_NeutNuc(Atm,lun=lun)
             else
                call Create_Table_AFP_NeutNuc(Atm)
             end if

       end select

       if (.not. err_sfac) SF_Initialized=.true.

       return
    End Subroutine Init_Calc_hkl_StrFactors

    !!----
    !!---- Subroutine Init_Calc_StrFactors(Reflex,Atm,Grp,Mode,lambda,lun)
    !!----    type(reflection_list_type),          intent(in) :: Reflex
    !!----    type(atom_list_type),                intent(in) :: Atm
    !!----    type(space_group_type),              intent(in) :: Grp
    !!----    character(len=*),          optional, intent(in) :: Mode
    !!----    real(kind=cp),             optional, intent(in) :: lambda
    !!----    integer,                   optional, intent(in) :: lun  !Logical unit for writing scatt-factors
    !!----
    !!----    Allocates and initializes arrays for Calc_StrFactors calculations.
    !!----    Calculations of fixed tables are performed. Should be called before using
    !!----    the subroutine Calc_StrFactor
    !!----
    !!---- Update: April - 2009
    !!
    Subroutine Init_Calc_StrFactors(Reflex,Atm,Grp,Mode,lambda,lun)
       !---Arguments ---!
       type(reflection_list_type),          intent(in) :: Reflex
       type(atom_list_type),                intent(in) :: Atm
       type(space_group_type),              intent(in) :: Grp
       character(len=*),          optional, intent(in) :: Mode
       real(kind=cp),             optional, intent(in) :: lambda
       integer,                   optional, intent(in) :: lun

       !--- Local variables ---!

       Call Init_Structure_Factors(Reflex,Atm,Grp,Mode,lambda,lun)
       !---- Table TH ----!
       Call Calc_Table_TH(Reflex,Atm)

       return
    End Subroutine Init_Calc_StrFactors

    !!----
    !!---- Subroutine Init_Structure_Factors(Reflex,Atm,Grp,Mode,lambda,lun)
    !!----    type(reflection_list_type),          intent(in) :: Reflex
    !!----    type(atom_list_type),                intent(in) :: Atm
    !!----    type(space_group_type),              intent(in) :: Grp
    !!----    character(len=*),          optional, intent(in) :: Mode
    !!----    real(kind=cp),             optional, intent(in) :: lambda
    !!----    integer,                   optional, intent(in) :: lun  !Logical unit for writing scatt-factors
    !!----
    !!----    Allocates and initializes arrays for Structure Factors calculations.
    !!----    A calculation of fixed tables is also performed.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Init_Structure_Factors(Reflex,Atm,Grp,Mode,lambda,lun)
       !---Arguments ---!
       type(reflection_list_type),          intent(in) :: Reflex
       type(atom_list_type),                intent(in) :: Atm
       type(space_group_type),              intent(in) :: Grp
       character(len=*),          optional, intent(in) :: Mode
       real(kind=cp),             optional, intent(in) :: lambda
       integer,                   optional, intent(in) :: lun

       !--- Local variables ---!

       integer :: Natm, Multr
       integer :: ierr

       err_sfac=.false.
       Natm = Atm%natoms
       Multr= Grp%Numops

       !---- Scattering factor tables ----!
       if (allocated(AF0)) deallocate(AF0)
       allocate(AF0(Natm,Reflex%Nref),stat=ierr)
       if (ierr /=0) then
          err_sfac=.true.
          ERR_SFac_Mess="Error on memory for AF0"
          return
       end if
       AF0=0.0

       !---- Anomalous Scattering factor tables ----!
       if (allocated(AFP)) deallocate(AFP)
       allocate(AFP(Natm),stat=ierr)
       if (ierr /=0) then
          err_sfac=.true.
          ERR_SFac_Mess="Error on memory for AFP"
          return
       end if
       AFP=0.0

       if (allocated(AFPP)) deallocate(AFPP)
       allocate(AFPP(Natm),stat=ierr)
       if (ierr /=0) then
          err_sfac=.true.
          ERR_SFac_Mess="Error on memory for AFPP"
          return
       end if
       AFPP=0.0

       !---- HR Table ----!
       if (allocated(HR)) deallocate(HR)
       allocate(HR(Multr,Reflex%Nref),stat=ierr)
       if (ierr /=0) then
          err_sfac=.true.
          ERR_SFac_Mess="Error on memory for HR"
          return
       end if
       HR=HR_Type(0)

       !---- HT Table ----!
       if (allocated(HT)) deallocate(HT)
       allocate(HT(Multr,Reflex%Nref),stat=ierr)
       if (ierr /=0) then
          err_sfac=.true.
          ERR_SFac_Mess="Error on memory for HTR"
          return
       end if
       HT=0.0

       if (allocated(TH)) deallocate(TH)
       allocate(TH(Natm,Reflex%Nref),stat=ierr)
       if (ierr /=0) then
          err_sfac=.true.
          ERR_SFac_Mess="Error on memory for HTR"
          return
       end if
       TH=0.0

       if (allocated(Ajh)) deallocate(Ajh)
       allocate(Ajh(Natm,Reflex%Nref), stat=ierr)
       if (ierr /=0) then
          err_sfac=.true.
          ERR_SFac_Mess="Error in Memory for Aj(h)"
          return
       end if
       Ajh=0.0

       if (allocated(Bjh)) deallocate(Bjh)
       allocate(Bjh(Natm,Reflex%Nref), stat=ierr)
       if (ierr /=0) then
          err_sfac=.true.
          ERR_SFac_Mess="Error in Memory for Bj(h)"
          return
       end if
       Bjh=0.0

       if (present(mode)) then
          if (present(lambda)) then
             if (present(lun)) then
                call Set_Fixed_Tables(Reflex,Atm,Grp,Mode,lambda,lun)
             else
                call Set_Fixed_Tables(Reflex,Atm,Grp,Mode,lambda)
             end if
          else
             if (present(lun)) then
                call Set_Fixed_Tables(Reflex,Atm,Grp,Mode,lun=lun)
             else
                call Set_Fixed_Tables(Reflex,Atm,Grp,Mode)
             end if
          end if
       else
          if (present(lambda)) then
             if (present(lun)) then
                call Set_Fixed_Tables(Reflex,Atm,Grp,lambda=lambda,lun=lun)
             else
                call Set_Fixed_Tables(Reflex,Atm,Grp,lambda=lambda)
             end if
          else
             if (present(lun)) then
                call Set_Fixed_Tables(Reflex,Atm,Grp,lun=lun)
             else
                call Set_Fixed_Tables(Reflex,Atm,Grp)
             end if
          end if
       end if

       if (.not. err_sfac) SF_Initialized=.true.

       return
    End Subroutine Init_Structure_Factors

    !!----
    !!---- Subroutine Magnetic_Structure_Factors(Cell,Atm,Grp,maxs,Reflex,Stf)
    !!----    type(Crystal_Cell_type),            intent(in)     :: Cell
    !!----    type(atom_list_type),               intent(in)     :: Atm
    !!----    type(magnetic_space_group_type),    intent(in)     :: Grp
    !!----    Real(kind=cp)                       intent(in)     :: maxs (maximum sinTheta/Lambda)
    !!----    type(reflect_list_type),            intent(out)    :: Reflex
    !!----    type(Strf_list_type),               intent(out)    :: Stf
    !!----
    !!----    Calculation of Structure Factors (nuclear and magnetic) when
    !!----    the crystal and magnetic structure is described by a Shubnikov
    !!----    group
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Magnetic_Structure_Factors(Cell,Atm,Grp,maxs,Reflex,Stf,lun)
       type(Crystal_Cell_type),            intent(in)     :: Cell
       type(atom_list_type),               intent(in out) :: Atm
       type(magnetic_space_group_type),    intent(in)     :: Grp
       Real(kind=cp),                      intent(in)     :: maxs !(maximum sinTheta/Lambda)
       type(reflect_list_type),            intent(out)    :: Reflex
       type(Strf_list_type),               intent(out)    :: Stf
       Integer, optional,                  intent(in)     :: lun

       !---- Local variables ----!
       character(len=50) :: mess
       integer           :: i
       logical           :: ok
       Type(Scattering_Species_Type) :: Scf
       type(reflect_type),dimension(:),allocatable :: rf

       call Hkl_Gen_Shub(Cell,Grp,maxs,Reflex%NRef,rf)
       Stf%Nref=Reflex%NRef
       allocate(Stf%Strf(Stf%Nref),Reflex%Ref(Reflex%NRef))
       if(present(lun)) then
         call Set_Form_Factors(Atm,Scf,ok,mess,lun=lun,mag=.true.)
       else
         call Set_Form_Factors(Atm,Scf,ok,mess,mag=.true.)
       end if
       if(.not. ok) then
         write(unit=*,fmt="(a)") " => "//trim(mess)
         return
       end if
       do i=1,Reflex%Nref
         Reflex%Ref(i)=rf(i)
         call Calc_Mag_Structure_Factor("P",Reflex%Ref(i),Cell,Grp,Atm,Scf,Stf%Strf(i))
       end do
       return
    End Subroutine Magnetic_Structure_Factors


    !!----
    !!---- Subroutine Modify_SF(Reflex,Atm,Grp,List,Nlist,Mode)
    !!----    type(reflection_list_type),         intent(in out) :: Reflex
    !!----    type(atom_list_type),              intent(in)     :: Atm
    !!----    type(space_group_type),             intent(in)     :: Grp
    !!----    integer,dimension(:),               intent(in)     :: List
    !!----    integer,                            intent(in)     :: Nlist
    !!----    character(len=*),optional,          intent(in)     :: Mode
    !!----
    !!----    Recalculation of Structure Factors because a list of Atoms
    !!----    parameters were modified. List variable
    !!----    contains the number of atoms to be changed.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Modify_SF(Reflex,Atm,Grp,List,Nlist,partyp,Mode)
       !---- Arguments ----!
       type(reflection_list_type),   intent(in out) :: Reflex
       type(atom_list_type),         intent(in)     :: Atm
       type(space_group_type),       intent(in)     :: Grp
       integer,dimension(:),         intent(in)     :: List
       integer,                      intent(in)     :: NList
       character(len=*),optional,    intent(in)     :: partyp
       character(len=*),optional,    intent(in)     :: Mode

       !---- Local variables ----!
       character(len=2) :: typ
       integer          :: i,j,k,ii
       real(kind=cp)    :: arg,b,s

       typ="CO"
       if (present(partyp)) typ=adjustl(partyp)
       typ=U_Case(typ)

       select case (typ)

          case ("CO") ! by coordinates

            if(Grp%Centred == 2) then

               do j=1,Reflex%Nref
                  do ii=1,Nlist
                     i=list(ii)
                     Ajh(i,j)=0.0
                     arg=0.0
                     do k=1,grp%NumOps
                        arg=tpi*(dot_product(hr(k,j)%h,Atm%atom(i)%x)+ht(k,j))
                        Ajh(i,j)=Ajh(i,j)+cos(arg)
                     end do ! symmetry
                  end do ! NList
               end do ! Reflections

            else

               do j=1,Reflex%Nref
                  do ii=1,Nlist
                     i=list(ii)
                     arg=0.0
                     Ajh(i,j)=0.0
                     Bjh(i,j)=0.0
                     do k=1,grp%NumOps
                        arg=tpi*(dot_product(hr(k,j)%h,Atm%atom(i)%x)+ht(k,j))
                        Ajh(i,j)=Ajh(i,j)+cos(arg)
                        Bjh(i,j)=Bjh(i,j)+sin(arg)
                     end do ! symmetry
                  end do ! NList
               end do ! Reflections

            end if

          case ("TH") ! by thermal parameter or occupation number

             do j=1,Reflex%Nref
                s=reflex%ref(j)%s
                do ii=1,Nlist
                   i=list(ii)
                   b=atm%atom(i)%biso
                   th(i,j)=atm%atom(i)%occ*exp(-b*s*s)
                end do ! NList
             end do ! Reflections

       end select

       !---- Recalculation of SF ----!
       if(present(mode)) then
         if(mode == "XRA" .or. mode == "ELE") then
            call Sum_AB(Reflex,Atm%Natoms,Grp%Centred)
         else if(mode == "NUC") then
            call Sum_AB_NeutNuc(Reflex,Atm%Natoms,Grp%Centred)
         end if
       else
         call Sum_AB(Reflex,Atm%Natoms,Grp%Centred)
       end if


       return
    End Subroutine Modify_SF

    !!--++
    !!--++ Subroutine Set_Fixed_Tables(Reflex,Atm,Grp,mode,lambda,lun)
    !!--++    type(reflection_list_type),         intent(in) :: Reflex
    !!--++    type(atom_list_type),              intent(in) :: Atm
    !!--++    type(space_group_type),             intent(in) :: Grp
    !!--++    character(len=*), optional,         intent(in) :: Mode
    !!--++    real(kind=cp), optional,            intent(in) :: lambda
    !!--++    integer, optional,                  intent(in) :: lun
    !!--++
    !!--++    (Private)
    !!--++    Calculates arrays that are fixed during all further
    !!--++    calculations
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Set_Fixed_Tables(Reflex,Atm,Grp,Mode,lambda,lun)
       !---- Arguments ----!
       type(reflection_list_type),         intent(in) :: Reflex
       type(atom_list_type),               intent(in) :: Atm
       type(space_group_type),             intent(in) :: Grp
       character(len=*), optional,         intent(in) :: Mode
       real(kind=cp), optional,            intent(in) :: lambda
       integer, optional,                  intent(in) :: lun

       !---- Local variables ----!
       character(len=3) :: tipo

       tipo="XRA"
       if (present(mode)) tipo=adjustl(mode)
       tipo=U_Case(tipo)

       !---- Table HR - HT ----!
       call Create_Table_HR_HT(Reflex,Grp)

       !---- Table AF0 ----!
       select case (tipo)

          case ("XRA")
             if (present(lambda)) then
                if (present(lun)) then
                   call Create_Table_AF0_Xray(Reflex,Atm,lambda,lun)
                else
                   call Create_Table_AF0_Xray(Reflex,Atm,lambda)
                end if
             else
                if (present(lun)) then
                   call Create_Table_AF0_Xray(Reflex,Atm,lun=lun)
                else
                   call Create_Table_AF0_Xray(Reflex,Atm)
                end if
             end if

             !---- Modify the scattering factor tables to include the
             !---- multipliers factors concerning centre of symmetry and
             !---- centred translations
             if (Grp%Centred == 2) then
                af0=2.0*af0
                afpp=2.0*afpp
             end if

             if (Grp%NumLat  > 1) then
                af0=Grp%NumLat*af0
                afpp=Grp%NumLat*afpp
             end if

          case ("ELE")

             if (present(lun)) then
                call Create_Table_AF0_Electrons(Reflex,Atm,lun=lun)
             else
                call Create_Table_AF0_Electrons(Reflex,Atm)
             end if

             !---- Modify the scattering factor tables to include the
             !---- multipliers factors concerning centre of symmetry and
             !---- centred translations
             if (Grp%Centred == 2) then
                af0=2.0*af0
                afpp=0.0
             end if

             if (Grp%NumLat  > 1) then
                af0=Grp%NumLat*af0
                afpp=0.0
             end if

          case ("NUC","NEU")
             if (present(lun)) then
                call Create_Table_AFP_NeutNuc(Atm,lun=lun)
             else
                call Create_Table_AFP_NeutNuc(Atm)
             end if
             if (Grp%Centred == 2) afp=2.0*afp
             if (Grp%NumLat  > 1) afp=Grp%NumLat*afp

       end select

       return
    End Subroutine Set_Fixed_Tables

    !!----
    !!---- Subroutine Set_Form_Factors(Atm,Scf,ok,mess,lambda,lun,Add_Scatt)
    !!----   type(Atom_List_Type),                             intent(in out):: Atm
    !!----   Type(Scattering_Species_Type),                    intent(out)   :: Scf
    !!----   logical,                                          intent(out)   :: ok
    !!----   character(len=*),                                 intent(out)   :: mess
    !!----   real(kind=cp),                optional,           intent(in)    :: lambda
    !!----   integer,                      optional,           intent(in)    :: lun
    !!----   Type(Scattering_Species_Type),optional,           intent(in)    :: Add_Scatt
    !!----
    !!----  Constructor subroutine of object Scf of Scattering_Species_Type, by reading
    !!----  the database contained in module CFML_Scattering_Chemical_Tables and, in the
    !!----  appropriate case, the object Add_Scatt.
    !!----
    !!----  Created (JRC): October-2015
    !!----
    Subroutine Set_Form_Factors(Atm,Scf,ok,mess,lambda,lun,Add_Scatt,mag)
      type(Atom_List_Type),                    intent(in out):: Atm
      Type(Scattering_Species_Type),           intent(out)   :: Scf
      logical,                                 intent(out)   :: ok
      character(len=*),                        intent(out)   :: mess
      real(kind=cp),                optional,  intent(in)    :: lambda
      integer,                      optional,  intent(in)    :: lun
      Type(Scattering_Species_Type),optional,  intent(in)    :: Add_Scatt
      logical,                      optional,  intent(in)    :: mag

      !---- Local variables ----!
      character(len=4)                        :: symbcar
      integer                                 :: i,j,k,n,m,L
      character(len=4), dimension(atm%natoms) :: symb
      character(len=4), dimension(atm%natoms) :: elem
      real(kind=cp),    dimension(atm%natoms) :: bs
      integer,          dimension(atm%natoms) :: ix,jx
      real(kind=cp)                           :: b,dmin,d
      character(len=*), parameter             :: digpm="0123456789+-"
      logical                                 :: found

      call set_chem_info()
       !---- Getting Fermi Lengths of atoms ----!
       symb="    "; Elem="    "
       bs=0.0
       n=0
       ok=.true.
       do i=1,atm%natoms
          symbcar=u_case(atm%atom(i)%chemsymb)
          call Get_ChemSymb(symbcar, atm%atom(i)%chemsymb, atm%atom(i)%Z) !Getting the atomic number
          call Get_Fermi_Length(symbcar,b)                                !equal to the charge of the nuclei
          if (abs(b) < 0.0001) then
             ok=.false.
             Mess="The Fermi Length of Species "//symbcar//" was not found"
             return
          else
             if(any(Elem == symbcar)) cycle
             n=n+1
             bs(n) = b
             Elem(n)=u_case(atm%atom(i)%chemsymb)
             symb(n)=atm%atom(i)%SfacSymb
          end if
       end do
       call Remove_chem_info()

       Call Allocate_Scattering_Species(n,Scf)

       Scf%br=bs(1:n)

       if(present(lambda)) then
         call Set_Delta_Fp_Fpp()

          !---- Select wavelength (by default is CuKalpha1: k=5 in the list) ----!
          dmin=1000.0
          do i=1,5
             d=abs(lambda-Xray_Wavelengths(i)%Kalfa(1))
             if (d < dmin) then
                dmin=d
                k=i        !Selection of the index for fp and fpp lists
             end if
          end do

          !---- Found Species on Anomalous_ScFac ----!
          do i=1,Scf%Num_species
             symbcar=l_case(Elem(i))
             do j=1,Num_Delta_Fp
                if (symbcar /= Anomalous_ScFac(j)%Symb) cycle
                Scf%delta_fp(i)=Anomalous_ScFac(j)%fp(k)
                Scf%delta_fpp(i)=Anomalous_ScFac(j)%fpp(k)
                exit
             end do
          end do
         call Remove_Delta_Fp_Fpp()
       end if

       call Set_Xray_Form()

       !Look for X-ray scattering coefficients
       ix=0
       do i=1,Scf%Num_species
          symbcar=l_case(symb(i)) !Scattering factors
          k=index(symbcar,"+")
          j=index(symbcar,"-")
          if(k == 0 .and. j == 0) then     !Simple element or magnetic form Factor
            if(len_trim(symbcar) > 2) then !Magnetic form Factor -> use the chemical symbol
               symbcar=l_case(Elem(i))
               ix(i) = i  !Magnetic atom
            end if
          end if
          found=.false.
          do j=1,Num_Xray_Form
             if (symbcar /= Xray_form(j)%Symb) cycle
             Scf%xcoef(i)=Xray_form(j)
             found=.true.
             exit
          end do
          if(.not. found) then
            ok=.false.
            mess="Error: X-ray scattering form factor coefficients not found for "//symbcar
            return
          end if
       end do
       call Remove_Xray_Form()

       m=0
       do i=1,atm%natoms
         symbcar=u_case(atm%atom(i)%chemsymb)
         if(u_case(trim(atm%atom(i)%SfacSymb)) == "MPOL") m=m+1
         do j=1,Scf%Num_species
           if(symbcar == Elem(j)) then
             atm%atom(i)%ind(1)=j
             Scf%Xcoef(j)%Z=atm%atom(i)%Z
             Scf%Symb(j)=atm%atom(i)%SfacSymb
             exit
           end if
         end do
       end do

       if(present(Add_Scatt)) then
         if(Add_Scatt%Num_Species > 0) then
           do i=1,Add_Scatt%Num_Species
             do j=1,Scf%Num_species
                if(Scf%Symb(j) == Add_Scatt%Symb(i)) then
                  if(abs(Add_Scatt%br(i))> 0.00001)  Scf%br(j)=Add_Scatt%br(i)
                  if(abs(Add_Scatt%bi(i))> 0.00001)  Scf%bi(j)=Add_Scatt%bi(i)
                  if(abs(Add_Scatt%delta_fp(i))> 0.00001)  Scf%delta_fp(j) =Add_Scatt%delta_fp(i)
                  if(abs(Add_Scatt%delta_fpp(i))> 0.00001) Scf%delta_fpp(j)=Add_Scatt%delta_fpp(i)
                  if(abs(sum(Add_Scatt%Xcoef(i)%a))> 0.00001) then
                    Scf%Xcoef(j)%a=Add_Scatt%Xcoef(i)%a
                    Scf%Xcoef(j)%b=Add_Scatt%Xcoef(i)%b
                    Scf%Xcoef(j)%c=Add_Scatt%Xcoef(i)%c
                  end if
                end if
             end do
           end do
         end if
       end if

       if(present(mag)) then
         !---- Load form factor values for Magnetic Scattering ----!
         call Set_Magnetic_Form()
         !---- Find Species in Magnetic_Form ----!
         ix=0
         jx=0
         n=0
         do i=1,atm%natoms
            symbcar=atm%atom(i)%SfacSymb
            if(symbcar(1:1) /= "M" .and. symbcar(1:1) /= "J") cycle
            do j=1,num_mag_form
               if (symbcar /= Magnetic_Form(j)%Symb) cycle
               if(any(jx == j) ) exit
               n=n+1
               jx(n)=j
               ix(n)=i
               exit
            end do
         end do
         Scf%Num_magspc=n
         allocate(Scf%Mcoef(n),Scf%Symb_mag(n))
         do k=1,Scf%Num_magspc
            j = jx(k)
            i = ix(k)
            Scf%Mcoef(k)=Magnetic_Form(j)
            Scf%Symb_mag(k)= atm%atom(i)%SfacSymb
         end do
         do i=1,atm%natoms
           symbcar=u_case(atm%atom(i)%SfacSymb)
           do j=1,Scf%Num_magspc
             if(symbcar == Scf%Symb_mag(j)) then
               atm%atom(i)%ind(2)=j
               exit
             end if
           end do
         end do
       end if

       !---- Printing Information ----!
       if (present(lun)) then
          if(present(lambda)) then
            write(unit=lun,fmt="(/,a,f10.6,a)")  "  WAVELENGTH: ",lambda," Angstroms"
          else
            write(unit=lun,fmt="(/,a)")  "  WAVELENGTH NOT PROVIDED! "
          end if
          write(unit=lun,fmt="(/,a)")  "  INFORMATION FROM TABULATED NEUTRON SCATTERING FACTORS"
          write(unit=lun,fmt="(a,/)")  "  ==================================================="
          write(unit=lun,fmt="(a)")    "  FERMI LENGTHS "
          write(unit=lun,fmt="(a,i3)") "   Number of chemically different species: ",Scf%Num_Species
          write(unit=lun,fmt="(/,a)")  "   Atom     Fermi Length (Br,Bi)[10^(-12) cm]      Atomic Number"
          do k=1,Scf%Num_species
             write(unit=lun,fmt="(a,2F10.6,tr20,i8)")  "     "//Scf%Symb(k), Scf%br(k), Scf%bi(k), Scf%Xcoef(k)%Z
          end do
          if(.not. present(mag)) then
             write(unit=lun,fmt="(/,/)")
             write(unit=lun,fmt="(/,a)")  "  INFORMATION FROM TABULATED X-RAY SCATTERING FACTORS"
             write(unit=lun,fmt="(a,/)")  "  ==================================================="
             write(unit=lun,fmt="(/,a,/)")    "   ATOMIC SCATTERING FACTOR COEFFICIENTS: {A(i),B(i),I=1,4},C  Dfp  Dfpp "
             write(unit=lun,fmt="(a,i3)")     "   Number of chemically different species: ",Scf%Num_Species
             write(unit=lun,fmt="(/,a)") &
                  " Atom(ScF)  Atom(ScFX)           a1       b1       a2       b2       a3       b3       a4       b4        c      Dfp     Dfpp"
             do k=1,Scf%Num_species
                write(unit=lun,fmt="(a,11F9.5)")    &
                               "  "//Scf%Symb(k)//"        "//Scf%Xcoef(k)%Symb//"        ", &
                              (Scf%xcoef(k)%a(L),Scf%xcoef(k)%b(L), L=1,4), Scf%xcoef(k)%c, &
                              Scf%delta_fp(k), Scf%delta_fpp(k)
             end do
             write(unit=lun,fmt="(/,/)")

          else

            write(unit=lun,fmt="(/,a)")  "  INFORMATION FROM TABULATED MAGNETIC FORM FACTORS"
            write(unit=lun,fmt="(a,/)")  "  ================================================"
            write(unit=lun,fmt="(/,a,/)")    "   MAGNETIC FORM FACTOR COEFFICIENTS: {A(i),B(i),I=1,3},C  "
            write(unit=lun,fmt="(a,i3)")     "   Number of chemically different species: ",Scf%Num_magspc
            write(unit=lun,fmt="(/,a)") &
                 "     Atom     a1       b1       a2       b2       a3       b3       c      "
            do k=1,Scf%Num_magspc
               write(unit=lun,fmt="(a,7F9.5)")    &
                             "     "//Scf%Symb_mag(k), (Scf%Mcoef(k)%SctM(L), L=1,7)
            end do
            write(unit=lun,fmt="(/,/)")

          end if
       end if

    End Subroutine Set_Form_Factors
    !!--++
    !!--++ Subroutine Create_Table_mFR_Shub(Reflex,Atm,Lun)
    !!--++    type(Reflect_List_Type),intent(in) :: Reflex
    !!--++    type(Matom_list_type),  intent(in) :: Atm
    !!--++    integer, optional,      intent(in) :: lun
    !!--++
    !!--++    (Private)
    !!--++    Calculate a Table of Magnetic form Factors
    !!--..     mFR(Natoms,Nref)
    !!--++
    !!--++ Update: April - 2005
    !!
    !Subroutine Create_Table_mFR_Shub(Reflex,Atm,lun)
    !   !---- Arguments ----!
    !   type(Reflect_List_Type), intent(in) :: Reflex
    !   type(Atom_list_type),    intent(in) :: Atm
    !   integer, optional,       intent(in) :: lun
    !
    !   !---- Local Variables ----!
    !   character(len=4)               :: symbcar
    !   integer                        :: i,j, k,n,L
    !   integer, dimension(atm%natoms) :: ix,jx,ia
    !
    !   !---- Init ----!
    !   err_msfac=.false.
    !
    !   !---- Load form factor values for Magnetic Scattering ----!
    !   call Set_Magnetic_Form()
    !
    !   !---- Find Species in Magnetic_Form ----!
    !   ix=0
    !   jx=0
    !   n=0
    !   do i=1,atm%natoms
    !      symbcar=u_case(atm%atom(i)%SfacSymb)
    !      do j=1,num_mag_form
    !         if (symbcar /= Magnetic_Form(j)%Symb) cycle
    !         ix(i)=j
    !         if(any(jx == j) ) exit
    !         n=n+1
    !         jx(n)=j
    !         ia(n)=i
    !         exit
    !      end do
    !   end do
    !
    !   if (present(lun)) then
    !      write(unit=lun,fmt="(/,a)") "  INFORMATION FROM TABULATED MAGNETIC FORM FACTORS"
    !      write(unit=lun,fmt="(a,/)") "  ================================================"
    !   End if
    !
    !   if (any(ix==0)) then
    !      err_msfac=.true.
    !      err_msfac_mess="The Species "//symbcar//" was not found"
    !      return
    !   else
    !      !---- Fill mFR Table ----!
    !      do j=1,reflex%nref
    !         do i=1,atm%natoms
    !            mFR(i,j)=mfj(reflex%h(j)%s,Magnetic_Form(ix(i))%SctM)
    !         end do
    !      end do
    !   end if
    !
    !   !---- Printing Information ----!
    !   if (present(lun)) then
    !      write(unit=lun,fmt="(/,a,/)")    "   MAGNETIC FORM FACTOR COEFFICIENTS: {A(i),B(i),I=1,3},C  "
    !      write(unit=lun,fmt="(a,i3)")     "   Number of chemically different species: ",n
    !      write(unit=lun,fmt="(/,a)") &
    !           "     Atom     a1       b1       a2       b2       a3       b3       c      "
    !      do k=1,n
    !         j = jx(k)
    !         i = ia(k)
    !         write(unit=lun,fmt="(a,7F9.5)")    &
    !                       "     "//atm%atom(i)%SfacSymb, (Magnetic_Form(j)%SctM(L), L=1,7)
    !      end do
    !      write(unit=lun,fmt="(/,/)")
    !   end if
    !
    !   call Remove_Magnetic_Form()
    !
    !   return
    !End Subroutine Create_Table_mFR_Shub



    !!----
    !!---- Subroutine Structure_Factors(Atm,Grp,Reflex,Mode,lambda)
    !!----    type(atom_list_type),               intent(in)     :: Atm    !List of atoms
    !!----    type(space_group_type),             intent(in)     :: Grp    !Space group
    !!----    type(reflection_list_type),         intent(in out) :: Reflex !It is completed on output
    !!----    character(len=*), optional,         intent(in)     :: Mode   !"NUC","ELE" for neutrons, electrons else: XRays
    !!----    real(kind=cp), optional,            intent(in)     :: lambda !Needed for Xrays
    !!----
    !!----    Calculate the Structure Factors from a list of Atoms
    !!----    and a set of reflections. A call to Init_Structure_Factors
    !!----    is a pre-requisite for using this subroutine. In any case
    !!----    the subroutine calls Init_Structure_Factors if SF_initialized=.false.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Structure_Factors(Atm,Grp,Reflex,Mode,lambda)
       !---- Arguments ----!
       type(atom_list_type),               intent(in)     :: Atm
       type(space_group_type),             intent(in)     :: Grp
       type(reflection_list_type),         intent(in out) :: Reflex
       character(len=*), optional,         intent(in)     :: Mode
       real(kind=cp), optional,            intent(in)     :: lambda

       !Provisional items
       ! integer::i,j
       !---------------
       if(present(Mode)) then
          if(present(lambda)) then
            if(.not. SF_Initialized) call Init_Structure_Factors(Reflex,Atm,Grp,Mode,lambda)
          else
            if(.not. SF_Initialized) call Init_Structure_Factors(Reflex,Atm,Grp,Mode)
          end if
       else
          if(present(lambda)) then
            if(.not. SF_Initialized) call Init_Structure_Factors(Reflex,Atm,Grp,Lambda=lambda)
          else
            if(.not. SF_Initialized) call Init_Structure_Factors(Reflex,Atm,Grp)
          end if
       end if

       !---- Table TH ----!
       Call Calc_Table_TH(Reflex,Atm)

       !---- Table AB ----!
       call Calc_Table_AB(Reflex%Nref,Atm,Grp)

       !Provisional items
       !open(unit=111,file="stfac.inf",status="replace",action="write")
       !do j=1,Nref
       !  write(111,"(a,3i4)") " Reflection:  ",hkl(j)%h
       !
       !  write(111,"(a)") " Atom              F0         occ*W         Ajh         Bjh"
       !  do i=1,Atm%natoms
       !       write(111,"(a,4f12.4)") "  "//atm%atom(i)%lab, af0(i,j),th(i,j),Ajh(i,j),Bjh(i,j)
       !  end do ! Atoms
       !end do ! Reflections
       !close(unit=111)
       !End Provisional items

       !---- Final Calculation ----!
       if(present(mode)) then
         if(mode == "XRA" .or. mode == "ELE" ) then
            call Sum_AB(Reflex,Atm%Natoms,Grp%Centred)
         else if(mode == "NUC") then
            call Sum_AB_NeutNuc(Reflex,Atm%Natoms,Grp%Centred)
         else if(mode == "MAG") then
         end if
       else
         call Sum_AB(Reflex,Atm%Natoms,Grp%Centred)
       end if
       return
    End Subroutine Structure_Factors

    !!--++
    !!--++ Subroutine Sum_AB(Reflex,Natm,icent)
    !!--++    type(reflection_list_type), intent(in out) :: Reflex
    !!--++    integer,                    intent(in)     :: Natm
    !!--++    integer,                    intent(in)     :: icent
    !!--++
    !!--++    (Private)
    !!--++    Calculate the Final Sum for Structure Factors calculations
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Sum_AB(Reflex,Natm,icent)
       !---- Arguments ----!
       type(reflection_list_type), intent(in out)  :: Reflex
       integer,                    intent(in)      :: Natm
       integer,                    intent(in)      :: icent

       !---- Local Variables ----!
       integer                                     :: i,j
       real(kind=cp)                               :: a,b, ph
       real(kind=cp), dimension(natm,reflex%nref)  :: aa,bb,cc,dd


       ! A(h)=SIG(i){(f0+Deltaf')*OCC*Tiso*Ag}    asfa=a-d
       ! C(h)=SIG(i){    Deltaf" *OCC*Tiso*Ag}    bsfa=b+c

       ! B(h)=SIG(i){(f0+Deltaf')*OCC*Tiso*Bg}
       ! D(h)=SIG(i){    Deltaf" *OCC*Tiso*Bg}

       !---- Fj(h)*Aj(h) ----!

       aa=af0*th*ajh

       if (icent == 2) then    !Calculation for centrosymmetric structures
          do j=1,reflex%nref
             cc(:,j)= afpp(:)*th(:,j)*ajh(:,j)
          end do

          !---- Final Sum ----!
          do i=1,reflex%Nref
             a=sum(aa(:,i))
             b=sum(cc(:,i))
             reflex%ref(i)%Fc=sqrt(a*a+b*b)
             ph = atan2d(b,a)
             if (ph < 0.0) ph=ph+360.0
             reflex%ref(i)%Phase = ph
             reflex%ref(i)%A=a
             reflex%ref(i)%B=b
          end do

       else       !Calculation for non-centrosymmetric structures
          !---- Fj(h)*Bj(h) ----!
          bb=af0*th*bjh

          do j=1,reflex%nref
             cc(:,j)= afpp(:)*th(:,j)*ajh(:,j)
             dd(:,j)= afpp(:)*th(:,j)*bjh(:,j)
          end do

          !---- Final Sum ----!
          do i=1,reflex%Nref
             a=sum(aa(:,i)-dd(:,i))
             b=sum(bb(:,i)+cc(:,i))
             reflex%ref(i)%Fc=sqrt(a*a+b*b)
             ph = atan2d(b,a)
             if (ph < 0.0) ph=ph+360.0
             reflex%ref(i)%Phase = ph
             reflex%ref(i)%A=a
             reflex%ref(i)%B=b
          end do
       end if

       return
    End Subroutine Sum_AB

    !!--++
    !!--++ Subroutine Sum_AB_NeutNuc(Reflex,Natm,icent)
    !!--++    type(reflection_list_type),         intent(in out) :: Reflex
    !!--++    integer,                            intent(in)     :: Natm
    !!--++    integer,                            intent(in)     :: icent
    !!--++
    !!--++    (Private)
    !!--++    Calculate the Final Sum for Structure Factors calculations
    !!--++    Adapted for Neutron Nuclear Scattering (real scattering lengths)
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Sum_AB_NeutNuc(Reflex,Natm,icent)
       !---- Arguments ----!
       type(reflection_list_type),   intent(in out) :: Reflex
       integer,                      intent(in)     :: Natm
       integer,                      intent(in)     :: icent

       !---- Local Variables ----!
       integer                                     :: i,j
       real(kind=cp)                               :: a,b, ph
       real(kind=cp), dimension(natm,reflex%nref)  :: aa,bb

       if (icent == 2) then    !Calculation for centrosymmetric structures

          !---- Fj(h)*Aj(h) ----!
          do j=1,reflex%nref
             aa(:,j)= afp(:)*th(:,j)*ajh(:,j)
          end do

          !---- Final Sum ----!
          do i=1,reflex%Nref
             a=sum(aa(:,i))
             reflex%ref(i)%Fc=abs(a)
             reflex%ref(i)%Phase = 90.0_cp - 90.0_cp * sign(1.0_cp,a)
             reflex%ref(i)%A=a
             reflex%ref(i)%B=0.0
          end do

       else       !Calculation for non-centrosymmetric structures
          !---- Fj(h)*Bj(h) ----!
          !---- Fj(h)*Aj(h) ----!
          do j=1,reflex%nref
             aa(:,j)= afp(:)*th(:,j)*ajh(:,j)
             bb(:,j)= afp(:)*th(:,j)*bjh(:,j)
          end do

          !---- Final Sum ----!
          do i=1,reflex%Nref
             a=sum(aa(:,i))
             b=sum(bb(:,i))
             reflex%ref(i)%Fc=sqrt(a*a+b*b)
             ph = atan2d(b,a)
             if (ph < 0.0) ph=ph+360.0
             reflex%ref(i)%Phase = ph
             reflex%ref(i)%A=a
             reflex%ref(i)%B=b
          end do
       end if

       return
    End Subroutine Sum_AB_NeutNuc

    !!--++
    !!--++ Subroutine Write_Structure_Factors(lun,Reflex,Mode)
    !!--++    integer,                            intent(in) :: lun
    !!--++    type(reflection_list_type),         intent(in) :: Reflex
    !!--++    Character(len=*), optional,         intent(in) :: Mode
    !!--++
    !!--++    Writes in logical unit=lun the list of structure factors
    !!--++    contained in the array hkl
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Write_Structure_Factors_Crys(lun,Reflex,Mode)
       !---- Argument ----!
       integer,                            intent(in) :: lun
       type(reflection_list_type),         intent(in) :: Reflex
       Character(len=*), optional,         intent(in) :: Mode
       !---- Local Variables ----!
       integer :: i

       If(present(mode)) then
         Select Case (mode(1:3))
           Case("NUC","nuc")
             write(unit=lun,fmt="(/,/,a)") "    LIST OF REFLECTIONS AND STRUCTURE FACTORS(NEUTRONS)"
             write(unit=lun,fmt="(a)")     "    ==================================================="
           Case("ELE","ele")
             write(unit=lun,fmt="(/,/,a)") "    LIST OF REFLECTIONS AND STRUCTURE FACTORS(ELECTRONS)"
             write(unit=lun,fmt="(a)")     "    ===================================================="
           Case default
             write(unit=lun,fmt="(/,/,a)") "    LIST OF REFLECTIONS AND STRUCTURE FACTORS(X-RAYS)"
             write(unit=lun,fmt="(a)")     "    ================================================="
         End Select
       else
         write(unit=lun,fmt="(a)")   "    LIST OF REFLECTIONS AND STRUCTURE FACTORS(X-RAYS)"
         write(unit=lun,fmt="(a)")   "    ================================================="
       end if

       write(unit=lun,fmt="(/,a,/)") &
 "   H   K   L   Mult    SinTh/Lda       dspc          |Fc|         Phase          F-Real        F-Imag       |Fc|^2      Num"
       do i=1,reflex%Nref
             write(unit=lun,fmt="(3i4,i5,6f14.5,f14.3,i8)") reflex%ref(i)%h, reflex%ref(i)%mult, &
                                 reflex%ref(i)%S,0.5/reflex%ref(i)%S, reflex%ref(i)%Fc, reflex%ref(i)%Phase,   &
                                 reflex%ref(i)%a, reflex%ref(i)%b, reflex%ref(i)%Fc*reflex%ref(i)%Fc,i
       end do
       return
    End Subroutine Write_Structure_Factors_Crys

    !!--++
    !!--++ Subroutine Write_Structure_Factors_Mag(lun,Reflex,stf,full)
    !!--++    integer,                   intent(in) :: lun
    !!--++    type(reflect_list_type),intent(in) :: Reflex
    !!--++    type(Strf_List_Type),intent(in) :: Reflex
    !!--++    logical, optional,         intent(in) :: full
    !!--++
    !!--++    Writes in logical unit=lun the list of structure factors
    !!--++    calculated for a structure described using Shubnikov groups.
    !!--++    If "full" is present, full information contained in stf is output
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Write_Structure_Factors_Mag(lun,Reflex,stf,full)
       !---- Argument ----!
       integer,                 intent(in) :: lun
       type(reflect_list_type), intent(in) :: Reflex
       type(Strf_List_Type),    intent(in) :: stf
       logical, optional,       intent(in) :: full
       !---- Local Variables ----!
       integer :: i
       Character(len=*),dimension(0:2),parameter :: rtyp=(/"  Nuc  ","  Mag  ","Nuc+Mag"/)

       write(unit=lun,fmt="(/,/,a)") "    LIST OF REFLECTIONS AND STRUCTURE FACTORS(NUCLEAR and MAGNETIC)"
       write(unit=lun,fmt="(a)")     "    ==============================================================="

       write(unit=lun,fmt="(/,a,/)") "   H   K   L Mult SinTh/Lda d-spacing  ref-type       sqNuc       sqMiV  NumRef"
       do i=1,reflex%Nref
           write(unit=lun,fmt="(3i4,i5,2f10.5,tr3,a,2f12.5,i8)") reflex%ref(i)%h, reflex%ref(i)%mult, &
                               reflex%ref(i)%S,0.5/reflex%ref(i)%S,rtyp(reflex%ref(i)%imag), stf%strf(i)%sqNuc, stf%strf(i)%sqMiV,i
       end do

       if(present(full)) then
         write(unit=lun,fmt="(/,a,/)") &
"   H   K   L    Nuc-Real    Nuc-Imag  (MsFx-Real  MsFx-Imag)  (MsFy-Real  MsFy-Imag)  (MsFz-Real  MsFz-Imag)  (MiVx-Real  MiVx-Imag)  (MiVy-Real  MiVy-Imag)  (MiVz-Real  MiVz-Imag)  NumRef"
         do i=1,reflex%Nref
             write(unit=lun,fmt="(3i4,14f12.5,i8)") reflex%ref(i)%h,stf%strf(i)%NsF, stf%strf(i)%MsF, stf%strf(i)%MiV,i
         end do
       else
       end if
       return
    End Subroutine Write_Structure_Factors_Mag

 End Module CFML_Structure_Factors
