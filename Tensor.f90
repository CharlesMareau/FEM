! **********************************************************************
! *** module tensor                                                  ***
! *** contains the subroutines for tensor calculus                   ***
! **********************************************************************
module tensor
 use utils
 use iso_c_binding
 use OMP_LIB
 implicit none
#ifdef USE_MPI
 include "mpif.h"
#endif /* USE_MPI */
 private
 logical, private :: is_init_module_Tensor=.false.
 public :: init_module_Tensor
#ifdef USE_MPI
 integer, public :: MPI_TENSOR1
 integer, public :: MPI_TENSOR2
 integer, public :: MPI_TENSOR2S
 integer, public :: MPI_TENSOR2W
 integer, public :: MPI_TENSOR3
 integer, public :: MPI_TENSOR4S
#endif /* USE_MPI */
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! defines first order tensors                                       !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 type, public :: Tensor1
  real(8), public :: val(3) ! (x,y,z)
  contains
!   procedure, public :: rotate => rotate_Tensor1
   procedure, public :: writea => writea_Tensor1
   procedure, public :: writeb => writeb_Tensor1
   procedure, public :: readb  => readb_Tensor1
#ifndef COMPILER_LACKS_DERIVED_TYPE_IO
   generic :: write(FORMATTED) => writea
   generic :: write(UNFORMATTED) => writeb
   generic :: read (UNFORMATTED) => readb
#endif /* COMPILER_LACKS_DERIVED_TYPE_IO */
 end type Tensor1
 type, public :: VTensor1
  type(Tensor1), allocatable :: T(:)
 end type VTensor1
 type, public :: MTensor1
  type(Tensor1), allocatable :: T(:,:)
 end type MTensor1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! defines second order tensors                                      !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 type, public ::  Tensor2
  real(8), public :: val(3,3)
  contains
!!!   procedure, public :: I1     => I1_Tensor2
!!!   procedure, public :: J2     => J2_Tensor2
!!!   procedure, public :: vprint => writea_Tensor2
   procedure, public :: writea => writea_Tensor2
   procedure, public :: writeb => writeb_Tensor2
   procedure, public :: readb  => readb_Tensor2
#ifndef COMPILER_LACKS_DERIVED_TYPE_IO
   generic :: write(FORMATTED) => writea
   generic :: write(UNFORMATTED) => writeb
   generic :: read (UNFORMATTED) => readb
#endif /* COMPILER_LACKS_DERIVED_TYPE_IO */
 end type Tensor2
 type, public :: VTensor2
  type(Tensor2), allocatable :: T(:)
 end type VTensor2
 type, public :: MTensor2
  type(Tensor2), allocatable :: T(:,:)
 end type MTensor2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! defines symmetric second order tensors                            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 type, public ::  Tensor2s
  real(8), public :: val(6) ! VOIGT (11 22 33 23 31 12)
  contains
!!!   procedure, public  :: I1 => I1_Tensor2s
!!!   procedure, public  :: J2 => J2_Tensor2s
   procedure, public :: writea => writea_Tensor2s
   procedure, public :: writeb => writeb_Tensor2s
   procedure, public :: readb  => readb_Tensor2s
#ifndef COMPILER_LACKS_DERIVED_TYPE_IO
   generic :: write(FORMATTED) => writea
   generic :: write(UNFORMATTED) => writeb
   generic :: read (UNFORMATTED) => readb
#endif /* COMPILER_LACKS_DERIVED_TYPE_IO */
 end type Tensor2s
 type, public :: VTensor2s
  type(Tensor2s), allocatable :: T(:)
 end type VTensor2s
 type, public :: MTensor2s
  type(Tensor2s), allocatable :: T(:,:)
 end type MTensor2s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! defines skew-symmetric second order tensors                       !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 type, public ::  Tensor2w
  real(8), public :: val(3) ! (23 31 12)
  contains
!!!   procedure, public  :: I1 => I1_Tensor2w
!!!   procedure, public  :: J2 => J2_Tensor2w
   procedure, public :: writea => writea_Tensor2w
   procedure, public :: writeb => writeb_Tensor2w
   procedure, public :: readb  => readb_Tensor2w
#ifndef COMPILER_LACKS_DERIVED_TYPE_IO
   generic :: write(FORMATTED) => writea
   generic :: write(UNFORMATTED) => writeb
   generic :: read (UNFORMATTED) => readb
#endif /* COMPILER_LACKS_DERIVED_TYPE_IO */
 end type Tensor2w
 type, public :: VTensor2w
  type(Tensor2w), allocatable :: T(:)
 end type VTensor2w
 type, public :: MTensor2w
  type(Tensor2w), allocatable :: T(:,:)
 end type MTensor2w
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! defines third order tensors                                       !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 type, public ::  Tensor3
  real(8), public :: val(3,3,3)
  contains
   procedure, public :: writea => writea_Tensor3
   procedure, public :: writeb => writeb_Tensor3
   procedure, public :: readb  => readb_Tensor3
#ifndef COMPILER_LACKS_DERIVED_TYPE_IO
   generic :: write(FORMATTED) => writea
   generic :: write(UNFORMATTED) => writeb
   generic :: read (UNFORMATTED) => readb
#endif /* COMPILER_LACKS_DERIVED_TYPE_IO */
 end type Tensor3
 type, public :: VTensor3
  type(Tensor3), allocatable :: T(:)
 end type VTensor3
 type, public :: MTensor3
  type(Tensor3), allocatable :: T(:,:)
 end type MTensor3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! defines fourth order tensors                                      !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 type, public ::  Tensor4 !!! CM Définition du type tenseur d'ordre 4 général
  real(8), public :: val(3,3,3,3)
  contains
!!!   procedure, public :: rotate => rotate_Tensor4
   procedure, public :: writea => writea_Tensor4
   procedure, public :: writeb => writeb_Tensor4
   procedure, public :: readb  => readb_Tensor4
#ifndef COMPILER_LACKS_DERIVED_TYPE_IO
!!!   generic :: write(FORMATTED) => writea
!!!   generic :: write(UNFORMATTED) => writeb
!!!   generic :: read (UNFORMATTED) => readb
#endif /* COMPILER_LACKS_DERIVED_TYPE_IO */
 end type Tensor4
 type, public :: VTensor4
  type(Tensor4), allocatable :: T(:)
 end type VTensor4
 type, public :: MTensor4
  type(Tensor4), allocatable :: T(:,:)
 end type MTensor4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! defines symmetric fourth order tensors                            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 type, public ::  Tensor4ss !!! CM Redéfinition du type tenseur d'ordre 4 sym-sym
  real(8), public :: val(6,6)
  contains
!!!   procedure, public :: rotate => rotate_Tensor4ss
   procedure, public :: writea => writea_Tensor4ss
   procedure, public :: writeb => writeb_Tensor4ss
   procedure, public :: readb  => readb_Tensor4ss
#ifndef COMPILER_LACKS_DERIVED_TYPE_IO
   generic :: write(FORMATTED) => writea
   generic :: write(UNFORMATTED) => writeb
   generic :: read (UNFORMATTED) => readb
#endif /* COMPILER_LACKS_DERIVED_TYPE_IO */
 end type Tensor4ss
 type, public :: VTensor4ss
  type(Tensor4ss), allocatable :: T(:)
 end type VTensor4ss
 type, public :: MTensor4ss
  type(Tensor4ss), allocatable :: T(:,:)
 end type MTensor4ss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! defines interfaces for determination of eigen values and eigen    !!!
!!! vectors                                                           !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 public :: eigen
 interface eigen
  module procedure eigen2
  module procedure eigen2s
 end interface eigen
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! defines interfaces for tensor assignments                         !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 public :: assignment(=)
 interface assignment(=)
  module procedure equal_T1T1
  module procedure equal_T2T2
  module procedure equal_T2T2s
  module procedure equal_T2sT2
  module procedure equal_T2sT2s
  module procedure equal_T2T2w
  module procedure equal_T2wT2
  module procedure equal_T2wT2w
  module procedure equal_T3T3
  module procedure equal_T4T4
  module procedure equal_T4ssT4
  module procedure equal_T4T4ss
  module procedure equal_T4ssT4ss
 end interface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! defines interfaces for tensor conversion                          !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 public :: full
 interface full
  module procedure full_Tensor2s
  module procedure full_Tensor2w
  module procedure full_Tensor4ss
 end interface full
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! defines interfaces for tensor norm calculation                    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 public :: nrm
 interface nrm
  module procedure nrm_Tensor1
  module procedure nrm_Tensor2
  module procedure nrm_Tensor2s
  module procedure nrm_Tensor2w
  module procedure nrm_Tensor3
  module procedure nrm_Tensor4
  module procedure nrm_Tensor4ss
 end interface nrm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! defines interfaces for extracting symmetric, skew-symmetric       !!!
!!! deviatoric and spherical parts of second order tensors            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 public :: sym
 interface sym
  module procedure sym_Tensor2
 end interface sym
 public :: skw
 interface skw
  module procedure skw_Tensor2
 end interface skw
 public :: dev
 interface dev
  module procedure dev_Tensor2,dev_Tensor2s !!! CM dev_Tensor2w supprimé car le deviateur de Tensor2w est lui meme
 end interface dev
 public :: sph
 interface sph
  module procedure sph_Tensor2,sph_Tensor2s !!! CM ajout de la partie spherique
 end interface sph
 public :: trc
 interface trc
  module procedure trc_Tensor2,trc_Tensor2s !!! CM ajout de la trace
 end interface trc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! defines interfaces for extracting symmetric and skew-symmetric    !!!
!!! of fourth order tensors                                           !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 public :: symsym
 interface symsym
  module procedure symsym_Tensor4
 end interface symsym
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! defines interfaces for tensor determinant calculation             !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 public :: det
 interface det
  module procedure det_Tensor2
  module procedure det_Tensor2s !,inv_Tensor2s,inv_Tensor2w,inv_Tensor4ss
 end interface det
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! defines interfaces for tensor inversion                           !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 public :: inv
 interface inv
  module procedure inv_Tensor2
  module procedure inv_Tensor2s
  module procedure inv_Tensor4
  module procedure inv_Tensor4ss!,inv_Tensor2s,inv_Tensor2w,inv_Tensor4ss
 end interface inv
 public :: psdinv
 interface psdinv
  module procedure psdinv_Tensor2
  module procedure psdinv_Tensor4
 end interface psdinv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! defines interfaces for tensor transposition                       !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 public :: trsp
 interface trsp
  module procedure transpose_Tensor2
!!!  module procedure transpose_Tensor2s !!! CM à enlever car inutile pour un tenseur symétrique
  module procedure transpose_Tensor2w
  module procedure transpose_Tensor4
  module procedure transpose_Tensor4ss
 end interface
 public :: ltrsp
 interface ltrsp
  module procedure lefttranspose_Tensor4
 end interface
 public :: rtrsp
 interface rtrsp
  module procedure righttranspose_Tensor4
 end interface
 public :: mtrsp
 interface mtrsp
  module procedure middletranspose_Tensor4
  module procedure middletranspose_Tensor4ss
 end interface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! defines interfaces for tensor addition                            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 public :: operator(+)
 interface operator(+)
  module procedure plus_T1T1
  module procedure plus_T2T2
  module procedure plus_T2T2s
  module procedure plus_T2sT2
  module procedure plus_T2sT2s
  module procedure plus_T2T2w
  module procedure plus_T2wT2
  module procedure plus_T2wT2w
  module procedure plus_T2wT2s
  module procedure plus_T2sT2w
  module procedure plus_T3T3
  module procedure plus_T4T4
  module procedure plus_T4ssT4ss
  module procedure plus_T4T4ss
  module procedure plus_T4ssT4
 end interface
!  public :: add
!  interface add !!! CM Ajout interface add
!   module procedure plus_T1T1
!   module procedure plus_T2T2
!   module procedure plus_T2T2s
!   module procedure plus_T2sT2
!   module procedure plus_T2sT2s
!   module procedure plus_T2T2w
!   module procedure plus_T2wT2
!   module procedure plus_T2wT2w
!   module procedure plus_T2wT2s
!   module procedure plus_T2sT2w
!   module procedure plus_T3T3
!   module procedure plus_T4T4
!   module procedure plus_T4ssT4ss
!   module procedure plus_T4T4ss
!   module procedure plus_T4ssT4
!  end interface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! defines interfaces for tensor subtraction                         !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 public :: operator(-)
 interface operator(-)
  module procedure minus_T1T1
  module procedure minus_T2T2
  module procedure minus_T2T2s
  module procedure minus_T2sT2
  module procedure minus_T2sT2s
  module procedure minus_T2T2w
  module procedure minus_T2wT2
  module procedure minus_T2wT2w
  module procedure minus_T2wT2s
  module procedure minus_T2sT2w
  module procedure minus_T3T3
  module procedure minus_T4T4
  module procedure minus_T4ssT4ss
  module procedure minus_T4T4ss
  module procedure minus_T4ssT4
 end interface
!  public :: subtract !!! CM Ajout interface subtract
!  interface subtract
!   module procedure minus_T1T1
!   module procedure minus_T2T2
!   module procedure minus_T2T2s
!   module procedure minus_T2sT2
!   module procedure minus_T2sT2s
!   module procedure minus_T2T2w
!   module procedure minus_T2wT2
!   module procedure minus_T2wT2w
!   module procedure minus_T2wT2s
!   module procedure minus_T2sT2w
!   module procedure minus_T3T3
!   module procedure minus_T4T4
!   module procedure minus_T4ssT4ss
!   module procedure minus_T4T4ss
!   module procedure minus_T4ssT4
!  end interface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! defines interfaces for tensor scaling                             !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 public :: operator(*)
 interface operator(*)
  module procedure product_T0T1
  module procedure product_T1T0
  module procedure product_T0T2
  module procedure product_T2T0
  module procedure product_T0T2s
  module procedure product_T2sT0
  module procedure product_T0T2w
  module procedure product_T2wT0
  module procedure product_T0T3
  module procedure product_T3T0
  module procedure product_T0T4
  module procedure product_T4T0
  module procedure product_T0T4ss
  module procedure product_T4ssT0
  module procedure tensor_product_T1T1
  module procedure tensor_product_T1T2
  module procedure tensor_product_T2T1
  module procedure tensor_product_T2T2
  module procedure tensor_product_T2sT2s
 end interface
 public :: operator(.times.)
 interface operator(.times.)
  module procedure times_T2T2
 end interface
 public :: operator(.bartimes.)
 interface operator(.bartimes.)
  module procedure bartimes_T2T2
  module procedure bartimes_T2T2s
  module procedure bartimes_T2sT2
  module procedure bartimes_T2sT2s
 end interface
 public :: operator(.timesbar.)
 interface operator(.timesbar.)
  module procedure timesbar_T2T2
 end interface
 public :: operator(.bartimesbar.)
 interface operator(.bartimesbar.)
  module procedure bartimesbar_T2T2
 end interface
 
!  public :: multiply  !!! CM Ajout interface multiply
!  interface multiply
!   module procedure product_T0T1
!   module procedure product_T1T0
!   module procedure product_T0T2
!   module procedure product_T2T0
!   module procedure product_T0T2s
!   module procedure product_T2sT0
!   module procedure product_T0T2w
!   module procedure product_T2wT0
!   module procedure product_T0T3
!   module procedure product_T3T0
!   module procedure product_T0T4
!   module procedure product_T4T0
!   module procedure product_T0T4ss
!   module procedure product_T4ssT0
!  end interface
 public :: operator(/)
 interface operator(/)
  module procedure divide_T1T0
  module procedure divide_T2T0
  module procedure divide_T2sT0
  module procedure divide_T2wT0
  module procedure divide_T3T0
  module procedure divide_T4T0
  module procedure divide_T4ssT0
 end interface
!  public :: divide  !!! CM Ajout interface divide
!  interface divide
!   module procedure divide_T1T0
!   module procedure divide_T2T0
!   module procedure divide_T2sT0
!   module procedure divide_T2wT0
!   module procedure divide_T3T0
!   module procedure divide_T4T0
!   module procedure divide_T4ssT0
!  end interface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! defines interfaces for tensor, dot and cross products             !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  public :: tensor_product  !!! CM Ajout interface tensor_product
!  interface tensor_product
!   module procedure tensor_product_T1T1
!   module procedure tensor_product_T1T2
!   module procedure tensor_product_T2T1
!   module procedure tensor_product_T2T2
!   module procedure tensor_product_T2sT2s
!  end interface
!  public :: dot_product
!  interface dot_product
!   module procedure dot_T1T1
!   module procedure dot_T2T1
!   module procedure dot_T1T2
!   module procedure dot_T2sT1
!   module procedure dot_T1T2s
!   module procedure dot_T2wT1
!   module procedure dot_T1T2w
!   module procedure dot_T2T2
!   module procedure dot_T2T2s
!   module procedure dot_T2sT2
!   module procedure dot_T2sT2s
!   module procedure dot_T2T2w
!   module procedure dot_T2wT2
!   module procedure dot_T2wT2w
!   module procedure dot_T2sT2w
!   module procedure dot_T2wT2s
!   module procedure dot_T2T4
!   module procedure dot_T4T2
!   module procedure dot_T2T4ss
!   module procedure dot_T4ssT2
!  end interface dot_product
 public :: operator(.dot.)
 interface operator(.dot.)
  module procedure dot_T1T1
  module procedure dot_T2T1
  module procedure dot_T1T2
  module procedure dot_T2sT1
  module procedure dot_T1T2s
  module procedure dot_T2wT1
  module procedure dot_T1T2w
  module procedure dot_T2T2
  module procedure dot_T2T2s
  module procedure dot_T2sT2
  module procedure dot_T2sT2s
  module procedure dot_T2T2w
  module procedure dot_T2wT2
  module procedure dot_T2wT2w
  module procedure dot_T2sT2w
  module procedure dot_T2wT2s
  module procedure dot_T2T4
  module procedure dot_T4T2
  module procedure dot_T2sT4ss
  module procedure dot_T4ssT2s
  module procedure dot_T2T4ss
  module procedure dot_T4ssT2
 end interface
!  public :: ddot_product
!  interface ddot_product
!   module procedure ddot_T2T2
!   module procedure ddot_T2T2s
!   module procedure ddot_T2sT2
!   module procedure ddot_T2T2w
!   module procedure ddot_T2wT2
!   module procedure ddot_T2sT2s
!   module procedure ddot_T2wT2w
!   module procedure ddot_T2sT2w
!   module procedure ddot_T2wT2s
!   module procedure ddot_T4T2
!   module procedure ddot_T2T4
!   module procedure ddot_T4T2s
!   module procedure ddot_T2sT4
!   module procedure ddot_T4ssT2s
!   module procedure ddot_T2sT4ss
!   module procedure ddot_T4T2w
!   module procedure ddot_T2wT4
!   module procedure ddot_T4T4
!   module procedure ddot_T4T4ss
!   module procedure ddot_T4ssT4
!   module procedure ddot_T4ssT4ss
!  end interface ddot_product
 public :: operator(.ddot.)
 interface operator(.ddot.)
  module procedure ddot_T2T2
  module procedure ddot_T2T2s
  module procedure ddot_T2sT2
  module procedure ddot_T2T2w
  module procedure ddot_T2wT2
  module procedure ddot_T2sT2s
  module procedure ddot_T2wT2w
  module procedure ddot_T2sT2w
  module procedure ddot_T2wT2s
  module procedure ddot_T4T2
  module procedure ddot_T2T4
  module procedure ddot_T4T2s
  module procedure ddot_T2sT4
  module procedure ddot_T4ssT2s
  module procedure ddot_T2sT4ss
  module procedure ddot_T4T2w
  module procedure ddot_T2wT4
  module procedure ddot_T4T4
  module procedure ddot_T4T4ss
  module procedure ddot_T4ssT4
  module procedure ddot_T4ssT4ss
 end interface
 public :: operator(.ddddot.)
 interface operator(.ddddot.)
  module procedure ddddot_T4T4
  module procedure ddddot_T4ssT4ss
 end interface
 public :: cross_product
 interface cross_product
  module procedure cross_product_T1T1
!!!  module procedure cross_product_Tensor2_Tensor2
 end interface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! defines interfaces for tensor tranformations                      !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 public :: rotate
 interface rotate
  module procedure rotate_Tensor1
  module procedure rotate_Tensor2
  module procedure rotate_Tensor2s
  module procedure rotate_Tensor2w
  module procedure rotate_Tensor3
  module procedure rotate_Tensor4
  module procedure rotate_Tensor4ss
 end interface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! defines interfaces for tensor weighted average                    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 public :: average
 interface average
  module procedure average_Tensor1
  module procedure average_Tensor2
  module procedure average_Tensor2s
  module procedure average_Tensor2w
  module procedure average_Tensor3
  module procedure average_Tensor4
  module procedure average_Tensor4ss
 end interface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! defines interfaces for isotropic functions of tensors             !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 public :: exp
 interface exp
  module procedure exp_Tensor2
  module procedure exp_Tensor2s
  module procedure exp_Tensor2w
 end interface
 public :: log
 interface log
  module procedure log_Tensor2
  module procedure log_Tensor2s
 end interface
 public :: sqrt
 interface sqrt
  module procedure sqrt_Tensor2s
  module procedure sqrt_Tensor4ss
 end interface
 public :: pos
 interface pos
  module procedure pos_Tensor2s
 end interface
 public :: plrdcmp
 public :: isotropic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! defines specific parameter arrays for tensor operations           !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 integer, public, parameter :: ijndx(9,2) =reshape([1,2,3,1,2,3,1,2,3,1,1,1,2,2,2,3,3,3],(/9,2/))
 integer, public, parameter :: ijndxs(6,2) =reshape([1,2,3,2,3,1,1,2,3,3,1,2],(/6,2/)) 
 integer, public, parameter :: ijndxw(3,2) =reshape([2,3,1,3,1,2],(/3,2/))
 integer, public, parameter :: indx(3,3) =reshape([1,2,3,4,5,6,7,8,9],(/3,3/))
 integer, public, parameter :: indxs(3,3) =reshape([1, 6, 5, 6, 2, 4, 5, 4, 3],(/3,3/))
 integer, public, parameter :: indxw(3,3) =reshape([1, 3, 2, 3, 2, 1, 2, 1, 3],(/3,3/))
 real(8), public, parameter :: sfac(6)=[1.d+0, 1.d+0, 1.d+0, 2.d+0, 2.d+0, 2.d+0]
 real(8), public, parameter :: wfac(3)=[2.d+0, 2.d+0, 2.d+0]
 real(8), parameter, public :: pfac(3,3) =reshape([0.d0,-1.d0, 1.d0, 1.d0, 0.d0,-1.d0,-1.d0, 1.d0, 0.d0],(/3,3/))
 real(8), parameter, public :: invfac(6,6)=reshape([1.0d+0, 1.0d+0, 1.0d+0, 5.0d-1, 5.0d-1, 5.0d-1, &
                                                    1.0d+0, 1.0d+0, 1.0d+0, 5.0d-1, 5.0d-1, 5.0d-1, &
                                                    1.0d+0, 1.0d+0, 1.0d+0, 5.0d-1, 5.0d-1, 5.0d-1, &
                                                    5.0d-1, 5.0d-1, 5.0d-1, 2.5d-1, 2.5d-1, 2.5d-1, &
                                                    5.0d-1, 5.0d-1, 5.0d-1, 2.5d-1, 2.5d-1, 2.5d-1, &
                                                    5.0d-1, 5.0d-1, 5.0d-1, 2.5d-1, 2.5d-1, 2.5d-1],(/6,6/))
 real(8), parameter, public :: perm(3,3,3)=reshape([0.d+0, 0.d+0, 0.d+0, 0.d0, 0.d0,-1.d0, 0.d0, 1.d0, 0.d0, &
                                                    0.d+0, 0.d+0, 1.d+0, 0.d0, 0.d0, 0.d0,-1.d0, 0.d0, 0.d0, &
                                                    0.d+0,-1.d+0, 0.d+0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0],(/3,3,3/))
 type(tensor1), public, parameter :: null1=tensor1(val=[0.d0,0.d0,0.d0])                                                   
!  type(tensor2), public, parameter :: null2=tensor2(val=reshape([0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0],(/3,3/)))    
 type(tensor2), public, parameter :: null2=tensor2(val=[0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0])
!  type(tensor2), public, parameter :: id2=tensor2(val=reshape([1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0],(/3,3/)))
 type(tensor2), public, parameter :: id2=tensor2(val=[1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0])
 type(tensor2), public, parameter :: one2=tensor2(val=[1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0])
 type(tensor2s), public, parameter :: null2s=tensor2s(val=[0.d0,0.d0,0.d0,0.d0,0.d0,0.d0])
 type(tensor2s), public, parameter :: id2s=tensor2s(val=[1.d0,1.d0,1.d0,0.d0,0.d0,0.d0])
 type(tensor2w), public, parameter :: null2w=tensor2w(val=[0.d0,0.d0,0.d0])
!  type(tensor3), public, parameter :: null3=tensor3(val=reshape([0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
!                                                                ,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
!                                                                ,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0],(/3,3,3/)))
 type(tensor3), public, parameter :: null3=tensor3(val=[0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
                                                       ,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
                                                       ,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0])                                                               
!  type(tensor4), public, parameter :: null4=tensor4(val=reshape([0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
!                                                                ,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
!                                                                ,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
!                                                                ,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
!                                                                ,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
!                                                                ,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
!                                                                ,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
!                                                                ,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
!                                                                ,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0],(/3,3,3,3/)))
 type(tensor4), public, parameter :: null4=tensor4(val=[0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
                                                       ,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
                                                       ,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
                                                       ,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
                                                       ,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
                                                       ,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
                                                       ,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
                                                       ,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
                                                       ,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0])                                                               
!  type(tensor4), public, parameter :: id4=tensor4(val=reshape([1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
!                                                              ,0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
!                                                              ,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
!                                                              ,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
!                                                              ,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,0.d0 &
!                                                              ,0.d0,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0 &
!                                                              ,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0 &
!                                                              ,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0 &
!                                                              ,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,1.d0],(/3,3,3,3/)))
 type(tensor4), public, parameter :: id4=tensor4(val=[1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
                                                     ,0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
                                                     ,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
                                                     ,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0 &
                                                     ,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,0.d0 &
                                                     ,0.d0,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0 &
                                                     ,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0 &
                                                     ,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0 &
                                                     ,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,1.d0])                                                             
!  type(tensor4ss), public, parameter :: null4ss=tensor4ss(val=reshape([0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0 &
!                                                                      ,0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0 &
!                                                                      ,0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0 &
!                                                                      ,0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0 &
!                                                                      ,0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0 &
!                                                                      ,0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0],(/6,6/)))
 type(tensor4ss), public, parameter :: null4ss=tensor4ss(val=[0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0 &
                                                             ,0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0 &
                                                             ,0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0 &
                                                             ,0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0 &
                                                             ,0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0 &
                                                             ,0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0])
!  type(tensor4ss), public, parameter :: id4ss=tensor4ss(val=reshape([1.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0 &
!                                                                    ,0.0d+0, 1.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0 &
!                                                                    ,0.0d+0, 0.0d+0, 1.0d+0, 0.0d+0, 0.0d+0, 0.0d+0 &
!                                                                    ,0.0d+0, 0.0d+0, 0.0d+0, 0.5d+0, 0.0d+0, 0.0d+0 &
!                                                                    ,0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.5d+0, 0.0d+0 &
!                                                                    ,0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.5d+0],(/6,6/)))
 type(tensor4ss), public, parameter :: id4ss=tensor4ss(val=[1.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0 &
                                                           ,0.0d+0, 1.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0 &
                                                           ,0.0d+0, 0.0d+0, 1.0d+0, 0.0d+0, 0.0d+0, 0.0d+0 &
                                                           ,0.0d+0, 0.0d+0, 0.0d+0, 0.5d+0, 0.0d+0, 0.0d+0 &
                                                           ,0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.5d+0, 0.0d+0 &
                                                           ,0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.5d+0])                                                                   
!  type(tensor4ss), public, parameter :: pd4ss=tensor4ss(val=reshape([2.d0/3.d0,-1.d0/3.d0,-1.d0/3.d0, 0.0d+0, 0.0d+0, 0.0d+0 &
!                                                                   ,-1.d0/3.d0, 2.d0/3.d0,-1.d0/3.d0, 0.0d+0, 0.0d+0, 0.0d+0 &
!                                                                   ,-1.d0/3.d0,-1.d0/3.d0, 2.d0/3.d0, 0.0d+0, 0.0d+0, 0.0d+0 &
!                                                                   ,0.0d+0, 0.0d+0, 0.0d+0, 0.5d+0, 0.0d+0, 0.0d+0 &
!                                                                   ,0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.5d+0, 0.0d+0 &
!                                                                   ,0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.5d+0],(/6,6/)))   
 type(tensor4ss), public, parameter :: pd4ss=tensor4ss(val=[2.d0/3.d0,-1.d0/3.d0,-1.d0/3.d0, 0.0d+0, 0.0d+0, 0.0d+0 &
                                                           ,-1.d0/3.d0, 2.d0/3.d0,-1.d0/3.d0, 0.0d+0, 0.0d+0, 0.0d+0 &
                                                           ,-1.d0/3.d0,-1.d0/3.d0, 2.d0/3.d0, 0.0d+0, 0.0d+0, 0.0d+0 &
                                                           ,0.0d+0, 0.0d+0, 0.0d+0, 0.5d+0, 0.0d+0, 0.0d+0 &
                                                           ,0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.5d+0, 0.0d+0 &
                                                           ,0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.5d+0])             
!  type(tensor4ss), public, parameter :: ps4ss=tensor4ss(val=reshape([1.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0, 0.0d+0, 0.0d+0, 0.0d+0 &
!                                                                    ,1.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0, 0.0d+0, 0.0d+0, 0.0d+0 &
!                                                                    ,1.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0, 0.0d+0, 0.0d+0, 0.0d+0 &
!                                                                    ,0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0 &
!                                                                    ,0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0 &
!                                                                    ,0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0],(/6,6/)))  
 type(tensor4ss), public, parameter :: ps4ss=tensor4ss(val=[1.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0, 0.0d+0, 0.0d+0, 0.0d+0 &
                                                           ,1.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0, 0.0d+0, 0.0d+0, 0.0d+0 &
                                                           ,1.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0, 0.0d+0, 0.0d+0, 0.0d+0 &
                                                           ,0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0 &
                                                           ,0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0 &
                                                           ,0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0])                                                                     
 contains
! **********************************************************************
! *** subroutine init_module_Tensor                                  ***
! *** does what ??? !!! CM                                           ***
! ***                                                                ***
! *** input :                                                        ***
! ***                                                                ***
! *** output :                                                       ***
! **********************************************************************
subroutine init_module_Tensor()
 implicit none
#ifdef USE_MPI
 integer :: code
#endif /* USE_MPI */

#ifdef USE_MPI
 if (.not.is_init_module_Tensor) then
  call MPI_TYPE_CONTIGUOUS(3,MPI_DOUBLE,MPI_TENSOR1,code)
  call MPI_TYPE_COMMIT(MPI_TENSOR1,code)
  call MPI_TYPE_CONTIGUOUS(9,MPI_DOUBLE,MPI_TENSOR2,code)
  call MPI_TYPE_COMMIT(MPI_TENSOR2,code)
  call MPI_TYPE_CONTIGUOUS(6,MPI_DOUBLE,MPI_TENSOR2S,code)
  call MPI_TYPE_COMMIT(MPI_TENSOR2S,code)
  call MPI_TYPE_CONTIGUOUS(3,MPI_DOUBLE,MPI_TENSOR2W,code)
  call MPI_TYPE_COMMIT(MPI_TENSOR2W,code)
  call MPI_TYPE_CONTIGUOUS(27,MPI_DOUBLE,MPI_TENSOR3,code)
  call MPI_TYPE_COMMIT(MPI_TENSOR3,code)
  call MPI_TYPE_CONTIGUOUS(36,MPI_DOUBLE,MPI_TENSOR4S,code)
  call MPI_TYPE_COMMIT(MPI_TENSOR4S,code)
 end if
#endif /* USE_MPI */
 is_init_module_Tensor=.true.
end  subroutine init_module_Tensor
! **********************************************************************
! *** subroutine writea_Tensor1                                      ***
! *** writes a first order tensor in a text file                     ***
! ***                                                                ***
! *** input : T, first order tensor to be written                    ***
! ***         u, unit to write in                                    ***
! ***         io, format of the output                               ***
! ***         list, !!! CM non utilisé ???                           ***
! ***         eta, !!! utilité ???                                   ***
! ***                                                                ***
! *** in/output : msg, message to be passed                          ***
! ***                                                                ***
! *** output :                                                       ***
! **********************************************************************
subroutine writea_Tensor1(T,u,io,list,eta,msg)
 implicit none
 class(Tensor1), intent(in) :: T
 integer, intent(in) :: u
 character(len=*), intent(in) :: io
 integer, dimension(:), intent(in) :: list
 integer, intent(out) :: eta
 character(len=*), intent(inout) :: msg
 integer :: i

 select case(io)
  case ("DTVERBOSE")
   do i=1,3
    write(u,'(a,1e12.4,a,a)',advance="no") '|', T%val(i),'|',C_NEW_LINE
   end do
  case default
  write(u,'(3e12.4)') T%val !!! CM, Attention à l'ordre d'affichage
 end select
 eta=0

end subroutine writea_Tensor1
! **********************************************************************
! *** subroutine writeb_Tensor1                                      ***
! *** writes a first order tensor in a binary file                   ***
! ***                                                                ***
! *** input : T, first order tensor to be written                    ***
! ***         u, unit to write in                                    ***
! ***         eta, !!! utilité ???                                   ***
! ***                                                                ***
! *** in/output : msg, message to be passed                          ***
! ***                                                                ***
! *** output :                                                       ***
! **********************************************************************
subroutine writeb_Tensor1(T,u,eta,msg)
 implicit none
 class(Tensor1), intent(in) :: T
 integer, intent(in) :: u
 integer, intent(out) :: eta
 character(len=*), intent(inout) :: msg
 integer :: i
 write(u) T%val(1:3)
 eta=0
end subroutine writeb_Tensor1
! **********************************************************************
! *** subroutine readb_Tensor1                                       ***
! *** reads a first order tensor from a binary file                  ***
! ***                                                                ***
! *** input : u, unit to read in                                     ***
! ***         eta, !!! utilité ???                                   ***
! ***                                                                ***
! *** in/output : msg, message to be passed                          ***
! ***                                                                ***
! *** output : T, first order tensor to be read                      ***
! **********************************************************************
subroutine readb_Tensor1(T,u,eta,msg)
 implicit none
 class(Tensor1), intent(inout) :: T !!! CM pourquoi inout et pas out
 integer, intent(in) :: u
 integer, intent(out) :: eta
 character(len=*), intent(inout) :: msg
 integer :: i
 read(unit=u,iostat=eta,iomsg=msg) T%val
 eta=0
end subroutine readb_Tensor1
! **********************************************************************
! *** subroutine writea_Tensor2                                      ***
! *** writes a second order tensor in a text file                    ***
! ***                                                                ***
! *** input : T, second order tensor to be written                   ***
! ***         u, unit to write in                                    ***
! ***         io, format of the output                               ***
! ***         list, !!! CM non utilisé ???                           ***
! ***         eta, !!! utilité ???                                   ***
! ***                                                                ***
! *** in/output : msg, message to be passed                          ***
! ***                                                                ***
! *** output :                                                       ***
! **********************************************************************
subroutine writea_Tensor2(T,u,io,list,eta,msg)
 implicit none
 class(Tensor2), intent(in) :: T
 integer, intent(in) :: u
 character(len=*), intent(in) :: io
 integer, dimension(:), intent(in) :: list
 integer, intent(out) :: eta
 character(len=*), intent(inout) :: msg
 integer :: i,j

 select case(io)
  case ("DTVERBOSE")
   do i=1,3
    write(u,'(a,3e12.4,a,a)',advance="no") '|', T%val(i,1:3),'|',C_NEW_LINE
   end do
  case default
   write(u,'(9e12.4)') ((T%val(i,j),j=1,3),i=1,3)
  end select
  eta=0

end subroutine writea_Tensor2
! **********************************************************************
! *** subroutine writeb_Tensor2                                      ***
! *** writes a second order tensor in a binary file                  ***
! ***                                                                ***
! *** input : T, second order tensor to be written                   ***
! ***         u, unit to write in                                    ***
! ***         eta, !!! utilité ???                                   ***
! ***                                                                ***
! *** in/output : msg, message to be passed                          ***
! ***                                                                ***
! *** output :                                                       ***
! **********************************************************************
subroutine writeb_Tensor2(T,u,eta,msg)
 implicit none
 class(Tensor2), intent(in) :: T
 integer, intent(in) :: u
 integer, intent(out) :: eta
 character(len=*), intent(inout) :: msg
 integer :: i !!! CM to be removed ??
 
 write(u) T%val
 eta=0
 
end subroutine writeb_Tensor2
! **********************************************************************
! *** subroutine readb_Tensor2                                       ***
! *** reads a second order tensor from a binary file                 ***
! ***                                                                ***
! *** input : u, unit to read in                                     ***
! ***         eta, !!! utilité ???                                   ***
! ***                                                                ***
! *** in/output : msg, message to be passed                          ***
! ***                                                                ***
! *** output : T, second order tensor to be read                     ***
! **********************************************************************
subroutine readb_Tensor2(T,u,eta,msg)
 implicit none
 class(Tensor2), intent(inout) :: T
 integer, intent(in) :: u
 integer, intent(out) :: eta
 character(len=*), intent(inout) :: msg
 integer :: i !!! CM to be removed ??
 
 read(unit=u,iostat=eta,iomsg=msg) T%val
 eta=0
 
end subroutine readb_Tensor2
! **********************************************************************
! *** subroutine writea_Tensor2s                                     ***
! *** writes a symmetric second order tensor in a text file          ***
! ***                                                                ***
! *** input : T, second order tensor to be written                   ***
! ***         u, unit to write in                                    ***
! ***         io, format of the output                               ***
! ***         list, !!! CM non utilisé ???                           ***
! ***         eta, !!! utilité ???                                   ***
! ***                                                                ***
! *** in/output : msg, message to be passed                          ***
! ***                                                                ***
! *** output :                                                       ***
! **********************************************************************
subroutine writea_Tensor2s(T,u,io,list,eta,msg)
 implicit none
 class(Tensor2s), intent(in) :: T
 integer, intent(in) :: u
 character(len=*), intent(in) :: io
 integer, dimension(:), intent(in) :: list
 integer, intent(out) :: eta
 character(len=*), intent(inout) :: msg
 type(Tensor2) :: T2

 select case(io)
  case ("DTVERBOSE")
   T2=full(T)
   call T2%writea(u,io,list,eta,msg)
  case ("DTFULL")
   T2=full(T)
   call T2%writea(u,io,list,eta,msg)
  case default
   write(u,'(6e12.4)') T%val(1:6)
  end select
  eta=0

end subroutine writea_Tensor2s
! **********************************************************************
! *** subroutine writeb_Tensor2s                                     ***
! *** writes a symmetric second order tensor in a binary file        ***
! ***                                                                ***
! *** input : T, second order tensor to be written                   ***
! ***         u, unit to write in                                    ***
! ***         eta, !!! utilité ???                                   ***
! ***                                                                ***
! *** in/output : msg, message to be passed                          ***
! ***                                                                ***
! *** output :                                                       ***
! **********************************************************************
subroutine writeb_Tensor2s(T,u,eta,msg)
 implicit none
 class(Tensor2s), intent(in) :: T
 integer, intent(in) :: u
 integer, intent(out) :: eta
 character(len=*), intent(inout) :: msg
 integer :: i
 
 write(u) T%val
 eta=0

end subroutine writeb_Tensor2s
! **********************************************************************
! *** subroutine readb_Tensor2s                                      ***
! *** reads a symmetric second order tensor from a binary file       ***
! ***                                                                ***
! *** input : u, unit to read in                                     ***
! ***         eta, !!! utilité ???                                   ***
! ***                                                                ***
! *** in/output : msg, message to be passed                          ***
! ***                                                                ***
! *** output : T, second order tensor to be read                     ***
! **********************************************************************
subroutine readb_Tensor2s(T,u,eta,msg)
 implicit none
 class(Tensor2s), intent(inout) :: T
 integer, intent(in) :: u
 integer, intent(out) :: eta
 character(len=*), intent(inout) :: msg
 integer :: i
 
 read(unit=u,iostat=eta,iomsg=msg) T%val
 eta=0
 
end subroutine readb_Tensor2s
! **********************************************************************
! *** subroutine writea_Tensor2w                                     ***
! *** writes a skew-symmetric second order tensor in a text file     ***
! ***                                                                ***
! *** input : T, second order tensor to be written                   ***
! ***         u, unit to write in                                    ***
! ***         io, format of the output                               ***
! ***         list, !!! CM non utilisé ???                           ***
! ***         eta, !!! utilité ???                                   ***
! ***                                                                ***
! *** in/output : msg, message to be passed                          ***
! ***                                                                ***
! *** output :                                                       ***
! **********************************************************************
subroutine writea_Tensor2w(T,u,io,list,eta,msg)
 implicit none
 class(Tensor2w), intent(in) :: T
 integer, intent(in) :: u
 character(len=*), intent(in) :: io
 integer, dimension(:), intent(in) :: list
 integer, intent(out) :: eta
 character(len=*), intent(inout) :: msg
 type(Tensor2) :: T2

 select case(io)
  case ("DTVERBOSE")
   T2=full(T)
   call T2%writea(u,io,list,eta,msg)
  case ("DTFULL")
   T2=full(T)
   call T2%writea(u,io,list,eta,msg)
  case default
   write(u,'(3e12.4)') T%val(1:3)
 end select
 eta=0
 
end subroutine writea_Tensor2w
! **********************************************************************
! *** subroutine writeb_Tensor2w                                     ***
! *** writes a skew-symmetric second order tensor in a binary file   ***
! ***                                                                ***
! *** input : T, second order tensor to be written                   ***
! ***         u, unit to write in                                    ***
! ***         eta, !!! utilité ???                                   ***
! ***                                                                ***
! *** in/output : msg, message to be passed                          ***
! ***                                                                ***
! *** output :                                                       ***
! **********************************************************************
subroutine writeb_Tensor2w(T,u,eta,msg)
 implicit none
 class(Tensor2w), intent(in) :: T
 integer, intent(in) :: u
 integer, intent(out) :: eta
 character(len=*), intent(inout) :: msg
 integer :: i

 write(u) T%val
 eta=0

end subroutine writeb_Tensor2w
! **********************************************************************
! *** subroutine readb_Tensor2w                                      ***
! *** reads a skew-symmetric second order tensor from a binary file  ***
! ***                                                                ***
! *** input : u, unit to read in                                     ***
! ***         eta, !!! utilité ???                                   ***
! ***                                                                ***
! *** in/output : msg, message to be passed                          ***
! ***                                                                ***
! *** output : T, second order tensor to be read                     ***
! **********************************************************************
subroutine readb_Tensor2w(T,u,eta,msg)
 implicit none
 class(Tensor2w), intent(inout) :: T
 integer, intent(in) :: u
 integer, intent(out) :: eta
 character(len=*), intent(inout) :: msg
 integer :: i

 read(unit=u,iostat=eta,iomsg=msg) T%val
 eta=0

end subroutine readb_Tensor2w
! **********************************************************************
! *** subroutine writea_Tensor3                                      ***
! *** writes a third order tensor in a text file                     ***
! ***                                                                ***
! *** input : T, third order tensor to be written                    ***
! ***         u, unit to write in                                    ***
! ***         io, format of the output                               ***
! ***         list, !!! CM non utilisé ???                           ***
! ***         eta, !!! utilité ???                                   ***
! ***                                                                ***
! *** in/output : msg, message to be passed                          ***
! ***                                                                ***
! *** output :                                                       ***
! **********************************************************************
subroutine writea_Tensor3(T,u,io,list,eta,msg)
 implicit none
 class(Tensor3), intent(in) :: T
 integer, intent(in) :: u
 character(len=*), intent(in) :: io
 integer, dimension(:), intent(in) :: list
 integer, intent(out) :: eta
 character(len=*), intent(inout) :: msg
 integer :: i

!   select case(io)
!    case ("DTVERBOSE")
!     do i=1,3
!      write(u,'(a,3e12.4,a,a)',advance="no") '|', T%val(i,1:3),'|',C_NEW_LINE
!     end do
!    case default
 write(u,'(27e12.4)') T%val
!   end select
 eta=0

end subroutine writea_Tensor3
! **********************************************************************
! *** subroutine writeb_Tensor3                                      ***
! *** writes a third order tensor in a binary file                   ***
! ***                                                                ***
! *** input : T, third order tensor to be written                    ***
! ***         u, unit to write in                                    ***
! ***         eta, !!! utilité ???                                   ***
! ***                                                                ***
! *** in/output : msg, message to be passed                          ***
! ***                                                                ***
! *** output :                                                       ***
! **********************************************************************
subroutine writeb_Tensor3(T,u,eta,msg)
 implicit none
 class(Tensor3), intent(in) :: T
 integer, intent(in) :: u
 integer, intent(out) :: eta
 character(len=*), intent(inout) :: msg
 integer :: i

 write(u) T%val
 eta=0
 
end subroutine writeb_Tensor3
! **********************************************************************
! *** subroutine readb_Tensor3                                       ***
! *** reads a third order tensor from a binary file                  ***
! ***                                                                ***
! *** input : u, unit to read in                                     ***
! ***         eta, !!! utilité ???                                   ***
! ***                                                                ***
! *** in/output : msg, message to be passed                          ***
! ***                                                                ***
! *** output : T, third order tensor to be read                      ***
! **********************************************************************
subroutine readb_Tensor3(T,u,eta,msg)
 implicit none
 class(Tensor3), intent(inout) :: T
 integer, intent(in) :: u
 integer, intent(out) :: eta
 character(len=*), intent(inout) :: msg
 integer :: i
 
 read(unit=u,iostat=eta,iomsg=msg) T%val
 eta=0

end subroutine readb_Tensor3
! **********************************************************************
! *** subroutine writea_Tensor4                                      ***
! *** writes a fourth order tensor in a text file                    ***
! ***                                                                ***
! *** input : T, fourth order tensor to be written                   ***
! ***         u, unit to write in                                    ***
! ***         io, format of the output                               ***
! ***         list, !!! CM non utilisé ???                           ***
! ***         eta, !!! utilité ???                                   ***
! ***                                                                ***
! *** in/output : msg, message to be passed                          ***
! ***                                                                ***
! *** output :                                                       ***
! **********************************************************************
subroutine writea_Tensor4(T,u,io,list,eta,msg)
 implicit none
 class(Tensor4), intent(in) :: T
 integer, intent(in) :: u
 character(len=*), intent(in) :: io
 integer, dimension(:), intent(in) :: list
 integer, intent(out) :: eta
 character(len=*), intent(inout) :: msg
 integer :: i

!   select case(io)
!    case ("DTVERBOSE")
!     do i=1,3
!      write(u,'(a,3e12.4,a,a)',advance="no") '|', T%val(i,1:3),'|',C_NEW_LINE
!     end do
!    case default
 write(u,'(81e12.4)') T%val
!   end select
 eta=0

end subroutine writea_Tensor4
! **********************************************************************
! *** subroutine writeb_Tensor3                                      ***
! *** writes a third order tensor in a binary file                   ***
! ***                                                                ***
! *** input : T, fourth order tensor to be written                   ***
! ***         u, unit to write in                                    ***
! ***         eta, !!! utilité ???                                   ***
! ***                                                                ***
! *** in/output : msg, message to be passed                          ***
! ***                                                                ***
! *** output :                                                       ***
! **********************************************************************
subroutine writeb_Tensor4(T,u,eta,msg)
 implicit none
 class(Tensor4), intent(in) :: T
 integer, intent(in) :: u
 integer, intent(out) :: eta
 character(len=*), intent(inout) :: msg
 integer :: i

 write(u) T%val
 eta=0
 
end subroutine writeb_Tensor4
! **********************************************************************
! *** subroutine readb_Tensor4                                       ***
! *** reads a fourth order tensor from a binary file                 ***
! ***                                                                ***
! *** input : u, unit to read in                                     ***
! ***         eta, !!! utilité ???                                   ***
! ***                                                                ***
! *** in/output : msg, message to be passed                          ***
! ***                                                                ***
! *** output : T, fourth order tensor to be read                     ***
! **********************************************************************
subroutine readb_Tensor4(T,u,eta,msg)
 implicit none
 class(Tensor4), intent(inout) :: T
 integer, intent(in) :: u
 integer, intent(out) :: eta
 character(len=*), intent(inout) :: msg
 integer :: i
 
 read(unit=u,iostat=eta,iomsg=msg) T%val
 eta=0

end subroutine readb_Tensor4
! **********************************************************************
! *** subroutine writea_Tensor4ss                                    ***
! *** writes a symmetric fourth order tensor in a text file          ***
! ***                                                                ***
! *** input : T, fourth order tensor to be written                   ***
! ***         u, unit to write in                                    ***
! ***         io, format of the output                               ***
! ***         list, !!! CM non utilisé ???                           ***
! ***         eta, !!! utilité ???                                   ***
! ***                                                                ***
! *** in/output : msg, message to be passed                          ***
! ***                                                                ***
! *** output :                                                       ***
! **********************************************************************
subroutine writea_Tensor4ss(T,u,io,list,eta,msg)
 implicit none
 class(Tensor4ss), intent(in) :: T
 integer, intent(in) :: u
 character(len=*), intent(in) :: io
 integer, dimension(:), intent(in) :: list
 integer, intent(out) :: eta
 character(len=*), intent(inout) :: msg
 integer :: i

 select case(io)
  case ("DTVERBOSE")
   do i=1,6
    write(u,'(a,6e12.4,a,a)',advance="no") '|', T%val(i,1:6),'|',C_NEW_LINE
   end do
  case default
   write(u,'(36e12.4)') T%val
 end select
 eta=0

end subroutine writea_Tensor4ss
! **********************************************************************
! *** subroutine writeb_Tensor4ss                                    ***
! *** writes a symmetric fourth order tensor in a binary file        ***
! ***                                                                ***
! *** input : T, fourth order tensor to be written                   ***
! ***         u, unit to write in                                    ***
! ***         eta, !!! utilité ???                                   ***
! ***                                                                ***
! *** in/output : msg, message to be passed                          ***
! ***                                                                ***
! *** output :                                                       ***
! **********************************************************************
subroutine writeb_Tensor4ss(T,u,eta,msg)
 implicit none
 class(Tensor4ss), intent(in) :: T
 integer, intent(in) :: u
 integer, intent(out) :: eta
 character(len=*), intent(inout) :: msg
 integer :: i

 write(u) T%val
 eta=0
 
end subroutine writeb_Tensor4ss
! **********************************************************************
! *** subroutine readb_Tensor4ss                                     ***
! *** reads a symmetric fourth order tensor from a binary file       ***
! ***                                                                ***
! *** input : u, unit to read in                                     ***
! ***         eta, !!! utilité ???                                   ***
! ***                                                                ***
! *** in/output : msg, message to be passed                          ***
! ***                                                                ***
! *** output : T, fourth order tensor to be read                     ***
! **********************************************************************
subroutine readb_Tensor4ss(T,u,eta,msg)
 implicit none
 class(Tensor4ss), intent(inout) :: T
 integer, intent(in) :: u
 integer, intent(out) :: eta
 character(len=*), intent(inout) :: msg
 integer :: i
 
 read(unit=u,iostat=eta,iomsg=msg) T%val
 eta=0

end subroutine readb_Tensor4ss
! **********************************************************************
! *** function full_Tensor2w                                         ***
! *** converts a skew-symmetric second order tensor to a general     ***
! *** second order tensor                                            ***
! ***                                                                ***
! *** input : T, skew-symmetric second order tensor to be converted  ***
! **********************************************************************
function full_Tensor2w(T)
 implicit none
 type(Tensor2w), intent(in) :: T
 type(Tensor2) :: full_Tensor2w

 full_Tensor2w%val(1,1)=0.d0
 full_Tensor2w%val(2,1)=-T%val(3)
 full_Tensor2w%val(3,1)=T%val(2)
 full_Tensor2w%val(1,2)=T%val(3)
 full_Tensor2w%val(2,2)=0.d0
 full_Tensor2w%val(3,2)=-T%val(1)
 full_Tensor2w%val(1,3)=-T%val(2)
 full_Tensor2w%val(2,3)=T%val(1)
 full_Tensor2w%val(3,3)=0.d0

end function full_Tensor2w
! **********************************************************************
! *** function full_Tensor2s                                         ***
! *** converts a symmetric second order tensor to a general second   ***
! *** order tensor                                                   ***
! ***                                                                ***
! *** input : T, symmetric second order tensor to be converted       ***
! **********************************************************************
function full_Tensor2s(T)
 implicit none
 type(Tensor2s), intent(in) :: T
 type(Tensor2) :: full_Tensor2s

 full_Tensor2s%val(1,1)=T%val(1)
 full_Tensor2s%val(2,1)=T%val(6)
 full_Tensor2s%val(3,1)=T%val(5)
 full_Tensor2s%val(1,2)=T%val(6)
 full_Tensor2s%val(2,2)=T%val(2)
 full_Tensor2s%val(3,2)=T%val(4)
 full_Tensor2s%val(1,3)=T%val(5)
 full_Tensor2s%val(2,3)=T%val(4)
 full_Tensor2s%val(3,3)=T%val(3)
 
end function full_Tensor2s
! **********************************************************************
! *** function full_Tensor4ss                                        ***
! *** converts a symmetric fourth order tensor to a general fourth   ***
! *** order tensor                                                   ***
! ***                                                                ***
! *** input : T, symmetric fourth order tensor to be converted       ***
! **********************************************************************
function full_Tensor4ss(T)
 implicit none
 type(Tensor4ss), intent(in) :: T
 type(Tensor4) :: full_Tensor4ss
 integer :: i,j,k,l,ij,kl
 
 do l=1,3
  do k=1,3
   kl=indxs(k,l)
   do j=1,3
    do i=1,3
     ij=indxs(i,j)
     full_Tensor4ss%val(i,j,k,l)=T%val(ij,kl)
    enddo
   enddo
  enddo
 enddo
    
 
end function full_Tensor4ss
! **********************************************************************
! *** function nrm_Tensor1                                           ***
! *** returns the norm of a first order tensor                       ***
! ***                                                                ***
! *** input : T, first order tensor for which the norm needs to be   ***
! ***            calculated                                          ***
! **********************************************************************
function nrm_Tensor1(T)
 implicit none
 type(Tensor1), intent(in) :: T
 real(8) :: nrm_Tensor1

 nrm_Tensor1=sqrt(T%val(1)*T%val(1)+T%val(2)*T%val(2)+T%val(3)*T%val(3)) 
 
end function nrm_Tensor1
! **********************************************************************
! *** function nrm_Tensor2                                           ***
! *** returns the norm of a second order tensor                      ***
! ***                                                                ***
! *** input : T, second order tensor for which the norm needs to be  ***
! ***            calculated                                          ***
! **********************************************************************
function nrm_Tensor2(T)
 implicit none
 type(Tensor2), intent(in) :: T
 real(8) :: nrm_Tensor2

 nrm_Tensor2=sqrt(T%val(1,1)*T%val(1,1)+T%val(2,1)*T%val(2,1)+T%val(3,1)*T%val(3,1) &
                 +T%val(1,2)*T%val(1,2)+T%val(2,2)*T%val(2,2)+T%val(3,2)*T%val(3,2) &
                 +T%val(1,3)*T%val(1,3)+T%val(2,3)*T%val(2,3)+T%val(3,3)*T%val(3,3))
 
end function nrm_Tensor2
! **********************************************************************
! *** function nrm_Tensor2s                                          ***
! *** returns the norm of a symmetric second order tensor            ***
! ***                                                                ***
! *** input : T, second order tensor for which the norm needs to be  ***
! ***            calculated                                          ***
! **********************************************************************
function nrm_Tensor2s(T)
 implicit none
 type(Tensor2s), intent(in) :: T
 real(8) :: nrm_Tensor2s

 nrm_Tensor2s=sqrt(T%val(1)*T%val(1)+T%val(2)*T%val(2)+T%val(3)*T%val(3) &
             +2.d0*T%val(4)*T%val(4)+2.d0*T%val(5)*T%val(5)+2.d0*T%val(6)*T%val(6))
 
end function nrm_Tensor2s
! **********************************************************************
! *** function nrm_Tensor2w                                          ***
! *** returns the norm of a skew-symmetric second order tensor       ***
! ***                                                                ***
! *** input : T, second order tensor for which the norm needs to be  ***
! ***            calculated                                          ***
! **********************************************************************
function nrm_Tensor2w(T)
 implicit none
 type(Tensor2w), intent(in) :: T
 real(8) :: nrm_Tensor2w

 nrm_Tensor2w=sqrt(2.d0*T%val(1)*T%val(1)+2.d0*T%val(2)*T%val(2)+2.d0*T%val(3)*T%val(3))
 
end function nrm_Tensor2w
! **********************************************************************
! *** function nrm_Tensor3                                           ***
! *** returns the norm of a third order tensor                       ***
! ***                                                                ***
! *** input : T, third order tensor for which the norm needs to be   ***
! ***            calculated                                          ***
! **********************************************************************
function nrm_Tensor3(T)
 implicit none
 type(Tensor3), intent(in) :: T
 real(8) :: nrm_Tensor3

 nrm_Tensor3=sqrt(T%val(1,1,1)*T%val(1,1,1)+T%val(2,1,1)*T%val(2,1,1)+T%val(3,1,1)*T%val(3,1,1) &
                 +T%val(1,2,1)*T%val(1,2,1)+T%val(2,2,1)*T%val(2,2,1)+T%val(3,2,1)*T%val(3,2,1) &
                 +T%val(1,3,1)*T%val(1,3,1)+T%val(2,3,1)*T%val(2,3,1)+T%val(3,3,1)*T%val(3,3,1) &
                 +T%val(1,1,2)*T%val(1,1,2)+T%val(2,1,2)*T%val(2,1,2)+T%val(3,1,2)*T%val(3,1,2) &
                 +T%val(1,2,2)*T%val(1,2,2)+T%val(2,2,2)*T%val(2,2,2)+T%val(3,2,2)*T%val(3,2,2) &
                 +T%val(1,3,2)*T%val(1,3,2)+T%val(2,3,2)*T%val(2,3,2)+T%val(3,3,2)*T%val(3,3,2) &
                 +T%val(1,1,3)*T%val(1,1,3)+T%val(2,1,3)*T%val(2,1,3)+T%val(3,1,3)*T%val(3,1,3) &
                 +T%val(1,2,3)*T%val(1,2,3)+T%val(2,2,3)*T%val(2,2,3)+T%val(3,2,3)*T%val(3,2,3) &
                 +T%val(1,3,3)*T%val(1,3,3)+T%val(2,3,3)*T%val(2,3,3)+T%val(3,3,3)*T%val(3,3,3))
 
end function nrm_Tensor3
! **********************************************************************
! *** function nrm_Tensor4                                           ***
! *** returns the norm of a fourth order tensor                      ***
! ***                                                                ***
! *** input : T, fourth order tensor for which the norm needs to be  ***
! ***            calculated                                          ***
! **********************************************************************
function nrm_Tensor4(T)
 implicit none
 type(Tensor4), intent(in) :: T
 real(8) :: nrm_Tensor4,norm2
 integer :: i,j,k,l
 
 norm2=0.d0
 do l=1,3
  do k=1,3
   do j=1,3
    do i=1,3
     norm2=norm2+T%val(i,j,k,l)*T%val(i,j,k,l)
    enddo
   enddo
  enddo
 enddo
 nrm_Tensor4=sqrt(norm2)
 
end function nrm_Tensor4
! **********************************************************************
! *** function nrm_Tensor4ss                                         ***
! *** returns the norm of a symmetric fourth order tensor            ***
! ***                                                                ***
! *** input : T, third order tensor for which the norm needs to be   ***
! ***            calculated                                          ***
! **********************************************************************
function nrm_Tensor4ss(T)
 implicit none
 type(Tensor4ss), intent(in) :: T
 real(8) :: nrm_Tensor4ss

 nrm_Tensor4ss=sqrt(1.d0*T%val(1,1)*T%val(1,1)+1.d0*T%val(2,1)*T%val(2,1)+1.d0*T%val(3,1)*T%val(3,1) &
                   +2.d0*T%val(4,1)*T%val(4,1)+2.d0*T%val(5,1)*T%val(5,1)+2.d0*T%val(6,1)*T%val(6,1) &
                   +1.d0*T%val(1,2)*T%val(1,2)+1.d0*T%val(2,2)*T%val(2,2)+1.d0*T%val(3,2)*T%val(3,2) &
                   +2.d0*T%val(4,2)*T%val(4,2)+2.d0*T%val(5,2)*T%val(5,2)+2.d0*T%val(6,2)*T%val(6,2) &
                   +1.d0*T%val(1,3)*T%val(1,3)+1.d0*T%val(2,3)*T%val(2,3)+1.d0*T%val(3,3)*T%val(3,3) &
                   +2.d0*T%val(4,3)*T%val(4,3)+2.d0*T%val(5,3)*T%val(5,3)+2.d0*T%val(6,3)*T%val(6,3) &
                   +2.d0*T%val(1,4)*T%val(1,4)+2.d0*T%val(2,4)*T%val(2,4)+2.d0*T%val(3,4)*T%val(3,4) &
                   +4.d0*T%val(4,4)*T%val(4,4)+4.d0*T%val(5,4)*T%val(5,4)+4.d0*T%val(6,4)*T%val(6,4) &
                   +2.d0*T%val(1,5)*T%val(1,5)+2.d0*T%val(2,5)*T%val(2,5)+2.d0*T%val(3,5)*T%val(3,5) &
                   +4.d0*T%val(4,5)*T%val(4,5)+4.d0*T%val(5,5)*T%val(5,5)+4.d0*T%val(6,5)*T%val(6,5) &
                   +2.d0*T%val(1,6)*T%val(1,6)+2.d0*T%val(2,6)*T%val(2,6)+2.d0*T%val(3,6)*T%val(3,6) &
                   +4.d0*T%val(4,6)*T%val(4,6)+4.d0*T%val(5,6)*T%val(5,6)+4.d0*T%val(6,6)*T%val(6,6))
 
end function nrm_Tensor4ss
! **********************************************************************
! *** function symsym_Tensor4                                        ***
! *** returns the minor symmetric part of a fourth order tensor      ***
! ***                                                                ***
! *** input : T4, fourth order tensor                                ***
! **********************************************************************
function symsym_Tensor4(T4)
 implicit none
 type(Tensor4), intent(in) :: T4
 type(Tensor4ss) :: symsym_Tensor4

 
 symsym_Tensor4%val(1,1)=T4%val(1,1,1,1)
 symsym_Tensor4%val(2,1)=T4%val(2,2,1,1)
 symsym_Tensor4%val(3,1)=T4%val(3,3,1,1)
 symsym_Tensor4%val(4,1)=(T4%val(2,3,1,1)+T4%val(3,2,1,1))*0.5d0
 symsym_Tensor4%val(5,1)=(T4%val(3,1,1,1)+T4%val(1,3,1,1))*0.5d0
 symsym_Tensor4%val(6,1)=(T4%val(1,2,1,1)+T4%val(2,1,1,1))*0.5d0
 
 symsym_Tensor4%val(1,2)=T4%val(1,1,2,2)
 symsym_Tensor4%val(2,2)=T4%val(2,2,2,2)
 symsym_Tensor4%val(3,2)=T4%val(3,3,2,2)
 symsym_Tensor4%val(4,2)=(T4%val(2,3,2,2)+T4%val(3,2,2,2))*0.5d0
 symsym_Tensor4%val(5,2)=(T4%val(3,1,2,2)+T4%val(1,3,2,2))*0.5d0
 symsym_Tensor4%val(6,2)=(T4%val(1,2,2,2)+T4%val(2,1,2,2))*0.5d0
 
 symsym_Tensor4%val(1,3)=T4%val(1,1,3,3)
 symsym_Tensor4%val(2,3)=T4%val(2,2,3,3)
 symsym_Tensor4%val(3,3)=T4%val(3,3,3,3)
 symsym_Tensor4%val(4,3)=(T4%val(2,3,3,3)+T4%val(3,2,3,3))*0.5d0
 symsym_Tensor4%val(5,3)=(T4%val(3,1,3,3)+T4%val(1,3,3,3))*0.5d0
 symsym_Tensor4%val(6,3)=(T4%val(1,2,3,3)+T4%val(2,1,3,3))*0.5d0
 
 symsym_Tensor4%val(1,4)=(T4%val(1,1,2,3)+T4%val(1,1,3,2))*0.5d0
 symsym_Tensor4%val(2,4)=(T4%val(2,2,2,3)+T4%val(2,2,3,2))*0.5d0
 symsym_Tensor4%val(3,4)=(T4%val(3,3,2,3)+T4%val(3,3,3,2))*0.5d0
 symsym_Tensor4%val(4,4)=(T4%val(2,3,2,3)+T4%val(3,2,2,3)+T4%val(3,2,3,2)+T4%val(2,3,3,2))*0.25d0
 symsym_Tensor4%val(5,4)=(T4%val(3,1,2,3)+T4%val(1,3,2,3)+T4%val(1,3,3,2)+T4%val(3,1,3,2))*0.25d0
 symsym_Tensor4%val(6,4)=(T4%val(1,2,2,3)+T4%val(2,1,2,3)+T4%val(2,1,3,2)+T4%val(1,2,3,2))*0.25d0
 
 symsym_Tensor4%val(1,5)=(T4%val(1,1,3,1)+T4%val(1,1,1,3))*0.5d0
 symsym_Tensor4%val(2,5)=(T4%val(2,2,3,1)+T4%val(2,2,1,3))*0.5d0
 symsym_Tensor4%val(3,5)=(T4%val(3,3,3,1)+T4%val(3,3,1,3))*0.5d0
 symsym_Tensor4%val(4,5)=(T4%val(2,3,3,1)+T4%val(3,2,3,1)+T4%val(3,2,1,3)+T4%val(2,3,1,3))*0.25d0
 symsym_Tensor4%val(5,5)=(T4%val(3,1,3,1)+T4%val(1,3,3,1)+T4%val(1,3,1,3)+T4%val(3,1,1,3))*0.25d0
 symsym_Tensor4%val(6,5)=(T4%val(1,2,3,1)+T4%val(2,1,3,1)+T4%val(2,1,1,3)+T4%val(1,2,1,3))*0.25d0
 
 symsym_Tensor4%val(1,6)=(T4%val(1,1,1,2)+T4%val(1,1,2,1))*0.5d0
 symsym_Tensor4%val(2,6)=(T4%val(2,2,1,2)+T4%val(2,2,2,1))*0.5d0
 symsym_Tensor4%val(3,6)=(T4%val(3,3,1,2)+T4%val(3,3,2,1))*0.5d0
 symsym_Tensor4%val(4,6)=(T4%val(2,3,1,2)+T4%val(3,2,1,2)+T4%val(3,2,2,1)+T4%val(2,3,2,1))*0.25d0
 symsym_Tensor4%val(5,6)=(T4%val(3,1,1,2)+T4%val(1,3,1,2)+T4%val(1,3,2,1)+T4%val(3,1,2,1))*0.25d0
 symsym_Tensor4%val(6,6)=(T4%val(1,2,1,2)+T4%val(2,1,1,2)+T4%val(2,1,2,1)+T4%val(1,2,2,1))*0.25d0

end function symsym_Tensor4
! **********************************************************************
! *** function sym_Tensor2                                           ***
! *** returns the symmetric part of a second order tensor            ***
! ***                                                                ***
! *** input : T2, second order tensor                                ***
! **********************************************************************
function sym_Tensor2(T2)
 implicit none
 type(Tensor2), intent(in) :: T2
 type(Tensor2s) :: sym_Tensor2
 
 sym_Tensor2%val(1)=T2%val(1,1)
 sym_Tensor2%val(2)=T2%val(2,2)
 sym_Tensor2%val(3)=T2%val(3,3)
 sym_Tensor2%val(4)=0.5d0*(T2%val(2,3)+T2%val(3,2))
 sym_Tensor2%val(5)=0.5d0*(T2%val(3,1)+T2%val(1,3))
 sym_Tensor2%val(6)=0.5d0*(T2%val(1,2)+T2%val(2,1))

end function sym_Tensor2
! **********************************************************************
! *** function skw_Tensor2                                           ***
! *** returns the skew-symmetric part of a second order tensor       ***
! ***                                                                ***
! *** input : T2, second order tensor                                ***
! **********************************************************************
function skw_Tensor2(T2)
 implicit none
 type(Tensor2), intent(in) :: T2
 type(Tensor2w) :: skw_Tensor2
 
 skw_Tensor2%val(1)=0.5d0*(T2%val(2,3)-T2%val(3,2))
 skw_Tensor2%val(2)=0.5d0*(T2%val(3,1)-T2%val(1,3))
 skw_Tensor2%val(3)=0.5d0*(T2%val(1,2)-T2%val(2,1))

end function skw_Tensor2
! **********************************************************************
! *** function trc_Tensor2                                           ***
! *** returns the trace of a general second order tensor             ***
! ***                                                                ***
! *** input : T2, second order tensor                                ***
! **********************************************************************
function trc_Tensor2(T2)
 implicit none
 type(Tensor2), intent(in) :: T2
 real(8) :: trc_Tensor2
 
 trc_Tensor2=T2%val(1,1)+T2%val(2,2)+T2%val(3,3)
 
end function trc_Tensor2
! **********************************************************************
! *** function trc_Tensor2s                                          ***
! *** returns the trace of a symmetric second order tensor           ***
! ***                                                                ***
! *** input : T2, second order tensor                                ***
! **********************************************************************
function trc_Tensor2s(T2s)
 implicit none
 type(Tensor2s), intent(in) :: T2s
 real(8) :: trc_Tensor2s
 
 trc_Tensor2s=T2s%val(1)+T2s%val(2)+T2s%val(3)
 
end function trc_Tensor2s
! **********************************************************************
! *** function dev_Tensor2                                           ***
! *** returns the deviatoric part of a general second order tensor   ***
! ***                                                                ***
! *** input : T2, second order tensor                                ***
! **********************************************************************
function dev_Tensor2(T2)
 implicit none
 type(Tensor2), intent(in) :: T2
 type(Tensor2) :: dev_Tensor2
 real(8) :: hyd
 
 hyd=trc_Tensor2(T2)/3.d0
 
 dev_Tensor2%val(1,1)=T2%val(1,1)-hyd
 dev_Tensor2%val(2,1)=T2%val(2,1)
 dev_Tensor2%val(3,1)=T2%val(3,1)
 dev_Tensor2%val(1,2)=T2%val(1,2)
 dev_Tensor2%val(2,2)=T2%val(2,2)-hyd
 dev_Tensor2%val(3,2)=T2%val(3,2)
 dev_Tensor2%val(1,3)=T2%val(1,3)
 dev_Tensor2%val(2,3)=T2%val(2,3)
 dev_Tensor2%val(3,3)=T2%val(3,3)-hyd

end function dev_Tensor2
! **********************************************************************
! *** function dev_Tensor2s                                          ***
! *** returns the deviatoric part of a symmetric second order tensor ***
! ***                                                                ***
! *** input : T2s, second order tensor                               ***
! **********************************************************************
function dev_Tensor2s(T2s)
 implicit none
 type(Tensor2s), intent(in) :: T2s
 type(Tensor2s) :: dev_Tensor2s
 real(8) :: hyd
 
 hyd=trc_Tensor2s(T2s)/3.d0
 
 dev_Tensor2s%val(1)=T2s%val(1)-hyd
 dev_Tensor2s%val(2)=T2s%val(2)-hyd
 dev_Tensor2s%val(3)=T2s%val(3)-hyd
 dev_Tensor2s%val(4)=T2s%val(4)
 dev_Tensor2s%val(5)=T2s%val(5)
 dev_Tensor2s%val(6)=T2s%val(6)

end function dev_Tensor2s
! **********************************************************************
! *** function sph_Tensor2                                           ***
! *** returns the spherical part of a general second order tensor    ***
! ***                                                                ***
! *** input : T2, second order tensor                                ***
! **********************************************************************
function sph_Tensor2(T2)
 implicit none
 type(Tensor2), intent(in) :: T2
 type(Tensor2) :: sph_Tensor2
 real(8) :: hyd
 
 hyd=trc_Tensor2(T2)/3.d0
 
 sph_Tensor2%val(1,1)=hyd
 sph_Tensor2%val(2,1)=0.d0
 sph_Tensor2%val(3,1)=0.d0
 sph_Tensor2%val(1,2)=0.d0
 sph_Tensor2%val(2,2)=hyd
 sph_Tensor2%val(3,2)=0.d0
 sph_Tensor2%val(1,3)=0.d0
 sph_Tensor2%val(2,3)=0.d0
 sph_Tensor2%val(3,3)=hyd

end function sph_Tensor2
! **********************************************************************
! *** function sph_Tensor2s                                          ***
! *** returns the spherical part of a symmetric second order tensor  ***
! ***                                                                ***
! *** input : T2s, second order tensor                               ***
! **********************************************************************
function sph_Tensor2s(T2s)
 implicit none
 type(Tensor2s), intent(in) :: T2s
 type(Tensor2s) :: sph_Tensor2s
 real(8) :: hyd
 
 hyd=trc_Tensor2s(T2s)/3.d0
 
 sph_Tensor2s%val(1)=hyd
 sph_Tensor2s%val(2)=hyd
 sph_Tensor2s%val(3)=hyd
 sph_Tensor2s%val(4)=0.d0
 sph_Tensor2s%val(5)=0.d0
 sph_Tensor2s%val(6)=0.d0

end function sph_Tensor2s
! **********************************************************************
! *** function det_Tensor2                                           ***
! *** returns the determinant of a general second order tensor       ***
! ***                                                                ***
! *** input : T2, second order tensor                                ***
! **********************************************************************
function det_Tensor2(T2)
 implicit none
 type(Tensor2), intent(in) :: T2
 real(8) :: det_Tensor2

 det_Tensor2=T2%val(1,1)*(T2%val(3,3)*T2%val(2,2)-T2%val(3,2)*T2%val(2,3)) &
            -T2%val(2,1)*(T2%val(3,3)*T2%val(1,2)-T2%val(3,2)*T2%val(1,3)) &
            +T2%val(3,1)*(T2%val(2,3)*T2%val(1,2)-T2%val(2,2)*T2%val(1,3))
  
end function det_Tensor2
! **********************************************************************
! *** function det_Tensor2s                                          ***
! *** returns the determinant of a symmetric second order tensor     ***
! ***                                                                ***
! *** input : T2s, second order tensor                               ***
! **********************************************************************
function det_Tensor2s(T2s)
 implicit none
 type(Tensor2s), intent(in) :: T2s
 real(8) :: det_Tensor2s

 det_Tensor2s=T2s%val(1)*(T2s%val(3)*T2s%val(2)-T2s%val(4)*T2s%val(4)) &
             -T2s%val(6)*(T2s%val(3)*T2s%val(6)-T2s%val(4)*T2s%val(5)) &
             +T2s%val(5)*(T2s%val(4)*T2s%val(6)-T2s%val(2)*T2s%val(5))
  
end function det_Tensor2s
! **********************************************************************
! *** function inv_Tensor2                                           ***
! *** returns the inverse of a general second order tensor           ***
! ***                                                                ***
! *** input : T2, second order tensor to be inverted                 ***
! **********************************************************************
function inv_Tensor2(T2)
 implicit none
 type(Tensor2), intent(in) :: T2
 type(Tensor2) :: inv_Tensor2
 real(8) :: det

 det=det_Tensor2(T2)
 
 inv_Tensor2%val(1,1)=(T2%val(3,3)*T2%val(2,2)-T2%val(3,2)*T2%val(2,3))/det
 inv_Tensor2%val(2,1)=(T2%val(3,1)*T2%val(2,3)-T2%val(3,3)*T2%val(2,1))/det
 inv_Tensor2%val(3,1)=(T2%val(3,2)*T2%val(2,1)-T2%val(3,1)*T2%val(2,2))/det 
 inv_Tensor2%val(1,2)=(T2%val(3,2)*T2%val(1,3)-T2%val(3,3)*T2%val(1,2))/det 
 inv_Tensor2%val(2,2)=(T2%val(3,3)*T2%val(1,1)-T2%val(3,1)*T2%val(1,3))/det 
 inv_Tensor2%val(3,2)=(T2%val(3,1)*T2%val(1,2)-T2%val(3,2)*T2%val(1,1))/det
 inv_Tensor2%val(1,3)=(T2%val(2,3)*T2%val(1,2)-T2%val(2,2)*T2%val(1,3))/det
 inv_Tensor2%val(2,3)=(T2%val(2,1)*T2%val(1,3)-T2%val(2,3)*T2%val(1,1))/det
 inv_Tensor2%val(3,3)=(T2%val(2,2)*T2%val(1,1)-T2%val(2,1)*T2%val(1,2))/det
  
end function inv_Tensor2
! **********************************************************************
! *** function psdinv_Tensor2                                        ***
! *** returns the pseudo-inverse of a general second order tensor    ***
! ***                                                                ***
! *** input : T4, fourth order tensor to be inverted                 ***
! **********************************************************************
function psdinv_Tensor2(T2)
 implicit none
 type(Tensor2), intent(in) :: T2
 type(Tensor2) :: psdinv_Tensor2

 psdinv_Tensor2%val(:,:)=psdinvmat(T2%val(:,:),3,3,3,3)

end function psdinv_Tensor2
! **********************************************************************
! *** function inv_Tensor2s                                          ***
! *** returns the inverse of a symmetric second order tensor         ***
! ***                                                                ***
! *** input : T2s, second order tensor to be inverted                ***
! **********************************************************************
function inv_Tensor2s(T2s)
 implicit none
 type(Tensor2s), intent(in) :: T2s
 type(Tensor2s) :: inv_Tensor2s
 real(8) :: det

 det=det_Tensor2s(T2s)
 
 inv_Tensor2s%val(1)=(T2s%val(3)*T2s%val(2)-T2s%val(4)*T2s%val(4))/det
 inv_Tensor2s%val(2)=(T2s%val(3)*T2s%val(1)-T2s%val(5)*T2s%val(5))/det 
 inv_Tensor2s%val(3)=(T2s%val(2)*T2s%val(1)-T2s%val(6)*T2s%val(6))/det
 inv_Tensor2s%val(4)=(T2s%val(6)*T2s%val(5)-T2s%val(4)*T2s%val(1))/det
 inv_Tensor2s%val(5)=(T2s%val(4)*T2s%val(6)-T2s%val(5)*T2s%val(2))/det 
 inv_Tensor2s%val(6)=(T2s%val(5)*T2s%val(4)-T2s%val(6)*T2s%val(3))/det 
   
end function inv_Tensor2s
! **********************************************************************
! *** function inv_Tensor4                                           ***
! *** returns the inverse of a general fourth order tensor           ***
! ***                                                                ***
! *** input : T4, fourth order tensor to be inverted                 ***
! **********************************************************************
function inv_Tensor4(T4)
 implicit none
 type(Tensor4), intent(in) :: T4
 type(Tensor4) :: inv_Tensor4
 
 inv_Tensor4%val(:,:,:,:)=reshape(inverse('G',reshape(T4%val(:,:,:,:),(/9,9/)),9),(/3,3,3,3/))
   
end function inv_Tensor4
! **********************************************************************
! *** function inv_Tensor4ss                                         ***
! *** returns the inverse of a symmetric fourth order tensor         ***
! ***                                                                ***
! *** input : T4ss, fourth order tensor to be inverted               ***
! **********************************************************************
function inv_Tensor4ss(T4ss)
 implicit none
 type(Tensor4ss), intent(in) :: T4ss
 type(Tensor4ss) :: inv_Tensor4ss
 
 inv_Tensor4ss%val(:,:)=inverse('G',T4ss%val(:,:),6)*invfac(:,:)
 
end function inv_Tensor4ss
! **********************************************************************
! *** function psdinv_Tensor4                                        ***
! *** returns the pseudo-inverse of a general fourth order tensor    ***
! ***                                                                ***
! *** input : T4, fourth order tensor to be inverted                 ***
! **********************************************************************
function psdinv_Tensor4(T4)
 implicit none
 type(Tensor4), intent(in) :: T4
 type(Tensor4) :: psdinv_Tensor4
 
 psdinv_Tensor4%val(:,:,:,:)=reshape(psdinvmat(reshape(T4%val(:,:,:,:),(/9,9/)),9,9,9,9),(/3,3,3,3/))
   
end function psdinv_Tensor4
! **********************************************************************
! *** function transpose_Tensor2                                     ***
! *** returns the transpose of a general second order tensor         ***
! ***                                                                ***
! *** input : T2, second order tensor to be transposed               ***
! **********************************************************************
function transpose_Tensor2(T2)
 implicit none
 type(Tensor2), intent(in) :: T2
 type(Tensor2) :: transpose_Tensor2

 transpose_Tensor2%val(1,1)=T2%val(1,1)
 transpose_Tensor2%val(2,1)=T2%val(1,2)
 transpose_Tensor2%val(3,1)=T2%val(1,3)
 transpose_Tensor2%val(1,2)=T2%val(2,1)
 transpose_Tensor2%val(2,2)=T2%val(2,2)
 transpose_Tensor2%val(3,2)=T2%val(2,3)
 transpose_Tensor2%val(1,3)=T2%val(3,1)
 transpose_Tensor2%val(2,3)=T2%val(3,2)
 transpose_Tensor2%val(3,3)=T2%val(3,3)
 
end function transpose_Tensor2
! **********************************************************************
! *** function transpose_Tensor2w                                    ***
! *** returns the transpose of a skew-symmetric second order tensor  ***
! ***                                                                ***
! *** input : T2w, second order tensor to be transposed              ***
! **********************************************************************
function transpose_Tensor2w(T2w)
 implicit none
 type(Tensor2w), intent(in) :: T2w
 type(Tensor2w) :: transpose_Tensor2w

 transpose_Tensor2w%val(1)=-T2w%val(1)
 transpose_Tensor2w%val(2)=-T2w%val(2)
 transpose_Tensor2w%val(3)=-T2w%val(3)
 
end function transpose_Tensor2w
! **********************************************************************
! *** function transpose_Tensor4                                     ***
! *** returns the (major) transpose of a general fourth order tensor ***
! ***                                                                ***
! *** input : T4, fourth order tensor to be transposed               ***
! **********************************************************************
function transpose_Tensor4(T4)
 implicit none
 type(Tensor4), intent(in) :: T4
 type(Tensor4) :: transpose_Tensor4
 integer :: i,j,k,l
 
 do l=1,3
  do k=1,3
   do j=1,3
    do i=1,3
     transpose_Tensor4%val(i,j,k,l)=T4%val(k,l,i,j)
    enddo
   enddo
  enddo
 enddo
 
end function transpose_Tensor4
! **********************************************************************
! *** function transpose_Tensor4ss                                   ***
! *** returns the (major) transpose of a symmetric fourth order      ***
! *** tensor                                                         ***
! ***                                                                ***
! *** input : T4ss, fourth order tensor to be transposed             ***
! **********************************************************************
function transpose_Tensor4ss(T4ss)
 implicit none
 type(Tensor4ss), intent(in) :: T4ss
 type(Tensor4ss) :: transpose_Tensor4ss
 integer :: i,j
 
 do j=1,6
  do i=1,6
   transpose_Tensor4ss%val(i,j)=T4ss%val(j,i)
  enddo
 enddo
 
end function transpose_Tensor4ss
! **********************************************************************
! *** function righttranspose_Tensor4                                ***
! *** returns the right minor transpose of a general fourth order    ***
! *** tensor                                                         ***
! ***                                                                ***
! *** input : T4, fourth order tensor to be transposed               ***
! **********************************************************************
function righttranspose_Tensor4(T4)
 implicit none
 type(Tensor4), intent(in) :: T4
 type(Tensor4) :: righttranspose_Tensor4
 integer :: i,j,k,l
 
 do l=1,3
  do k=1,3
   do j=1,3
    do i=1,3
     righttranspose_Tensor4%val(i,j,k,l)=T4%val(i,j,l,k)
    enddo
   enddo
  enddo
 enddo
 
end function righttranspose_Tensor4
! **********************************************************************
! *** function lefttranspose_Tensor4                                 ***
! *** returns the left minor transpose of a general fourth order     ***
! *** tensor                                                         ***
! ***                                                                ***
! *** input : T4, fourth order tensor to be transposed               ***
! **********************************************************************
function lefttranspose_Tensor4(T4)
 implicit none
 type(Tensor4), intent(in) :: T4
 type(Tensor4) :: lefttranspose_Tensor4
 integer :: i,j,k,l
 
 do l=1,3
  do k=1,3
   do j=1,3
    do i=1,3
     lefttranspose_Tensor4%val(i,j,k,l)=T4%val(j,i,k,l)
    enddo
   enddo
  enddo
 enddo
 
end function lefttranspose_Tensor4
! **********************************************************************
! *** function middletranspose_Tensor4                               ***
! *** returns the middle minor transpose of a general fourth order   ***
! *** tensor                                                         ***
! ***                                                                ***
! *** input : T4, fourth order tensor to be transposed               ***
! **********************************************************************
function middletranspose_Tensor4(T4)
 implicit none
 type(Tensor4), intent(in) :: T4
 type(Tensor4) :: middletranspose_Tensor4
 integer :: i,j,k,l
 
 do l=1,3
  do k=1,3
   do j=1,3
    do i=1,3
     middletranspose_Tensor4%val(i,j,k,l)=T4%val(i,k,j,l)
    enddo
   enddo
  enddo
 enddo
 
end function middletranspose_Tensor4
! **********************************************************************
! *** function middletranspose_Tensor4ss                             ***
! *** returns the middle minor transpose of a general fourth order   ***
! *** tensor                                                         ***
! ***                                                                ***
! *** input : T4ss, fourth order tensor to be transposed             ***
! **********************************************************************
function middletranspose_Tensor4ss(T4ss)
 implicit none
 type(Tensor4ss), intent(in) :: T4ss
 type(Tensor4) :: middletranspose_Tensor4ss
 integer :: i,j,k,l,ij,kl,ik,jl
 
 do l=1,3
  do k=1,3
   do j=1,3
    jl=indxs(j,l)
    do i=1,3
     ik=indxs(i,k)
     middletranspose_Tensor4ss%val(i,j,k,l)=T4ss%val(ik,jl)
    enddo
   enddo
  enddo
 enddo
 
end function middletranspose_Tensor4ss
! **********************************************************************
! *** function plus_T1T1                                             ***
! *** returns the sum of two first order tensors                     ***
! ***                                                                ***
! *** input : T1A, first order tensor to be added                    ***
! ***         T1B, first order tensor to be added                    ***
! **********************************************************************
function plus_T1T1(T1A,T1B)
 implicit none
 type(Tensor1), intent(in) :: T1A,T1B
 type(Tensor1) :: plus_T1T1

 plus_T1T1%val(1)=T1A%val(1)+T1B%val(1)
 plus_T1T1%val(2)=T1A%val(2)+T1B%val(2)
 plus_T1T1%val(3)=T1A%val(3)+T1B%val(3)
 
end function plus_T1T1
! **********************************************************************
! *** function plus_T2T2                                             ***
! *** returns the sum of two second order tensors                    ***
! ***                                                                ***
! *** input : T2A, second order tensor to be added                   ***
! ***         T2B, second order tensor to be added                   ***
! **********************************************************************
function plus_T2T2(T2A,T2B)
 implicit none
 type(Tensor2), intent(in) :: T2A,T2B
 type(Tensor2) :: plus_T2T2

 plus_T2T2%val(1,1)=T2A%val(1,1)+T2B%val(1,1)
 plus_T2T2%val(2,1)=T2A%val(2,1)+T2B%val(2,1)
 plus_T2T2%val(3,1)=T2A%val(3,1)+T2B%val(3,1)
 plus_T2T2%val(1,2)=T2A%val(1,2)+T2B%val(1,2)
 plus_T2T2%val(2,2)=T2A%val(2,2)+T2B%val(2,2)
 plus_T2T2%val(3,2)=T2A%val(3,2)+T2B%val(3,2)
 plus_T2T2%val(1,3)=T2A%val(1,3)+T2B%val(1,3)
 plus_T2T2%val(2,3)=T2A%val(2,3)+T2B%val(2,3)
 plus_T2T2%val(3,3)=T2A%val(3,3)+T2B%val(3,3)
 
end function plus_T2T2
! **********************************************************************
! *** function plus_T2sT2                                            ***
! *** returns the sum of a symmetric and a general second order      ***
! *** tensors                                                        ***
! ***                                                                ***
! *** input : T2sA, second order tensor to be added                  ***
! ***         T2B, second order tensor to be added                   ***
! **********************************************************************
function plus_T2sT2(T2sA,T2B)
 implicit none
 type(Tensor2s), intent(in) :: T2sA
 type(Tensor2), intent(in) :: T2B
 type(Tensor2) :: plus_T2sT2

 plus_T2sT2%val(1,1)=T2sA%val(1)+T2B%val(1,1)
 plus_T2sT2%val(2,1)=T2sA%val(6)+T2B%val(2,1)
 plus_T2sT2%val(3,1)=T2sA%val(5)+T2B%val(3,1)
 plus_T2sT2%val(1,2)=T2sA%val(6)+T2B%val(1,2)
 plus_T2sT2%val(2,2)=T2sA%val(2)+T2B%val(2,2)
 plus_T2sT2%val(3,2)=T2sA%val(4)+T2B%val(3,2)
 plus_T2sT2%val(1,3)=T2sA%val(5)+T2B%val(1,3)
 plus_T2sT2%val(2,3)=T2sA%val(4)+T2B%val(2,3)
 plus_T2sT2%val(3,3)=T2sA%val(3)+T2B%val(3,3)
 
end function plus_T2sT2
! **********************************************************************
! *** function plus_T2T2s                                            ***
! *** returns the sum of a general and a symmetric second order      ***
! *** tensors                                                        ***
! ***                                                                ***
! *** input : T2A, second order tensor to be added                   ***
! ***         T2sB, second order tensor to be added                  ***
! **********************************************************************
function plus_T2T2s(T2A,T2sB)
 implicit none
 type(Tensor2), intent(in) :: T2A
 type(Tensor2s), intent(in) :: T2sB
 type(Tensor2) :: plus_T2T2s

 plus_T2T2s%val(1,1)=T2A%val(1,1)+T2sB%val(1)
 plus_T2T2s%val(2,1)=T2A%val(2,1)+T2sB%val(6)
 plus_T2T2s%val(3,1)=T2A%val(3,1)+T2sB%val(5)
 plus_T2T2s%val(1,2)=T2A%val(1,2)+T2sB%val(6)
 plus_T2T2s%val(2,2)=T2A%val(2,2)+T2sB%val(2)
 plus_T2T2s%val(3,2)=T2A%val(3,2)+T2sB%val(4)
 plus_T2T2s%val(1,3)=T2A%val(1,3)+T2sB%val(5)
 plus_T2T2s%val(2,3)=T2A%val(2,3)+T2sB%val(4)
 plus_T2T2s%val(3,3)=T2A%val(3,3)+T2sB%val(3)
 
end function plus_T2T2s
! **********************************************************************
! *** function plus_T2sT2s                                           ***
! *** returns the sum of two symmetric second order tensors          ***
! ***                                                                ***
! *** input : T2sA, second order tensor to be added                  ***
! ***         T2sB, second order tensor to be added                  ***
! **********************************************************************
function plus_T2sT2s(T2sA,T2sB)
 implicit none
 type(Tensor2s), intent(in) :: T2sA,T2sB
 type(Tensor2s) :: plus_T2sT2s

 plus_T2sT2s%val(1)=T2sA%val(1)+T2sB%val(1)
 plus_T2sT2s%val(2)=T2sA%val(2)+T2sB%val(2)
 plus_T2sT2s%val(3)=T2sA%val(3)+T2sB%val(3)
 plus_T2sT2s%val(4)=T2sA%val(4)+T2sB%val(4)
 plus_T2sT2s%val(5)=T2sA%val(5)+T2sB%val(5)
 plus_T2sT2s%val(6)=T2sA%val(6)+T2sB%val(6)
 
end function plus_T2sT2s
! **********************************************************************
! *** function plus_T2wT2                                            ***
! *** returns the sum of a skew-symmetric and a general second order ***
! *** tensors                                                        ***
! ***                                                                ***
! *** input : T2wA, second order tensor to be added                  ***
! ***         T2B, second order tensor to be added                   ***
! **********************************************************************
function plus_T2wT2(T2wA,T2B)
 implicit none
 type(Tensor2w), intent(in) :: T2wA
 type(Tensor2), intent(in) :: T2B
 type(Tensor2) :: plus_T2wT2

 plus_T2wT2%val(1,1)=T2B%val(1,1)
 plus_T2wT2%val(2,1)=-T2wA%val(3)+T2B%val(2,1)
 plus_T2wT2%val(3,1)=T2wA%val(2)+T2B%val(3,1)
 plus_T2wT2%val(1,2)=T2wA%val(3)+T2B%val(1,2)
 plus_T2wT2%val(2,2)=T2B%val(2,2)
 plus_T2wT2%val(3,2)=-T2wA%val(1)+T2B%val(3,2)
 plus_T2wT2%val(1,3)=-T2wA%val(2)+T2B%val(1,3)
 plus_T2wT2%val(2,3)=T2wA%val(1)+T2B%val(2,3)
 plus_T2wT2%val(3,3)=T2B%val(3,3)
 
end function plus_T2wT2
! **********************************************************************
! *** function plus_T2T2w                                            ***
! *** returns the sum of a general and a skew-symmetric second order ***
! *** tensors                                                        ***
! ***                                                                ***
! *** input : T2A, second order tensor to be added                   ***
! ***         T2wB, second order tensor to be added                  ***
! **********************************************************************
function plus_T2T2w(T2A,T2wB)
 implicit none
 type(Tensor2), intent(in) :: T2A
 type(Tensor2w), intent(in) :: T2wB
 type(Tensor2) :: plus_T2T2w

 plus_T2T2w%val(1,1)=T2A%val(1,1)
 plus_T2T2w%val(2,1)=T2A%val(2,1)-T2wB%val(3)
 plus_T2T2w%val(3,1)=T2A%val(3,1)+T2wB%val(2)
 plus_T2T2w%val(1,2)=T2A%val(1,2)+T2wB%val(3)
 plus_T2T2w%val(2,2)=T2A%val(2,2)
 plus_T2T2w%val(3,2)=T2A%val(3,2)-T2wB%val(1)
 plus_T2T2w%val(1,3)=T2A%val(1,3)-T2wB%val(2)
 plus_T2T2w%val(2,3)=T2A%val(2,3)+T2wB%val(1)
 plus_T2T2w%val(3,3)=T2A%val(3,3)
 
end function plus_T2T2w
! **********************************************************************
! *** function plus_T2wT2w                                           ***
! *** returns the sum of two skew-symmetric second order tensors     ***
! ***                                                                ***
! *** input : T2wA, second order tensor to be added                  ***
! ***         T2wB, second order tensor to be added                  ***
! **********************************************************************
function plus_T2wT2w(T2wA,T2wB)
 implicit none
 type(Tensor2w), intent(in) :: T2wA,T2wB
 type(Tensor2w) :: plus_T2wT2w

 plus_T2wT2w%val(1)=T2wA%val(1)+T2wB%val(1)
 plus_T2wT2w%val(2)=T2wA%val(2)+T2wB%val(2)
 plus_T2wT2w%val(3)=T2wA%val(3)+T2wB%val(3)
 
end function plus_T2wT2w
! **********************************************************************
! *** function plus_T2sT2w                                           ***
! *** returns the sum of one symmetric and one skew-symmetric second ***
! *** order tensors                                                  ***
! ***                                                                ***
! *** input : T2sA, second order tensor to be added                  ***
! ***         T2wB, second order tensor to be added                  ***
! **********************************************************************
function plus_T2sT2w(T2sA,T2wB)
 implicit none
 type(Tensor2s), intent(in) :: T2sA
 type(Tensor2w), intent(in) :: T2wB
 type(Tensor2) :: plus_T2sT2w

 plus_T2sT2w%val(1,1)=T2sA%val(1)
 plus_T2sT2w%val(2,1)=T2sA%val(6)-T2wB%val(3)
 plus_T2sT2w%val(3,1)=T2sA%val(5)+T2wB%val(2)
 plus_T2sT2w%val(1,2)=T2sA%val(6)+T2wB%val(3)
 plus_T2sT2w%val(2,2)=T2sA%val(2)
 plus_T2sT2w%val(3,2)=T2sA%val(4)-T2wB%val(1)
 plus_T2sT2w%val(1,3)=T2sA%val(5)-T2wB%val(2)
 plus_T2sT2w%val(2,3)=T2sA%val(4)+T2wB%val(1)
 plus_T2sT2w%val(3,3)=T2sA%val(3)
 
end function plus_T2sT2w
! **********************************************************************
! *** function plus_T2wT2s                                           ***
! *** returns the sum of one skew-symmetric and one symmetric second ***
! *** order tensors                                                  ***
! ***                                                                ***
! *** input : T2wA, second order tensor to be added                  ***
! ***         T2sB, second order tensor to be added                  ***
! **********************************************************************
function plus_T2wT2s(T2wA,T2sB)
 implicit none
 type(Tensor2w), intent(in) :: T2wA
 type(Tensor2s), intent(in) :: T2sB
 type(Tensor2) :: plus_T2wT2s

 plus_T2wT2s%val(1,1)=T2sB%val(1)
 plus_T2wT2s%val(2,1)=-T2wA%val(3)+T2sB%val(6)
 plus_T2wT2s%val(3,1)=T2wA%val(2)+T2sB%val(5)
 plus_T2wT2s%val(1,2)=T2wA%val(3)+T2sB%val(6)
 plus_T2wT2s%val(2,2)=T2sB%val(2)
 plus_T2wT2s%val(3,2)=-T2wA%val(1)+T2sB%val(4)
 plus_T2wT2s%val(1,3)=-T2wA%val(2)+T2sB%val(5)
 plus_T2wT2s%val(2,3)=T2wA%val(1)+T2sB%val(4)
 plus_T2wT2s%val(3,3)=T2sB%val(3)
 
end function plus_T2wT2s
! **********************************************************************
! *** function plus_T3T3                                             ***
! *** returns the sum of two third order tensors                     ***
! ***                                                                ***
! *** input : T3A, third order tensor to be added                    ***
! ***         T3B, third order tensor to be added                    ***
! **********************************************************************
function plus_T3T3(T3A,T3B)
 implicit none
 type(Tensor3), intent(in) :: T3A,T3B
 type(Tensor3) :: plus_T3T3
 integer :: i,j,k
 
 do k=1,3
  do j=1,3
   do i=1,3
    plus_T3T3%val(i,j,k)=T3A%val(i,j,k)+T3B%val(i,j,k)
   enddo
  enddo
 enddo
 
end function plus_T3T3
! **********************************************************************
! *** function plus_T4T4                                             ***
! *** returns the sum of two fourth order tensors                    ***
! ***                                                                ***
! *** input : T4A, fourth order tensor to be added                   ***
! ***         T4B, fourth order tensor to be added                   ***
! **********************************************************************
function plus_T4T4(T4A,T4B)
 implicit none
 type(Tensor4), intent(in) :: T4A,T4B
 type(Tensor4) :: plus_T4T4
 integer :: i,j,k,l
 
 do l=1,3
  do k=1,3
   do j=1,3
    do i=1,3
     plus_T4T4%val(i,j,k,l)=T4A%val(i,j,k,l)+T4B%val(i,j,k,l)
    enddo
   enddo
  enddo
 enddo
 
end function plus_T4T4
! **********************************************************************
! *** function plus_T4ssT4ss                                         ***
! *** returns the sum of two symmetric fourth order tensors          ***
! ***                                                                ***
! *** input : T4ssA, fourth order tensor to be added                 ***
! ***         T4ssB, fourth order tensor to be added                 ***
! **********************************************************************
function plus_T4ssT4ss(T4ssA,T4ssB)
 implicit none
 type(Tensor4ss), intent(in) :: T4ssA,T4ssB
 type(Tensor4ss) :: plus_T4ssT4ss
 integer :: i,j
 
 do i=1,6
  do j=1,6
   plus_T4ssT4ss%val(i,j)=T4ssA%val(i,j)+T4ssB%val(i,j)
  enddo
 enddo

end function plus_T4ssT4ss
! **********************************************************************
! *** function plus_T4T4ss                                           ***
! *** returns the sum of one general and one symmetric fourth order  ***
! *** tensors                                                        ***
! ***                                                                ***
! *** input : T4A, fourth order tensor to be added                   ***
! ***         T4ssB, fourth order tensor to be added                 ***
! **********************************************************************
function plus_T4T4ss(T4A,T4ssB)
 implicit none
 type(Tensor4), intent(in) :: T4A
 type(Tensor4ss), intent(in) :: T4ssB
 type(Tensor4) :: plus_T4T4ss
 integer :: i,j,k,l,ij,kl
 
 do l=1,3
  do k=1,3
   kl=indxs(k,l)
   do j=1,3
    do i=1,3
     ij=indxs(i,j)
     plus_T4T4ss%val(i,j,k,l)=T4A%val(i,j,k,l)+T4ssB%val(ij,kl)
    enddo
   enddo
  enddo
 enddo
 
end function plus_T4T4ss
! **********************************************************************
! *** function plus_T4ssT4                                           ***
! *** returns the sum of one symmetric and one general fourth order  ***
! *** tensors                                                        ***
! ***                                                                ***
! *** input : T4ssA, fourth order tensor to be added                 ***
! ***         T4B, fourth order tensor to be added                   ***
! **********************************************************************
function plus_T4ssT4(T4ssA,T4B)
 implicit none
 type(Tensor4ss), intent(in) :: T4ssA
 type(Tensor4), intent(in) :: T4B
 type(Tensor4) :: plus_T4ssT4
 integer :: i,j,k,l,ij,kl
 
 do l=1,3
  do k=1,3
   kl=indxs(k,l)
   do j=1,3
    do i=1,3
     ij=indxs(i,j)
     plus_T4ssT4%val(i,j,k,l)=T4ssA%val(ij,kl)+T4B%val(i,j,k,l)
    enddo
   enddo
  enddo
 enddo
 
end function plus_T4ssT4
! **********************************************************************
! *** function minus_T1T1                                            ***
! *** returns the difference between two first order tensors         ***
! ***                                                                ***
! *** input : T1A, minuend first order tensor                        ***
! ***         T1B, subtrahend first order tensor                     ***
! **********************************************************************
function minus_T1T1(T1A,T1B)
 implicit none
 type(Tensor1), intent(in) :: T1A,T1B
 type(Tensor1) :: minus_T1T1

 minus_T1T1%val(1)=T1A%val(1)-T1B%val(1)
 minus_T1T1%val(2)=T1A%val(2)-T1B%val(2)
 minus_T1T1%val(3)=T1A%val(3)-T1B%val(3)
 
end function minus_T1T1
! **********************************************************************
! *** function minus_T2T2                                            ***
! *** returns the difference between two second order tensors        ***
! ***                                                                ***
! *** input : T2A, minuend second order tensor                       ***
! ***         T2B, subtrahend second order tensor                    ***
! **********************************************************************
function minus_T2T2(T2A,T2B)
 implicit none
 type(Tensor2), intent(in) :: T2A,T2B
 type(Tensor2) :: minus_T2T2

 minus_T2T2%val(1,1)=T2A%val(1,1)-T2B%val(1,1)
 minus_T2T2%val(2,1)=T2A%val(2,1)-T2B%val(2,1)
 minus_T2T2%val(3,1)=T2A%val(3,1)-T2B%val(3,1)
 minus_T2T2%val(1,2)=T2A%val(1,2)-T2B%val(1,2)
 minus_T2T2%val(2,2)=T2A%val(2,2)-T2B%val(2,2)
 minus_T2T2%val(3,2)=T2A%val(3,2)-T2B%val(3,2)
 minus_T2T2%val(1,3)=T2A%val(1,3)-T2B%val(1,3)
 minus_T2T2%val(2,3)=T2A%val(2,3)-T2B%val(2,3)
 minus_T2T2%val(3,3)=T2A%val(3,3)-T2B%val(3,3)
 
end function minus_T2T2
! **********************************************************************
! *** function minus_T2sT2                                           ***
! *** returns the difference between a symmetric and a general       ***
! *** second order tensors                                           ***
! ***                                                                ***
! *** input : T2sA, minuend second order tensor                      ***
! ***         T2B, subtrahend second order tensor                    ***
! **********************************************************************
function minus_T2sT2(T2sA,T2B)
 implicit none
 type(Tensor2s), intent(in) :: T2sA
 type(Tensor2), intent(in) :: T2B
 type(Tensor2) :: minus_T2sT2

 minus_T2sT2%val(1,1)=T2sA%val(1)-T2B%val(1,1)
 minus_T2sT2%val(2,1)=T2sA%val(6)-T2B%val(2,1)
 minus_T2sT2%val(3,1)=T2sA%val(5)-T2B%val(3,1)
 minus_T2sT2%val(1,2)=T2sA%val(6)-T2B%val(1,2)
 minus_T2sT2%val(2,2)=T2sA%val(2)-T2B%val(2,2)
 minus_T2sT2%val(3,2)=T2sA%val(4)-T2B%val(3,2)
 minus_T2sT2%val(1,3)=T2sA%val(5)-T2B%val(1,3)
 minus_T2sT2%val(2,3)=T2sA%val(4)-T2B%val(2,3)
 minus_T2sT2%val(3,3)=T2sA%val(3)-T2B%val(3,3)
 
end function minus_T2sT2
! **********************************************************************
! *** function minus_T2T2s                                           ***
! *** returns the difference between a general and a symmetric       ***
! *** second order tensors                                           ***
! ***                                                                ***
! *** input : T2A, minuend second order tensor                       ***
! ***         T2sB, subtrahend second order tensor                   ***
! **********************************************************************
function minus_T2T2s(T2A,T2sB)
 implicit none
 type(Tensor2), intent(in) :: T2A
 type(Tensor2s), intent(in) :: T2sB
 type(Tensor2) :: minus_T2T2s

 minus_T2T2s%val(1,1)=T2A%val(1,1)-T2sB%val(1)
 minus_T2T2s%val(2,1)=T2A%val(2,1)-T2sB%val(6)
 minus_T2T2s%val(3,1)=T2A%val(3,1)-T2sB%val(5)
 minus_T2T2s%val(1,2)=T2A%val(1,2)-T2sB%val(6)
 minus_T2T2s%val(2,2)=T2A%val(2,2)-T2sB%val(2)
 minus_T2T2s%val(3,2)=T2A%val(3,2)-T2sB%val(4)
 minus_T2T2s%val(1,3)=T2A%val(1,3)-T2sB%val(5)
 minus_T2T2s%val(2,3)=T2A%val(2,3)-T2sB%val(4)
 minus_T2T2s%val(3,3)=T2A%val(3,3)-T2sB%val(3)
 
end function minus_T2T2s
! **********************************************************************
! *** function minus_T2sT2s                                          ***
! *** returns the difference between two symmetric second order      ***
! *** tensors                                                        ***
! ***                                                                ***
! *** input : T2sA, minuend second order tensor                      ***
! ***         T2sB, subtrahend second order tensor                   ***
! **********************************************************************
function minus_T2sT2s(T2sA,T2sB)
 implicit none
 type(Tensor2s), intent(in) :: T2sA,T2sB
 type(Tensor2s) :: minus_T2sT2s

 minus_T2sT2s%val(1)=T2sA%val(1)-T2sB%val(1)
 minus_T2sT2s%val(2)=T2sA%val(2)-T2sB%val(2)
 minus_T2sT2s%val(3)=T2sA%val(3)-T2sB%val(3)
 minus_T2sT2s%val(4)=T2sA%val(4)-T2sB%val(4)
 minus_T2sT2s%val(5)=T2sA%val(5)-T2sB%val(5)
 minus_T2sT2s%val(6)=T2sA%val(6)-T2sB%val(6)
 
end function minus_T2sT2s
! **********************************************************************
! *** function minus_T2wT2                                           ***
! *** returns the difference between a skew-symmetric and a general  ***
! *** second order tensors                                           ***
! ***                                                                ***
! *** input : T2wA, minuend second order tensor                      ***
! ***         T2B, subtrahend second order tensor                    ***
! **********************************************************************
function minus_T2wT2(T2wA,T2B)
 implicit none
 type(Tensor2w), intent(in) :: T2wA
 type(Tensor2), intent(in) :: T2B
 type(Tensor2) :: minus_T2wT2

 minus_T2wT2%val(1,1)=-T2B%val(1,1)
 minus_T2wT2%val(2,1)=-T2wA%val(3)-T2B%val(2,1)
 minus_T2wT2%val(3,1)=T2wA%val(2)-T2B%val(3,1)
 minus_T2wT2%val(1,2)=T2wA%val(3)-T2B%val(1,2)
 minus_T2wT2%val(2,2)=-T2B%val(2,2)
 minus_T2wT2%val(3,2)=-T2wA%val(1)-T2B%val(3,2)
 minus_T2wT2%val(1,3)=-T2wA%val(2)-T2B%val(1,3)
 minus_T2wT2%val(2,3)=T2wA%val(1)-T2B%val(2,3)
 minus_T2wT2%val(3,3)=-T2B%val(3,3)
 
end function minus_T2wT2
! **********************************************************************
! *** function minus_T2T2w                                           ***
! *** returns the difference between a general and a skew-symmetric  ***
! *** second order tensors                                           ***
! ***                                                                ***
! *** input : T2A, minuend second order tensor                       ***
! ***         T2wB, subtrahend second order tensor                   ***
! **********************************************************************
function minus_T2T2w(T2A,T2wB)
 implicit none
 type(Tensor2), intent(in) :: T2A
 type(Tensor2w), intent(in) :: T2wB
 type(Tensor2) :: minus_T2T2w

 minus_T2T2w%val(1,1)=T2A%val(1,1)
 minus_T2T2w%val(2,1)=T2A%val(2,1)+T2wB%val(3)
 minus_T2T2w%val(3,1)=T2A%val(3,1)-T2wB%val(2)
 minus_T2T2w%val(1,2)=T2A%val(1,2)-T2wB%val(3)
 minus_T2T2w%val(2,2)=T2A%val(2,2)
 minus_T2T2w%val(3,2)=T2A%val(3,2)+T2wB%val(1)
 minus_T2T2w%val(1,3)=T2A%val(1,3)+T2wB%val(2)
 minus_T2T2w%val(2,3)=T2A%val(2,3)-T2wB%val(1)
 minus_T2T2w%val(3,3)=T2A%val(3,3)
 
end function minus_T2T2w
! **********************************************************************
! *** function minus_T2wT2w                                          ***
! *** returns the difference between two skew-symmetric second order ***
! *** tensors                                                        ***
! ***                                                                ***
! *** input : T2wA, second order tensor                              ***
! ***         T2wB, second order tensor                              ***
! **********************************************************************
function minus_T2wT2w(T2wA,T2wB)
 implicit none
 type(Tensor2w), intent(in) :: T2wA,T2wB
 type(Tensor2w) :: minus_T2wT2w

 minus_T2wT2w%val(1)=T2wA%val(1)-T2wB%val(1)
 minus_T2wT2w%val(2)=T2wA%val(2)-T2wB%val(2)
 minus_T2wT2w%val(3)=T2wA%val(3)-T2wB%val(3)
 
end function minus_T2wT2w
! **********************************************************************
! *** function minus_T2sT2w                                          ***
! *** returns the difference between one symmetric and one           ***
! *** skew-symmetric second order tensors                            ***
! ***                                                                ***
! *** input : T2sA, minuend second order tensor                      ***
! ***         T2wB, subtrahend second order tensor                   ***
! **********************************************************************
function minus_T2sT2w(T2sA,T2wB)
 implicit none
 type(Tensor2s), intent(in) :: T2sA
 type(Tensor2w), intent(in) :: T2wB
 type(Tensor2) :: minus_T2sT2w

 minus_T2sT2w%val(1,1)=T2sA%val(1)
 minus_T2sT2w%val(2,1)=T2sA%val(6)+T2wB%val(3)
 minus_T2sT2w%val(3,1)=T2sA%val(5)-T2wB%val(2)
 minus_T2sT2w%val(1,2)=T2sA%val(6)-T2wB%val(3)
 minus_T2sT2w%val(2,2)=T2sA%val(2)
 minus_T2sT2w%val(3,2)=T2sA%val(4)+T2wB%val(1)
 minus_T2sT2w%val(1,3)=T2sA%val(5)+T2wB%val(2)
 minus_T2sT2w%val(2,3)=T2sA%val(4)-T2wB%val(1)
 minus_T2sT2w%val(3,3)=T2sA%val(3)
 
end function minus_T2sT2w
! **********************************************************************
! *** function minus_T2wT2s                                          ***
! *** returns the difference between one skew-symmetric and one      ***
! *** symmetric second order tensors                                 ***
! ***                                                                ***
! *** input : T2wA, minuend second order tensor                      ***
! ***         T2sB, subtrahend second order tensor                   ***
! **********************************************************************
function minus_T2wT2s(T2wA,T2sB)
 implicit none
 type(Tensor2w), intent(in) :: T2wA
 type(Tensor2s), intent(in) :: T2sB
 type(Tensor2) :: minus_T2wT2s

 minus_T2wT2s%val(1,1)=-T2sB%val(1)
 minus_T2wT2s%val(2,1)=-T2wA%val(3)-T2sB%val(6)
 minus_T2wT2s%val(3,1)=T2wA%val(2)-T2sB%val(5)
 minus_T2wT2s%val(1,2)=T2wA%val(3)-T2sB%val(6)
 minus_T2wT2s%val(2,2)=-T2sB%val(2)
 minus_T2wT2s%val(3,2)=-T2wA%val(1)-T2sB%val(4)
 minus_T2wT2s%val(1,3)=-T2wA%val(2)-T2sB%val(5)
 minus_T2wT2s%val(2,3)=T2wA%val(1)-T2sB%val(4)
 minus_T2wT2s%val(3,3)=-T2sB%val(3)
 
end function minus_T2wT2s
! **********************************************************************
! *** function minus_T3T3                                            ***
! *** returns the difference between two third order tensors         ***
! ***                                                                ***
! *** input : T3A, minuend third order tensor                        ***
! ***         T3B, subtrahend third order tensor                     ***
! **********************************************************************
function minus_T3T3(T3A,T3B)
 implicit none
 type(Tensor3), intent(in) :: T3A,T3B
 type(Tensor3) :: minus_T3T3
 integer :: i,j,k
 
 do k=1,3
  do j=1,3
   do i=1,3
    minus_T3T3%val(i,j,k)=T3A%val(i,j,k)-T3B%val(i,j,k)
   enddo
  enddo
 enddo
 
end function minus_T3T3
! **********************************************************************
! *** function minus_T4T4                                            ***
! *** returns the difference between two fourth order tensors        ***
! ***                                                                ***
! *** input : T4A, minuend fourth order tensor                       ***
! ***         T4B, subtrahend fourth order tensor                    ***
! **********************************************************************
function minus_T4T4(T4A,T4B)
 implicit none
 type(Tensor4), intent(in) :: T4A,T4B
 type(Tensor4) :: minus_T4T4
 integer :: i,j,k,l
 
 do l=1,3
  do k=1,3
   do j=1,3
    do i=1,3
     minus_T4T4%val(i,j,k,l)=T4A%val(i,j,k,l)-T4B%val(i,j,k,l)
    enddo
   enddo
  enddo
 enddo
 
end function minus_T4T4
! **********************************************************************
! *** function minus_T4ssT4ss                                        ***
! *** returns the difference between two symmetric fourth order      ***
! *** tensors                                                        ***
! ***                                                                ***
! *** input : T4ssA, minuend fourth order tensor                     ***
! ***         T4ssB, subtrahend fourth order tensor                  ***
! **********************************************************************
function minus_T4ssT4ss(T4ssA,T4ssB)
 implicit none
 type(Tensor4ss), intent(in) :: T4ssA,T4ssB
 type(Tensor4ss) :: minus_T4ssT4ss
 integer :: i,j
 
 do j=1,6
  do i=1,6
   minus_T4ssT4ss%val(i,j)=T4ssA%val(i,j)-T4ssB%val(i,j)
  enddo
 enddo
 
end function minus_T4ssT4ss
! **********************************************************************
! *** function minus_T4T4ss                                          ***
! *** returns the difference between one general and one symmetric   ***
! *** fourth order tensors                                           ***
! ***                                                                ***
! *** input : T4A, minuend fourth order tensor                       ***
! ***         T4ssB, subtrahend fourth order tensor                  ***
! **********************************************************************
function minus_T4T4ss(T4A,T4ssB)
 implicit none
 type(Tensor4), intent(in) :: T4A
 type(Tensor4ss), intent(in) :: T4ssB
 type(Tensor4) :: minus_T4T4ss
 integer :: i,j,k,l,ij,kl
 
 do l=1,3
  do k=1,3
   kl=indxs(k,l)
   do j=1,3
    do i=1,3
     ij=indxs(i,j)
     minus_T4T4ss%val(i,j,k,l)=T4A%val(i,j,k,l)-T4ssB%val(ij,kl)
    enddo
   enddo
  enddo
 enddo

end function minus_T4T4ss
! **********************************************************************
! *** function minus_T4ssT4                                          ***
! *** returns the difference between one symmetric and one general   ***
! *** fourth order tensors                                           ***
! ***                                                                ***
! *** input : T4ssA, minuend fourth order tensor                     ***
! ***         T4B, subtrahend fourth order tensor                    ***
! **********************************************************************
function minus_T4ssT4(T4ssA,T4B)
 implicit none
 type(Tensor4ss), intent(in) :: T4ssA
 type(Tensor4), intent(in) :: T4B
 type(Tensor4) :: minus_T4ssT4
 integer :: i,j,k,l,ij,kl
 
 do l=1,3
  do k=1,3
   kl=indxs(k,l)
   do j=1,3
    do i=1,3
     ij=indxs(i,j)
     minus_T4ssT4%val(i,j,k,l)=T4ssA%val(ij,kl)-T4B%val(i,j,k,l)
    enddo
   enddo
  enddo
 enddo

end function minus_T4ssT4
! **********************************************************************
! *** function product_T0T1                                          ***
! *** returns the product between one scalar value and one first     ***
! *** order tensor                                                   ***
! ***                                                                ***
! *** input : T0, scalar value                                       ***
! ***         T1, first order tensor                                 ***
! **********************************************************************
function product_T0T1(T0,T1)
 implicit none
 real(8), intent(in) :: T0
 type(Tensor1), intent(in) :: T1
 type(Tensor1) :: product_T0T1
 
 product_T0T1%val(1)=T1%val(1)*T0
 product_T0T1%val(2)=T1%val(2)*T0
 product_T0T1%val(3)=T1%val(3)*T0

end function product_T0T1
! **********************************************************************
! *** function product_T1T0                                          ***
! *** returns the product between one scalar value and one first     ***
! *** order tensor                                                   ***
! ***                                                                ***
! *** input : T1, first order tensor                                 ***
! ***         T0, scaling coefficient                                ***
! **********************************************************************
function product_T1T0(T1,T0)
 implicit none
 type(Tensor1), intent(in) :: T1
 real(8), intent(in) :: T0
 type(Tensor1) :: product_T1T0
 
 product_T1T0=product_T0T1(T0,T1)

end function product_T1T0
! **********************************************************************
! *** function product_T0T2                                          ***
! *** returns the product between one scalar value and one second    ***
! *** order tensor                                                   ***
! ***                                                                ***
! *** input : T0, scalar value                                       ***
! ***         T2, second order tensor                                ***
! **********************************************************************
function product_T0T2(T0,T2)
 implicit none
 real(8), intent(in) :: T0
 type(Tensor2), intent(in) :: T2
 type(Tensor2) :: product_T0T2
 
 product_T0T2%val(1,1)=T2%val(1,1)*T0
 product_T0T2%val(2,1)=T2%val(2,1)*T0
 product_T0T2%val(3,1)=T2%val(3,1)*T0
 product_T0T2%val(1,2)=T2%val(1,2)*T0
 product_T0T2%val(2,2)=T2%val(2,2)*T0
 product_T0T2%val(3,2)=T2%val(3,2)*T0
 product_T0T2%val(1,3)=T2%val(1,3)*T0
 product_T0T2%val(2,3)=T2%val(2,3)*T0
 product_T0T2%val(3,3)=T2%val(3,3)*T0

end function product_T0T2
! **********************************************************************
! *** function product_T2T0                                          ***
! *** returns the product between one scalar value and one second    ***
! *** order tensor                                                   ***
! ***                                                                ***
! *** input : T2, second order tensor                                ***
! ***         T0, scaling coefficient                                ***
! **********************************************************************
function product_T2T0(T2,T0)
 implicit none
 type(Tensor2), intent(in) :: T2
 real(8), intent(in) :: T0
 type(Tensor2) :: product_T2T0
 
 product_T2T0=product_T0T2(T0,T2)

end function product_T2T0
! **********************************************************************
! *** function product_T0T2s                                         ***
! *** returns the product between one scalar value and one symmetric ***
! *** second order tensor                                            ***
! ***                                                                ***
! *** input : T0, scalar value                                       ***
! ***         T2s, second order tensor                               ***
! **********************************************************************
function product_T0T2s(T0,T2s)
 implicit none
 real(8), intent(in) :: T0
 type(Tensor2s), intent(in) :: T2s
 type(Tensor2s) :: product_T0T2s
 
 product_T0T2s%val(1)=T2s%val(1)*T0
 product_T0T2s%val(2)=T2s%val(2)*T0
 product_T0T2s%val(3)=T2s%val(3)*T0
 product_T0T2s%val(4)=T2s%val(4)*T0
 product_T0T2s%val(5)=T2s%val(5)*T0
 product_T0T2s%val(6)=T2s%val(6)*T0

end function product_T0T2s
! **********************************************************************
! *** function product_T2sT0                                         ***
! *** returns the product between one scalar value and one symmetric ***
! *** second order tensor                                            ***
! ***                                                                ***
! *** input : T2s, second order tensor                               ***
! ***         T0, scaling coefficient                                ***
! **********************************************************************
function product_T2sT0(T2s,T0)
 implicit none
 type(Tensor2s), intent(in) :: T2s
 real(8), intent(in) :: T0
 type(Tensor2s) :: product_T2sT0
 
 product_T2sT0=product_T0T2s(T0,T2s)

end function product_T2sT0
! **********************************************************************
! *** function product_T0T2w                                         ***
! *** returns the product between one scalar value and one           ***
! *** skew-symmetric second order tensor                             ***
! ***                                                                ***
! *** input : T0, scalar value                                       ***
! ***         T2w, second order tensor                               ***
! **********************************************************************
function product_T0T2w(T0,T2w)
 implicit none
 real(8), intent(in) :: T0
 type(Tensor2w), intent(in) :: T2w
 type(Tensor2w) :: product_T0T2w
 
 product_T0T2w%val(1)=T2w%val(1)*T0
 product_T0T2w%val(2)=T2w%val(2)*T0
 product_T0T2w%val(3)=T2w%val(3)*T0

end function product_T0T2w
! **********************************************************************
! *** function product_T2wT0                                         ***
! *** returns the product between one scalar value and one           ***
! *** skew-symmetric second order tensor                             ***
! ***                                                                ***
! *** input : T2w, second order tensor                               ***
! ***         T0, scaling coefficient                                ***
! **********************************************************************
function product_T2wT0(T2w,T0)
 implicit none
 type(Tensor2w), intent(in) :: T2w
 real(8), intent(in) :: T0
 type(Tensor2w) :: product_T2wT0
 
 product_T2wT0=product_T0T2w(T0,T2w)

end function product_T2wT0
! **********************************************************************
! *** function product_T0T3                                          ***
! *** returns the product between one scalar value and one third     ***
! *** order tensor                                                   ***
! ***                                                                ***
! *** input : T0, scalar value                                       ***
! ***         T3, third order tensor                                ***
! **********************************************************************
function product_T0T3(T0,T3)
 implicit none
 real(8), intent(in) :: T0
 type(Tensor3), intent(in) :: T3
 type(Tensor3) :: product_T0T3
 integer :: i,j,k
 
 do k=1,3
  do j=1,3
   do i=1,3 
    product_T0T3%val(i,j,k)=T3%val(i,j,k)*T0
   enddo
  enddo
 enddo

end function product_T0T3
! **********************************************************************
! *** function product_T3T0                                          ***
! *** returns the product between one scalar value and one third     ***
! *** order tensor                                                   ***
! ***                                                                ***
! *** input : T3, second order tensor                                ***
! ***         T0, scaling coefficient                                ***
! **********************************************************************
function product_T3T0(T3,T0)
 implicit none
 type(Tensor3), intent(in) :: T3
 real(8), intent(in) :: T0
 type(Tensor3) :: product_T3T0
 
 product_T3T0=product_T0T3(T0,T3)

end function product_T3T0
! **********************************************************************
! *** function product_T0T4                                          ***
! *** returns the product between one scalar value and one fourth    ***
! *** order tensor                                                   ***
! ***                                                                ***
! *** input : T0, scalar value                                       ***
! ***         T4, fourth order tensor                                ***
! **********************************************************************
function product_T0T4(T0,T4)
 implicit none
 real(8), intent(in) :: T0
 type(Tensor4), intent(in) :: T4
 type(Tensor4) :: product_T0T4
 integer :: i,j,k,l
 
 do l=1,3
  do k=1,3
   do j=1,3
    do i=1,3 
     product_T0T4%val(i,j,k,l)=T4%val(i,j,k,l)*T0
    enddo
   enddo
  enddo
 enddo
 
end function product_T0T4
! **********************************************************************
! *** function product_T4T0                                          ***
! *** returns the product between one scalar value and one fourth    ***
! *** order tensor                                                   ***
! ***                                                                ***
! *** input : T4, second order tensor                                ***
! ***         T0, scaling coefficient                                ***
! **********************************************************************
function product_T4T0(T4,T0)
 implicit none
 type(Tensor4), intent(in) :: T4
 real(8), intent(in) :: T0
 type(Tensor4) :: product_T4T0
 
 product_T4T0=product_T0T4(T0,T4)

end function product_T4T0
! **********************************************************************
! *** function product_T0T4ss                                        ***
! *** returns the product between one scalar value and one symmetric ***
! *** fourth order tensor                                            ***
! ***                                                                ***
! *** input : T0, scalar value                                       ***
! ***         T4ss, fourth order tensor                              ***
! **********************************************************************
function product_T0T4ss(T0,T4ss)
 implicit none
 real(8), intent(in) :: T0
 type(Tensor4ss), intent(in) :: T4ss
 type(Tensor4ss) :: product_T0T4ss
 integer :: i,j
 
 do i=1,6
  do j=1,6
   product_T0T4ss%val(i,j)=T4ss%val(i,j)*T0
  enddo
 enddo
 
end function product_T0T4ss
! **********************************************************************
! *** function product_T4ssT0                                        ***
! *** returns the product between one scalar value and one symmetric ***
! *** fourth order tensor                                            ***
! ***                                                                ***
! *** input : T4ss, second order tensor                              ***
! ***         T0, scaling coefficient                                ***
! **********************************************************************
function product_T4ssT0(T4ss,T0)
 implicit none
 type(Tensor4ss), intent(in) :: T4ss
 real(8), intent(in) :: T0
 type(Tensor4ss) :: product_T4ssT0
 
 product_T4ssT0=product_T0T4ss(T0,T4ss)

end function product_T4ssT0
! **********************************************************************
! *** function divide_T1T0                                           ***
! *** returns the division of one first order tensor with one scalar ***
! *** value                                                          ***
! ***                                                                ***
! *** input : T0, scalar value                                       ***
! ***         T1, first order tensor                                 ***
! **********************************************************************
function divide_T1T0(T1,T0)
 implicit none
 real(8), intent(in) :: T0
 type(Tensor1), intent(in) :: T1
 type(Tensor1) :: divide_T1T0
 real(8) :: T0inv
 
 T0inv=1.d0/T0
 divide_T1T0=product_T0T1(T0inv,T1)

end function divide_T1T0
! **********************************************************************
! *** function divide_T2T0                                           ***
! *** returns the division of one second order tensor with one scalar***
! *** value                                                          ***
! ***                                                                ***
! *** input : T0, scalar value                                       ***
! ***         T2, first order tensor                                 ***
! **********************************************************************
function divide_T2T0(T2,T0)
 implicit none
 real(8), intent(in) :: T0
 type(Tensor2), intent(in) :: T2
 type(Tensor2) :: divide_T2T0
 real(8) :: T0inv
 
 T0inv=1.d0/T0
 divide_T2T0=product_T0T2(T0inv,T2)

end function divide_T2T0
! **********************************************************************
! *** function divide_T2sT0                                          ***
! *** returns the division of one symmetric second order tensor with ***
! *** one scalar value                                               ***
! ***                                                                ***
! *** input : T0, scalar value                                       ***
! ***         T2s, second order tensor                               ***
! **********************************************************************
function divide_T2sT0(T2s,T0)
 implicit none
 real(8), intent(in) :: T0
 type(Tensor2s), intent(in) :: T2s
 type(Tensor2s) :: divide_T2sT0
 real(8) :: T0inv
 
 T0inv=1.d0/T0
 divide_T2sT0=product_T0T2s(T0inv,T2s)

end function divide_T2sT0
! **********************************************************************
! *** function divide_T2wT0                                          ***
! *** returns the division of one skew-symmetric second order tensor ***
! *** with one scalar value                                          ***
! ***                                                                ***
! *** input : T0, scalar value                                       ***
! ***         T2w, second order tensor                               ***
! **********************************************************************
function divide_T2wT0(T2w,T0)
 implicit none
 real(8), intent(in) :: T0
 type(Tensor2w), intent(in) :: T2w
 type(Tensor2w) :: divide_T2wT0
 real(8) :: T0inv
 
 T0inv=1.d0/T0
 divide_T2wT0=product_T0T2w(T0inv,T2w)

end function divide_T2wT0
! **********************************************************************
! *** function divide_T3T0                                           ***
! *** returns the division of one third order tensor with one scalar ***
! *** value                                                          ***
! ***                                                                ***
! *** input : T0, scalar value                                       ***
! ***         T3, third order tensor                                 ***
! **********************************************************************
function divide_T3T0(T3,T0)
 implicit none
 real(8), intent(in) :: T0
 type(Tensor3), intent(in) :: T3
 type(Tensor3) :: divide_T3T0
 real(8) :: T0inv
 
 T0inv=1.d0/T0
 divide_T3T0=product_T0T3(T0inv,T3)

end function divide_T3T0
! **********************************************************************
! *** function divide_T4T0                                           ***
! *** returns the division of one fourth order tensor with one scalar***
! *** value                                                          ***
! ***                                                                ***
! *** input : T0, scalar value                                       ***
! ***         T4, fourth order tensor                                ***
! **********************************************************************
function divide_T4T0(T4,T0)
 implicit none
 real(8), intent(in) :: T0
 type(Tensor4), intent(in) :: T4
 type(Tensor4) :: divide_T4T0
 real(8) :: T0inv
 
 T0inv=1.d0/T0
 divide_T4T0=product_T0T4(T0inv,T4)

end function divide_T4T0
! **********************************************************************
! *** function divide_T0T4ss                                         ***
! *** returns the division of one symmetric fourth order tensor with ***
! *** one scalar value                                               ***
! ***                                                                ***
! *** input : T0, scalar value                                       ***
! ***         T4ss, fourth order tensor                              ***
! **********************************************************************
function divide_T4ssT0(T4ss,T0)
 implicit none
 real(8), intent(in) :: T0
 type(Tensor4ss), intent(in) :: T4ss
 type(Tensor4ss) :: divide_T4ssT0
 real(8) :: T0inv
 
 T0inv=1.d0/T0
 divide_T4ssT0=product_T0T4ss(T0inv,T4ss)

end function divide_T4ssT0
! **********************************************************************
! *** function dot_T1T1                                              ***
! *** returns the dot product of two first order tensors             ***
! ***                                                                ***
! *** input : T1A, first order tensor                                ***
! ***         T1B, first order tensor                                ***
! **********************************************************************
function dot_T1T1(T1A,T1B)
 implicit none
 type(Tensor1), intent(in) :: T1A,T1B
 real(8) :: dot_T1T1

 dot_T1T1=T1A%val(1)*T1B%val(1)+T1A%val(2)*T1B%val(2)+T1A%val(3)*T1B%val(3)
 
end function dot_T1T1
! **********************************************************************
! *** function dot_T2T1                                              ***
! *** returns the dot product between one second order tensor and    ***
! *** one first order tensor                                         ***
! ***                                                                ***
! *** input : T2A, second order tensor                               ***
! ***         T1A, first order tensor                                ***
! **********************************************************************
function dot_T2T1(T2A,T1A)
 implicit none
 type(Tensor2), intent(in) :: T2A
 type(Tensor1), intent(in) :: T1A
 type(Tensor1) :: dot_T2T1

 dot_T2T1%val(1)=T2A%val(1,1)*T1A%val(1)+T2A%val(1,2)*T1A%val(2)+T2A%val(1,3)*T1A%val(3)
 dot_T2T1%val(2)=T2A%val(2,1)*T1A%val(1)+T2A%val(2,2)*T1A%val(2)+T2A%val(2,3)*T1A%val(3)
 dot_T2T1%val(3)=T2A%val(3,1)*T1A%val(1)+T2A%val(3,2)*T1A%val(2)+T2A%val(3,3)*T1A%val(3)
 
end function dot_T2T1
! **********************************************************************
! *** function dot_T1T2                                              ***
! *** returns the dot product between one first order tensor and     ***
! *** one second order tensor                                        ***
! ***                                                                ***
! *** input : T2A, second order tensor                               ***
! ***         T1A, first order tensor                                ***
! **********************************************************************
function dot_T1T2(T1A,T2A)
 implicit none
 type(Tensor1), intent(in) :: T1A
 type(Tensor2), intent(in) :: T2A
 type(Tensor1) :: dot_T1T2

 dot_T1T2%val(1)=T1A%val(1)*T2A%val(1,1)+T1A%val(2)*T2A%val(2,1)+T1A%val(3)*T2A%val(3,1)
 dot_T1T2%val(2)=T1A%val(1)*T2A%val(1,2)+T1A%val(2)*T2A%val(2,2)+T1A%val(3)*T2A%val(3,2)
 dot_T1T2%val(3)=T1A%val(1)*T2A%val(1,3)+T1A%val(2)*T2A%val(2,3)+T1A%val(3)*T2A%val(3,3)
 
end function dot_T1T2
! **********************************************************************
! *** function dot_T2sT1                                             ***
! *** returns the dot product between one symmetric second order     ***
! *** tensor and one first order tensor                              ***
! ***                                                                ***
! *** input : T2sA, second order tensor                              ***
! ***         T1A, first order tensor                                ***
! **********************************************************************
function dot_T2sT1(T2sA,T1A)
 implicit none
 type(Tensor2s), intent(in) :: T2sA
 type(Tensor1), intent(in) :: T1A
 type(Tensor1) :: dot_T2sT1

 dot_T2sT1%val(1)=T2sA%val(1)*T1A%val(1)+T2sA%val(6)*T1A%val(2)+T2sA%val(5)*T1A%val(3)
 dot_T2sT1%val(2)=T2sA%val(6)*T1A%val(1)+T2sA%val(2)*T1A%val(2)+T2sA%val(4)*T1A%val(3)
 dot_T2sT1%val(3)=T2sA%val(5)*T1A%val(1)+T2sA%val(4)*T1A%val(2)+T2sA%val(3)*T1A%val(3)
 
end function dot_T2sT1
! **********************************************************************
! *** function dot_T1T2s                                             ***
! *** returns the dot product between one first order tensor and     ***
! *** one symmetric second order tensor                              ***
! ***                                                                ***
! *** input : T2sA, second order tensor                              ***
! ***         T1A, first order tensor                                ***
! **********************************************************************
function dot_T1T2s(T1A,T2sA)
 implicit none
 type(Tensor1), intent(in) :: T1A
 type(Tensor2s), intent(in) :: T2sA
 type(Tensor1) :: dot_T1T2s

 dot_T1T2s%val(1)=T1A%val(1)*T2sA%val(1)+T1A%val(2)*T2sA%val(6)+T1A%val(3)*T2sA%val(5)
 dot_T1T2s%val(2)=T1A%val(1)*T2sA%val(6)+T1A%val(2)*T2sA%val(2)+T1A%val(3)*T2sA%val(4)
 dot_T1T2s%val(3)=T1A%val(1)*T2sA%val(5)+T1A%val(2)*T2sA%val(4)+T1A%val(3)*T2sA%val(3)
 
end function dot_T1T2s
! **********************************************************************
! *** function dot_T2wT1                                             ***
! *** returns the dot product between one skew-symmetric second order***
! *** tensor and one first order tensor                              ***
! ***                                                                ***
! *** input : T2wA, second order tensor                              ***
! ***         T1A, first order tensor                                ***
! **********************************************************************
function dot_T2wT1(T2wA,T1A)
 implicit none
 type(Tensor2w), intent(in) :: T2wA
 type(Tensor1), intent(in) :: T1A
 type(Tensor1) :: dot_T2wT1

 dot_T2wT1%val(1)=T2wA%val(3)*T1A%val(2)-T2wA%val(2)*T1A%val(3)
 dot_T2wT1%val(2)=T2wA%val(1)*T1A%val(3)-T2wA%val(3)*T1A%val(1)
 dot_T2wT1%val(3)=T2wA%val(2)*T1A%val(1)-T2wA%val(1)*T1A%val(2)
 
end function dot_T2wT1
! **********************************************************************
! *** function dot_T1T2w                                             ***
! *** returns the dot product between one first order tensor and     ***
! *** one skew-symmetric second order tensor                         ***
! ***                                                                ***
! *** input : T2wA, second order tensor                              ***
! ***         T1A, first order tensor                                ***
! **********************************************************************
function dot_T1T2w(T1A,T2wA)
 implicit none
 type(Tensor1), intent(in) :: T1A
 type(Tensor2w), intent(in) :: T2wA
 type(Tensor1) :: dot_T1T2w

 dot_T1T2w%val(1)=T1A%val(3)*T2wA%val(2)-T1A%val(2)*T2wA%val(3)
 dot_T1T2w%val(2)=T1A%val(1)*T2wA%val(3)-T1A%val(3)*T2wA%val(1)
 dot_T1T2w%val(3)=T1A%val(2)*T2wA%val(1)-T1A%val(1)*T2wA%val(2)
 
end function dot_T1T2w
! **********************************************************************
! *** function dot_T2T2                                              ***
! *** returns the dot product of two second order tensors            ***
! ***                                                                ***
! *** input : T2A, second order tensor                               ***
! ***         T2B, second order tensor                               ***
! **********************************************************************
function dot_T2T2(T2A,T2B)
 implicit none
 type(Tensor2), intent(in) :: T2A,T2B
 type(Tensor2) :: dot_T2T2

 dot_T2T2%val(1,1)=T2A%val(1,1)*T2B%val(1,1)+T2A%val(1,2)*T2B%val(2,1)+T2A%val(1,3)*T2B%val(3,1)
 dot_T2T2%val(2,1)=T2A%val(2,1)*T2B%val(1,1)+T2A%val(2,2)*T2B%val(2,1)+T2A%val(2,3)*T2B%val(3,1)
 dot_T2T2%val(3,1)=T2A%val(3,1)*T2B%val(1,1)+T2A%val(3,2)*T2B%val(2,1)+T2A%val(3,3)*T2B%val(3,1)
 dot_T2T2%val(1,2)=T2A%val(1,1)*T2B%val(1,2)+T2A%val(1,2)*T2B%val(2,2)+T2A%val(1,3)*T2B%val(3,2)
 dot_T2T2%val(2,2)=T2A%val(2,1)*T2B%val(1,2)+T2A%val(2,2)*T2B%val(2,2)+T2A%val(2,3)*T2B%val(3,2)
 dot_T2T2%val(3,2)=T2A%val(3,1)*T2B%val(1,2)+T2A%val(3,2)*T2B%val(2,2)+T2A%val(3,3)*T2B%val(3,2)
 dot_T2T2%val(1,3)=T2A%val(1,1)*T2B%val(1,3)+T2A%val(1,2)*T2B%val(2,3)+T2A%val(1,3)*T2B%val(3,3)
 dot_T2T2%val(2,3)=T2A%val(2,1)*T2B%val(1,3)+T2A%val(2,2)*T2B%val(2,3)+T2A%val(2,3)*T2B%val(3,3)
 dot_T2T2%val(3,3)=T2A%val(3,1)*T2B%val(1,3)+T2A%val(3,2)*T2B%val(2,3)+T2A%val(3,3)*T2B%val(3,3)
 
end function dot_T2T2
! **********************************************************************
! *** function dot_T2T2s                                             ***
! *** returns the dot product of one general second order tensor and ***
! *** one symmetric second order tensor                              ***
! ***                                                                ***
! *** input : T2A, second order tensor                               ***
! ***         T2sA, second order tensor                              ***
! **********************************************************************
function dot_T2T2s(T2A,T2sA)
 implicit none
 type(Tensor2), intent(in) :: T2A
 type(Tensor2s), intent(in) :: T2sA
 type(Tensor2) :: dot_T2T2s

 dot_T2T2s%val(1,1)=T2A%val(1,1)*T2sA%val(1)+T2A%val(1,2)*T2sA%val(6)+T2A%val(1,3)*T2sA%val(5)
 dot_T2T2s%val(2,1)=T2A%val(2,1)*T2sA%val(1)+T2A%val(2,2)*T2sA%val(6)+T2A%val(2,3)*T2sA%val(5)
 dot_T2T2s%val(3,1)=T2A%val(3,1)*T2sA%val(1)+T2A%val(3,2)*T2sA%val(6)+T2A%val(3,3)*T2sA%val(5)
 dot_T2T2s%val(1,2)=T2A%val(1,1)*T2sA%val(6)+T2A%val(1,2)*T2sA%val(2)+T2A%val(1,3)*T2sA%val(4)
 dot_T2T2s%val(2,2)=T2A%val(2,1)*T2sA%val(6)+T2A%val(2,2)*T2sA%val(2)+T2A%val(2,3)*T2sA%val(4)
 dot_T2T2s%val(3,2)=T2A%val(3,1)*T2sA%val(6)+T2A%val(3,2)*T2sA%val(2)+T2A%val(3,3)*T2sA%val(4)
 dot_T2T2s%val(1,3)=T2A%val(1,1)*T2sA%val(5)+T2A%val(1,2)*T2sA%val(4)+T2A%val(1,3)*T2sA%val(3)
 dot_T2T2s%val(2,3)=T2A%val(2,1)*T2sA%val(5)+T2A%val(2,2)*T2sA%val(4)+T2A%val(2,3)*T2sA%val(3)
 dot_T2T2s%val(3,3)=T2A%val(3,1)*T2sA%val(5)+T2A%val(3,2)*T2sA%val(4)+T2A%val(3,3)*T2sA%val(3)
 
end function dot_T2T2s
! **********************************************************************
! *** function dot_T2sT2                                             ***
! *** returns the dot product of one symmetric order tensor and one  ***
! *** general second order tensor                                    ***
! ***                                                                ***
! *** input : T2A, second order tensor                               ***
! ***         T2sA, second order tensor                              ***
! **********************************************************************
function dot_T2sT2(T2sA,T2A)
 implicit none
 type(Tensor2s), intent(in) :: T2sA
 type(Tensor2), intent(in) :: T2A
 type(Tensor2) :: dot_T2sT2

 dot_T2sT2%val(1,1)=T2sA%val(1)*T2A%val(1,1)+T2sA%val(6)*T2A%val(2,1)+T2sA%val(5)*T2A%val(3,1)
 dot_T2sT2%val(2,1)=T2sA%val(6)*T2A%val(1,1)+T2sA%val(2)*T2A%val(2,1)+T2sA%val(4)*T2A%val(3,1)
 dot_T2sT2%val(3,1)=T2sA%val(5)*T2A%val(1,1)+T2sA%val(4)*T2A%val(2,1)+T2sA%val(3)*T2A%val(3,1)
 dot_T2sT2%val(1,2)=T2sA%val(1)*T2A%val(1,2)+T2sA%val(6)*T2A%val(2,2)+T2sA%val(5)*T2A%val(3,2)
 dot_T2sT2%val(2,2)=T2sA%val(6)*T2A%val(1,2)+T2sA%val(2)*T2A%val(2,2)+T2sA%val(4)*T2A%val(3,2)
 dot_T2sT2%val(3,2)=T2sA%val(5)*T2A%val(1,2)+T2sA%val(4)*T2A%val(2,2)+T2sA%val(3)*T2A%val(3,2)
 dot_T2sT2%val(1,3)=T2sA%val(1)*T2A%val(1,3)+T2sA%val(6)*T2A%val(2,3)+T2sA%val(5)*T2A%val(3,3)
 dot_T2sT2%val(2,3)=T2sA%val(6)*T2A%val(1,3)+T2sA%val(2)*T2A%val(2,3)+T2sA%val(4)*T2A%val(3,3)
 dot_T2sT2%val(3,3)=T2sA%val(5)*T2A%val(1,3)+T2sA%val(4)*T2A%val(2,3)+T2sA%val(3)*T2A%val(3,3)
 
end function dot_T2sT2
! **********************************************************************
! *** function dot_T2sT2s                                            ***
! *** returns the dot product of two symmetric second order tensors  ***
! ***                                                                ***
! *** input : T2sA, second order tensor                              ***
! ***         T2sB, second order tensor                              ***
! **********************************************************************
function dot_T2sT2s(T2sA,T2sB)
 implicit none
 type(Tensor2s), intent(in) :: T2sA,T2sB
 type(Tensor2) :: dot_T2sT2s

 dot_T2sT2s%val(1,1)=T2sA%val(1)*T2sB%val(1)+T2sA%val(6)*T2sB%val(6)+T2sA%val(5)*T2sB%val(5)
 dot_T2sT2s%val(2,1)=T2sA%val(6)*T2sB%val(1)+T2sA%val(2)*T2sB%val(6)+T2sA%val(4)*T2sB%val(5)
 dot_T2sT2s%val(3,1)=T2sA%val(5)*T2sB%val(1)+T2sA%val(4)*T2sB%val(6)+T2sA%val(3)*T2sB%val(5)
 dot_T2sT2s%val(1,2)=T2sA%val(1)*T2sB%val(6)+T2sA%val(6)*T2sB%val(2)+T2sA%val(5)*T2sB%val(4)
 dot_T2sT2s%val(2,2)=T2sA%val(6)*T2sB%val(6)+T2sA%val(2)*T2sB%val(2)+T2sA%val(4)*T2sB%val(4)
 dot_T2sT2s%val(3,2)=T2sA%val(5)*T2sB%val(6)+T2sA%val(4)*T2sB%val(2)+T2sA%val(3)*T2sB%val(4)
 dot_T2sT2s%val(1,3)=T2sA%val(1)*T2sB%val(5)+T2sA%val(6)*T2sB%val(4)+T2sA%val(5)*T2sB%val(3)
 dot_T2sT2s%val(2,3)=T2sA%val(6)*T2sB%val(5)+T2sA%val(2)*T2sB%val(4)+T2sA%val(4)*T2sB%val(3)
 dot_T2sT2s%val(3,3)=T2sA%val(5)*T2sB%val(5)+T2sA%val(4)*T2sB%val(4)+T2sA%val(3)*T2sB%val(3)
 
end function dot_T2sT2s
! **********************************************************************
! *** function dot_T2T2w                                             ***
! *** returns the dot product of one general second order tensor and ***
! *** one skew-symmetric second order tensor                         ***
! ***                                                                ***
! *** input : T2A, second order tensor                               ***
! ***         T2wA, second order tensor                              ***
! **********************************************************************
function dot_T2T2w(T2A,T2wA)
 implicit none
 type(Tensor2), intent(in) :: T2A
 type(Tensor2w), intent(in) :: T2wA
 type(Tensor2) :: dot_T2T2w

 dot_T2T2w%val(1,1)=-T2A%val(1,2)*T2wA%val(3)+T2A%val(1,3)*T2wA%val(2)
 dot_T2T2w%val(2,1)=-T2A%val(2,2)*T2wA%val(3)+T2A%val(2,3)*T2wA%val(2)
 dot_T2T2w%val(3,1)=-T2A%val(3,2)*T2wA%val(3)+T2A%val(3,3)*T2wA%val(2)
 dot_T2T2w%val(1,2)=T2A%val(1,1)*T2wA%val(3)-T2A%val(1,3)*T2wA%val(1)
 dot_T2T2w%val(2,2)=T2A%val(2,1)*T2wA%val(3)-T2A%val(2,3)*T2wA%val(1)
 dot_T2T2w%val(3,2)=T2A%val(3,1)*T2wA%val(3)-T2A%val(3,3)*T2wA%val(1)
 dot_T2T2w%val(1,3)=-T2A%val(1,1)*T2wA%val(2)+T2A%val(1,2)*T2wA%val(1)
 dot_T2T2w%val(2,3)=-T2A%val(2,1)*T2wA%val(2)+T2A%val(2,2)*T2wA%val(1)
 dot_T2T2w%val(3,3)=-T2A%val(3,1)*T2wA%val(2)+T2A%val(3,2)*T2wA%val(1)
 
end function dot_T2T2w
! **********************************************************************
! *** function dot_T2wT2                                             ***
! *** returns the dot product of one skew-symmetric order tensor and ***
! *** one general second order tensor                                ***
! ***                                                                ***
! *** input : T2A, second order tensor                               ***
! ***         T2wA, second order tensor                              ***
! **********************************************************************
function dot_T2wT2(T2wA,T2A)
 implicit none
 type(Tensor2w), intent(in) :: T2wA
 type(Tensor2), intent(in) :: T2A
 type(Tensor2) :: dot_T2wT2

 dot_T2wT2%val(1,1)=T2wA%val(3)*T2A%val(2,1)-T2wA%val(2)*T2A%val(3,1)
 dot_T2wT2%val(2,1)=-T2wA%val(3)*T2A%val(1,1)+T2wA%val(1)*T2A%val(3,1)
 dot_T2wT2%val(3,1)=T2wA%val(2)*T2A%val(1,1)-T2wA%val(1)*T2A%val(2,1)
 dot_T2wT2%val(1,2)=T2wA%val(3)*T2A%val(2,2)-T2wA%val(2)*T2A%val(3,2)
 dot_T2wT2%val(2,2)=-T2wA%val(3)*T2A%val(1,2)+T2wA%val(1)*T2A%val(3,2)
 dot_T2wT2%val(3,2)=T2wA%val(2)*T2A%val(1,2)-T2wA%val(1)*T2A%val(2,2)
 dot_T2wT2%val(1,3)=T2wA%val(3)*T2A%val(2,3)-T2wA%val(2)*T2A%val(3,3)
 dot_T2wT2%val(2,3)=-T2wA%val(3)*T2A%val(1,3)+T2wA%val(1)*T2A%val(3,3)
 dot_T2wT2%val(3,3)=T2wA%val(2)*T2A%val(1,3)-T2wA%val(1)*T2A%val(2,3)
 
end function dot_T2wT2
! **********************************************************************
! *** function dot_T2wT2w                                            ***
! *** returns the dot product of two skew-symmetric second order     ***
! *** tensors                                                        ***
! ***                                                                ***
! *** input : T2wA, second order tensor                              ***
! ***         T2wB, second order tensor                              ***
! **********************************************************************
function dot_T2wT2w(T2wA,T2wB)
 implicit none
 type(Tensor2w), intent(in) :: T2wA,T2wB
 type(Tensor2) :: dot_T2wT2w

 dot_T2wT2w%val(1,1)=-T2wA%val(3)*T2wB%val(3)-T2wA%val(2)*T2wB%val(2)
 dot_T2wT2w%val(2,1)=T2wA%val(1)*T2wB%val(2)
 dot_T2wT2w%val(3,1)=T2wA%val(1)*T2wB%val(3)
 dot_T2wT2w%val(1,2)=T2wA%val(2)*T2wB%val(1)
 dot_T2wT2w%val(2,2)=-T2wA%val(3)*T2wB%val(3)-T2wA%val(1)*T2wB%val(1)
 dot_T2wT2w%val(3,2)=T2wA%val(2)*T2wB%val(3)
 dot_T2wT2w%val(1,3)=T2wA%val(3)*T2wB%val(1)
 dot_T2wT2w%val(2,3)=T2wA%val(3)*T2wB%val(2)
 dot_T2wT2w%val(3,3)=-T2wA%val(2)*T2wB%val(2)-T2wA%val(1)*T2wB%val(1)
 
end function dot_T2wT2w
! **********************************************************************
! *** function dot_T2T2w                                             ***
! *** returns the dot product of one symmetric second order tensor   ***
! *** and one skew-symmetric second order tensor                     ***
! ***                                                                ***
! *** input : T2sA, second order tensor                              ***
! ***         T2wA, second order tensor                              ***
! **********************************************************************
function dot_T2sT2w(T2sA,T2wA)
 implicit none
 type(Tensor2s), intent(in) :: T2sA
 type(Tensor2w), intent(in) :: T2wA
 type(Tensor2) :: dot_T2sT2w

 dot_T2sT2w%val(1,1)=-T2sA%val(6)*T2wA%val(3)+T2sA%val(5)*T2wA%val(2)
 dot_T2sT2w%val(2,1)=-T2sA%val(2)*T2wA%val(3)+T2sA%val(4)*T2wA%val(2)
 dot_T2sT2w%val(3,1)=-T2sA%val(4)*T2wA%val(3)+T2sA%val(3)*T2wA%val(2)
 dot_T2sT2w%val(1,2)=T2sA%val(1)*T2wA%val(3)-T2sA%val(5)*T2wA%val(1)
 dot_T2sT2w%val(2,2)=T2sA%val(6)*T2wA%val(3)-T2sA%val(4)*T2wA%val(1)
 dot_T2sT2w%val(3,2)=T2sA%val(5)*T2wA%val(3)-T2sA%val(3)*T2wA%val(1)
 dot_T2sT2w%val(1,3)=-T2sA%val(1)*T2wA%val(2)+T2sA%val(6)*T2wA%val(1)
 dot_T2sT2w%val(2,3)=-T2sA%val(6)*T2wA%val(2)+T2sA%val(2)*T2wA%val(1)
 dot_T2sT2w%val(3,3)=-T2sA%val(5)*T2wA%val(2)+T2sA%val(4)*T2wA%val(1)
 
end function dot_T2sT2w
! **********************************************************************
! *** function dot_T2wT2s                                            ***
! *** returns the dot product of one skew-symmetric order tensor and ***
! *** one symmetric second order tensor                              ***
! ***                                                                ***
! *** input : T2sA, second order tensor                              ***
! ***         T2wA, second order tensor                              ***
! **********************************************************************
function dot_T2wT2s(T2wA,T2sA)
 implicit none
 type(Tensor2w), intent(in) :: T2wA
 type(Tensor2s), intent(in) :: T2sA
 type(Tensor2) :: dot_T2wT2s

 dot_T2wT2s%val(1,1)=T2wA%val(3)*T2sA%val(6)-T2wA%val(2)*T2sA%val(5)
 dot_T2wT2s%val(2,1)=-T2wA%val(3)*T2sA%val(1)+T2wA%val(1)*T2sA%val(5)
 dot_T2wT2s%val(3,1)=T2wA%val(2)*T2sA%val(1)-T2wA%val(1)*T2sA%val(6)
 dot_T2wT2s%val(1,2)=T2wA%val(3)*T2sA%val(2)-T2wA%val(2)*T2sA%val(4)
 dot_T2wT2s%val(2,2)=-T2wA%val(3)*T2sA%val(6)+T2wA%val(1)*T2sA%val(4)
 dot_T2wT2s%val(3,2)=T2wA%val(2)*T2sA%val(6)-T2wA%val(1)*T2sA%val(2)
 dot_T2wT2s%val(1,3)=T2wA%val(3)*T2sA%val(4)-T2wA%val(2)*T2sA%val(3)
 dot_T2wT2s%val(2,3)=-T2wA%val(3)*T2sA%val(5)+T2wA%val(1)*T2sA%val(3)
 dot_T2wT2s%val(3,3)=T2wA%val(2)*T2sA%val(5)-T2wA%val(1)*T2sA%val(4)
 
end function dot_T2wT2s
! **********************************************************************
! *** function dot_T2T4                                              ***
! *** returns the dot product of one second order tensor and one     ***
! *** fourth order tensor                                            ***
! ***                                                                ***
! *** input : T2A, second order tensor                               ***
! ***         T4A, fourth order tensor                               ***
! **********************************************************************
function dot_T2T4(T2A,T4A)
 implicit none
 type(Tensor2), intent(in) :: T2A
 type(Tensor4), intent(in) :: T4A
 type(Tensor4) :: dot_T2T4
 integer :: i,j,k,l,m

 do l=1,3
  do k=1,3
   do j=1,3
    do i=1,3
     dot_T2T4%val(i,j,k,l)=0.d0
     do m=1,3  
      dot_T2T4%val(i,j,k,l)=dot_T2T4%val(i,j,k,l)+T2A%val(i,m)*T4A%val(m,j,k,l)
     enddo
    enddo
   enddo
  enddo
 enddo
 
end function dot_T2T4
! **********************************************************************
! *** function dot_T4T2                                              ***
! *** returns the dot product of one fourth order tensor and one     ***
! *** second order tensor                                            ***
! ***                                                                ***
! *** input : T2A, second order tensor                               ***
! ***         T4A, fourth order tensor                               ***
! **********************************************************************
function dot_T4T2(T4A,T2A)
 implicit none
 type(Tensor4), intent(in) :: T4A
 type(Tensor2), intent(in) :: T2A
 type(Tensor4) :: dot_T4T2
 integer :: i,j,k,l,m

 do l=1,3
  do k=1,3
   do j=1,3
    do i=1,3
     dot_T4T2%val(i,j,k,l)=0.d0
     do m=1,3  
      dot_T4T2%val(i,j,k,l)=dot_T4T2%val(i,j,k,l)+T4A%val(i,j,k,m)*T2A%val(m,l)
     enddo
    enddo
   enddo
  enddo
 enddo
 
end function dot_T4T2
! **********************************************************************
! *** function dot_T2sT4ss                                           ***
! *** returns the dot product of one second order tensor and one     ***
! *** symmetric fourth order tensor                                  ***
! ***                                                                ***
! *** input : T2sA, second order tensor                              ***
! ***         T4ssA, fourth order tensor                             ***
! **********************************************************************
function dot_T2sT4ss(T2sA,T4ssA)
 implicit none
 type(Tensor2s), intent(in) :: T2sA
 type(Tensor4ss), intent(in) :: T4ssA
 type(Tensor4) :: dot_T2sT4ss
 integer :: i,j,k,l,m,im,mj,kl

 do l=1,3
  do k=1,3
   kl=indxs(k,l)
   do j=1,3
    do i=1,3
     dot_T2sT4ss%val(i,j,k,l)=0.d0
     do m=1,3  
      mj=indxs(m,j)
      im=indxs(i,m)
      dot_T2sT4ss%val(i,j,k,l)=dot_T2sT4ss%val(i,j,k,l)+T2sA%val(im)*T4ssA%val(mj,kl)
     enddo
    enddo
   enddo
  enddo
 enddo
 
end function dot_T2sT4ss
! **********************************************************************
! *** function dot_T4ssT2s                                           ***
! *** returns the dot product of one symmetric fourth order tensor   ***
! *** and one second order tensor                                    ***
! ***                                                                ***
! *** input : T2A, second order tensor                               ***
! ***         T4ssA, fourth order tensor                             ***
! **********************************************************************
function dot_T4ssT2s(T4ssA,T2sA)
 implicit none
 type(Tensor4ss), intent(in) :: T4ssA
 type(Tensor2s), intent(in) :: T2sA
 type(Tensor4) :: dot_T4ssT2s
 integer :: i,j,k,l,m,ij,km,ml

 do l=1,3
  do k=1,3
   do j=1,3
    do i=1,3
     ij=indxs(i,j)
     dot_T4ssT2s%val(i,j,k,l)=0.d0
     do m=1,3  
      km=indxs(k,m)
      ml=indxs(m,l)
      dot_T4ssT2s%val(i,j,k,l)=dot_T4ssT2s%val(i,j,k,l)+T4ssA%val(ij,km)*T2sA%val(ml)
     enddo
    enddo
   enddo
  enddo
 enddo
 
end function dot_T4ssT2s
! **********************************************************************
! *** function dot_T2T4ss                                            ***
! *** returns the dot product of one second order tensor and one     ***
! *** symmetric fourth order tensor                                  ***
! ***                                                                ***
! *** input : T2A, second order tensor                               ***
! ***         T4ssA, fourth order tensor                             ***
! **********************************************************************
function dot_T2T4ss(T2A,T4ssA)
 implicit none
 type(Tensor2), intent(in) :: T2A
 type(Tensor4ss), intent(in) :: T4ssA
 type(Tensor4) :: dot_T2T4ss
 integer :: i,j,k,l,m,mj,kl

 do l=1,3
  do k=1,3
   kl=indxs(k,l)
   do j=1,3
    do i=1,3
     dot_T2T4ss%val(i,j,k,l)=0.d0
     do m=1,3  
      mj=indxs(m,j)
      dot_T2T4ss%val(i,j,k,l)=dot_T2T4ss%val(i,j,k,l)+T2A%val(i,m)*T4ssA%val(mj,kl)
     enddo
    enddo
   enddo
  enddo
 enddo
 
end function dot_T2T4ss
! **********************************************************************
! *** function dot_T4ssT2                                            ***
! *** returns the dot product of one symmetric fourth order tensor   ***
! *** and one second order tensor                                    ***
! ***                                                                ***
! *** input : T2A, second order tensor                               ***
! ***         T4ssA, fourth order tensor                             ***
! **********************************************************************
function dot_T4ssT2(T4ssA,T2A)
 implicit none
 type(Tensor4ss), intent(in) :: T4ssA
 type(Tensor2), intent(in) :: T2A
 type(Tensor4) :: dot_T4ssT2
 integer :: i,j,k,l,m,ij,km

 do l=1,3
  do k=1,3
   do j=1,3
    do i=1,3
     ij=indxs(i,j)
     dot_T4ssT2%val(i,j,k,l)=0.d0
     do m=1,3  
      km=indxs(k,m)
      dot_T4ssT2%val(i,j,k,l)=dot_T4ssT2%val(i,j,k,l)+T4ssA%val(ij,km)*T2A%val(m,l)
     enddo
    enddo
   enddo
  enddo
 enddo
 
end function dot_T4ssT2
! **********************************************************************
! *** function ddot_T2T2                                             ***
! *** returns the double dot product of two second order tensors     ***
! ***                                                                ***
! *** input : T2A, second order tensor                               ***
! ***         T2B, second order tensor                               ***
! **********************************************************************
function ddot_T2T2(T2A,T2B)
 implicit none
 type(Tensor2), intent(in) :: T2A,T2B
 real(8) :: ddot_T2T2

 ddot_T2T2=T2A%val(1,1)*T2B%val(1,1)+T2A%val(2,1)*T2B%val(2,1)+T2A%val(3,1)*T2B%val(3,1) &
          +T2A%val(1,2)*T2B%val(1,2)+T2A%val(2,2)*T2B%val(2,2)+T2A%val(3,2)*T2B%val(3,2) &
          +T2A%val(1,3)*T2B%val(1,3)+T2A%val(2,3)*T2B%val(2,3)+T2A%val(3,3)*T2B%val(3,3)
 
end function ddot_T2T2
! **********************************************************************
! *** function ddot_T2sT2                                            ***
! *** returns the double dot product of one symmetric and one general***
! *** second order tensors                                           ***
! ***                                                                ***
! *** input : T2sA, second order tensor                              ***
! ***         T2A, second order tensor                               ***
! **********************************************************************
function ddot_T2sT2(T2sA,T2A)
 implicit none
 type(Tensor2s), intent(in) :: T2sA
 type(Tensor2), intent(in) :: T2A
 real(8) :: ddot_T2sT2

 ddot_T2sT2=T2sA%val(1)*T2A%val(1,1)+T2sA%val(6)*T2A%val(2,1)+T2sA%val(5)*T2A%val(3,1) &
           +T2sA%val(6)*T2A%val(1,2)+T2sA%val(2)*T2A%val(2,2)+T2sA%val(4)*T2A%val(3,2) &
           +T2sA%val(5)*T2A%val(1,3)+T2sA%val(4)*T2A%val(2,3)+T2sA%val(3)*T2A%val(3,3)
 
end function ddot_T2sT2
! **********************************************************************
! *** function ddot_T2T2s                                            ***
! *** returns the double dot product of one general and one symmetric***
! *** second order tensors                                           ***
! ***                                                                ***
! *** input : T2A, second order tensor                              ***
! ***         T2sA, second order tensor                               ***
! **********************************************************************
function ddot_T2T2s(T2A,T2sA)
 implicit none
 type(Tensor2), intent(in) :: T2A
 type(Tensor2s), intent(in) :: T2sA
 real(8) :: ddot_T2T2s

 ddot_T2T2s=T2A%val(1,1)*T2sA%val(1)+T2A%val(2,1)*T2sA%val(6)+T2A%val(3,1)*T2sA%val(5) &
           +T2A%val(1,2)*T2sA%val(6)+T2A%val(2,2)*T2sA%val(2)+T2A%val(3,2)*T2sA%val(4) &
           +T2A%val(1,3)*T2sA%val(5)+T2A%val(2,3)*T2sA%val(4)+T2A%val(3,3)*T2sA%val(3)
          
end function ddot_T2T2s
! **********************************************************************
! *** function ddot_T2wT2                                            ***
! *** returns the double dot product of one skew-symmetric and one   ***
! *** general second order tensors                                   ***
! ***                                                                ***
! *** input : T2wA, second order tensor                              ***
! ***         T2A, second order tensor                               ***
! **********************************************************************
function ddot_T2wT2(T2wA,T2A)
 implicit none
 type(Tensor2w), intent(in) :: T2wA
 type(Tensor2), intent(in) :: T2A
 real(8) :: ddot_T2wT2

 ddot_T2wT2=-T2wA%val(3)*T2A%val(2,1)+T2wA%val(2)*T2A%val(3,1) &
            +T2wA%val(3)*T2A%val(1,2)-T2wA%val(1)*T2A%val(3,2) &
            -T2wA%val(2)*T2A%val(1,3)+T2wA%val(1)*T2A%val(2,3)
 
end function ddot_T2wT2
! **********************************************************************
! *** function ddot_T2T2w                                            ***
! *** returns the double dot product of one general and one          ***
! *** skew-symmetric second order tensors                            ***
! ***                                                                ***
! *** input : T2A, second order tensor                               ***
! ***         T2wA, second order tensor                              ***
! **********************************************************************
function ddot_T2T2w(T2A,T2wA)
 implicit none
 type(Tensor2), intent(in) :: T2A
 type(Tensor2w), intent(in) :: T2wA
 real(8) :: ddot_T2T2w

 ddot_T2T2w=-T2A%val(2,1)*T2wA%val(3)+T2A%val(3,1)*T2wA%val(2) &
            +T2A%val(1,2)*T2wA%val(3)-T2A%val(3,2)*T2wA%val(1) &
            -T2A%val(1,3)*T2wA%val(2)+T2A%val(2,3)*T2wA%val(1)
          
end function ddot_T2T2w
! **********************************************************************
! *** function ddot_T2sT2s                                           ***
! *** returns the double dot product of two symmetric second order   ***
! *** tensors                                                        ***
! ***                                                                ***
! *** input : T2sA, second order tensor                              ***
! ***         T2sB, second order tensor                              ***
! **********************************************************************
function ddot_T2sT2s(T2sA,T2sB)
 implicit none
 type(Tensor2s), intent(in) :: T2sA,T2sB
 real(8) :: ddot_T2sT2s

 ddot_T2sT2s=T2sA%val(1)*T2sB%val(1)+T2sA%val(2)*T2sB%val(2)+T2sA%val(3)*T2sB%val(3) &
           +(T2sA%val(4)*T2sB%val(4)+T2sA%val(5)*T2sB%val(5)+T2sA%val(6)*T2sB%val(6))*2.d0
 
end function ddot_T2sT2s
! **********************************************************************
! *** function ddot_T2wT2w                                           ***
! *** returns the double dot product of two skew-symmetric second    ***
! *** order tensors                                                  ***
! ***                                                                ***
! *** input : T2wA, second order tensor                              ***
! ***         T2wB, second order tensor                              ***
! **********************************************************************
function ddot_T2wT2w(T2wA,T2wB)
 implicit none
 type(Tensor2w), intent(in) :: T2wA,T2wB
 real(8) :: ddot_T2wT2w

 ddot_T2wT2w=(T2wA%val(1)*T2wB%val(1)+T2wA%val(2)*T2wB%val(2)+T2wA%val(3)*T2wB%val(3))*2.d0
 
end function ddot_T2wT2w
! **********************************************************************
! *** function ddot_T2wT2s                                           ***
! *** returns the double dot product of one skew-symmetric and one   ***
! *** symmetric second order tensors                                 ***
! ***                                                                ***
! *** input : T2wA, second order tensor                              ***
! ***         T2sA, second order tensor                              ***
! **********************************************************************
function ddot_T2wT2s(T2wA,T2sA)
 implicit none
 type(Tensor2w), intent(in) :: T2wA
 type(Tensor2s), intent(in) :: T2sA
 real(8) :: ddot_T2wT2s

 ddot_T2wT2s=-T2wA%val(3)*T2sA%val(6)+T2wA%val(2)*T2sA%val(5) &
             +T2wA%val(3)*T2sA%val(6)-T2wA%val(1)*T2sA%val(4) &
             -T2wA%val(2)*T2sA%val(5)+T2wA%val(1)*T2sA%val(4)
 
end function ddot_T2wT2s
! **********************************************************************
! *** function ddot_T2sT2w                                           ***
! *** returns the double dot product of one symmetric and one        ***
! *** skew-symmetric second order tensors                            ***
! ***                                                                ***
! *** input : T2sA, second order tensor                              ***
! ***         T2wA, second order tensor                              ***
! **********************************************************************
function ddot_T2sT2w(T2sA,T2wA)
 implicit none
 type(Tensor2s), intent(in) :: T2sA
 type(Tensor2w), intent(in) :: T2wA
 real(8) :: ddot_T2sT2w

 ddot_T2sT2w=-T2sA%val(6)*T2wA%val(3)+T2sA%val(5)*T2wA%val(2) &
             +T2sA%val(6)*T2wA%val(3)-T2sA%val(4)*T2wA%val(1) &
             -T2sA%val(5)*T2wA%val(2)+T2sA%val(4)*T2wA%val(1)
          
end function ddot_T2sT2w
! **********************************************************************
! *** function ddot_T4T2                                             ***
! *** returns the double dot product of one general fourth rank      ***
! *** tensor and one general second order tensor                     ***
! ***                                                                ***
! *** input : T4A, fourth order tensor                               ***
! ***         T2A, second order tensor                               ***
! **********************************************************************
function ddot_T4T2(T4A,T2A)
 implicit none
 type(Tensor4), intent(in) :: T4A
 type(Tensor2), intent(in) :: T2A
 type(Tensor2) :: ddot_T4T2

 ddot_T4T2%val(1,1)=T4A%val(1,1,1,1)*T2A%val(1,1)+T4A%val(1,1,2,1)*T2A%val(2,1)+T4A%val(1,1,3,1)*T2A%val(3,1) &
                   +T4A%val(1,1,1,2)*T2A%val(1,2)+T4A%val(1,1,2,2)*T2A%val(2,2)+T4A%val(1,1,3,2)*T2A%val(3,2) &
                   +T4A%val(1,1,1,3)*T2A%val(1,3)+T4A%val(1,1,2,3)*T2A%val(2,3)+T4A%val(1,1,3,3)*T2A%val(3,3)
 ddot_T4T2%val(2,1)=T4A%val(2,1,1,1)*T2A%val(1,1)+T4A%val(2,1,2,1)*T2A%val(2,1)+T4A%val(2,1,3,1)*T2A%val(3,1) &
                   +T4A%val(2,1,1,2)*T2A%val(1,2)+T4A%val(2,1,2,2)*T2A%val(2,2)+T4A%val(2,1,3,2)*T2A%val(3,2) &
                   +T4A%val(2,1,1,3)*T2A%val(1,3)+T4A%val(2,1,2,3)*T2A%val(2,3)+T4A%val(2,1,3,3)*T2A%val(3,3)
 ddot_T4T2%val(3,1)=T4A%val(3,1,1,1)*T2A%val(1,1)+T4A%val(3,1,2,1)*T2A%val(2,1)+T4A%val(3,1,3,1)*T2A%val(3,1) &
                   +T4A%val(3,1,1,2)*T2A%val(1,2)+T4A%val(3,1,2,2)*T2A%val(2,2)+T4A%val(3,1,3,2)*T2A%val(3,2) &
                   +T4A%val(3,1,1,3)*T2A%val(1,3)+T4A%val(3,1,2,3)*T2A%val(2,3)+T4A%val(3,1,3,3)*T2A%val(3,3)
 ddot_T4T2%val(1,2)=T4A%val(1,2,1,1)*T2A%val(1,1)+T4A%val(1,2,2,1)*T2A%val(2,1)+T4A%val(1,2,3,1)*T2A%val(3,1) &
                   +T4A%val(1,2,1,2)*T2A%val(1,2)+T4A%val(1,2,2,2)*T2A%val(2,2)+T4A%val(1,2,3,2)*T2A%val(3,2) &
                   +T4A%val(1,2,1,3)*T2A%val(1,3)+T4A%val(1,2,2,3)*T2A%val(2,3)+T4A%val(1,2,3,3)*T2A%val(3,3)
 ddot_T4T2%val(2,2)=T4A%val(2,2,1,1)*T2A%val(1,1)+T4A%val(2,2,2,1)*T2A%val(2,1)+T4A%val(2,2,3,1)*T2A%val(3,1) &
                   +T4A%val(2,2,1,2)*T2A%val(1,2)+T4A%val(2,2,2,2)*T2A%val(2,2)+T4A%val(2,2,3,2)*T2A%val(3,2) &
                   +T4A%val(2,2,1,3)*T2A%val(1,3)+T4A%val(2,2,2,3)*T2A%val(2,3)+T4A%val(2,2,3,3)*T2A%val(3,3)
 ddot_T4T2%val(3,2)=T4A%val(3,2,1,1)*T2A%val(1,1)+T4A%val(3,2,2,1)*T2A%val(2,1)+T4A%val(3,2,3,1)*T2A%val(3,1) &
                   +T4A%val(3,2,1,2)*T2A%val(1,2)+T4A%val(3,2,2,2)*T2A%val(2,2)+T4A%val(3,2,3,2)*T2A%val(3,2) &
                   +T4A%val(3,2,1,3)*T2A%val(1,3)+T4A%val(3,2,2,3)*T2A%val(2,3)+T4A%val(3,2,3,3)*T2A%val(3,3)
 ddot_T4T2%val(1,3)=T4A%val(1,3,1,1)*T2A%val(1,1)+T4A%val(1,3,2,1)*T2A%val(2,1)+T4A%val(1,3,3,1)*T2A%val(3,1) &
                   +T4A%val(1,3,1,2)*T2A%val(1,2)+T4A%val(1,3,2,2)*T2A%val(2,2)+T4A%val(1,3,3,2)*T2A%val(3,2) &
                   +T4A%val(1,3,1,3)*T2A%val(1,3)+T4A%val(1,3,2,3)*T2A%val(2,3)+T4A%val(1,3,3,3)*T2A%val(3,3)
 ddot_T4T2%val(2,3)=T4A%val(2,3,1,1)*T2A%val(1,1)+T4A%val(2,3,2,1)*T2A%val(2,1)+T4A%val(2,3,3,1)*T2A%val(3,1) &
                   +T4A%val(2,3,1,2)*T2A%val(1,2)+T4A%val(2,3,2,2)*T2A%val(2,2)+T4A%val(2,3,3,2)*T2A%val(3,2) &
                   +T4A%val(2,3,1,3)*T2A%val(1,3)+T4A%val(2,3,2,3)*T2A%val(2,3)+T4A%val(2,3,3,3)*T2A%val(3,3)
 ddot_T4T2%val(3,3)=T4A%val(3,3,1,1)*T2A%val(1,1)+T4A%val(3,3,2,1)*T2A%val(2,1)+T4A%val(3,3,3,1)*T2A%val(3,1) &
                   +T4A%val(3,3,1,2)*T2A%val(1,2)+T4A%val(3,3,2,2)*T2A%val(2,2)+T4A%val(3,3,3,2)*T2A%val(3,2) &
                   +T4A%val(3,3,1,3)*T2A%val(1,3)+T4A%val(3,3,2,3)*T2A%val(2,3)+T4A%val(3,3,3,3)*T2A%val(3,3)
                   
end function ddot_T4T2
! **********************************************************************
! *** function ddot_T2T4                                             ***
! *** returns the double dot product of one general second rank      ***
! *** tensor and one general fourth order tensor                     ***
! ***                                                                ***
! *** input : T4A, fourth order tensor                               ***
! ***         T2A, second order tensor                               ***
! **********************************************************************
function ddot_T2T4(T2A,T4A)
 implicit none
 type(Tensor2), intent(in) :: T2A
 type(Tensor4), intent(in) :: T4A
 type(Tensor2) :: ddot_T2T4

 ddot_T2T4%val(1,1)=T2A%val(1,1)*T4A%val(1,1,1,1)+T2A%val(2,1)*T4A%val(2,1,1,1)+T2A%val(3,1)*T4A%val(3,1,1,1) &
                   +T2A%val(1,2)*T4A%val(1,2,1,1)+T2A%val(2,2)*T4A%val(2,2,1,1)+T2A%val(3,2)*T4A%val(3,2,1,1) &
                   +T2A%val(1,3)*T4A%val(1,3,1,1)+T2A%val(2,3)*T4A%val(2,3,1,1)+T2A%val(3,3)*T4A%val(3,3,1,1)
 ddot_T2T4%val(2,1)=T2A%val(1,1)*T4A%val(1,1,2,1)+T2A%val(2,1)*T4A%val(2,1,2,1)+T2A%val(3,1)*T4A%val(3,1,2,1) &
                   +T2A%val(1,2)*T4A%val(1,2,2,1)+T2A%val(2,2)*T4A%val(2,2,2,1)+T2A%val(3,2)*T4A%val(3,2,2,1) &
                   +T2A%val(1,3)*T4A%val(1,3,2,1)+T2A%val(2,3)*T4A%val(2,3,2,1)+T2A%val(3,3)*T4A%val(3,3,2,1)
 ddot_T2T4%val(3,1)=T2A%val(1,1)*T4A%val(1,1,3,1)+T2A%val(2,1)*T4A%val(2,1,3,1)+T2A%val(3,1)*T4A%val(3,1,3,1) &
                   +T2A%val(1,2)*T4A%val(1,2,3,1)+T2A%val(2,2)*T4A%val(2,2,3,1)+T2A%val(3,2)*T4A%val(3,2,3,1) &
                   +T2A%val(1,3)*T4A%val(1,3,3,1)+T2A%val(2,3)*T4A%val(2,3,3,1)+T2A%val(3,3)*T4A%val(3,3,3,1)
 ddot_T2T4%val(1,2)=T2A%val(1,1)*T4A%val(1,1,1,2)+T2A%val(2,1)*T4A%val(2,1,1,2)+T2A%val(3,1)*T4A%val(3,1,1,2) &
                   +T2A%val(1,2)*T4A%val(1,2,1,2)+T2A%val(2,2)*T4A%val(2,2,1,2)+T2A%val(3,2)*T4A%val(3,2,1,2) &
                   +T2A%val(1,3)*T4A%val(1,3,1,2)+T2A%val(2,3)*T4A%val(2,3,1,2)+T2A%val(3,3)*T4A%val(3,3,1,2)
 ddot_T2T4%val(2,2)=T2A%val(1,1)*T4A%val(1,1,2,2)+T2A%val(2,1)*T4A%val(2,1,2,2)+T2A%val(3,1)*T4A%val(3,1,2,2) &
                   +T2A%val(1,2)*T4A%val(1,2,2,2)+T2A%val(2,2)*T4A%val(2,2,2,2)+T2A%val(3,2)*T4A%val(3,2,2,2) &
                   +T2A%val(1,3)*T4A%val(1,3,2,2)+T2A%val(2,3)*T4A%val(2,3,2,2)+T2A%val(3,3)*T4A%val(3,3,2,2)
 ddot_T2T4%val(3,2)=T2A%val(1,1)*T4A%val(1,1,3,2)+T2A%val(2,1)*T4A%val(2,1,3,2)+T2A%val(3,1)*T4A%val(3,1,3,2) &
                   +T2A%val(1,2)*T4A%val(1,2,3,2)+T2A%val(2,2)*T4A%val(2,2,3,2)+T2A%val(3,2)*T4A%val(3,2,3,2) &
                   +T2A%val(1,3)*T4A%val(1,3,3,2)+T2A%val(2,3)*T4A%val(2,3,3,2)+T2A%val(3,3)*T4A%val(3,3,3,2)
 ddot_T2T4%val(1,3)=T2A%val(1,1)*T4A%val(1,1,1,3)+T2A%val(2,1)*T4A%val(2,1,1,3)+T2A%val(3,1)*T4A%val(3,1,1,3) &
                   +T2A%val(1,2)*T4A%val(1,2,1,3)+T2A%val(2,2)*T4A%val(2,2,1,3)+T2A%val(3,2)*T4A%val(3,2,1,3) &
                   +T2A%val(1,3)*T4A%val(1,3,1,3)+T2A%val(2,3)*T4A%val(2,3,1,3)+T2A%val(3,3)*T4A%val(3,3,1,3)
 ddot_T2T4%val(2,3)=T2A%val(1,1)*T4A%val(1,1,2,3)+T2A%val(2,1)*T4A%val(2,1,2,3)+T2A%val(3,1)*T4A%val(3,1,2,3) &
                   +T2A%val(1,2)*T4A%val(1,2,2,3)+T2A%val(2,2)*T4A%val(2,2,2,3)+T2A%val(3,2)*T4A%val(3,2,2,3) &
                   +T2A%val(1,3)*T4A%val(1,3,2,3)+T2A%val(2,3)*T4A%val(2,3,2,3)+T2A%val(3,3)*T4A%val(3,3,2,3)
 ddot_T2T4%val(3,3)=T2A%val(1,1)*T4A%val(1,1,3,3)+T2A%val(2,1)*T4A%val(2,1,3,3)+T2A%val(3,1)*T4A%val(3,1,3,3) &
                   +T2A%val(1,2)*T4A%val(1,2,3,3)+T2A%val(2,2)*T4A%val(2,2,3,3)+T2A%val(3,2)*T4A%val(3,2,3,3) &
                   +T2A%val(1,3)*T4A%val(1,3,3,3)+T2A%val(2,3)*T4A%val(2,3,3,3)+T2A%val(3,3)*T4A%val(3,3,3,3)
                   
end function ddot_T2T4
! **********************************************************************
! *** function ddot_T4ssT2s                                          ***
! *** returns the double dot product of one symmetric fourth rank    ***
! *** tensor and one symmetric second order tensor                   ***
! ***                                                                ***
! *** input : T4ssA, fourth order tensor                             ***
! ***         T2sA, second order tensor                              ***
! **********************************************************************
function ddot_T4ssT2s(T4ssA,T2sA)
 implicit none
 type(Tensor4ss), intent(in) :: T4ssA
 type(Tensor2s), intent(in) :: T2sA
 type(Tensor2s) :: ddot_T4ssT2s

 ddot_T4ssT2s%val(1)=T4ssA%val(1,1)*T2sA%val(1)+T4ssA%val(1,2)*T2sA%val(2)+T4ssA%val(1,3)*T2sA%val(3) &
                   +(T4ssA%val(1,4)*T2sA%val(4)+T4ssA%val(1,5)*T2sA%val(5)+T4ssA%val(1,6)*T2sA%val(6))*2.d0
                   
 ddot_T4ssT2s%val(2)=T4ssA%val(2,1)*T2sA%val(1)+T4ssA%val(2,2)*T2sA%val(2)+T4ssA%val(2,3)*T2sA%val(3) &
                   +(T4ssA%val(2,4)*T2sA%val(4)+T4ssA%val(2,5)*T2sA%val(5)+T4ssA%val(2,6)*T2sA%val(6))*2.d0
                   
 ddot_T4ssT2s%val(3)=T4ssA%val(3,1)*T2sA%val(1)+T4ssA%val(3,2)*T2sA%val(2)+T4ssA%val(3,3)*T2sA%val(3) &
                   +(T4ssA%val(3,4)*T2sA%val(4)+T4ssA%val(3,5)*T2sA%val(5)+T4ssA%val(3,6)*T2sA%val(6))*2.d0
                   
 ddot_T4ssT2s%val(4)=T4ssA%val(4,1)*T2sA%val(1)+T4ssA%val(4,2)*T2sA%val(2)+T4ssA%val(4,3)*T2sA%val(3) &
                   +(T4ssA%val(4,4)*T2sA%val(4)+T4ssA%val(4,5)*T2sA%val(5)+T4ssA%val(4,6)*T2sA%val(6))*2.d0
                   
 ddot_T4ssT2s%val(5)=T4ssA%val(5,1)*T2sA%val(1)+T4ssA%val(5,2)*T2sA%val(2)+T4ssA%val(5,3)*T2sA%val(3) &
                   +(T4ssA%val(5,4)*T2sA%val(4)+T4ssA%val(5,5)*T2sA%val(5)+T4ssA%val(5,6)*T2sA%val(6))*2.d0
                   
 ddot_T4ssT2s%val(6)=T4ssA%val(6,1)*T2sA%val(1)+T4ssA%val(6,2)*T2sA%val(2)+T4ssA%val(6,3)*T2sA%val(3) &
                   +(T4ssA%val(6,4)*T2sA%val(4)+T4ssA%val(6,5)*T2sA%val(5)+T4ssA%val(6,6)*T2sA%val(6))*2.d0
 
end function ddot_T4ssT2s
! **********************************************************************
! *** function ddot_T2sT4ss                                          ***
! *** returns the double dot product of one symmetric second rank    ***
! *** tensor and one symmetric fourth order tensor                   ***
! ***                                                                ***
! *** input : T4ssA, fourth order tensor                             ***
! ***         T2sA, second order tensor                              ***
! **********************************************************************
function ddot_T2sT4ss(T2sA,T4ssA)
 implicit none
 type(Tensor2s), intent(in) :: T2sA
 type(Tensor4ss), intent(in) :: T4ssA
 type(Tensor2s) :: ddot_T2sT4ss
 
 ddot_T2sT4ss%val(1)=T2sA%val(1)*T4ssA%val(1,1)+T2sA%val(2)*T4ssA%val(2,1)+T2sA%val(3)*T4ssA%val(3,1) &
                   +(T2sA%val(4)*T4ssA%val(4,1)+T2sA%val(5)*T4ssA%val(5,1)+T2sA%val(6)*T4ssA%val(6,1))*2.d0

 ddot_T2sT4ss%val(2)=T2sA%val(1)*T4ssA%val(1,2)+T2sA%val(2)*T4ssA%val(2,2)+T2sA%val(3)*T4ssA%val(3,2) &
                   +(T2sA%val(4)*T4ssA%val(4,2)+T2sA%val(5)*T4ssA%val(5,2)+T2sA%val(6)*T4ssA%val(6,2))*2.d0
                   
 ddot_T2sT4ss%val(3)=T2sA%val(1)*T4ssA%val(1,3)+T2sA%val(2)*T4ssA%val(2,3)+T2sA%val(3)*T4ssA%val(3,3) &
                   +(T2sA%val(4)*T4ssA%val(4,3)+T2sA%val(5)*T4ssA%val(5,3)+T2sA%val(6)*T4ssA%val(6,3))*2.d0
                   
 ddot_T2sT4ss%val(4)=T2sA%val(1)*T4ssA%val(1,4)+T2sA%val(2)*T4ssA%val(2,4)+T2sA%val(3)*T4ssA%val(3,4) &
                   +(T2sA%val(4)*T4ssA%val(4,4)+T2sA%val(5)*T4ssA%val(5,4)+T2sA%val(6)*T4ssA%val(6,4))*2.d0
                   
 ddot_T2sT4ss%val(5)=T2sA%val(1)*T4ssA%val(1,5)+T2sA%val(2)*T4ssA%val(2,5)+T2sA%val(3)*T4ssA%val(3,5) &
                   +(T2sA%val(4)*T4ssA%val(4,5)+T2sA%val(5)*T4ssA%val(5,5)+T2sA%val(6)*T4ssA%val(6,5))*2.d0
                   
 ddot_T2sT4ss%val(6)=T2sA%val(1)*T4ssA%val(1,6)+T2sA%val(2)*T4ssA%val(2,6)+T2sA%val(3)*T4ssA%val(3,6) &
                   +(T2sA%val(4)*T4ssA%val(4,6)+T2sA%val(5)*T4ssA%val(5,6)+T2sA%val(6)*T4ssA%val(6,6))*2.d0
                   
end function ddot_T2sT4ss
! **********************************************************************
! *** function ddot_T4T2s                                            ***
! *** returns the double dot product of one general fourth rank      ***
! *** tensor and one symmetric second order tensor                   ***
! ***                                                                ***
! *** input : T4A, fourth order tensor                               ***
! ***         T2sA, second order tensor                              ***
! **********************************************************************
function ddot_T4T2s(T4A,T2sA)
 implicit none
 type(Tensor4), intent(in) :: T4A
 type(Tensor2s), intent(in) :: T2sA
 type(Tensor2) :: ddot_T4T2s

 ddot_T4T2s%val(1,1)=T4A%val(1,1,1,1)*T2sA%val(1)+T4A%val(1,1,2,1)*T2sA%val(6)+T4A%val(1,1,3,1)*T2sA%val(5) &
                    +T4A%val(1,1,1,2)*T2sA%val(6)+T4A%val(1,1,2,2)*T2sA%val(2)+T4A%val(1,1,3,2)*T2sA%val(4) &
                    +T4A%val(1,1,1,3)*T2sA%val(5)+T4A%val(1,1,2,3)*T2sA%val(4)+T4A%val(1,1,3,3)*T2sA%val(3)
 ddot_T4T2s%val(2,1)=T4A%val(2,1,1,1)*T2sA%val(1)+T4A%val(2,1,2,1)*T2sA%val(6)+T4A%val(2,1,3,1)*T2sA%val(5) &
                    +T4A%val(2,1,1,2)*T2sA%val(6)+T4A%val(2,1,2,2)*T2sA%val(2)+T4A%val(2,1,3,2)*T2sA%val(4) &
                    +T4A%val(2,1,1,3)*T2sA%val(5)+T4A%val(2,1,2,3)*T2sA%val(4)+T4A%val(2,1,3,3)*T2sA%val(3)
 ddot_T4T2s%val(3,1)=T4A%val(3,1,1,1)*T2sA%val(1)+T4A%val(3,1,2,1)*T2sA%val(6)+T4A%val(3,1,3,1)*T2sA%val(5) &
                    +T4A%val(3,1,1,2)*T2sA%val(6)+T4A%val(3,1,2,2)*T2sA%val(2)+T4A%val(3,1,3,2)*T2sA%val(4) &
                    +T4A%val(3,1,1,3)*T2sA%val(5)+T4A%val(3,1,2,3)*T2sA%val(4)+T4A%val(3,1,3,3)*T2sA%val(3)
 ddot_T4T2s%val(1,2)=T4A%val(1,2,1,1)*T2sA%val(1)+T4A%val(1,2,2,1)*T2sA%val(6)+T4A%val(1,2,3,1)*T2sA%val(5) &
                    +T4A%val(1,2,1,2)*T2sA%val(6)+T4A%val(1,2,2,2)*T2sA%val(2)+T4A%val(1,2,3,2)*T2sA%val(4) &
                    +T4A%val(1,2,1,3)*T2sA%val(5)+T4A%val(1,2,2,3)*T2sA%val(4)+T4A%val(1,2,3,3)*T2sA%val(3)
 ddot_T4T2s%val(2,2)=T4A%val(2,2,1,1)*T2sA%val(1)+T4A%val(2,2,2,1)*T2sA%val(6)+T4A%val(2,2,3,1)*T2sA%val(5) &
                    +T4A%val(2,2,1,2)*T2sA%val(6)+T4A%val(2,2,2,2)*T2sA%val(2)+T4A%val(2,2,3,2)*T2sA%val(4) &
                    +T4A%val(2,2,1,3)*T2sA%val(5)+T4A%val(2,2,2,3)*T2sA%val(4)+T4A%val(2,2,3,3)*T2sA%val(3)
 ddot_T4T2s%val(3,2)=T4A%val(3,2,1,1)*T2sA%val(1)+T4A%val(3,2,2,1)*T2sA%val(6)+T4A%val(3,2,3,1)*T2sA%val(5) &
                    +T4A%val(3,2,1,2)*T2sA%val(6)+T4A%val(3,2,2,2)*T2sA%val(2)+T4A%val(3,2,3,2)*T2sA%val(4) &
                    +T4A%val(3,2,1,3)*T2sA%val(5)+T4A%val(3,2,2,3)*T2sA%val(4)+T4A%val(3,2,3,3)*T2sA%val(3)
 ddot_T4T2s%val(1,3)=T4A%val(1,3,1,1)*T2sA%val(1)+T4A%val(1,3,2,1)*T2sA%val(6)+T4A%val(1,3,3,1)*T2sA%val(5) &
                    +T4A%val(1,3,1,2)*T2sA%val(6)+T4A%val(1,3,2,2)*T2sA%val(2)+T4A%val(1,3,3,2)*T2sA%val(4) &
                    +T4A%val(1,3,1,3)*T2sA%val(5)+T4A%val(1,3,2,3)*T2sA%val(4)+T4A%val(1,3,3,3)*T2sA%val(3)
 ddot_T4T2s%val(2,3)=T4A%val(2,3,1,1)*T2sA%val(1)+T4A%val(2,3,2,1)*T2sA%val(6)+T4A%val(2,3,3,1)*T2sA%val(5) &
                    +T4A%val(2,3,1,2)*T2sA%val(6)+T4A%val(2,3,2,2)*T2sA%val(2)+T4A%val(2,3,3,2)*T2sA%val(4) &
                    +T4A%val(2,3,1,3)*T2sA%val(5)+T4A%val(2,3,2,3)*T2sA%val(4)+T4A%val(2,3,3,3)*T2sA%val(3)
 ddot_T4T2s%val(3,3)=T4A%val(3,3,1,1)*T2sA%val(1)+T4A%val(3,3,2,1)*T2sA%val(6)+T4A%val(3,3,3,1)*T2sA%val(5) &
                    +T4A%val(3,3,1,2)*T2sA%val(6)+T4A%val(3,3,2,2)*T2sA%val(2)+T4A%val(3,3,3,2)*T2sA%val(4) &
                    +T4A%val(3,3,1,3)*T2sA%val(5)+T4A%val(3,3,2,3)*T2sA%val(4)+T4A%val(3,3,3,3)*T2sA%val(3)
                   
end function ddot_T4T2s
! **********************************************************************
! *** function ddot_T2sT4                                            ***
! *** returns the double dot product of one symmetric second rank    ***
! *** tensor and one general fourth order tensor                     ***
! ***                                                                ***
! *** input : T4A, fourth order tensor                               ***
! ***         T2sA, second order tensor                              ***
! **********************************************************************
function ddot_T2sT4(T2sA,T4A)
 implicit none
 type(Tensor2s), intent(in) :: T2sA
 type(Tensor4), intent(in) :: T4A
 type(Tensor2) :: ddot_T2sT4

 ddot_T2sT4%val(1,1)=T2sA%val(1)*T4A%val(1,1,1,1)+T2sA%val(6)*T4A%val(2,1,1,1)+T2sA%val(5)*T4A%val(3,1,1,1) &
                    +T2sA%val(6)*T4A%val(1,2,1,1)+T2sA%val(2)*T4A%val(2,2,1,1)+T2sA%val(4)*T4A%val(3,2,1,1) &
                    +T2sA%val(5)*T4A%val(1,3,1,1)+T2sA%val(4)*T4A%val(2,3,1,1)+T2sA%val(3)*T4A%val(3,3,1,1)
 ddot_T2sT4%val(2,1)=T2sA%val(1)*T4A%val(1,1,2,1)+T2sA%val(6)*T4A%val(2,1,2,1)+T2sA%val(5)*T4A%val(3,1,2,1) &
                    +T2sA%val(6)*T4A%val(1,2,2,1)+T2sA%val(2)*T4A%val(2,2,2,1)+T2sA%val(4)*T4A%val(3,2,2,1) &
                    +T2sA%val(5)*T4A%val(1,3,2,1)+T2sA%val(4)*T4A%val(2,3,2,1)+T2sA%val(3)*T4A%val(3,3,2,1)
 ddot_T2sT4%val(3,1)=T2sA%val(1)*T4A%val(1,1,3,1)+T2sA%val(6)*T4A%val(2,1,3,1)+T2sA%val(5)*T4A%val(3,1,3,1) &
                    +T2sA%val(6)*T4A%val(1,2,3,1)+T2sA%val(2)*T4A%val(2,2,3,1)+T2sA%val(4)*T4A%val(3,2,3,1) &
                    +T2sA%val(5)*T4A%val(1,3,3,1)+T2sA%val(4)*T4A%val(2,3,3,1)+T2sA%val(3)*T4A%val(3,3,3,1)
 ddot_T2sT4%val(1,2)=T2sA%val(1)*T4A%val(1,1,1,2)+T2sA%val(6)*T4A%val(2,1,1,2)+T2sA%val(5)*T4A%val(3,1,1,2) &
                    +T2sA%val(6)*T4A%val(1,2,1,2)+T2sA%val(2)*T4A%val(2,2,1,2)+T2sA%val(4)*T4A%val(3,2,1,2) &
                    +T2sA%val(5)*T4A%val(1,3,1,2)+T2sA%val(4)*T4A%val(2,3,1,2)+T2sA%val(3)*T4A%val(3,3,1,2)
 ddot_T2sT4%val(2,2)=T2sA%val(1)*T4A%val(1,1,2,2)+T2sA%val(6)*T4A%val(2,1,2,2)+T2sA%val(5)*T4A%val(3,1,2,2) &
                    +T2sA%val(6)*T4A%val(1,2,2,2)+T2sA%val(2)*T4A%val(2,2,2,2)+T2sA%val(4)*T4A%val(3,2,2,2) &
                    +T2sA%val(5)*T4A%val(1,3,2,2)+T2sA%val(4)*T4A%val(2,3,2,2)+T2sA%val(3)*T4A%val(3,3,2,2)
 ddot_T2sT4%val(3,2)=T2sA%val(1)*T4A%val(1,1,3,2)+T2sA%val(6)*T4A%val(2,1,3,2)+T2sA%val(5)*T4A%val(3,1,3,2) &
                    +T2sA%val(6)*T4A%val(1,2,3,2)+T2sA%val(2)*T4A%val(2,2,3,2)+T2sA%val(4)*T4A%val(3,2,3,2) &
                    +T2sA%val(5)*T4A%val(1,3,3,2)+T2sA%val(4)*T4A%val(2,3,3,2)+T2sA%val(3)*T4A%val(3,3,3,2)
 ddot_T2sT4%val(1,3)=T2sA%val(1)*T4A%val(1,1,1,3)+T2sA%val(6)*T4A%val(2,1,1,3)+T2sA%val(5)*T4A%val(3,1,1,3) &
                    +T2sA%val(6)*T4A%val(1,2,1,3)+T2sA%val(2)*T4A%val(2,2,1,3)+T2sA%val(4)*T4A%val(3,2,1,3) &
                    +T2sA%val(5)*T4A%val(1,3,1,3)+T2sA%val(4)*T4A%val(2,3,1,3)+T2sA%val(3)*T4A%val(3,3,1,3)
 ddot_T2sT4%val(2,3)=T2sA%val(1)*T4A%val(1,1,2,3)+T2sA%val(6)*T4A%val(2,1,2,3)+T2sA%val(5)*T4A%val(3,1,2,3) &
                    +T2sA%val(6)*T4A%val(1,2,2,3)+T2sA%val(2)*T4A%val(2,2,2,3)+T2sA%val(4)*T4A%val(3,2,2,3) &
                    +T2sA%val(5)*T4A%val(1,3,2,3)+T2sA%val(4)*T4A%val(2,3,2,3)+T2sA%val(3)*T4A%val(3,3,2,3)
 ddot_T2sT4%val(3,3)=T2sA%val(1)*T4A%val(1,1,3,3)+T2sA%val(6)*T4A%val(2,1,3,3)+T2sA%val(5)*T4A%val(3,1,3,3) &
                    +T2sA%val(6)*T4A%val(1,2,3,3)+T2sA%val(2)*T4A%val(2,2,3,3)+T2sA%val(4)*T4A%val(3,2,3,3) &
                    +T2sA%val(5)*T4A%val(1,3,3,3)+T2sA%val(4)*T4A%val(2,3,3,3)+T2sA%val(3)*T4A%val(3,3,3,3)
                   
end function ddot_T2sT4
! **********************************************************************
! *** function ddot_T4T2w                                            ***
! *** returns the double dot product of one general fourth rank      ***
! *** tensor and one skew-symmetric second order tensor              ***
! ***                                                                ***
! *** input : T4A, fourth order tensor                               ***
! ***         T2wA, second order tensor                              ***
! **********************************************************************
function ddot_T4T2w(T4A,T2wA)
 implicit none
 type(Tensor4), intent(in) :: T4A
 type(Tensor2w), intent(in) :: T2wA
 type(Tensor2) :: ddot_T4T2w

 ddot_T4T2w%val(1,1)=-T4A%val(1,1,2,1)*T2wA%val(3)+T4A%val(1,1,3,1)*T2wA%val(2) &
                     +T4A%val(1,1,1,2)*T2wA%val(3)-T4A%val(1,1,3,2)*T2wA%val(1) &
                     -T4A%val(1,1,1,3)*T2wA%val(2)+T4A%val(1,1,2,3)*T2wA%val(1)
 ddot_T4T2w%val(2,1)=-T4A%val(2,1,2,1)*T2wA%val(3)+T4A%val(2,1,3,1)*T2wA%val(2) &
                     +T4A%val(2,1,1,2)*T2wA%val(3)-T4A%val(2,1,3,2)*T2wA%val(1) &
                     -T4A%val(2,1,1,3)*T2wA%val(2)+T4A%val(2,1,2,3)*T2wA%val(1)
 ddot_T4T2w%val(3,1)=-T4A%val(3,1,2,1)*T2wA%val(3)+T4A%val(3,1,3,1)*T2wA%val(2) &
                     +T4A%val(3,1,1,2)*T2wA%val(3)-T4A%val(3,1,3,2)*T2wA%val(1) &
                     -T4A%val(3,1,1,3)*T2wA%val(2)+T4A%val(3,1,2,3)*T2wA%val(1)
 ddot_T4T2w%val(1,2)=-T4A%val(1,2,2,1)*T2wA%val(3)+T4A%val(1,2,3,1)*T2wA%val(2) &
                     +T4A%val(1,2,1,2)*T2wA%val(3)-T4A%val(1,2,3,2)*T2wA%val(1) &
                     -T4A%val(1,2,1,3)*T2wA%val(2)+T4A%val(1,2,2,3)*T2wA%val(1)
 ddot_T4T2w%val(2,2)=-T4A%val(2,2,2,1)*T2wA%val(3)+T4A%val(2,2,3,1)*T2wA%val(2) &
                     +T4A%val(2,2,1,2)*T2wA%val(3)-T4A%val(2,2,3,2)*T2wA%val(1) &
                     -T4A%val(2,2,1,3)*T2wA%val(2)+T4A%val(2,2,2,3)*T2wA%val(1)
 ddot_T4T2w%val(3,2)=-T4A%val(3,2,2,1)*T2wA%val(3)+T4A%val(3,2,3,1)*T2wA%val(2) &
                     +T4A%val(3,2,1,2)*T2wA%val(3)-T4A%val(3,2,3,2)*T2wA%val(1) &
                     -T4A%val(3,2,1,3)*T2wA%val(2)+T4A%val(3,2,2,3)*T2wA%val(1)
 ddot_T4T2w%val(1,3)=-T4A%val(1,3,2,1)*T2wA%val(3)+T4A%val(1,3,3,1)*T2wA%val(2) &
                     +T4A%val(1,3,1,2)*T2wA%val(3)-T4A%val(1,3,3,2)*T2wA%val(1) &
                     -T4A%val(1,3,1,3)*T2wA%val(2)+T4A%val(1,3,2,3)*T2wA%val(1)
 ddot_T4T2w%val(2,3)=-T4A%val(2,3,2,1)*T2wA%val(3)+T4A%val(2,3,3,1)*T2wA%val(2) &
                     +T4A%val(2,3,1,2)*T2wA%val(3)-T4A%val(2,3,3,2)*T2wA%val(1) &
                     -T4A%val(2,3,1,3)*T2wA%val(2)+T4A%val(2,3,2,3)*T2wA%val(1)
 ddot_T4T2w%val(3,3)=-T4A%val(3,3,2,1)*T2wA%val(3)+T4A%val(3,3,3,1)*T2wA%val(2) &
                     +T4A%val(3,3,1,2)*T2wA%val(3)-T4A%val(3,3,3,2)*T2wA%val(1) &
                     -T4A%val(3,3,1,3)*T2wA%val(2)+T4A%val(3,3,2,3)*T2wA%val(1)

                   
end function ddot_T4T2w
! **********************************************************************
! *** function ddot_T2wT4                                            ***
! *** returns the double dot product of one skew-symmetric second    ***
! *** order tensor and one general fourth order tensor               ***
! ***                                                                ***
! *** input : T4A, fourth order tensor                               ***
! ***         T2wA, second order tensor                              ***
! **********************************************************************
function ddot_T2wT4(T2wA,T4A)
 implicit none
 type(Tensor2w), intent(in) :: T2wA
 type(Tensor4), intent(in) :: T4A
 type(Tensor2) :: ddot_T2wT4

 ddot_T2wT4%val(1,1)=-T2wA%val(3)*T4A%val(2,1,1,1)+T2wA%val(2)*T4A%val(3,1,1,1) &
                     +T2wA%val(3)*T4A%val(1,2,1,1)-T2wA%val(1)*T4A%val(3,2,1,1) &
                     -T2wA%val(2)*T4A%val(1,3,1,1)+T2wA%val(1)*T4A%val(2,3,1,1)
 ddot_T2wT4%val(2,1)=-T2wA%val(3)*T4A%val(2,1,2,1)+T2wA%val(2)*T4A%val(3,1,2,1) &
                     +T2wA%val(3)*T4A%val(1,2,2,1)-T2wA%val(1)*T4A%val(3,2,2,1) &
                     -T2wA%val(2)*T4A%val(1,3,2,1)+T2wA%val(1)*T4A%val(2,3,2,1)
 ddot_T2wT4%val(3,1)=-T2wA%val(3)*T4A%val(2,1,3,1)+T2wA%val(2)*T4A%val(3,1,3,1) &
                     +T2wA%val(3)*T4A%val(1,2,3,1)-T2wA%val(1)*T4A%val(3,2,3,1) &
                     -T2wA%val(2)*T4A%val(1,3,3,1)+T2wA%val(1)*T4A%val(2,3,3,1)
 ddot_T2wT4%val(1,2)=-T2wA%val(3)*T4A%val(2,1,1,2)+T2wA%val(2)*T4A%val(3,1,1,2) &
                     +T2wA%val(3)*T4A%val(1,2,1,2)-T2wA%val(1)*T4A%val(3,2,1,2) &
                     -T2wA%val(2)*T4A%val(1,3,1,2)+T2wA%val(1)*T4A%val(2,3,1,2)
 ddot_T2wT4%val(2,2)=-T2wA%val(3)*T4A%val(2,1,2,2)+T2wA%val(2)*T4A%val(3,1,2,2) &
                     +T2wA%val(3)*T4A%val(1,2,2,2)-T2wA%val(1)*T4A%val(3,2,2,2) &
                     -T2wA%val(2)*T4A%val(1,3,2,2)+T2wA%val(1)*T4A%val(2,3,2,2)
 ddot_T2wT4%val(3,2)=-T2wA%val(3)*T4A%val(2,1,3,2)+T2wA%val(2)*T4A%val(3,1,3,2) &
                     +T2wA%val(3)*T4A%val(1,2,3,2)-T2wA%val(1)*T4A%val(3,2,3,2) &
                     -T2wA%val(2)*T4A%val(1,3,3,2)+T2wA%val(1)*T4A%val(2,3,3,2)
 ddot_T2wT4%val(1,3)=-T2wA%val(3)*T4A%val(2,1,1,3)+T2wA%val(2)*T4A%val(3,1,1,3) &
                     +T2wA%val(3)*T4A%val(1,2,1,3)-T2wA%val(1)*T4A%val(3,2,1,3) &
                     -T2wA%val(2)*T4A%val(1,3,1,3)+T2wA%val(1)*T4A%val(2,3,1,3)
 ddot_T2wT4%val(2,3)=-T2wA%val(3)*T4A%val(2,1,2,3)+T2wA%val(2)*T4A%val(3,1,2,3) &
                     +T2wA%val(3)*T4A%val(1,2,2,3)-T2wA%val(1)*T4A%val(3,2,2,3) &
                     -T2wA%val(2)*T4A%val(1,3,2,3)+T2wA%val(1)*T4A%val(2,3,2,3)
 ddot_T2wT4%val(3,3)=-T2wA%val(3)*T4A%val(2,1,3,3)+T2wA%val(2)*T4A%val(3,1,3,3) &
                     +T2wA%val(3)*T4A%val(1,2,3,3)-T2wA%val(1)*T4A%val(3,2,3,3) &
                     -T2wA%val(2)*T4A%val(1,3,3,3)+T2wA%val(1)*T4A%val(2,3,3,3)
                   
end function ddot_T2wT4
! **********************************************************************
! *** function ddot_T4 T4                                            ***
! *** returns the double dot product of two general fourth order     ***
! *** tensors                                                        ***
! ***                                                                ***
! *** input : T4A, fourth order tensor                               ***
! ***         T4B, fourth order tensor                               ***
! **********************************************************************
function ddot_T4T4(T4A,T4B)
 implicit none
 type(Tensor4), intent(in) :: T4A,T4B
 type(Tensor4) :: ddot_T4T4
 integer :: i,j,k,l,m,n

 do l=1,3
  do k=1,3
   do j=1,3
    do i=1,3
     ddot_T4T4%val(i,j,k,l)=0.d0
     do n=1,3
      do m=1,3
       ddot_T4T4%val(i,j,k,l)=ddot_T4T4%val(i,j,k,l)+T4A%val(i,j,m,n)*T4B%val(m,n,k,l)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
                   
end function ddot_T4T4
! **********************************************************************
! *** function ddot_T4T4ss                                           ***
! *** returns the double dot product of one general and one symmetric***
! *** fourth order tensors                                           ***
! ***                                                                ***
! *** input : T4A, fourth order tensor                               ***
! ***         T4ssB, fourth order tensor                             ***
! **********************************************************************
function ddot_T4T4ss(T4A,T4ssB)
 implicit none
 type(Tensor4), intent(in) :: T4A
 type(Tensor4ss), intent(in) :: T4ssB
 type(Tensor4) :: ddot_T4T4ss
 integer :: i,j,k,l,m,n,kl,mn

 do l=1,3
  do k=1,3
   kl=indxs(k,l)
   do j=1,3
    do i=1,3
     ddot_T4T4ss%val(i,j,k,l)=0.d0
     do n=1,3
      do m=1,3
       mn=indxs(m,n)
       ddot_T4T4ss%val(i,j,k,l)=ddot_T4T4ss%val(i,j,k,l)+T4A%val(i,j,m,n)*T4ssB%val(mn,kl)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
                   
end function ddot_T4T4ss
! **********************************************************************
! *** function ddot_T4ssT4                                           ***
! *** returns the double dot product of one symmetric and one general***
! *** fourth order tensors                                           ***
! ***                                                                ***
! *** input : T4ssA, fourth order tensor                             ***
! ***         T4B, fourth order tensor                               ***
! **********************************************************************
function ddot_T4ssT4(T4ssA,T4B)
 implicit none
 type(Tensor4ss), intent(in) :: T4ssA
 type(Tensor4), intent(in) :: T4B
 type(Tensor4) :: ddot_T4ssT4
 integer :: i,j,k,l,m,n,ij,mn

 do l=1,3
  do k=1,3
   do j=1,3
    do i=1,3
     ij=indxs(i,j)
     ddot_T4ssT4%val(i,j,k,l)=0.d0
     do n=1,3
      do m=1,3
       mn=indxs(m,n)
       ddot_T4ssT4%val(i,j,k,l)=ddot_T4ssT4%val(i,j,k,l)+T4ssA%val(ij,mn)*T4B%val(m,n,k,l)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
                   
end function ddot_T4ssT4
! **********************************************************************
! *** function ddot_T4ssT4ss                                         ***
! *** returns the double dot product of two symmetric fourth order   ***
! *** tensors                                                        ***
! ***                                                                ***
! *** input : T4ssA, fourth order tensor                             ***
! ***         T4ssB, fourth order tensor                             ***
! **********************************************************************
function ddot_T4ssT4ss(T4ssA,T4ssB)
 implicit none
 type(Tensor4ss), intent(in) :: T4ssA,T4ssB
 type(Tensor4ss) :: ddot_T4ssT4ss
 integer :: i,j,k

 do j=1,6
  do i=1,6
   ddot_T4ssT4ss%val(i,j)=0.d0
   do k=1,6
    ddot_T4ssT4ss%val(i,j)=ddot_T4ssT4ss%val(i,j)+T4ssA%val(i,k)*T4ssB%val(k,j)*sfac(k)
   enddo
  enddo
 enddo
                   
end function ddot_T4ssT4ss
! **********************************************************************
! *** function ddddot_T4T4                                           ***
! *** returns the quadruple dot product of two fourth order tensors  ***
! ***                                                                ***
! *** input : T4A, fourth order tensor                               ***
! ***         T4B, fourth order tensor                               ***
! **********************************************************************
function ddddot_T4T4(T4A,T4B)
 implicit none
 type(Tensor4), intent(in) :: T4A,T4B
 real(8) :: ddddot_T4T4
 integer :: i,j,k,l
 
 ddddot_T4T4=0.d0
 do l=1,3
  do k=1,3
   do j=1,3
    do i=1,3
     ddddot_T4T4=ddddot_T4T4+T4A%val(i,j,k,l)*T4B%val(i,j,k,l)
    enddo
   enddo
  enddo
 enddo
                   
end function ddddot_T4T4
! **********************************************************************
! *** function ddddot_T4ssT4ss                                       ***
! *** returns the quadruple dot product of two symmetric fourth order***
! *** tensors                                                        ***
! ***                                                                ***
! *** input : T4ssA, fourth order tensor                             ***
! ***         T4ssB, fourth order tensor                             ***
! **********************************************************************
function ddddot_T4ssT4ss(T4ssA,T4ssB)
 implicit none
 type(Tensor4ss), intent(in) :: T4ssA,T4ssB
 real(8) :: ddddot_T4ssT4ss
 integer :: i,j
 
 ddddot_T4ssT4ss=0.d0
 do j=1,6
  do i=1,6
   ddddot_T4ssT4ss=ddddot_T4ssT4ss+T4ssA%val(i,j)*T4ssB%val(i,j)*sfac(i)*sfac(j)
  enddo
 enddo
                   
end function ddddot_T4ssT4ss
! **********************************************************************
! *** function tensor_product_T1T1                                   ***
! *** returns the tensor product of two first order tensors          ***
! ***                                                                ***
! *** input : T1A, first order tensor                                ***
! ***         T1B, first order tensor                                ***
! **********************************************************************
function tensor_product_T1T1(T1A,T1B)
 implicit none
 type(Tensor1), intent(in) :: T1A,T1B
 type(Tensor2) :: tensor_product_T1T1

 tensor_product_T1T1%val(1,1)=T1A%val(1)*T1B%val(1)
 tensor_product_T1T1%val(2,1)=T1A%val(2)*T1B%val(1)
 tensor_product_T1T1%val(3,1)=T1A%val(3)*T1B%val(1)
 
 tensor_product_T1T1%val(1,2)=T1A%val(1)*T1B%val(2)
 tensor_product_T1T1%val(2,2)=T1A%val(2)*T1B%val(2)
 tensor_product_T1T1%val(3,2)=T1A%val(3)*T1B%val(2)
 
 tensor_product_T1T1%val(1,3)=T1A%val(1)*T1B%val(3)
 tensor_product_T1T1%val(2,3)=T1A%val(2)*T1B%val(3)
 tensor_product_T1T1%val(3,3)=T1A%val(3)*T1B%val(3)
                   
end function tensor_product_T1T1
! **********************************************************************
! *** function tensor_product_T2T1                                   ***
! *** returns the tensor product of one second and one first order   ***
! *** tensors                                                        ***
! ***                                                                ***
! *** input : T2A, second order tensor                               ***
! ***         T1B, first order tensor                                ***
! **********************************************************************
function tensor_product_T2T1(T2A,T1B)
 implicit none
 type(Tensor2), intent(in) :: T2A
 type(Tensor1), intent(in) :: T1B
 type(Tensor3) :: tensor_product_T2T1
 integer :: i,j,k
 
 do k=1,3
  do j=1,3
   do i=1,3
    tensor_product_T2T1%val(i,j,k)=T2A%val(i,j)*T1B%val(k)
   enddo
  enddo
 enddo
 
end function tensor_product_T2T1
! **********************************************************************
! *** function tensor_product_T1T2                                   ***
! *** returns the tensor product of one first and one second order   ***
! *** tensors                                                        ***
! ***                                                                ***
! *** input : T1A, first order tensor                                ***
! ***         T2B, second order tensor                               ***
! **********************************************************************
function tensor_product_T1T2(T1A,T2B)
 implicit none
 type(Tensor1), intent(in) :: T1A
 type(Tensor2), intent(in) :: T2B
 type(Tensor3) :: tensor_product_T1T2
 integer :: i,j,k
 
 do k=1,3
  do j=1,3
   do i=1,3
    tensor_product_T1T2%val(i,j,k)=T1A%val(i)*T2B%val(j,k)
   enddo
  enddo
 enddo
 
end function tensor_product_T1T2
! **********************************************************************
! *** function tensor_product_T2T2                                   ***
! *** returns the tensor product of two general second order tensors ***
! ***                                                                ***
! *** input : T2A, second order tensor                               ***
! ***         T2B, second order tensor                               ***
! **********************************************************************
function tensor_product_T2T2(T2A,T2B)
 implicit none
 type(Tensor2), intent(in) :: T2A,T2B
 type(Tensor4) :: tensor_product_T2T2
 integer :: i,j,k,l
 
 do l=1,3
  do k=1,3
   do j=1,3
    do i=1,3
     tensor_product_T2T2%val(i,j,k,l)=T2A%val(i,j)*T2B%val(k,l)
    enddo
   enddo
  enddo
 enddo
 
end function tensor_product_T2T2
! **********************************************************************
! *** function tensor_product_T2sT2s                                 ***
! *** returns the tensor product of two symmetric second order       ***
! *** tensors                                                        ***
! ***                                                                ***
! *** input : T2A, second order tensor                               ***
! ***         T2B, second order tensor                               ***
! **********************************************************************
function tensor_product_T2sT2s(T2sA,T2sB)
 implicit none
 type(Tensor2s), intent(in) :: T2sA,T2sB
 type(Tensor4ss) :: tensor_product_T2sT2s
 integer :: i,j
 
 do j=1,6
  do i=1,6
   tensor_product_T2sT2s%val(i,j)=T2sA%val(i)*T2sB%val(j)
  enddo
 enddo
 
end function tensor_product_T2sT2s
! **********************************************************************
! *** function times_T2T2                                            ***
! *** returns the tensor product of two general second order tensors ***
! ***                                                                ***
! *** input : T2A, second order tensor                               ***
! ***         T2B, second order tensor                               ***
! **********************************************************************
function times_T2T2(T2A,T2B)
 implicit none
 type(Tensor2), intent(in) :: T2A,T2B
 type(Tensor4) :: times_T2T2
 integer :: i,j,k,l
 
 do l=1,3
  do k=1,3
   do j=1,3
    do i=1,3
     times_T2T2%val(i,j,k,l)=T2A%val(i,j)*T2B%val(k,l)
    enddo
   enddo
  enddo
 enddo
 
end function times_T2T2
! **********************************************************************
! *** function bartimes_T2T2                                         ***
! *** returns the tensor product of two general second order tensors ***
! ***                                                                ***
! *** input : T2A, second order tensor                               ***
! ***         T2B, second order tensor                               ***
! **********************************************************************
function bartimes_T2T2(T2A,T2B)
 implicit none
 type(Tensor2), intent(in) :: T2A,T2B
 type(Tensor4) :: bartimes_T2T2
 integer :: i,j,k,l
 
 do l=1,3
  do k=1,3
   do j=1,3
    do i=1,3
     bartimes_T2T2%val(i,j,k,l)=T2A%val(i,k)*T2B%val(j,l)
    enddo
   enddo
  enddo
 enddo
 
end function bartimes_T2T2
! **********************************************************************
! *** function bartimes_T2T2s                                        ***
! *** returns the tensor product of one general second order tensor  ***
! *** and one symmetric second order tensor                          ***
! ***                                                                ***
! *** input : T2A, second order tensor                               ***
! ***         T2sB, symmetric second order tensor                    ***
! **********************************************************************
function bartimes_T2T2s(T2A,T2sB)
 implicit none
 type(Tensor2), intent(in) :: T2A
 type(Tensor2s), intent(in) :: T2sB
 type(Tensor4) :: bartimes_T2T2s
 integer :: i,j,k,l,jl
 
 do l=1,3
  do k=1,3
   do j=1,3
    jl=indxs(j,l)
    do i=1,3
     bartimes_T2T2s%val(i,j,k,l)=T2A%val(i,k)*T2sB%val(jl)
    enddo
   enddo
  enddo
 enddo
 
end function bartimes_T2T2s
! **********************************************************************
! *** function bartimes_T2sT2                                        ***
! *** returns the tensor product of one symmetric second order tensor***
! *** and one general second order tensor                            ***
! ***                                                                ***
! *** input : T2sA, symmetrix second order tensor                    ***
! ***         T2B, second order tensor                               ***
! **********************************************************************
function bartimes_T2sT2(T2sA,T2B)
 implicit none
 type(Tensor2s), intent(in) :: T2sA
 type(Tensor2), intent(in) :: T2B
 type(Tensor4) :: bartimes_T2sT2
 integer :: i,j,k,l,ik
 
 do l=1,3
  do k=1,3
   do j=1,3
    do i=1,3
     ik=indxs(i,k)
     bartimes_T2sT2%val(i,j,k,l)=T2sA%val(ik)*T2B%val(j,l)
    enddo
   enddo
  enddo
 enddo
 
end function bartimes_T2sT2

! **********************************************************************
! *** function bartimes_T2sT2s                                       ***
! *** returns the tensor product of two symmetric second order       ***
! *** tensors                                                        ***
! ***                                                                ***
! *** input : T2sA, symmetric second order tensor                    ***
! ***         T2sB, symmetric second order tensor                    ***
! **********************************************************************
function bartimes_T2sT2s(T2sA,T2sB)
 implicit none
 type(Tensor2s), intent(in) :: T2sA
 type(Tensor2s), intent(in) :: T2sB
 type(Tensor4) :: bartimes_T2sT2s
 integer :: i,j,k,l,ik,jl

 do l=1,3
  do k=1,3
   do j=1,3
    do i=1,3
     ik=indxs(i,k)
     jl=indxs(j,l)
     bartimes_T2sT2s%val(i,j,k,l)=T2sA%val(ik)*T2sB%val(jl)
    enddo
   enddo
  enddo
 enddo

end function bartimes_T2sT2s
! **********************************************************************
! *** function timesbar_T2T2                                         ***
! *** returns the tensor product of two general second order tensors ***
! ***                                                                ***
! *** input : T2A, second order tensor                               ***
! ***         T2B, second order tensor                               ***
! **********************************************************************
function timesbar_T2T2(T2A,T2B)
 implicit none
 type(Tensor2), intent(in) :: T2A,T2B
 type(Tensor4) :: timesbar_T2T2
 integer :: i,j,k,l
 
 do l=1,3
  do k=1,3
   do j=1,3
    do i=1,3
     timesbar_T2T2%val(i,j,k,l)=T2A%val(i,l)*T2B%val(j,k)
    enddo
   enddo
  enddo
 enddo
 
end function timesbar_T2T2
! **********************************************************************
! *** function bartimesbar_T2T2                                      ***
! *** returns the tensor product of two general second order tensors ***
! ***                                                                ***
! *** input : T2A, second order tensor                               ***
! ***         T2B, second order tensor                               ***
! **********************************************************************
function bartimesbar_T2T2(T2A,T2B)
 implicit none
 type(Tensor2), intent(in) :: T2A,T2B
 type(Tensor4) :: bartimesbar_T2T2
 integer :: i,j,k,l
 
 do l=1,3
  do k=1,3
   do j=1,3
    do i=1,3
     bartimesbar_T2T2%val(i,j,k,l)=0.5d0*(T2A%val(i,l)*T2B%val(j,k)+T2A%val(j,l)*T2B%val(i,k))
    enddo
   enddo
  enddo
 enddo
 
end function bartimesbar_T2T2
! **********************************************************************
! *** function cross_product_T1T1                                    ***
! *** returns the cross product of two first order tensors           ***
! ***                                                                ***
! *** input : T1A, first order tensor                                ***
! ***         T1B, first order tensor                                ***
! **********************************************************************
function cross_product_T1T1(T1A,T1B)
 implicit none
 type(Tensor1), intent(in) :: T1A,T1B
 type(Tensor1) :: cross_product_T1T1

 cross_product_T1T1%val(1)=T1A%val(2)*T1B%val(3)-T1A%val(3)*T1B%val(2)
 cross_product_T1T1%val(2)=T1A%val(3)*T1B%val(1)-T1A%val(1)*T1B%val(3)
 cross_product_T1T1%val(3)=T1A%val(1)*T1B%val(2)-T1A%val(2)*T1B%val(1)
                   
end function cross_product_T1T1
! **********************************************************************
! *** function roteate_Tensor1                                        ***
! *** rotates a first order tensor according to a specified rotation ***
! *** matrix                                                         ***
! ***                                                                ***
! *** input : rot, rotation matrix (eAi=rot.eBj)                     ***
! ***         T1A, first order tensor to be rotated                  ***
! **********************************************************************
function rotate_Tensor1(T1A,rot)
 implicit none
 type(Tensor1), intent(in) :: T1A
 real(8), intent(in) :: rot(3,3)
 type(Tensor1) :: rotate_Tensor1
 
 rotate_Tensor1%val(1)=T1A%val(1)*rot(1,1)+T1A%val(2)*rot(2,1)+T1A%val(3)*rot(3,1)
 rotate_Tensor1%val(2)=T1A%val(1)*rot(1,2)+T1A%val(2)*rot(2,2)+T1A%val(3)*rot(3,2)
 rotate_Tensor1%val(3)=T1A%val(1)*rot(1,3)+T1A%val(2)*rot(2,3)+T1A%val(3)*rot(3,3)
 
end function rotate_Tensor1
! **********************************************************************
! *** function rotate_Tensor2                                        ***
! *** rotates a second order tensor according to a specified rotation***
! *** matrix                                                         ***
! ***                                                                ***
! *** input : rot, rotation matrix (eAi=rot.eBI)                     ***
! ***         T2A, second order tensor to be rotated                 ***
! **********************************************************************
function rotate_Tensor2(T2A,rot)
 implicit none
 type(Tensor2), intent(in) :: T2A
 real(8), intent(in) :: rot(3,3)
 type(Tensor2) :: rotate_Tensor2
 
 rotate_Tensor2%val(1,1)=T2A%val(1,1)*rot(1,1)*rot(1,1)+T2A%val(2,1)*rot(2,1)*rot(1,1)+T2A%val(3,1)*rot(3,1)*rot(1,1) &
                        +T2A%val(1,2)*rot(1,1)*rot(2,1)+T2A%val(2,2)*rot(2,1)*rot(2,1)+T2A%val(3,2)*rot(3,1)*rot(2,1) &
                        +T2A%val(1,3)*rot(1,1)*rot(3,1)+T2A%val(2,3)*rot(2,1)*rot(3,1)+T2A%val(3,3)*rot(3,1)*rot(3,1)
 rotate_Tensor2%val(2,1)=T2A%val(1,1)*rot(1,2)*rot(1,1)+T2A%val(2,1)*rot(2,2)*rot(1,1)+T2A%val(3,1)*rot(3,2)*rot(1,1) &
                        +T2A%val(1,2)*rot(1,2)*rot(2,1)+T2A%val(2,2)*rot(2,2)*rot(2,1)+T2A%val(3,2)*rot(3,2)*rot(2,1) &
                        +T2A%val(1,3)*rot(1,2)*rot(3,1)+T2A%val(2,3)*rot(2,2)*rot(3,1)+T2A%val(3,3)*rot(3,2)*rot(3,1)
 rotate_Tensor2%val(3,1)=T2A%val(1,1)*rot(1,3)*rot(1,1)+T2A%val(2,1)*rot(2,3)*rot(1,1)+T2A%val(3,1)*rot(3,3)*rot(1,1) &
                        +T2A%val(1,2)*rot(1,3)*rot(2,1)+T2A%val(2,2)*rot(2,3)*rot(2,1)+T2A%val(3,2)*rot(3,3)*rot(2,1) &
                        +T2A%val(1,3)*rot(1,3)*rot(3,1)+T2A%val(2,3)*rot(2,3)*rot(3,1)+T2A%val(3,3)*rot(3,3)*rot(3,1)
             
 rotate_Tensor2%val(1,2)=T2A%val(1,1)*rot(1,1)*rot(1,2)+T2A%val(2,1)*rot(2,1)*rot(1,2)+T2A%val(3,1)*rot(3,1)*rot(1,2) &
                        +T2A%val(1,2)*rot(1,1)*rot(2,2)+T2A%val(2,2)*rot(2,1)*rot(2,2)+T2A%val(3,2)*rot(3,1)*rot(2,2) &
                        +T2A%val(1,3)*rot(1,1)*rot(3,2)+T2A%val(2,3)*rot(2,1)*rot(3,2)+T2A%val(3,3)*rot(3,1)*rot(3,2)
 rotate_Tensor2%val(2,2)=T2A%val(1,1)*rot(1,2)*rot(1,2)+T2A%val(2,1)*rot(2,2)*rot(1,2)+T2A%val(3,1)*rot(3,2)*rot(1,2) &
                        +T2A%val(1,2)*rot(1,2)*rot(2,2)+T2A%val(2,2)*rot(2,2)*rot(2,2)+T2A%val(3,2)*rot(3,2)*rot(2,2) &
                        +T2A%val(1,3)*rot(1,2)*rot(3,2)+T2A%val(2,3)*rot(2,2)*rot(3,2)+T2A%val(3,3)*rot(3,2)*rot(3,2)
 rotate_Tensor2%val(3,2)=T2A%val(1,1)*rot(1,3)*rot(1,2)+T2A%val(2,1)*rot(2,3)*rot(1,2)+T2A%val(3,1)*rot(3,3)*rot(1,2) &
                        +T2A%val(1,2)*rot(1,3)*rot(2,2)+T2A%val(2,2)*rot(2,3)*rot(2,2)+T2A%val(3,2)*rot(3,3)*rot(2,2) &
                        +T2A%val(1,3)*rot(1,3)*rot(3,2)+T2A%val(2,3)*rot(2,3)*rot(3,2)+T2A%val(3,3)*rot(3,3)*rot(3,2)

 rotate_Tensor2%val(1,3)=T2A%val(1,1)*rot(1,1)*rot(1,3)+T2A%val(2,1)*rot(2,1)*rot(1,3)+T2A%val(3,1)*rot(3,1)*rot(1,3) &
                        +T2A%val(1,2)*rot(1,1)*rot(2,3)+T2A%val(2,2)*rot(2,1)*rot(2,3)+T2A%val(3,2)*rot(3,1)*rot(2,3) &
                        +T2A%val(1,3)*rot(1,1)*rot(3,3)+T2A%val(2,3)*rot(2,1)*rot(3,3)+T2A%val(3,3)*rot(3,1)*rot(3,3)
 rotate_Tensor2%val(2,3)=T2A%val(1,1)*rot(1,2)*rot(1,3)+T2A%val(2,1)*rot(2,2)*rot(1,3)+T2A%val(3,1)*rot(3,2)*rot(1,3) &
                        +T2A%val(1,2)*rot(1,2)*rot(2,3)+T2A%val(2,2)*rot(2,2)*rot(2,3)+T2A%val(3,2)*rot(3,2)*rot(2,3) &
                        +T2A%val(1,3)*rot(1,2)*rot(3,3)+T2A%val(2,3)*rot(2,2)*rot(3,3)+T2A%val(3,3)*rot(3,2)*rot(3,3)
 rotate_Tensor2%val(3,3)=T2A%val(1,1)*rot(1,3)*rot(1,3)+T2A%val(2,1)*rot(2,3)*rot(1,3)+T2A%val(3,1)*rot(3,3)*rot(1,3) &
                        +T2A%val(1,2)*rot(1,3)*rot(2,3)+T2A%val(2,2)*rot(2,3)*rot(2,3)+T2A%val(3,2)*rot(3,3)*rot(2,3) &
                        +T2A%val(1,3)*rot(1,3)*rot(3,3)+T2A%val(2,3)*rot(2,3)*rot(3,3)+T2A%val(3,3)*rot(3,3)*rot(3,3)
  
end function rotate_Tensor2
! **********************************************************************
! *** function rotate_Tensor2s                                       ***
! *** rotates a symmetric second order tensor according to a         ***
! *** specified rotation matrix                                      ***
! ***                                                                ***
! *** input : rot, rotation matrix (eAi=rot.eBI)                     ***
! ***         T2sA, second order tensor to be rotated                ***
! **********************************************************************
function rotate_Tensor2s(T2sA,rot)
 implicit none
 type(Tensor2s), intent(in) :: T2sA
 real(8), intent(in) :: rot(3,3)
 type(Tensor2s) :: rotate_Tensor2s
 
 rotate_Tensor2s%val(1)=T2sA%val(1)*rot(1,1)*rot(1,1)+T2sA%val(6)*rot(2,1)*rot(1,1)+T2sA%val(5)*rot(3,1)*rot(1,1) &
                       +T2sA%val(6)*rot(1,1)*rot(2,1)+T2sA%val(2)*rot(2,1)*rot(2,1)+T2sA%val(4)*rot(3,1)*rot(2,1) &
                       +T2sA%val(5)*rot(1,1)*rot(3,1)+T2sA%val(4)*rot(2,1)*rot(3,1)+T2sA%val(3)*rot(3,1)*rot(3,1)
            
 rotate_Tensor2s%val(2)=T2sA%val(1)*rot(1,2)*rot(1,2)+T2sA%val(6)*rot(2,2)*rot(1,2)+T2sA%val(5)*rot(3,2)*rot(1,2) &
                       +T2sA%val(6)*rot(1,2)*rot(2,2)+T2sA%val(2)*rot(2,2)*rot(2,2)+T2sA%val(4)*rot(3,2)*rot(2,2) &
                       +T2sA%val(5)*rot(1,2)*rot(3,2)+T2sA%val(4)*rot(2,2)*rot(3,2)+T2sA%val(3)*rot(3,2)*rot(3,2)
            
 rotate_Tensor2s%val(3)=T2sA%val(1)*rot(1,3)*rot(1,3)+T2sA%val(6)*rot(2,3)*rot(1,3)+T2sA%val(5)*rot(3,3)*rot(1,3) &
                       +T2sA%val(6)*rot(1,3)*rot(2,3)+T2sA%val(2)*rot(2,3)*rot(2,3)+T2sA%val(4)*rot(3,3)*rot(2,3) &
                       +T2sA%val(5)*rot(1,3)*rot(3,3)+T2sA%val(4)*rot(2,3)*rot(3,3)+T2sA%val(3)*rot(3,3)*rot(3,3)
            
 rotate_Tensor2s%val(4)=T2sA%val(1)*rot(1,3)*rot(1,2)+T2sA%val(6)*rot(2,3)*rot(1,2)+T2sA%val(5)*rot(3,3)*rot(1,2) &
                       +T2sA%val(6)*rot(1,3)*rot(2,2)+T2sA%val(2)*rot(2,3)*rot(2,2)+T2sA%val(4)*rot(3,3)*rot(2,2) &
                       +T2sA%val(5)*rot(1,3)*rot(3,2)+T2sA%val(4)*rot(2,3)*rot(3,2)+T2sA%val(3)*rot(3,3)*rot(3,2)
 
 rotate_Tensor2s%val(5)=T2sA%val(1)*rot(1,1)*rot(1,3)+T2sA%val(6)*rot(2,1)*rot(1,3)+T2sA%val(5)*rot(3,1)*rot(1,3) &
                       +T2sA%val(6)*rot(1,1)*rot(2,3)+T2sA%val(2)*rot(2,1)*rot(2,3)+T2sA%val(4)*rot(3,1)*rot(2,3) &
                       +T2sA%val(5)*rot(1,1)*rot(3,3)+T2sA%val(4)*rot(2,1)*rot(3,3)+T2sA%val(3)*rot(3,1)*rot(3,3)
            
 rotate_Tensor2s%val(6)=T2sA%val(1)*rot(1,2)*rot(1,1)+T2sA%val(6)*rot(2,2)*rot(1,1)+T2sA%val(5)*rot(3,2)*rot(1,1) &
                       +T2sA%val(6)*rot(1,2)*rot(2,1)+T2sA%val(2)*rot(2,2)*rot(2,1)+T2sA%val(4)*rot(3,2)*rot(2,1) &
                       +T2sA%val(5)*rot(1,2)*rot(3,1)+T2sA%val(4)*rot(2,2)*rot(3,1)+T2sA%val(3)*rot(3,2)*rot(3,1)
             
end function rotate_Tensor2s
! **********************************************************************
! *** function rotate_Tensor2w                                       ***
! *** rotates a skew-symmetric second order tensor according to a    ***
! *** specified rotation matrix                                      ***
! ***                                                                ***
! *** input : rot, rotation matrix (eAi=rot.eBI)                     ***
! ***         T2A, second order tensor to be rotated                 ***
! **********************************************************************
function rotate_Tensor2w(T2wA,rot)
 implicit none
 type(Tensor2w), intent(in) :: T2wA
 real(8), intent(in) :: rot(3,3)
 type(Tensor2w) :: rotate_Tensor2w


 rotate_Tensor2w%val(1)=-T2wA%val(3)*rot(2,2)*rot(1,3)+T2wA%val(2)*rot(3,2)*rot(1,3) &
                        +T2wA%val(3)*rot(1,2)*rot(2,3)-T2wA%val(1)*rot(3,2)*rot(2,3) &
                        -T2wA%val(2)*rot(1,2)*rot(3,3)+T2wA%val(1)*rot(2,2)*rot(3,3)

 rotate_Tensor2w%val(2)=-T2wA%val(3)*rot(2,3)*rot(1,1)+T2wA%val(2)*rot(3,3)*rot(1,1) &
                        +T2wA%val(3)*rot(1,3)*rot(2,1)-T2wA%val(1)*rot(3,3)*rot(2,1) &
                        -T2wA%val(2)*rot(1,3)*rot(3,1)+T2wA%val(1)*rot(2,3)*rot(3,1)
             
 rotate_Tensor2w%val(3)=-T2wA%val(3)*rot(2,1)*rot(1,2)+T2wA%val(2)*rot(3,1)*rot(1,2) &
                        +T2wA%val(3)*rot(1,1)*rot(2,2)-T2wA%val(1)*rot(3,1)*rot(2,2) &
                        -T2wA%val(2)*rot(1,1)*rot(3,2)+T2wA%val(1)*rot(2,1)*rot(3,2)

  
end function rotate_Tensor2w
! **********************************************************************
! *** function rotate_Tensor3                                        ***
! *** rotates a third order tensor according to a specified rotation ***
! *** matrix                                                         ***
! ***                                                                ***
! *** input : rot, rotation matrix (eAi=rot.eBI)                     ***
! ***         T3A, third order tensor to be rotated                  ***
! **********************************************************************
function rotate_Tensor3(T3A,rot)
 implicit none
 type(Tensor3), intent(in) :: T3A
 real(8), intent(in) :: rot(3,3)
 type(Tensor3) :: rotate_Tensor3
 integer :: i,j,k,ii,jj,kk
 
 do i=1,3
  do j=1,3
   do k=1,3
    rotate_Tensor3%val(i,j,k)=0.d0
    do ii=1,3
     do jj=1,3
      do kk=1,3
       rotate_Tensor3%val(i,j,k)=rotate_Tensor3%val(i,j,k)+T3A%val(ii,jj,kk)*rot(ii,i)*rot(jj,j)*rot(kk,k)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 
end function rotate_Tensor3
! **********************************************************************
! *** function rotate_Tensor4                                        ***
! *** rotates a fourth order tensor according to a specified rotation***
! *** matrix                                                         ***
! ***                                                                ***
! *** input : rot, rotation matrix (eAi=rot.eBI)                     ***
! ***         T4A, fourth order tensor to be rotated                 ***
! **********************************************************************
function rotate_Tensor4(T4A,rot)
 implicit none
 type(Tensor4), intent(in) :: T4A
 real(8), intent(in) :: rot(3,3)
 type(Tensor4) :: rotate_Tensor4
 integer :: i,j,k,l,ii,jj,kk,ll
 
 do i=1,3
  do j=1,3
   do k=1,3
    do l=1,3
     rotate_Tensor4%val(i,j,k,l)=0.d0
     do ii=1,3
      do jj=1,3
       do kk=1,3
        do ll=1,3
         rotate_Tensor4%val(i,j,k,l)=rotate_Tensor4%val(i,j,k,l)+T4A%val(ii,jj,kk,ll)*rot(ii,i)*rot(jj,j)*rot(kk,k)*rot(ll,l)
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 
end function rotate_Tensor4
! **********************************************************************
! *** function rotate_Tensor4ss                                      ***
! *** rotates a symmetric fourth order tensor according to a         ***
! *** specified rotation matrix                                      ***
! ***                                                                ***
! *** input : rot, rotation matrix (eAi=rot.eBI)                     ***
! ***         T4ssA, fourth order tensor to be rotated               ***
! **********************************************************************
function rotate_Tensor4ss(T4ssA,rot)
 implicit none
 type(Tensor4ss), intent(in) :: T4ssA
 real(8), intent(in) :: rot(3,3)
 type(Tensor4ss) :: rotate_Tensor4ss
 real(8) :: rot2(6,6)
 integer :: i,j,ii,jj
 
! a11=r11*r11*a'11  =  A1=R11*A'1
 rot2(1,1)=rot(1,1)*rot(1,1)
! a11=r21*r21*a'22  =  A1=R21*A'2
 rot2(2,1)=rot(2,1)*rot(2,1)
! a11=r31*r31*a'33  =  A1=R31*A'3
 rot2(3,1)=rot(3,1)*rot(3,1)
! a11=r11*r21*a'12+r21*r11*a'21  =  A1=R61*A'6
 rot2(6,1)=rot(1,1)*rot(2,1)+rot(2,1)*rot(1,1)
! a11=r21*r31*a'23+r31*r21*a'32  =  A1=R41*A'4
 rot2(4,1)=rot(2,1)*rot(3,1)+rot(3,1)*rot(2,1)
! a11=r31*r11*a'31+r11*r31*a'13  =  A1=R51*A'5
 rot2(5,1)=rot(3,1)*rot(1,1)+rot(1,1)*rot(3,1)
 rot2(1,2)=rot(1,2)*rot(1,2)
 rot2(2,2)=rot(2,2)*rot(2,2)
 rot2(3,2)=rot(3,2)*rot(3,2)
 rot2(6,2)=rot(1,2)*rot(2,2)+rot(2,2)*rot(1,2)
 rot2(4,2)=rot(2,2)*rot(3,2)+rot(3,2)*rot(2,2)
 rot2(5,2)=rot(3,2)*rot(1,2)+rot(1,2)*rot(3,2)
 rot2(1,3)=rot(1,3)*rot(1,3)
 rot2(2,3)=rot(2,3)*rot(2,3)
 rot2(3,3)=rot(3,3)*rot(3,3)
 rot2(6,3)=rot(1,3)*rot(2,3)+rot(2,3)*rot(1,3)
 rot2(4,3)=rot(2,3)*rot(3,3)+rot(3,3)*rot(2,3)
 rot2(5,3)=rot(3,3)*rot(1,3)+rot(1,3)*rot(3,3)
! a12=r11*r12*a'11  =  A6=R16*A'1
 rot2(1,6)=rot(1,1)*rot(1,2)
! a12=r21*r22*a'22  =  A6=R26*A'2
 rot2(2,6)=rot(2,1)*rot(2,2)
! a12=r31*r32*a'33  =  A6=R36*A'3
 rot2(3,6)=rot(3,1)*rot(3,2)
! a12=r11*r22*a'12+r21*r12*a'21  =  A6=R66*A'6
 rot2(6,6)=rot(1,1)*rot(2,2)+rot(2,1)*rot(1,2)
! a12=r21*r32*a'23+r31*r22*a'32  =  A6=R46*A'4
 rot2(4,6)=rot(2,1)*rot(3,2)+rot(3,1)*rot(2,2)
! a12=r31*r12*a'31+r11*r32*a'13  =  A6=R56*A'5
 rot2(5,6)=rot(3,1)*rot(1,2)+rot(1,1)*rot(3,2)
 rot2(1,4)=rot(1,2)*rot(1,3)
 rot2(2,4)=rot(2,2)*rot(2,3)
 rot2(3,4)=rot(3,2)*rot(3,3)
 rot2(6,4)=rot(1,2)*rot(2,3)+rot(2,2)*rot(1,3)
 rot2(4,4)=rot(2,2)*rot(3,3)+rot(3,2)*rot(2,3)
 rot2(5,4)=rot(3,2)*rot(1,3)+rot(1,2)*rot(3,3)
 rot2(1,5)=rot(1,3)*rot(1,1)
 rot2(2,5)=rot(2,3)*rot(2,1)
 rot2(3,5)=rot(3,3)*rot(3,1)
 rot2(6,5)=rot(1,3)*rot(2,1)+rot(2,3)*rot(1,1)
 rot2(4,5)=rot(2,3)*rot(3,1)+rot(3,3)*rot(2,1)
 rot2(5,5)=rot(3,3)*rot(1,1)+rot(1,3)*rot(3,1)
 
 do i=1,6
  do j=1,6
   rotate_Tensor4ss%val(i,j)=0.d0
   do ii=1,6
    do jj=1,6
     rotate_Tensor4ss%val(i,j)=rotate_Tensor4ss%val(i,j)+T4ssA%val(ii,jj)*rot2(ii,i)*rot2(jj,j)
    enddo
   enddo
  enddo
 enddo
 
end function rotate_Tensor4ss
! **********************************************************************
! *** function average_Tensor1                                       ***
! *** calculates the average first order tensor associated with a    ***
! *** weighted list of first order tensors                           ***
! ***                                                                ***
! *** input : nT1, number of elements                                ***
! ***         wT1, list of weights                                   ***
! ***         T1, list of first order tensors                        ***
! **********************************************************************
function average_Tensor1(nT1,wT1,T1)
 implicit none
 integer, intent(in) :: nT1
 type(Tensor1), intent(in) :: T1(nT1)
 real(8), intent(in) :: wT1(nT1)
 type(Tensor1) :: average_Tensor1
 integer :: iT1
 
 average_Tensor1=null1
 do iT1=1,nT1
  average_Tensor1=average_Tensor1+(wT1(iT1)*T1(iT1))
 enddo 
 
end function average_Tensor1
! **********************************************************************
! *** function average_Tensor2                                       ***
! *** calculates the average second order tensor associated with a   ***
! *** weighted list of second order tensors                          ***
! ***                                                                ***
! *** input : nT2, number of elements                                ***
! ***         wT2, list of weights                                   ***
! ***         T2, list of second order tensors                       ***
! **********************************************************************
function average_Tensor2(nT2,wT2,T2)
 implicit none
 integer, intent(in) :: nT2
 type(Tensor2), intent(in) :: T2(nT2)
 real(8), intent(in) :: wT2(nT2)
 type(Tensor2) :: average_Tensor2
 integer :: iT2
 
 average_Tensor2=null2
 do iT2=1,nT2
   average_Tensor2=average_Tensor2+(wT2(iT2)*T2(iT2))
 enddo
 
end function average_Tensor2
! **********************************************************************
! *** function average_Tensor2s                                      ***
! *** calculates the average symmetric second order tensor associated***
! *** with a weighted list of symmetric second order tensors         ***
! ***                                                                ***
! *** input : nT2s, number of elements                               ***
! ***         wT2s, list of weights                                  ***
! ***         T2s, list of second order tensors                      ***
! **********************************************************************
function average_Tensor2s(nT2s,wT2s,T2s)
 implicit none
 integer, intent(in) :: nT2s
 type(Tensor2s), intent(in) :: T2s(nT2s)
 real(8), intent(in) :: wT2s(nT2s)
 type(Tensor2s) :: average_Tensor2s
 integer :: iT2s
 
 average_Tensor2s=null2s
 do iT2s=1,nT2s
   average_Tensor2s=average_Tensor2s+(wT2s(iT2s)*T2s(iT2s))
 enddo
 
end function average_Tensor2s
! **********************************************************************
! *** function average_Tensor2w                                      ***
! *** calculates the average skew-symmetric second order tensor      ***
! *** associated with a weighted list of skew-symmetric second order ***
! *** tensors                                                        ***
! ***                                                                ***
! *** input : nT2w, number of elements                               ***
! ***         wT2w, list of weights                                  ***
! ***         T2w, list of second order tensors                      ***
! **********************************************************************
function average_Tensor2w(nT2w,wT2w,T2w)
 implicit none
 integer, intent(in) :: nT2w
 type(Tensor2w), intent(in) :: T2w(nT2w)
 real(8), intent(in) :: wT2w(nT2w)
 type(Tensor2w) :: average_Tensor2w
 integer :: iT2w
 
 average_Tensor2w=null2w
 do iT2w=1,nT2w
   average_Tensor2w=average_Tensor2w+(wT2w(iT2w)*T2w(iT2w))
 enddo
 
end function average_Tensor2w
! **********************************************************************
! *** function average_Tensor3                                       ***
! *** calculates the average third order tensor associated with a    ***
! *** weighted list of third order tensors                           ***
! ***                                                                ***
! *** input : nT3, number of elements                                ***
! ***         wT3, list of weights                                   ***
! ***         T3, list of third order tensors                        ***
! **********************************************************************
function average_Tensor3(nT3,wT3,T3)
 implicit none
 integer, intent(in) :: nT3
 type(Tensor3), intent(in) :: T3(nT3)
 real(8), intent(in) :: wT3(nT3)
 type(Tensor3) :: average_Tensor3
 integer :: iT3
 
 average_Tensor3=null3
 do iT3=1,nT3
   average_Tensor3=average_Tensor3+(wT3(iT3)*T3(iT3))
 enddo
 
end function average_Tensor3
! **********************************************************************
! *** function average_Tensor4                                       ***
! *** calculates the average fourth order tensor associated with a   ***
! *** weighted list of fourth order tensors                          ***
! ***                                                                ***
! *** input : nT4, number of elements                                ***
! ***         wT4, list of weights                                   ***
! ***         T4, list of fourth order tensors                        ***
! **********************************************************************
function average_Tensor4(nT4,wT4,T4)
 implicit none
 integer, intent(in) :: nT4
 type(Tensor4), intent(in) :: T4(nT4)
 real(8), intent(in) :: wT4(nT4)
 type(Tensor4) :: average_Tensor4
 integer :: iT4
 
 average_Tensor4=null4
 do iT4=1,nT4
   average_Tensor4=average_Tensor4+(wT4(iT4)*T4(iT4))
 enddo
 
end function average_Tensor4
! **********************************************************************
! *** function average_Tensor4ss                                     ***
! *** calculates the average symmetric fourth order tensor associated***
! *** with a weighted list of symmetric fourth order tensors         ***
! ***                                                                ***
! *** input : nT4ss, number of elements                              ***
! ***         wT4ss, list of weights                                 ***
! ***         T4ss, list of fourth order tensors                     ***
! **********************************************************************
function average_Tensor4ss(nT4ss,wT4ss,T4ss)
 implicit none
 integer, intent(in) :: nT4ss
 type(Tensor4ss), intent(in) :: T4ss(nT4ss)
 real(8), intent(in) :: wT4ss(nT4ss)
 type(Tensor4ss) :: average_Tensor4ss
 integer :: iT4ss
 
 average_Tensor4ss=null4ss
 do iT4ss=1,nT4ss
   average_Tensor4ss=average_Tensor4ss+(wT4ss(iT4ss)*T4ss(iT4ss))
 enddo
 
end function average_Tensor4ss
! **********************************************************************
! *** subroutine equal_T1T1                                          ***
! *** assigns a first order tensor                                   ***
! ***                                                                ***
! *** input : T1B, first order tensor to be copied                   ***
! ***                                                                ***
! *** in/output :                                                    ***
! ***                                                                ***
! *** output : T1A, first order tensor to be assigned                ***
! **********************************************************************
subroutine equal_T1T1(T1A,T1B)
 implicit none
 type(tensor1), intent(out) :: T1A
 type(tensor1), intent(in) :: T1B

 T1A%val(:)=T1B%val(:)
 
end subroutine
! **********************************************************************
! *** subroutine equal_T2T2                                          ***
! *** assigns a second order tensor                                  ***
! ***                                                                ***
! *** input : T2B, second order tensor to be copied                  ***
! ***                                                                ***
! *** in/output :                                                    ***
! ***                                                                ***
! *** output : T2A, second order tensor to be assigned               ***
! **********************************************************************
subroutine equal_T2T2(T2A,T2B)
 implicit none
 type(tensor2), intent(out) :: T2A
 type(tensor2), intent(in) :: T2B

 T2A%val(:,:)=T2B%val(:,:)
 
end subroutine
! **********************************************************************
! *** subroutine equal_T2T2s                                         ***
! *** assigns a second order tensor                                  ***
! ***                                                                ***
! *** input : T2sB, symmetric second order tensor to be copied       ***
! ***                                                                ***
! *** in/output :                                                    ***
! ***                                                                ***
! *** output : T2A, second order tensor to be assigned               ***
! **********************************************************************
subroutine equal_T2T2s(T2A,T2sB)
 implicit none
 type(tensor2), intent(out) :: T2A
 type(tensor2s), intent(in) :: T2sB
 
 T2A%val(1,1)=T2sB%val(1)
 T2A%val(2,1)=T2sB%val(6)
 T2A%val(3,1)=T2sB%val(5)
 T2A%val(1,2)=T2sB%val(6)
 T2A%val(2,2)=T2sB%val(2)
 T2A%val(3,2)=T2sB%val(4)
 T2A%val(1,3)=T2sB%val(5)
 T2A%val(2,3)=T2sB%val(4)
 T2A%val(3,3)=T2sB%val(3)
 
end subroutine
! **********************************************************************
! *** subroutine equal_T2sT2                                         ***
! *** THIS SUBROUTINE SHOULD BE USED WITH CAUTION. IT WILL ONLY      ***
! *** ASSIGN THE SYMMETRIC PART OF THE TENSOR.                       ***
! *** assigns a second order tensor                                  ***
! ***                                                                ***
! *** input : T2B, second order tensor to be copied                  ***
! ***                                                                ***
! *** in/output :                                                    ***
! ***                                                                ***
! *** output : T2sA, symmetric second order tensor to be assigned    ***
! **********************************************************************
subroutine equal_T2sT2(T2sA,T2B)
 implicit none
 type(tensor2s), intent(out) :: T2sA
 type(tensor2), intent(in) :: T2B
 
 T2sA%val(1)=T2B%val(1,1)
 T2sA%val(2)=T2B%val(2,2)
 T2sA%val(3)=T2B%val(3,3)
 T2sA%val(4)=(T2B%val(3,2)+T2B%val(2,3))*0.5d0
 T2sA%val(5)=(T2B%val(3,1)+T2B%val(1,3))*0.5d0
 T2sA%val(6)=(T2B%val(2,1)+T2B%val(1,2))*0.5d0
 
end subroutine
! **********************************************************************
! *** subroutine equal_T2T2w                                         ***
! *** assigns a second order tensor                                  ***
! ***                                                                ***
! *** input : T2wB, skew-symmetric second order tensor to be copied  ***
! ***                                                                ***
! *** in/output :                                                    ***
! ***                                                                ***
! *** output : T2A, second order tensor to be assigned               ***
! **********************************************************************
subroutine equal_T2T2w(T2A,T2wB)
 implicit none
 type(tensor2), intent(out) :: T2A
 type(tensor2w), intent(in) :: T2wB
 
 T2A%val(1,1)=0.d0
 T2A%val(2,1)=-T2wB%val(3)
 T2A%val(3,1)=T2wB%val(2)
 T2A%val(1,2)=T2wB%val(3)
 T2A%val(2,2)=0.d0
 T2A%val(3,2)=-T2wB%val(1)
 T2A%val(1,3)=-T2wB%val(2)
 T2A%val(2,3)=T2wB%val(1)
 T2A%val(3,3)=0.0
 
end subroutine
! **********************************************************************
! *** subroutine equal_T2wT2                                         ***
! *** assigns a second order tensor                                  ***
! *** THIS SUBROUTINE SHOULD BE USED WITH CAUTION. IT WILL ONLY      ***
! *** ASSIGN THE SKEW-SYMMETRIC PART OF THE TENSOR.                  ***
! ***                                                                ***
! *** input : T2B, second order tensor to be copied                  ***
! ***                                                                ***
! *** in/output :                                                    ***
! ***                                                                ***
! *** output : T2wA, skew-symmetric second order tensor to be        ***
! ***          assigned                                              ***
! **********************************************************************
subroutine equal_T2wT2(T2wA,T2B)
 implicit none
 type(tensor2w), intent(out) :: T2wA
 type(tensor2), intent(in) :: T2B
 
 T2wA%val(1)=0.5d0*(T2B%val(2,3)-T2B%val(3,2))
 T2wA%val(2)=0.5d0*(T2B%val(3,1)-T2B%val(1,3))
 T2wA%val(3)=0.5d0*(T2B%val(1,2)-T2B%val(2,1))
 
end subroutine
! **********************************************************************
! *** subroutine equal_T2sT2s                                        ***
! *** assigns a symmetric second order tensor                        ***
! ***                                                                ***
! *** input : T2sB, second order tensor to be copied                 ***
! ***                                                                ***
! *** in/output :                                                    ***
! ***                                                                ***
! *** output : T2sA, second order tensor to be assigned              ***
! **********************************************************************
subroutine equal_T2sT2s(T2sA,T2sB)
 implicit none
 type(tensor2s), intent(out) :: T2sA
 type(tensor2s), intent(in) :: T2sB

 T2sA%val(:)=T2sB%val(:)
 
end subroutine
! **********************************************************************
! *** subroutine equal_T2w                                           ***
! *** assigns a skew-symmetric second order tensor                   ***
! ***                                                                ***
! *** input : T2wB, second order tensor to be copied                 ***
! ***                                                                ***
! *** in/output :                                                    ***
! ***                                                                ***
! *** output : T2wA, second order tensor to be assigned              ***
! **********************************************************************
subroutine equal_T2wT2w(T2wA,T2wB)
 implicit none
 type(tensor2w), intent(out) :: T2wA
 type(tensor2w), intent(in) :: T2wB

 T2wA%val(:)=T2wB%val(:)
 
end subroutine
! **********************************************************************
! *** subroutine equal_T3T3                                          ***
! *** assigns a third order tensor                                   ***
! ***                                                                ***
! *** input : T3B, third order tensor to be copied                   ***
! ***                                                                ***
! *** in/output :                                                    ***
! ***                                                                ***
! *** output : T3A, third order tensor to be assigned                ***
! **********************************************************************
subroutine equal_T3T3(T3A,T3B)
 implicit none
 type(tensor3), intent(out) :: T3A
 type(tensor3), intent(in) :: T3B

 T3A%val(:,:,:)=T3B%val(:,:,:)
 
end subroutine
! **********************************************************************
! *** subroutine equal_T4T4                                          ***
! *** assigns a fourth order tensor                                  ***
! ***                                                                ***
! *** input : T4B, fourth order tensor to be copied                  ***
! ***                                                                ***
! *** in/output :                                                    ***
! ***                                                                ***
! *** output : T4A, fourth order tensor to be assigned               ***
! **********************************************************************
subroutine equal_T4T4(T4A,T4B)
 implicit none
 type(tensor4), intent(out) :: T4A
 type(tensor4), intent(in) :: T4B

 T4A%val(:,:,:,:)=T4B%val(:,:,:,:)
 
end subroutine
! **********************************************************************
! *** subroutine equal_T4ssT4ss                                      ***
! *** assigns a symmetric fourth order tensor                        ***
! ***                                                                ***
! *** input : T4ssB, fourth order tensor to be copied                ***
! ***                                                                ***
! *** in/output :                                                    ***
! ***                                                                ***
! *** output : T4ssA, fourth order tensor to be assigned             ***
! **********************************************************************
subroutine equal_T4ssT4ss(T4ssA,T4ssB)
 implicit none
 type(tensor4ss), intent(out) :: T4ssA
 type(tensor4ss), intent(in) :: T4ssB

 T4ssA%val(:,:)=T4ssB%val(:,:)
 
end subroutine
! **********************************************************************
! *** subroutine equal_T4T4ss                                        ***
! *** assigns a symmetric fourth order tensor                        ***
! ***                                                                ***
! *** input : T4ssB, fourth order tensor to be copied                ***
! ***                                                                ***
! *** in/output :                                                    ***
! ***                                                                ***
! *** output : T4A, fourth order tensor to be assigned               ***
! **********************************************************************
subroutine equal_T4T4ss(T4A,T4ssB)
 type(tensor4), intent(out) :: T4A
 type(tensor4ss), intent(in) :: T4ssB
 integer :: i,j,k,l,ij,kl
 
 do l=1,3
  do k=1,3
   kl=indxs(k,l)
   do j=1,3
    do i=1,3
     ij=indxs(i,j)
     T4A%val(i,j,k,l)=T4ssB%val(ij,kl)
    enddo
   enddo
  enddo
 enddo

end subroutine
! **********************************************************************
! *** subroutine equal_T4ssT4                                        ***
! *** assigns a symmetric fourth order tensor                        ***
! ***                                                                ***
! *** input : T4B, fourth order tensor to be copied                  ***
! ***                                                                ***
! *** in/output :                                                    ***
! ***                                                                ***
! *** output : T4ssA, fourth order tensor to be assigned             ***
! **********************************************************************
subroutine equal_T4ssT4(T4ssA,T4B)
 type(tensor4ss), intent(out) :: T4ssA
 type(tensor4), intent(in) :: T4B
 integer :: i,j,k,l,ij,kl
 
 do ij=1,6
  i=ijndxs(ij,1)
  j=ijndxs(ij,2)
  do kl=1,6
   k=ijndxs(kl,1)
   l=ijndxs(kl,2)
   T4ssA%val(ij,kl)=T4B%val(i,j,k,l)
  enddo
 enddo

end subroutine
! **********************************************************************
! *** subroutine eigen2                                              ***
! *** computes all eigenvalues and eigenvectors of a general         ***
! *** second order tensor                                            ***
! ***                                                                ***
! *** input :  T2, second order tensor                               ***
! ***                                                                ***
! *** in/output :                                                    ***
! ***                                                                ***
! *** output : val, eigenvalues                                      ***
! ***          vec, eigenvectors                                     ***
! **********************************************************************
subroutine eigen2(T2,val,vec)
 implicit none
 type(Tensor2), intent(in) :: T2
 real(8), intent(out) :: val(3)
 type(Tensor1), intent(out) :: vec(3)
 real(8) :: wr(3),wi(3),vl(3,3),vr(3,3),a(3,3)
 integer :: info
 integer, parameter :: lwork=24
 real(8) :: work(lwork)
 integer :: i,j
 
 a(1,1)=T2%val(1,1)
 a(2,1)=T2%val(2,1)
 a(3,1)=T2%val(3,1)
 a(1,2)=T2%val(1,2)
 a(2,2)=T2%val(2,2)
 a(3,2)=T2%val(3,2)
 a(1,3)=T2%val(1,3)
 a(2,3)=T2%val(2,3)
 a(3,3)=T2%val(3,3)
 
 call dgeev('N','V',3,a,3,wr,wi,vl,3,vr,3,work,lwork,info)
 if (info.ne.0) then
  write(*,'("The algorithm failed to compute eigenvalues.")')
  stop
 else
  do i=1,3
   val(i)=wr(i)
   vec(i)%val(:)=vr(:,i)
  enddo
 endif
 
end subroutine eigen2
! **********************************************************************
! *** subroutine eigen2s                                             ***
! *** computes all eigenvalues and eigenvectors of a symmetric       ***
! *** second order tensor                                            ***
! ***                                                                ***
! *** input :  T2s, second order tensor                              ***
! ***                                                                ***
! *** in/output :                                                    ***
! ***                                                                ***
! *** output : val, eigenvalues                                      ***
! ***          vec, eigenvectors                                     ***
! **********************************************************************
subroutine eigen2s(T2s,val,vec)
 implicit none
 type(Tensor2s), intent(in) :: T2s
 real(8), intent(out) :: val(3)
 type(Tensor1), intent(out) :: vec(3)
 integer, parameter :: n=3, lwork=2*n*n+6*n+1,liwork=5*n+3
 integer :: iwork(liwork),info
 real(8) :: a(n,n),work(lwork),w(n),b(n,n)

 a(1,1)=T2s%val(1)
 a(2,2)=T2s%val(2)
 a(3,3)=T2s%val(3)
 a(2,3)=T2s%val(4)
 a(1,3)=T2s%val(5)
 a(1,2)=T2s%val(6)
 call dsyevd('V','U',n,a,n,val,work,lwork,iwork,liwork,info)
 if (info.eq.0) then
  vec(1)%val(:)=a(:,1)
  vec(2)%val(:)=a(:,2)
  vec(3)%val(:)=a(:,3)
 else
  write(*,'("CANNOT FIND EIGENVALUES AND EIGENVECTORS. PROGRAM WILL STOP.")')
  stop
 endif
end subroutine eigen2s
! subroutine eigen2s(T2s,val,vec)
!  implicit none
!  type(Tensor2s), intent(in) :: T2s
!  real(8), intent(out) :: val(3)
!  type(Tensor1), intent(out) :: vec(3)
!  integer, parameter :: nmax=500
!  integer :: i,ip,iq,j,nrot
!  real(8) c,g,h,s,sm,t,tau,theta,tresh,b(nmax),z(nmax),aa(3,3)
!
!  aa(1,1)=T2s%val(1)
!  aa(2,1)=T2s%val(6)
!  aa(3,1)=T2s%val(5)
!  aa(1,2)=T2s%val(6)
!  aa(2,2)=T2s%val(2)
!  aa(3,2)=T2s%val(4)
!  aa(1,3)=T2s%val(5)
!  aa(2,3)=T2s%val(4)
!  aa(3,3)=T2s%val(3)
!
!  do 12 ip=1,3
!   do 11 iq=1,3
!    vec(iq)%val(ip)=0.d+0
!   11 continue
!   vec(ip)%val(ip)=1.d+0
!  12 continue
!  do 13 ip=1,3
!   b(ip)=aa(ip,ip)
!   val(ip)=b(ip)
!   z(ip)=0.d+0
!  13 continue
!  nrot=0
!  do 24 i=1,50
!   sm=0.d+0
!   do 15 ip=1,3-1
!    do 14 iq=ip+1,3
!     sm=sm+abs(aa(ip,iq))
!    14 continue
!   15 continue
!   if(sm.eq.0.d+0)return
!   if(i.lt.4)then
!    tresh=0.2d+0*sm/3**2.d+0
!   else
!    tresh=0.d+0
!   endif
!   do 22 ip=1,3-1
!    do 21 iq=ip+1,3
!     g=100.d+0*abs(aa(ip,iq))
!     if((i.gt.4).and.(abs(val(ip))+g.eq.abs(val(ip))).and.(abs(val(iq))+g.eq.abs(val(iq))))then
!      aa(ip,iq)=0.d+0
!     else if(abs(aa(ip,iq)).gt.tresh)then
!      h=val(iq)-val(ip)
!      if(abs(h)+g.eq.abs(h))then
!       t=aa(ip,iq)/h
!      else
!       theta=0.5d+0*h/aa(ip,iq)
!       t=1.d+0/(abs(theta)+sqrt(1.d+0+theta**2.d+0))
!       if(theta.lt.0.d+0)t=-t
!      endif
!      c=1.d+0/sqrt(1.d+0+t**2.d+0)
!      s=t*c
!      tau=s/(1.d+0+c)
!      h=t*aa(ip,iq)
!      z(ip)=z(ip)-h
!      z(iq)=z(iq)+h
!      val(ip)=val(ip)-h
!      val(iq)=val(iq)+h
!      aa(ip,iq)=0.d+0
!      do 16 j=1,ip-1
!       g=aa(j,ip)
!       h=aa(j,iq)
!       aa(j,ip)=g-s*(h+g*tau)
!       aa(j,iq)=h+s*(g-h*tau)
!      16 continue
!      do 17 j=ip+1,iq-1
!       g=aa(ip,j)
!       h=aa(j,iq)
!       aa(ip,j)=g-s*(h+g*tau)
!       aa(j,iq)=h+s*(g-h*tau)
!      17 continue
!      do 18 j=iq+1,3
!       g=aa(ip,j)
!       h=aa(iq,j)
!       aa(ip,j)=g-s*(h+g*tau)
!       aa(iq,j)=h+s*(g-h*tau)
!      18 continue
!      do 19 j=1,3
!       g=vec(ip)%val(j)
!       h=vec(iq)%val(j)
!       vec(ip)%val(j)=g-s*(h+g*tau)
!       vec(iq)%val(j)=h+s*(g-h*tau)
!      19 continue
!      nrot=nrot+1
!     endif
!    21 continue
!   22 continue
!   do 23 ip=1,3
!    b(ip)=b(ip)+z(ip)
!    val(ip)=b(ip)
!    z(ip)=0.d+0
!   23 continue
!  24 continue
!  write(*,'("TOO MANY ITERATIONS IN EIGEN2S. PROGRAM WILL STOP")')
!  stop
!  return
! end subroutine eigen2s
! **********************************************************************
! *** function exp_Tensor2                                           ***
! *** computes the exponential of a general second order tensor      ***
! ***                                                                ***
! *** input : T2sA, second order tensor                              ***
! **********************************************************************
function exp_Tensor2(T2A)
 implicit none
 type(Tensor2), intent(in) :: T2A
 type(Tensor2) :: exp_Tensor2
 real(8) :: wr(3),wi(3),vl(3,3),vr(3,3),a(3,3)
 integer :: info
 integer, parameter :: lwork=24
 real(8) :: work(lwork)
 integer :: i,j,k
 
 a(1,1)=T2A%val(1,1)
 a(2,1)=T2A%val(2,1)
 a(3,1)=T2A%val(3,1)
 a(1,2)=T2A%val(1,2)
 a(2,2)=T2A%val(2,2)
 a(3,2)=T2A%val(3,2)
 a(1,3)=T2A%val(1,3)
 a(2,3)=T2A%val(2,3)
 a(3,3)=T2A%val(3,3)
 
 call dgeev('N','V',3,a,3,wr,wi,vl,3,vr,3,work,lwork,info)
!  call dgeev('N','V',3,a,3,wr wi,vl,3,vr,3,work,lwork,info)
 if (info.ne.0) then
  write(*,'("The algorithm failed to compute eigenvalues.")')
  stop
 else
!   write(*,'("EIGEN VECTORS",3es12.4)') vr
!   write(*,'("EIGEN VALUES",3es12.4)') wr
!   vl=invmat(vr,3,3)
  vl=inverse('G',vr,3)
!   write(*,'("EIGEN VECTORS",3es12.4)') vr
  do i=1,3
   do j=1,3
    exp_Tensor2%val(i,j)=0.d0
    do k=1,3
     exp_Tensor2%val(i,j)=exp_Tensor2%val(i,j)+vr(i,k)*exp(wr(k))*vl(k,j)
    enddo
   enddo
  enddo
 endif 
             
end function exp_Tensor2
! **********************************************************************
! *** function exp_Tensor2s                                          ***
! *** computes the exponential of a symmetric second order tensor    ***
! ***                                                                ***
! *** input : T2sA, second order tensor                              ***
! **********************************************************************
function exp_Tensor2s(T2sA)
 implicit none
 type(Tensor2s), intent(in) :: T2sA
 type(Tensor1) :: evec(3)
 real(8) :: eval(3)
 type(Tensor2s) :: exp_Tensor2s
 
 call eigen2s(T2sA,eval,evec)
 exp_Tensor2s%val(1)=exp(eval(1))*evec(1)%val(1)*evec(1)%val(1) &
                    +exp(eval(2))*evec(2)%val(1)*evec(2)%val(1) &
                    +exp(eval(3))*evec(3)%val(1)*evec(3)%val(1)
 exp_Tensor2s%val(2)=exp(eval(1))*evec(1)%val(2)*evec(1)%val(2) &
                    +exp(eval(2))*evec(2)%val(2)*evec(2)%val(2) &
                    +exp(eval(3))*evec(3)%val(2)*evec(3)%val(2)
 exp_Tensor2s%val(3)=exp(eval(1))*evec(1)%val(3)*evec(1)%val(3) &
                    +exp(eval(2))*evec(2)%val(3)*evec(2)%val(3) &
                    +exp(eval(3))*evec(3)%val(3)*evec(3)%val(3)
 exp_Tensor2s%val(4)=exp(eval(1))*evec(1)%val(2)*evec(1)%val(3) &
                    +exp(eval(2))*evec(2)%val(2)*evec(2)%val(3) &
                    +exp(eval(3))*evec(3)%val(2)*evec(3)%val(3)
 exp_Tensor2s%val(5)=exp(eval(1))*evec(1)%val(3)*evec(1)%val(1) &
                    +exp(eval(2))*evec(2)%val(3)*evec(2)%val(1) &
                    +exp(eval(3))*evec(3)%val(3)*evec(3)%val(1)
 exp_Tensor2s%val(6)=exp(eval(1))*evec(1)%val(1)*evec(1)%val(2) &
                    +exp(eval(2))*evec(2)%val(1)*evec(2)%val(2) &
                    +exp(eval(3))*evec(3)%val(1)*evec(3)%val(2)
             
end function exp_Tensor2s
! **********************************************************************
! *** function exp_Tensor2w                                          ***
! *** computes the exponential of a skew-symmetric second order      ***
! *** tensor                                                         ***
! ***                                                                ***
! *** input : T2wA, second order tensor                              ***
! **********************************************************************
function exp_Tensor2w(T2wA)
 implicit none
 type(Tensor2w), intent(in) :: T2wA
 type(Tensor2w) :: k
 type(Tensor2) :: exp_Tensor2w,kk
 real(8) :: th
 real(8), parameter :: small=1.d-10

 th=nrm(T2wA)/sqrt(2.d0)
 if (abs(th).gt.small) then
  k=T2wA/th
  kk=k.dot.k
  exp_Tensor2w=id2+sin(th)*k+(1.d0-cos(th))*kk
 else
  exp_Tensor2w=id2
 endif

end function
! **********************************************************************
! *** function log_Tensor2                                           ***
! *** computes the natural logarithm of a general second order tensor***
! ***                                                                ***
! *** input : T2A, second order tensor                               ***
! **********************************************************************
function log_Tensor2(T2A)
 implicit none
 type(Tensor2), intent(in) :: T2A
 type(Tensor2) :: log_Tensor2
 real(8) :: wr(3),wi(3),vl(3,3),vr(3,3),a(3,3)
 integer :: info
 integer, parameter :: lwork=24
 real(8) :: work(lwork)
 integer :: i,j,k
 
 a(1,1)=T2A%val(1,1)
 a(2,1)=T2A%val(2,1)
 a(3,1)=T2A%val(3,1)
 a(1,2)=T2A%val(1,2)
 a(2,2)=T2A%val(2,2)
 a(3,2)=T2A%val(3,2)
 a(1,3)=T2A%val(1,3)
 a(2,3)=T2A%val(2,3)
 a(3,3)=T2A%val(3,3)
 
 call dgeev('N','V',3,a,3,wr,wi,vl,3,vr,3,work,lwork,info)
 if (info.ne.0) then
  write(*,'("The algorithm failed to compute eigenvalues.")')
  stop
 else
  do k=1,3
   if (wr(k).lt.0.d0) then
    write(*,'("Negative eigenvalue. Tensor log cannot be computed.")')
   endif
  enddo
  vl=inverse('G',vr,3)
  do i=1,3
   do j=1,3
    log_Tensor2%val(i,j)=0.d0
    do k=1,3
     log_Tensor2%val(i,j)=log_Tensor2%val(i,j)+vr(i,k)*log(wr(k))*vl(k,j)
    enddo
   enddo
  enddo
 endif 
             
end function log_Tensor2
! **********************************************************************
! *** function log_Tensor2s                                          ***
! *** computes the natural logarithm of a symmetric second order     ***
! *** tensor                                                         ***
! ***                                                                ***
! *** input : T2sA, second order tensor                              ***
! **********************************************************************
function log_Tensor2s(T2sA)
 implicit none
 type(Tensor2s), intent(in) :: T2sA
 type(Tensor1) :: evec(3)
 real(8) :: eval(3)
 type(Tensor2s) :: log_Tensor2s
 
 call eigen2s(T2sA,eval,evec)
 if (eval(1).le.0.d0.or.eval(2).le.0.d0.or.eval(3).le.0.d0) then
  write(*,'("NON POSITIVE TENSOR IN LOG_TENSOR2S. PROGRAM WILL STOP")') 
  stop
 endif
 log_Tensor2s%val(1)=log(eval(1))*evec(1)%val(1)*evec(1)%val(1) &
                    +log(eval(2))*evec(2)%val(1)*evec(2)%val(1) &
                    +log(eval(3))*evec(3)%val(1)*evec(3)%val(1)
 log_Tensor2s%val(2)=log(eval(1))*evec(1)%val(2)*evec(1)%val(2) &
                    +log(eval(2))*evec(2)%val(2)*evec(2)%val(2) &
                    +log(eval(3))*evec(3)%val(2)*evec(3)%val(2)
 log_Tensor2s%val(3)=log(eval(1))*evec(1)%val(3)*evec(1)%val(3) &
                    +log(eval(2))*evec(2)%val(3)*evec(2)%val(3) &
                    +log(eval(3))*evec(3)%val(3)*evec(3)%val(3)
 log_Tensor2s%val(4)=log(eval(1))*evec(1)%val(2)*evec(1)%val(3) &
                    +log(eval(2))*evec(2)%val(2)*evec(2)%val(3) &
                    +log(eval(3))*evec(3)%val(2)*evec(3)%val(3)
 log_Tensor2s%val(5)=log(eval(1))*evec(1)%val(3)*evec(1)%val(1) &
                    +log(eval(2))*evec(2)%val(3)*evec(2)%val(1) &
                    +log(eval(3))*evec(3)%val(3)*evec(3)%val(1)
 log_Tensor2s%val(6)=log(eval(1))*evec(1)%val(1)*evec(1)%val(2) &
                    +log(eval(2))*evec(2)%val(1)*evec(2)%val(2) &
                    +log(eval(3))*evec(3)%val(1)*evec(3)%val(2)

end function log_Tensor2s
! **********************************************************************
! *** function sqrt_Tensor2s                                         ***
! *** computes the square root of a symmetric second order tensor    ***
! ***                                                                ***
! *** input : T2sA, second order tensor                              ***
! **********************************************************************
function sqrt_Tensor2s(T2sA)
 implicit none
 type(Tensor2s), intent(in) :: T2sA
 type(Tensor1) :: evec(3)
 real(8) :: eval(3)
 type(Tensor2s) :: sqrt_Tensor2s
 
 call eigen2s(T2sA,eval,evec)
 if (eval(1).lt.0.d0.or.eval(2).lt.0.d0.or.eval(3).lt.0.d0) then
  write(*,'("NON POSITIVE TENSOR IN SQRT_TENSOR2S. PROGRAM WILL STOP")') 
  stop
 endif
 sqrt_Tensor2s%val(1)=sqrt(eval(1))*evec(1)%val(1)*evec(1)%val(1) &
                     +sqrt(eval(2))*evec(2)%val(1)*evec(2)%val(1) &
                     +sqrt(eval(3))*evec(3)%val(1)*evec(3)%val(1)
 sqrt_Tensor2s%val(2)=sqrt(eval(1))*evec(1)%val(2)*evec(1)%val(2) &
                     +sqrt(eval(2))*evec(2)%val(2)*evec(2)%val(2) &
                     +sqrt(eval(3))*evec(3)%val(2)*evec(3)%val(2)
 sqrt_Tensor2s%val(3)=sqrt(eval(1))*evec(1)%val(3)*evec(1)%val(3) &
                     +sqrt(eval(2))*evec(2)%val(3)*evec(2)%val(3) &
                     +sqrt(eval(3))*evec(3)%val(3)*evec(3)%val(3)
 sqrt_Tensor2s%val(4)=sqrt(eval(1))*evec(1)%val(2)*evec(1)%val(3) &
                     +sqrt(eval(2))*evec(2)%val(2)*evec(2)%val(3) &
                     +sqrt(eval(3))*evec(3)%val(2)*evec(3)%val(3)
 sqrt_Tensor2s%val(5)=sqrt(eval(1))*evec(1)%val(3)*evec(1)%val(1) &
                     +sqrt(eval(2))*evec(2)%val(3)*evec(2)%val(1) &
                     +sqrt(eval(3))*evec(3)%val(3)*evec(3)%val(1)
 sqrt_Tensor2s%val(6)=sqrt(eval(1))*evec(1)%val(1)*evec(1)%val(2) &
                     +sqrt(eval(2))*evec(2)%val(1)*evec(2)%val(2) &
                     +sqrt(eval(3))*evec(3)%val(1)*evec(3)%val(2)

end function sqrt_Tensor2s
! **********************************************************************
! *** function sqrt_Tensor4ss                                        ***
! *** computes the square root of a symmetric fourth order tensor    ***
! ***                                                                ***
! *** input : T4ssA, second order tensor                              ***
! **********************************************************************
! function sqrt_Tensor4ss(T4ssA)
!  implicit none
!  type(Tensor4ss), intent(in) :: T4ssA
!  type(Tensor4ss) :: y,z,a,b,q,sqrt_Tensor4ss
!  real(8), parameter :: toler=1.e-8
!  integer :: iter
!  integer, parameter :: niter=1000
!  real(8) :: error
!
!  y=T4ssA
!  z=id4ss
!  error=toler+1.d0
!  iter=0
!  do while (error.gt.toler)
!   sqrt_Tensor4ss=y
!   iter=iter+1
!   a=inv(y)
!   b=inv(z)
!   y=0.5d0*(y+b)
!   z=0.5d0*(z+a)
!   error=nrm(y-sqrt_Tensor4ss)/nrm(y)
!   if (iter.gt.niter) then
!    write(*,'("ERROR IN sqrt_Tensor4ss. PROGRAM WILL STOP")')
!    stop
!   endif
!  enddo
!
! end function sqrt_Tensor4ss
function sqrt_Tensor4ss(T4ssA)
 implicit none
 type(Tensor4ss), intent(in) :: T4ssA
 type(Tensor4ss) :: sqrt_Tensor4ss
 real(8) :: wr(6),wi(6),vl(6,6),vr(6,6),a(6,6)
 integer :: info
 integer, parameter :: lwork=24
 real(8) :: work(lwork)
 integer :: i,j,k
 real(8), parameter :: sqr=1.41421356237d0
! *** switch to Mandel notation
 a(1,1)=T4ssA%val(1,1)
 a(2,1)=T4ssA%val(2,1)
 a(3,1)=T4ssA%val(3,1)
 a(4,1)=T4ssA%val(4,1)*sqr
 a(5,1)=T4ssA%val(5,1)*sqr
 a(6,1)=T4ssA%val(6,1)*sqr
 a(1,2)=T4ssA%val(1,2)
 a(2,2)=T4ssA%val(2,2)
 a(3,2)=T4ssA%val(3,2)
 a(4,2)=T4ssA%val(4,2)*sqr
 a(5,2)=T4ssA%val(5,2)*sqr
 a(6,2)=T4ssA%val(6,2)*sqr
 a(1,3)=T4ssA%val(1,3)
 a(2,3)=T4ssA%val(2,3)
 a(3,3)=T4ssA%val(3,3)
 a(4,3)=T4ssA%val(4,3)*sqr
 a(5,3)=T4ssA%val(5,3)*sqr
 a(6,3)=T4ssA%val(6,3)*sqr
 a(1,4)=T4ssA%val(1,4)*sqr
 a(2,4)=T4ssA%val(2,4)*sqr
 a(3,4)=T4ssA%val(3,4)*sqr
 a(4,4)=T4ssA%val(4,4)*2.d0
 a(5,4)=T4ssA%val(5,4)*2.d0
 a(6,4)=T4ssA%val(6,4)*2.d0
 a(1,5)=T4ssA%val(1,5)*sqr
 a(2,5)=T4ssA%val(2,5)*sqr
 a(3,5)=T4ssA%val(3,5)*sqr
 a(4,5)=T4ssA%val(4,5)*2.d0
 a(5,5)=T4ssA%val(5,5)*2.d0
 a(6,5)=T4ssA%val(6,5)*2.d0
 a(1,6)=T4ssA%val(1,6)*sqr
 a(2,6)=T4ssA%val(2,6)*sqr
 a(3,6)=T4ssA%val(3,6)*sqr
 a(4,6)=T4ssA%val(4,6)*2.d0
 a(5,6)=T4ssA%val(5,6)*2.d0
 a(6,6)=T4ssA%val(6,6)*2.d0
! *** calculate eigenvalues and eigenvectors
 call dgeev('N','V',6,a,6,wr,wi,vl,6,vr,6,work,lwork,info)
 if (info.ne.0) then
  write(*,'("The algorithm failed to compute eigenvalues.")')
  stop
 else
  vl=inverse('G',vr,6)
  do i=1,6
   do j=1,6
    a(i,j)=0.d0
    do k=1,6
     a(i,j)=a(i,j)+vr(i,k)*sqrt(wr(k))*vl(k,j)
    enddo
   enddo
  enddo
 endif
! *** switch to standard notation
 sqrt_Tensor4ss%val(1,1)=a(1,1)
 sqrt_Tensor4ss%val(2,1)=a(2,1)
 sqrt_Tensor4ss%val(3,1)=a(3,1)
 sqrt_Tensor4ss%val(4,1)=a(4,1)/sqr
 sqrt_Tensor4ss%val(5,1)=a(5,1)/sqr
 sqrt_Tensor4ss%val(6,1)=a(6,1)/sqr
 sqrt_Tensor4ss%val(1,2)=a(1,2)
 sqrt_Tensor4ss%val(2,2)=a(2,2)
 sqrt_Tensor4ss%val(3,2)=a(3,2)
 sqrt_Tensor4ss%val(4,2)=a(4,2)/sqr
 sqrt_Tensor4ss%val(5,2)=a(5,2)/sqr
 sqrt_Tensor4ss%val(6,2)=a(6,2)/sqr
 sqrt_Tensor4ss%val(1,3)=a(1,3)
 sqrt_Tensor4ss%val(2,3)=a(2,3)
 sqrt_Tensor4ss%val(3,3)=a(3,3)
 sqrt_Tensor4ss%val(4,3)=a(4,3)/sqr
 sqrt_Tensor4ss%val(5,3)=a(5,3)/sqr
 sqrt_Tensor4ss%val(6,3)=a(6,3)/sqr
 sqrt_Tensor4ss%val(1,4)=a(1,4)/sqr
 sqrt_Tensor4ss%val(2,4)=a(2,4)/sqr
 sqrt_Tensor4ss%val(3,4)=a(3,4)/sqr
 sqrt_Tensor4ss%val(4,4)=a(4,4)/2.d0
 sqrt_Tensor4ss%val(5,4)=a(5,4)/2.d0
 sqrt_Tensor4ss%val(6,4)=a(6,4)/2.d0
 sqrt_Tensor4ss%val(1,5)=a(1,5)/sqr
 sqrt_Tensor4ss%val(2,5)=a(2,5)/sqr
 sqrt_Tensor4ss%val(3,5)=a(3,5)/sqr
 sqrt_Tensor4ss%val(4,5)=a(4,5)/2.d0
 sqrt_Tensor4ss%val(5,5)=a(5,5)/2.d0
 sqrt_Tensor4ss%val(6,5)=a(6,5)/2.d0
 sqrt_Tensor4ss%val(1,6)=a(1,6)/sqr
 sqrt_Tensor4ss%val(2,6)=a(2,6)/sqr
 sqrt_Tensor4ss%val(3,6)=a(3,6)/sqr
 sqrt_Tensor4ss%val(4,6)=a(4,6)/2.d0
 sqrt_Tensor4ss%val(5,6)=a(5,6)/2.d0
 sqrt_Tensor4ss%val(6,6)=a(6,6)/2.d0
end function sqrt_Tensor4ss
! **********************************************************************
! *** function pos_Tensor2s                                          ***
! *** computes the positive part of a symmetric second order tensor  ***
! ***                                                                ***
! *** input : T2sA, second order tensor                              ***
! **********************************************************************
function pos_Tensor2s(T2sA)
 implicit none
 type(Tensor2s), intent(in) :: T2sA
 type(Tensor1) :: evec(3)
 real(8) :: eval(3)
 type(Tensor2s) :: pos_Tensor2s
 
 call eigen2s(T2sA,eval,evec)

 if (eval(1).lt.0.d0) then
  eval(1)=0.d0
 endif
 if (eval(2).lt.0.d0) then
  eval(2)=0.d0
 endif
 if (eval(3).lt.0.d0) then
  eval(3)=0.d0
 endif
 
 pos_Tensor2s%val(1)=eval(1)*evec(1)%val(1)*evec(1)%val(1) &
                    +eval(2)*evec(2)%val(1)*evec(2)%val(1) &
                    +eval(3)*evec(3)%val(1)*evec(3)%val(1)
 pos_Tensor2s%val(2)=eval(1)*evec(1)%val(2)*evec(1)%val(2) &
                    +eval(2)*evec(2)%val(2)*evec(2)%val(2) &
                    +eval(3)*evec(3)%val(2)*evec(3)%val(2)
 pos_Tensor2s%val(3)=eval(1)*evec(1)%val(3)*evec(1)%val(3) &
                    +eval(2)*evec(2)%val(3)*evec(2)%val(3) &
                    +eval(3)*evec(3)%val(3)*evec(3)%val(3)
 pos_Tensor2s%val(4)=eval(1)*evec(1)%val(2)*evec(1)%val(3) &
                    +eval(2)*evec(2)%val(2)*evec(2)%val(3) &
                    +eval(3)*evec(3)%val(2)*evec(3)%val(3)
 pos_Tensor2s%val(5)=eval(1)*evec(1)%val(3)*evec(1)%val(1) &
                    +eval(2)*evec(2)%val(3)*evec(2)%val(1) &
                    +eval(3)*evec(3)%val(3)*evec(3)%val(1)
 pos_Tensor2s%val(6)=eval(1)*evec(1)%val(1)*evec(1)%val(2) &
                    +eval(2)*evec(2)%val(1)*evec(2)%val(2) &
                    +eval(3)*evec(3)%val(1)*evec(3)%val(2)

end function pos_Tensor2s
! **********************************************************************
! *** subroutine plrdcmp                                             ***
! *** decomposes a second order tensor according to the polar        ***
! *** decomposition                                                  ***
! ***                                                                ***
! *** input : f, second order tensor to be decomposed                ***
! ***                                                                ***
! *** in/output :                                                    ***
! ***                                                                ***
! *** output : r, rotation tensor                                    ***
! ***          u, right stretch tensor                               ***
! ***          v, left stretch tensor                                ***
! **********************************************************************
subroutine plrdcmp(f,r,u,v)
 implicit none
 type(Tensor2), intent(in) :: f
 type(Tensor2), intent(out) :: r
 type(Tensor2s), intent(out) :: u,v
 
 u=sqrt_Tensor2s(sym_Tensor2(dot_T2T2(transpose_Tensor2(f),f)))
 r=dot_T2T2s(f,inv_Tensor2s(u))
 v=sym_Tensor2(dot_T2T2(f,transpose_Tensor2(r)))
 
end subroutine plrdcmp
! **********************************************************************
! *** subroutine isotropic                                           ***
! ***                                                                ***
! *** input : func, string to determine what isotropic function to be***
! ***               used                                             ***
! ***         TA, symmetric second order tensor (mandatory)          ***
! ***         TB, symmetric second order tensor (optional)           ***
! ***                                                                ***
! *** in/output :                                                    ***
! ***                                                                ***
! *** output : TC=func(TA), symmetric second order tensor (optional) ***
! ***          dTCdTA=dfunc/dTA, symmetric fourth order tensor       ***
! ***                (optional)                                      ***
! ***          TBddTCddTA=TB:d²func/dTA², symmetric fourth order     ***
! ***                     tensor (optional)                          ***
! **********************************************************************
subroutine isotropic(func,TA,TB,TC,dTCdTA,TBddTCddTA)
 implicit none
 character(3), intent(in) :: func
 type(Tensor1) :: ncs(3)
 type(Tensor2s), intent(in) :: TA
 type(Tensor2s), intent(in), optional :: TB
 type(Tensor2s), intent(out), optional :: TC
 type(Tensor4ss), intent(out), optional :: dTCdTA,TBddTCddTA
 type(Tensor2s) :: mcs(3)
 type(Tensor2) :: bcs(3,3)
 type(Tensor4ss) :: gcs(3,3),hcs(3,3,3)
 real(8) :: t(3,3),d(3),s(3),z(3),x(3,3),f(3),e
 real(8), parameter :: small=1.d-4
 integer :: a,b,c,i,j,k,l,m,n,ij,kl,mn,ik,il,im,in,jk,jl,jm,jn,km,kn,lm,ln
 real(8) :: realTBddTCddTA(6,6),reallcs(6,6,6),realhcs(3,3,3,6,6,6)

! *** spectral decomposition of a second rank tensor
 call eigen(TA,s,ncs)
 do a=1,3
  mcs(a)=ncs(a)*ncs(a)
  select case(func)
  case("exp")
   z(a)=exp(s(a))
   d(a)=z(a)
   f(a)=z(a)
  case("log")
   z(a)=log(s(a))
   d(a)=1.d0/s(a)
   f(a)=-1.d0/(s(a)*s(a))
  case("sqr")
   z(a)=sqrt(s(a))
   d(a)=0.5d0/z(a)
   f(a)=-0.25d0/(z(a)*s(a))
  case("pos")
   if (s(a).ge.0.d0) then
    z(a)=s(a)
    d(a)=1.d0
    f(a)=0.d0
   else
    z(a)=0.d0
    d(a)=0.d0
    f(a)=0.d0
   endif
  case("neg")
   if (s(a).lt.0.d0) then
    z(a)=s(a)
    d(a)=1.d0
    f(a)=0.d0
   else
    z(a)=0.d0
    d(a)=0.d0
    f(a)=0.d0
   endif
  case default
   write(*,'("Unknown function in subroutine ISOTROPIC.")')
   stop
  end select
 enddo
 if (present(dTCdTA).or.(present(TB).and.present(TBddTCddTA))) then
  do a=1,3
   do b=1,3
    if (abs(s(a)-s(b)).gt.small) then
     t(a,b)=(z(a)-z(b))/(s(a)-s(b))
    else
     t(a,b)=0.5d0*(d(a)+d(b))
    endif
   enddo
  enddo
 endif
 if (present(TB).and.present(TBddTCddTA)) then
  do a=1,3
   do b=1,3
    if (abs(s(a)-s(b)).gt.small) then
     x(a,b)=0.25d0*(t(a,b)-d(b))/(s(a)-s(b))
    else
     x(a,b)=0.0625d0*(f(a)+f(b))
    endif
   enddo
  enddo
  if (abs(s(1)-s(2)).gt.small) then
   if (abs(s(2)-s(3)).gt.small) then
    e=0.25d0*z(1)/(s(1)-s(2))/(s(1)-s(3))+0.25d0*z(2)/(s(2)-s(3))/(s(2)-s(1))+0.25d0*z(3)/(s(3)-s(1))/(s(3)-s(2))  ! *** s1 is not s2 and s2 is not s3
   else
    e=0.5d0*(x(1,2)+x(1,3)) ! *** s1 is not s2 but s2 is s3
   endif
  else
   if (abs(s(2)-s(3)).gt.small) then
    e=0.5d0*(x(3,2)+x(3,1)) ! *** s1 is s2 but s2 is not s3
   else
    e=(f(1)+f(2)+f(3))/24.d0 ! *** s1 is s2 and s2 is s3
   endif
  endif
 endif
! *** isotropic tensor function
 if (present(TC)) then
  TC=null2s
  do a=1,3
   TC=TC+(z(a)*mcs(a))
  enddo
 endif
! *** first order derivative isotropic tensor function wrt tensor
 if (present(dTCdTA)) then
  do a=1,3
   do b=(a+1),3
    do ij=1,6
     i=ijndxs(ij,1)
     j=ijndxs(ij,2)
     do kl=1,6
      k=ijndxs(kl,1)
      l=ijndxs(kl,2)
      gcs(a,b)%val(ij,kl)=0.25d0*(ncs(a)%val(i)*ncs(b)%val(j)+ncs(b)%val(i)*ncs(a)%val(j))*(ncs(a)%val(k)*ncs(b)%val(l)+ncs(b)%val(k)*ncs(a)%val(l))
     enddo
    enddo
    gcs(b,a)=gcs(a,b)
   enddo
  enddo
  dTCdTA=null4ss
  do a=1,3
   dTCdTA=dTCdTA+(d(a)*(mcs(a)*mcs(a)))
  enddo
  do a=1,3
   do b=(a+1),3
    if (a.ne.b) then
     dTCdTA=dTCdTA+((2.d0*t(a,b))*gcs(a,b))
    endif
   enddo
  enddo
 endif
! *** second order derivative of isotropic tensor function wrt tensor
 if (present(TB).and.present(TBddTCddTA)) then
  do a=1,3
   do b=(a+1),3
    bcs(a,b)=(mcs(a).dot.(TB.dot.mcs(b)))
    bcs(b,a)=trsp(bcs(a,b))
   enddo
   bcs(a,a)=(mcs(a).dot.(TB.dot.mcs(a)))
  enddo
  do a=1,3
   do b=1,3
    do c=1,3
     do kl=1,6
      k=ijndxs(kl,1)
      l=ijndxs(kl,2)
      do mn=1,6
       m=ijndxs(mn,1)
       n=ijndxs(mn,2)
       km=indxs(k,m)
       kn=indxs(k,n)
       lm=indxs(l,m)
       ln=indxs(l,n)
       hcs(a,b,c)%val(kl,mn)=(bcs(a,b)%val(k,m)*mcs(c)%val(ln)+bcs(a,b)%val(k,n)*mcs(c)%val(lm)+bcs(a,b)%val(l,m)*mcs(c)%val(kn)+bcs(a,b)%val(l,n)*mcs(c)%val(km))*2.d0
      enddo
     enddo
    enddo
   enddo
  enddo
  TBddTCddTA=null4ss
  do a=1,3
   TBddTCddTA=TBddTCddTA+(f(a)*(sym(bcs(a,a))*mcs(a)))
  enddo
  do a=1,3
   do b=1,3
    if (a.ne.b) then
     TBddTCddTA=TBddTCddTA+(x(a,b)*(hcs(a,b,b)+hcs(b,a,b)+hcs(b,b,a)))
    endif
   enddo
  enddo
  do a=1,3
   do b=1,3
    do c=1,3
     if (a.ne.b.and.b.ne.c.and.c.ne.a) then
      TBddTCddTA=TBddTCddTA+(e*hcs(a,b,c))
     endif
    enddo
   enddo
  enddo
 endif

end subroutine
end module Tensor

