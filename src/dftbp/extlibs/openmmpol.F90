!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!
#:include "common.fypp"
#:include "error.fypp"

!> Proxy module for interfacing with the openmmpol library.
module dftbp_extlibs_openmmpol
    use dftbp_common_accuracy, only : dp, mc
    use dftbp_dftbplus_qdepextpotgen, only : TQDepExtPotGen
    use dftbp_io_message, only : error, warning
  #:if WITH_OPENMMPOL
    use ommp_interface, only : ommp_system, ommp_init_mmp, ommp_init_xyz, ommp_terminate,&
    ommp_print_summary, ommp_prepare_qm_ele_ene, ommp_set_external_field,&
    ommp_potential_mmpol2ext, ommp_qm_helper, ommp_init_qm_helper, ommp_terminate_qm_helper,&
    ommp_set_verbose, ommp_get_full_ele_energy, ommp_get_fixedelec_energy,&
    ommp_get_polelec_energy, ommp_get_full_bnd_energy, ommp_get_vdw_energy,&
    ommp_qm_helper_init_vdw_prm, ommp_qm_helper_vdw_energy, ommp_qm_helper_set_attype,&
    ommp_qm_helper_vdw_energy, ommp_print_summary_to_file, ommp_ignore_duplicated_angle_prm,&
    ommp_ignore_duplicated_opb_prm, OMMP_SOLVER_NONE, OMMP_FF_AMOEBA, OMMP_FF_WANG_AL,&
    OMMP_FF_WANG_DL, OMMP_VERBOSE_DEBUG, OMMP_VERBOSE_LOW, OMMP_VERBOSE_HIGH
  #:endif

    implicit none
    
    private
    public TOpenmmpolPotGen, TOpenmmpolPotGen_init, TOpenmmpolInput

    ! Constant for AMBER identification inside openmmpol;
    ! TODO: import from openmmpol
    integer, protected :: OMMP_FF_AMBER = 0
  
    !> Interface to control polarizable QM/MM
    !! through openmmpol external library
    type, extends(TQDepExtPotGen) :: TOpenmmpolPotGen
      private
    #:if WITH_OPENMMPOL
      !> Pointer to openmmpol system object
      type(ommp_system), pointer :: pSystem
      !> Pointer to openmmpol QM helper object
      type(ommp_qm_helper), pointer :: pQMHelper
    #:endif
      !> Potential to calculate total energy
      real(dp), allocatable :: pot(:)
      !> Effective potential to calculate Fock matrix elements;
      !! normally is equal to the ordinary potential,
      !! but in cases like AMOEBA it includes extra contributions
      !! from screened dipole-dipole interaction, which happen
      !! do not contribute to QM/MM interaction energy.
      real(dp), allocatable :: potFock(:)
      !> Total energy of all force field terms that
      !! do not depend on the density matrix
      real(dp) :: Eff
      
      contains
        !> method for initializing a potential generator object
        procedure :: init => TOpenmmpolPotGen_init
        !> method for evaluating external potential
        procedure :: getExternalPot => TOpenmmpolPotGen_getExternalPot
        !> method for evaluating screened Fock potential
        procedure :: getExternalPotFock => TOpenmmpolPotGen_getExternalPotFock
        !> method for evaluating gradient of external potential
        procedure :: getExternalPotGrad => TOpenmmpolPotGen_getExternalPotGrad
        !> method for adding internal energy contribution
        procedure :: getInternalEnergy => TOpenmmpolPotGen_getInternalEnergy
        !> method for writing openmmpol-specific data
        procedure :: writeOutput => TOpenmmpolPotGen_writeOutput
        ! TODO: add methods for controlling dynamics and excited states contributions
        !> method for updating coordinates in the QM zone
        ! procedure :: updateQMCoords => TOpenmmpolPotGen_updateQMCoords
    end type

    !> Data type for storing openmmpol-specific input parameters
    type :: TOpenmmpolInput
      !> Used input format ("Tinker" or "mmp")
      character(:), allocatable :: inputFormat
      !> Index of linear solver
      integer :: solver
      !> Path to MM geometry file
      character(:), allocatable :: mmGeomFilename
      !> Path to a separate parameter file, if present
      character(:), allocatable :: mmParamsFilename
      !> MM atom types for atoms in the QM zone;
      !! used for vdW parameter access
      integer, allocatable :: qmAtomTypes(:)
      !> Path to parameter file containing MM atom
      !! types for atoms in the QM zone
      character(:), allocatable :: qmParamsFilename
    end type

contains

  !> TODO: add descrption
  subroutine TOpenmmpolPotGen_init(this, input)
    class(TOpenmmpolPotGen), intent(inout) :: this
    type(TOpenmmpolInput), intent(in) ::  input

  #:if WITH_OPENMMPOL
    call notImplementedError
  #:else
    call notImplementedError
  #:endif
  end subroutine TOpenmmpolPotGen_init
  

  !> TODO: add descrption
  subroutine TOpenmmpolPotGen_getExternalPot(this, chargePerAtom, chargePerShell, extPotAtom, extPotShell)
    !> Instance.
    class(TOpenmmpolPotGen), intent(inout) :: this

    !> Net number of electrons on the atom with respect to its reference configuration, the
    !> neutral atom.  Shape: [nAtom].
    real(dp), intent(in) :: chargePerAtom(:)

    !> Shell-resolved net number of electrons. Shape: [mShell, nAtom].
    real(dp), intent(in) :: chargePerShell(:,:)

    !> Atom dependent external potential contribution. Shape: [nAtom]
    real(dp), intent(out) :: extPotAtom(:)

    !> Shell-resolved external potential contribution. Shape: [mShell, nAtom].
    real(dp), intent(out) :: extPotShell(:,:)

  #:if WITH_OPENMMPOL
    ! TODO: implement adding potential
    call notImplementedError
  #:else
    call notImplementedError
  #:endif

  end subroutine TOpenmmpolPotGen_getExternalPot


  !> TODO: add descrption
  subroutine TOpenmmpolPotGen_getExternalPotGrad(this, chargePerAtom, chargePerShell, extPotGrad)
    !> Class instance.
    class(TOpenmmpolPotGen), intent(inout) :: this

    !> Net number of electrons on the atom with respect to its reference configuration, the
    !> neutral atom.  Shape: [nAtom].
    real(dp), intent(in) :: chargePerAtom(:)

    !> Shell-resolved net number of electrons. Shape: [mShell, nAtom].
    real(dp), intent(in) :: chargePerShell(:,:)

    !> Gradient of the potential at each atom. Shape: [3, nAtom].
    real(dp), intent(out) :: extPotGrad(:,:)

  #:if WITH_OPENMMPOL
    call notImplementedError
  #:else
    call notImplementedError
  #:endif

  end subroutine TOpenmmpolPotGen_getExternalPotGrad

  !> TODO: add descrption
  subroutine TOpenmmpolPotGen_getExternalPotFock(this, chargePerAtom, chargePerShell, extPotAtom, extPotShell)
    !> Instance.
    class(TOpenmmpolPotGen), intent(inout) :: this

    !> Net number of electrons on the atom with respect to its reference configuration, the
    !> neutral atom.  Shape: [nAtom].
    real(dp), intent(in) :: chargePerAtom(:)

    !> Shell-resolved net number of electrons. Shape: [mShell, nAtom].
    real(dp), intent(in) :: chargePerShell(:,:)

    !> Atom dependent external potential contribution. Shape: [nAtom]
    real(dp), intent(out) :: extPotAtom(:)

    !> Shell-resolved external potential contribution. Shape: [mShell, nAtom].
    real(dp), intent(out) :: extPotShell(:,:)

  #:if WITH_OPENMMPOL
    ! TODO: implement adding potential
    call notImplementedError
  #:else
    call notImplementedError
  #:endif

  end subroutine TOpenmmpolPotGen_getExternalPotFock


  !> TODO: add descrption
  subroutine TOpenmmpolPotGen_getInternalEnergy(this, output)
    !> Class instance.
    class(TOpenmmpolPotGen), intent(inout) :: this

    !> External energy contribution
    real(dp), intent(out) :: output

  #:if WITH_OPENMMPOL
    output = this%Eff
  #:else
    call notImplementedError
  #:endif

  end subroutine TOpenmmpolPotGen_getInternalEnergy


  !> TODO: add descrption
  subroutine TOpenmmpolPotGen_writeOutput(this)
    class(TOpenmmpolPotGen) :: this

  #:if WITH_OPENMMPOL
    call ommp_print_summary_to_file(this%pSystem, "openmmpol.out")
  #:else
    call notImplementedError
  #:endif

  end subroutine TOpenmmpolPotGen_writeOutput

! TODO: enable later
! #:if not WITH_OPENMMPOL
  subroutine notImplementedError
    call error("DFTB+ compiled without support for the openmmpol library")
  end subroutine notImplementedError
! #:endif
    
end module dftbp_extlibs_openmmpol