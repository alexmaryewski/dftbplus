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
    public TTOpenmmpol, TOpenmmpolInput, TTOpenmmpol_init


    ! Constant for AMBER identification inside openmmpol
    integer, protected :: OMMP_FF_AMBER = 0
  

    !> Library interface handler
    type :: TTOpenmmpol
    #:if WITH_OPENMMPOL
      !> Pointer to openmmpol system object
      type(ommp_system), pointer :: pSystem

      !> Pointer to openmmpol QM helper object
      type(ommp_qm_helper), pointer :: pQmHelper
    #:endif

      !> Site-resolved potential
      real(dp), allocatable :: potential(:)

      !> Site-resolved potential due to the Lagrangian (if present)
      real(dp), allocatable :: potentialLagrangian(:)

      !> Total energy of all bonding terms in the force field
      real(dp) :: energyBonded

      !> Total energy of all non-bonded terms in the force field
      real(dp) :: energyNonbonded

      !> Linear solver
      integer :: solver
      
      contains
        procedure :: getPotential => TTOpenmmpol_getPotential
        procedure :: getPotentialLagrangian => TTOpenmmpol_getPotentialLagrangian
        procedure :: getPotentialGradient => TTOpenmmpol_getPotentialGradient
        procedure :: getPotentialGradientLagrangian => TTOpenmmpol_getPotentialGradientLagrangian
        procedure :: getInternalEnergy => TTOpenmmpol_getInternalEnergy
        procedure :: writeOutput => TTOpenmmpol_writeOutput
        procedure :: updateCharges => TTOpenmmpol_updateCharges
        procedure :: updateCoords => TTOpenmmpol_updateCoords
    end type


    !> Data type for storing openmmpol-specific input parameters
    type :: TOpenmmpolInput
      !> Used input format (either "Tinker" or "mmp")
      character(:), allocatable :: inputFormat

      !> Index of linear solver
      integer :: solver

      !> Path to MM geometry file
      character(:), allocatable :: mmGeomFilename

      !> Path to a separate parameter file, if present
      character(:), allocatable :: mmParamsFilename

      !> MM atom types for atoms in the QM zone
      integer, allocatable :: qmAtomTypes(:)

      !> Path to parameter file containing MM atom
      !! types for atoms in the QM zone
      character(:), allocatable :: qmParamsFilename
    end type

contains
  
  subroutine TTOpenmmpol_init(this, input, nAtom, species0, speciesNames, coords0)

    !> Instance of the library interface
    class(TTOpenmmpol), intent(out) :: this

    !> Input to construct the library interface from
    type(TOpenmmpolInput), intent(in) ::  input

    !> Nr. of atoms in the system
    integer, intent(in) :: nAtom

    !> Species of every atom in the unit cell
    integer, intent(in) :: species0(:)

    !> Atomic coordinates in the unit cell
    real(dp), intent(in) :: coords0(:,:)

    !> Symbols of the species
    character(len=*), intent(in) :: speciesNames(:)

    real(dp), allocatable, dimension(:) :: charges
    integer, allocatable :: Zvector(:)

  #:if WITH_OPENMMPOL
    this%solver = input%solver
    this%energyBonded = 0.0_dp
    this%energyNonbonded = 0.0_dp
    allocate(this%potential(nAtom), source=0.0_dp)

    ! Tinker compatibility by default
    call ommp_ignore_duplicated_opb_prm()
    call ommp_ignore_duplicated_angle_prm()

    select case (input%inputFormat)
    case ("tinker")
      call ommp_init_xyz(this%pSystem, input%mmGeomFilename, input%mmParamsFilename)
    case ("mmp")
      call ommp_init_mmp(this%pSystem, input%mmParamsFilename)
    case default
      call error("Bad openmmpol input format supplied to initializer!")
    end select

    allocate(charges(nAtom), source=0.0_dp)
    allocate(Zvector(nAtom))
    call getNuclChargeVector(Zvector, speciesNames, species0)
    call ommp_init_qm_helper(this%pQmHelper, nAtom, coords0, charges, Zvector)
    call ommp_qm_helper_set_attype(this%pQMHelper, input%qmAtomTypes)
    call ommp_qm_helper_init_vdw_prm(this%pQMHelper, input%qmParamsFilename)

    if (this%pSystem%amoeba) then
      allocate(this%potentialLagrangian(nAtom), source=0.0_dp)
      this%pQMHelper%V_pp2n_req = .true.
    end if
    
    this%energyBonded = ommp_get_full_bnd_energy(this%pSystem)
    this%energyNonbonded = ommp_get_vdw_energy(this%pSystem) + &
        & ommp_qm_helper_vdw_energy(this%pQMHelper, this%pSystem)

    call ommp_set_verbose(OMMP_VERBOSE_DEBUG)
    ! call ommp_set_verbose(OMMP_VERBOSE_LOW)

    deallocate(charges)
    deallocate(Zvector)

    call error("DEBUG: stopping")
  #:else
    call notImplementedError
  #:endif
  end subroutine TTOpenmmpol_init
  

  subroutine TTOpenmmpol_getPotential(this)
    !> Instance.
    class(TTOpenmmpol), intent(inout) :: this

  #:if WITH_OPENMMPOL
    call notImplementedError
  #:else
    call notImplementedError
  #:endif

  end subroutine TTOpenmmpol_getPotential


  subroutine TTOpenmmpol_getPotentialGradient(this)
    !> Class instance.
    class(TTOpenmmpol), intent(inout) :: this

  #:if WITH_OPENMMPOL
    call notImplementedError
  #:else
    call notImplementedError
  #:endif

  end subroutine TTOpenmmpol_getPotentialGradient


  subroutine TTOpenmmpol_getPotentialLagrangian(this)
    !> Instance.
    class(TTOpenmmpol), intent(inout) :: this


  #:if WITH_OPENMMPOL
    call notImplementedError
  #:else
    call notImplementedError
  #:endif

  end subroutine TTOpenmmpol_getPotentialLagrangian


  subroutine TTOpenmmpol_getPotentialGradientLagrangian(this)
    !> Instance.
    class(TTOpenmmpol), intent(inout) :: this

  #:if WITH_OPENMMPOL
    call notImplementedError
  #:else
    call notImplementedError
  #:endif

  end subroutine TTOpenmmpol_getPotentialGradientLagrangian


  subroutine TTOpenmmpol_getInternalEnergy(this, output)
    !> Class instance.
    class(TTOpenmmpol), intent(inout) :: this

    !> External energy contribution
    real(dp), intent(out) :: output

  #:if WITH_OPENMMPOL
    output = this%energyBonded
  #:else
    call notImplementedError
  #:endif

  end subroutine TTOpenmmpol_getInternalEnergy


  subroutine TTOpenmmpol_updateCharges(this)
    !> Class instance.
    class(TTOpenmmpol), intent(inout) :: this

  #:if WITH_OPENMMPOL
    call notImplementedError
  #:else
    call notImplementedError
  #:endif

  end subroutine TTOpenmmpol_updateCharges


  subroutine TTOpenmmpol_updateCoords(this)
    !> Class instance.
    class(TTOpenmmpol), intent(inout) :: this

  #:if WITH_OPENMMPOL
    call notImplementedError
  #:else
    call notImplementedError
  #:endif

  end subroutine TTOpenmmpol_updateCoords


  subroutine TTOpenmmpol_writeOutput(this)
    class(TTOpenmmpol) :: this

  #:if WITH_OPENMMPOL
    call ommp_print_summary_to_file(this%pSystem, "openmmpol.out")
  #:else
    call notImplementedError
  #:endif

  end subroutine TTOpenmmpol_writeOutput


  subroutine getNuclChargeVector(Zvector, speciesNames, species)
    !> Nuclear charge vector
    integer, intent(out) :: Zvector(:)

    !> Atom type names
    character(len=*), intent(in) :: speciesNames(:)

    !> Species in internal format
    integer, intent(in) :: species(:)

    integer :: i
    integer :: currentSpecies

    do i = 1, size(species)
      currentSpecies = species(i)
      select case (speciesNames(currentSpecies))
        case ('H')
          Zvector(i) = 1
        case ('He')
          Zvector(i) = 2
        case ('Li')
          Zvector(i) = 3
        case ('Be')
          Zvector(i) = 4
        case ('B')
          Zvector(i) = 5
        case ('C')
          Zvector(i) = 6
        case ('N')
          Zvector(i) = 7
        case ('O')
          Zvector(i) = 8
        case ('F')
          Zvector(i) = 9
        case ('Ne')
          Zvector(i) = 10
        case ('Na')
          Zvector(i) = 11
        case ('Mg')
          Zvector(i) = 12
        case ('Al')
          Zvector(i) = 13
        case ('Si')
          Zvector(i) = 14
        case ('P')
          Zvector(i) = 15
        case ('S')
          Zvector(i) = 16
        case ('Cl')
          Zvector(i) = 17
        case ('Ar')
          Zvector(i) = 18
        case ('K')
          Zvector(i) = 19
        case ('Ca')
          Zvector(i) = 20
        case ('Sc')
          Zvector(i) = 21
        case ('Ti')
          Zvector(i) = 22
        case ('V')
          Zvector(i) = 23
        case ('Cr')
          Zvector(i) = 24
        case ('Mn')
          Zvector(i) = 25
        case ('Fe')
          Zvector(i) = 26
        case ('Co')
          Zvector(i) = 27
        case ('Ni')
          Zvector(i) = 28
        case ('Cu')
          Zvector(i) = 29
        case ('Zn')
          Zvector(i) = 30
        case ('Ga')
          Zvector(i) = 31
        case ('Ge')
          Zvector(i) = 32
        case ('As')
          Zvector(i) = 33
        case ('Se')
          Zvector(i) = 34
        case ('Br')
          Zvector(i) = 35
        case ('Kr')
          Zvector(i) = 36
        case ('Rb')
          Zvector(i) = 37
        case ('Sr')
          Zvector(i) = 38
        case ('Y')
          Zvector(i) = 39
        case ('Zr')
          Zvector(i) = 40
        case ('Nb')
          Zvector(i) = 41
        case ('Mo')
          Zvector(i) = 42
        case ('Tc')
          Zvector(i) = 43
        case ('Ru')
          Zvector(i) = 44
        case ('Rh')
          Zvector(i) = 45
        case ('Pd')
          Zvector(i) = 46
        case ('Ag')
          Zvector(i) = 47
        case ('Cd')
          Zvector(i) = 48
        case ('In')
          Zvector(i) = 49
        case ('Sn')
          Zvector(i) = 50
        case ('Sb')
          Zvector(i) = 51
        case ('Te')
          Zvector(i) = 52
        case ('I')
          Zvector(i) = 53
        case ('Xe')
          Zvector(i) = 54
        case ('Cs')
          Zvector(i) = 55
        case ('Ba')
          Zvector(i) = 56
        case ('La')
          Zvector(i) = 57
        case ('Ce')
          Zvector(i) = 58
        case ('Pr')
          Zvector(i) = 59
        case ('Nd')
          Zvector(i) = 60
        case ('Pm')
          Zvector(i) = 61
        case ('Sm')
          Zvector(i) = 62
        case ('Eu')
          Zvector(i) = 63
        case ('Gd')
          Zvector(i) = 64
        case ('Tb')
          Zvector(i) = 65
        case ('Dy')
          Zvector(i) = 66
        case ('Ho')
          Zvector(i) = 67
        case ('Er')
          Zvector(i) = 68
        case ('Tm')
          Zvector(i) = 69
        case ('Yb')
          Zvector(i) = 70
        case ('Lu')
          Zvector(i) = 71
        case ('Hf')
          Zvector(i) = 72
        case ('Ta')
          Zvector(i) = 73
        case ('W')
          Zvector(i) = 74
        case ('Re')
          Zvector(i) = 75
        case ('Os')
          Zvector(i) = 76
        case ('Ir')
          Zvector(i) = 77
        case ('Pt')
          Zvector(i) = 78
        case ('Au')
          Zvector(i) = 79
        case ('Hg')
          Zvector(i) = 80
        case ('Tl')
          Zvector(i) = 81
        case ('Pb')
          Zvector(i) = 82
        case ('Bi')
          Zvector(i) = 83
        case ('Po')
          Zvector(i) = 84
        case ('At')
          Zvector(i) = 85
        case ('Rn')
          Zvector(i) = 86
        case ('Fr')
          Zvector(i) = 87
        case ('Ra')
          Zvector(i) = 88
        case default
          call error("Unrecognized atom name")
       end select
    end do

 end subroutine getNuclChargeVector

! TODO: enable conditional compilation later
! #:if not WITH_OPENMMPOL
  subroutine notImplementedError
    call error("DFTB+ compiled without support for the openmmpol library")
  end subroutine notImplementedError
! #:endif
    
end module dftbp_extlibs_openmmpol