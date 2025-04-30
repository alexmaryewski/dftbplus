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
    use dftbp_common_environment, only : TEnvironment
    use dftbp_common_status, only : TStatus
    use dftbp_dftb_charges, only : getSummedCharges
    use dftbp_dftb_periodic, only : TNeighbourList
    use dftbp_io_message, only : error, warning
    use dftbp_math_simplealgebra, only : determinant33
    use dftbp_solvation_solvation, only : TSolvation
    use dftbp_type_commontypes, only : TOrbitals
  #:if WITH_OPENMMPOL
    use ommp_interface, only : ommp_system, ommp_init_mmp, ommp_init_xyz, ommp_terminate,&
    ommp_print_summary, ommp_prepare_qm_ele_ene, ommp_set_external_field,&
    ommp_potential_mmpol2ext, ommp_qm_helper, ommp_init_qm_helper, ommp_terminate_qm_helper,&
    ommp_set_verbose, ommp_get_fixedelec_energy, ommp_get_polelec_energy,&
    ommp_get_full_bnd_energy, ommp_get_vdw_energy, ommp_qm_helper_init_vdw_prm,&
    ommp_qm_helper_vdw_energy, ommp_qm_helper_set_attype, ommp_qm_helper_vdw_energy,&
    ommp_print_summary_to_file, ommp_ignore_duplicated_angle_prm,&
    ommp_ignore_duplicated_opb_prm, OMMP_SOLVER_NONE, OMMP_FF_AMOEBA, OMMP_FF_WANG_AL,&
    OMMP_FF_WANG_DL, OMMP_VERBOSE_DEBUG, OMMP_VERBOSE_LOW, OMMP_VERBOSE_HIGH
  #:endif

    implicit none
    
    private
    public TOpenmmpol, TOpenmmpolInput, TOpenmmpol_init, TOpenmmpol_final
    public writeOpenmmpolInfo

    ! Constant for AMBER identification inside openmmpol
    integer, protected :: OMMP_FF_AMBER = 0
  
    !> Library interface handler
    type, extends(TSolvation) :: TOpenmmpol
      !> number of atoms
      integer :: nAtom = 0

      !> solvation free energy
      real(dp), allocatable :: energies(:)

      !> lattice vectors if periodic
      real(dp) :: latVecs(3, 3) = 0.0_dp

      !> Volume of the unit cell
      real(dp) :: volume = 0.0_dp

      !> stress tensor
      real(dp) :: stress(3, 3) = 0.0_dp

      !> is this periodic
      logical :: tPeriodic

      !> are the coordinates current?
      logical :: tCoordsUpdated = .false.

      !> are the charges current?
      logical :: tChargesUpdated = .false.

      !> QM/MM interaction cutoff
      integer :: rCutoff = 0.0_dp

    #:if WITH_OPENMMPOL
      !> Pointer to openmmpol system object
      type(ommp_system), pointer :: pSystem

      !> Pointer to openmmpol QM helper object
      type(ommp_qm_helper), pointer :: pQmHelper
    #:endif

      !> Site-resolved potential
      real(dp), allocatable :: potential(:)

      !> QM/MM electrostatic energy from constant multipoles
      !! (site-resolved)
      real(dp), allocatable :: energyQmElecStat(:)

      !> QM/MM electrostatic energy from polarizable multipoles
      !! (site-resolved)
      real(dp), allocatable :: energyQmElecPol(:)

      !> MM/MM electrostatic energy from constant multipoles
      real(dp) :: energyMmElecStat

      !> MM/MM electrostatic energy from polarizable multipoles
      real(dp) :: energyMmElecPol

      !> Total energy of all bonding terms in the force field
      real(dp) :: energyBonded

      !> Total energy of all non-bonded terms in the force field
      real(dp) :: energyNonbonded

      !> Linear solver
      integer :: solver
      
    contains
    
      !> update internal copy of coordinates
      procedure :: updateCoords

      !> update internal copy of lattice vectors
      procedure :: updateLatVecs
      
      !> get real space cutoff
      procedure :: getRCutoff

      !> get energy contributions
      procedure :: getEnergies

      !> get force contributions
      procedure :: addGradients

      !> get stress tensor contributions
      procedure :: getStress

      !> Updates with changed charges for the instance.
      procedure :: updateCharges

      !> Returns shifts per atom
      procedure :: getShifts

      !> Is the electrostic field modified by this solvent model?
      procedure :: isEFieldModified

      !> Relative dielectric constant for solvent
      procedure :: getEpsilon_r

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
  
  !> Initialize an openmmpol object
  subroutine TOpenmmpol_init(this, input, nAtom, species0, speciesNames, errStatus, coords, latVecs)

    !> Initialised instance at return
    type(TOpenmmpol), intent(out) :: this

    !> Specific input parameters for openmmpol
    type(TOpenmmpolInput), intent(in) :: input

    !> Nr. of atoms in the system
    integer, intent(in) :: nAtom

    !> Species of every atom in the unit cell
    integer, intent(in) :: species0(:)

    !> Symbols of the species
    character(len=*), intent(in) :: speciesNames(:)

    !> Error status
    type(TStatus), intent(out) :: errStatus

    !> Initial atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Lattice vectors, if the system is periodic
    real(dp), intent(in), optional :: latVecs(:,:)

    real(dp), allocatable, dimension(:) :: initCharges
    integer, allocatable :: zVector(:)

  #:if WITH_OPENMMPOL
    this%nAtom = nAtom
    this%solver = input%solver

    this%energyBonded = 0.0_dp
    this%energyNonbonded = 0.0_dp

    allocate(this%potential(nAtom), source=0.0_dp)
    allocate(this%energyQmElecStat(nAtom), source=0.0_dp)
    allocate(this%energyQmElecPol(nAtom), source=0.0_dp)

    ! Tinker compatibility is enforced by default
    call ommp_ignore_duplicated_opb_prm()
    call ommp_ignore_duplicated_angle_prm()

    ! Initialize from input
    select case (input%inputFormat)
    case ("tinker")
      call ommp_init_xyz(this%pSystem, input%mmGeomFilename, input%mmParamsFilename)
    case ("mmp")
      call ommp_init_mmp(this%pSystem, input%mmParamsFilename)
    case default
      call error("Openmmpol is unable to recognize the input file format.")
    end select

    ! Initialize the qmhelper object with nuclear charges
    allocate(initCharges(nAtom), source=0.0_dp)
    allocate(zVector(nAtom))

    ! Openmmpol requires nuclear charges for initialization
    call getNuclChargeVector(zVector, speciesNames, species0)
    call ommp_init_qm_helper(this%pQmHelper, nAtom, coords, initCharges, zVector)
    call ommp_qm_helper_set_attype(this%pQMHelper, input%qmAtomTypes)
    call ommp_qm_helper_init_vdw_prm(this%pQMHelper, input%qmParamsFilename)

    if (this%pSystem%amoeba) then
      this%pQMHelper%V_pp2n_req = .true.
    end if
    
    ! Evaluate bonded and non-bonded energy terms for the initial geometry
    this%energyBonded = ommp_get_full_bnd_energy(this%pSystem)
    this%energyNonbonded = ommp_get_vdw_energy(this%pSystem) + &
        & ommp_qm_helper_vdw_energy(this%pQMHelper, this%pSystem)

    ! Verbosity control
    ! TODO: initialize from input?
    call ommp_set_verbose(OMMP_VERBOSE_DEBUG)
    ! call ommp_set_verbose(OMMP_VERBOSE_LOW)

    deallocate(initCharges)
    deallocate(zVector)

    this%tCoordsUpdated = .false.
    this%tChargesUpdated = .false.

    ! TODO: remove debug lines
    call error("DEBUG: stopping")
  #:else
    call notImplementedError
  #:endif
  
  end subroutine TOpenmmpol_init


  !> Terminate an openmmpol instance
  subroutine TOpenmmpol_final(this)

    !> Instance.
    class(TOpenmmpol), intent(inout) :: this

  #:if WITH_OPENMMPOL
    call ommp_terminate_qm_helper(this%pQmHelper)
    call ommp_terminate(this%pSystem)
    deallocate(this%potential)
    deallocate(this%energyQmElecStat)
    deallocate(this%energyQmElecPol)
  #:else
    call notImplementedError
  #:endif

  end subroutine TOpenmmpol_final


  !> Returns shifts per atom
  subroutine getShifts(this, shiftPerAtom, shiftPerShell)

    !> Instance.
    class(TOpenmmpol), intent(inout) :: this

    !> Shift per atom
    real(dp), intent(out) :: shiftPerAtom(:)

    !> Shift per shell
    real(dp), intent(out) :: shiftPerShell(:,:)

    shiftPerAtom(:) = 0.0_dp
    shiftPerShell(:,:) = 0.0_dp

    shiftPerAtom(:) = shiftPerAtom + this%potential     

  end subroutine getShifts


  !> Get energy contributions
  subroutine getEnergies(this, energies)

    !> Instance.
    class(TOpenmmpol), intent(inout) :: this

    !> energy contributions for each atom
    real(dp), intent(out) :: energies(:)

    ! @:ASSERT(this%tCoordsUpdated)
    ! @:ASSERT(this%tChargesUpdated)
    ! @:ASSERT(size(energies) == this%nAtom)

    ! if (allocated(this%sasaCont)) then
      ! call this%sasaCont%getEnergies(energies)
    ! else
      ! energies(:) = 0.0_dp
    ! end if

    ! energies(:) = energies + 0.5_dp * (this%shift * this%chargesPerAtom) &
      !  & + this%param%freeEnergyShift / real(this%nAtom, dp)

  end subroutine getEnergies


  !> Get force contributions
  subroutine addGradients(this, env, neighList, species, coords, img2CentCell, gradients, errStatus)

    !> Data structure
    class(TOpenmmpol), intent(inout) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Neighbour list.
    type(TNeighbourList), intent(in) :: neighList

    !> Specie for each atom.
    integer, intent(in) :: species(:)

    !> Coordinate of each atom.
    real(dp), intent(in) :: coords(:,:)

    !> Mapping of atoms to cetnral cell.
    integer, intent(in) :: img2CentCell(:)

    !> Gradient contributions for each atom
    real(dp), intent(inout) :: gradients(:,:)

    !> Error status
    type(TStatus), intent(out) :: errStatus

    integer :: ii, iat, ig
    real(dp), allocatable :: fx(:,:), zeta(:), ef1(:,:), ef2(:,:)

  #:if WITH_OPENMMPOL
    call notImplementedError
  #:else
    call notImplementedError
  #:endif

  end subroutine addGradients


  subroutine getStress(this, stress)

    !> Class instance.
    class(TOpenmmpol), intent(inout) :: this

    !> Stress tensor contributions
    real(dp), intent(out) :: stress(:,:)

  #:if WITH_OPENMMPOL
    call notImplementedError
  #:else
    call notImplementedError
  #:endif

  end subroutine getStress


  !> Updates with changed charges for the instance.
  subroutine updateCharges(this, env, species, neighList, qq, q0, img2CentCell, orb, errStatus)

    !> Instance.
    class(TOpenmmpol), intent(inout) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Species, shape: [nAtom]
    integer, intent(in) :: species(:)

    !> Neighbour list.
    type(TNeighbourList), intent(in) :: neighList

    !> Orbital charges.
    real(dp), intent(in) :: qq(:,:,:)

    !> Reference orbital charges.
    real(dp), intent(in) :: q0(:,:,:)

    !> Mapping on atoms in central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> Error status
    type(TStatus), intent(out) :: errStatus

  #:if WITH_OPENMMPOL
    ! DFTB+ internally uses charges and potential that are opposite to the chemical
    ! convention, therefore we need a sign inversion when passing them to openmmpol
    call getSummedCharges(species, orb, qq, q0=q0, dQatom=this%pQMHelper%qqm)
    this%pQMHelper%qqm = -this%pQMHelper%qqm

    ! Compute electric field produced by QM part of the system on MM atoms
    this%pQMHelper%E_n2p_done = .false.

    ! Only computes E_n2p (and V_m2n/V_p2n at first call)
    call ommp_prepare_qm_ele_ene(this%pSystem, this%pQMHelper)

    if (this%pSystem%ff_type .ne. OMMP_FF_AMBER) then

      ! Set external field for MM, solve the polarization equations
      this%pSystem%eel%D2Mgg_done = .false.
      this%pSystem%eel%D2Dgg_done = .false.

      call ommp_set_external_field(this%pSystem, &
         this%pQMHelper%E_n2p, &
         this%solver, &
         OMMP_SOLVER_NONE, &
         .true.)

      ! Set flags for computing the MM-to-QM potential for IPDs
      this%pQMHelper%V_p2n_done = .false.
      if (this%pSystem%amoeba) then
         this%pQMHelper%V_pp2n_done = .false.
      end if

      ! Only computes V_p2n, after having updated the external field/IPDs
      call ommp_prepare_qm_ele_ene(this%pSystem, this%pQMHelper)

      ! Compute and store shifts
      if (this%pSystem%amoeba) then
        this%potential = -((this%pQMHelper%V_pp2n + this%pQMHelper%V_p2n) * 0.5_dp + &
            & this%pQMHelper%V_m2n)
      else
        this%potential = -(this%pQMHelper%V_p2n + this%pQMHelper%V_m2n)
      end if

    end if

    ! Compute and store electrostatic energies;
    ! QM/MM energies are site-resolved
    this%energyQmElecStat(:) = this%pQMHelper%V_m2n * this%pQMHelper%qqm
    this%energyMmElecStat = ommp_get_fixedelec_energy(this%pSystem)
    this%energyQmElecPol(:) = 0.5_dp * this%pQMHelper%V_p2n * this%pQMHelper%qqm
    this%energyMmElecPol = ommp_get_polelec_energy(this%pSystem)

    this%tChargesUpdated = .true.
  #:else
    call notImplementedError
  #:endif

  end subroutine updateCharges


  !> Update internal stored coordinates
  subroutine updateCoords(this, env, neighList, img2CentCell, coords, species0)

    !> Data structure
    class(TOpenmmpol), intent(inout) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> List of neighbours to atoms
    type(TNeighbourList), intent(in) :: neighList

    !> Image to central cell atom index
    integer, intent(in) :: img2CentCell(:)

    !> Atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Central cell chemical species
    integer, intent(in) :: species0(:)

    integer, allocatable :: nNeigh(:)

  #:if WITH_OPENMMPOL
    call notImplementedError
    
    !> TODO: do not forget to update the MM coordinates as well
    ! Evaluate bonded and non-bonded energy terms for the initial geometry
    this%energyBonded = ommp_get_full_bnd_energy(this%pSystem)
    this%energyNonbonded = ommp_get_vdw_energy(this%pSystem) + &
        & ommp_qm_helper_vdw_energy(this%pQMHelper, this%pSystem)

  #:else
    call notImplementedError
  #:endif

  end subroutine updateCoords


  !> Update internal copy of lattice vectors
  subroutine updateLatVecs(this, latVecs)

    !> Data structure
    class(TOpenmmpol), intent(inout) :: this

    !> Lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

    @:ASSERT(this%tPeriodic)
    @:ASSERT(all(shape(latvecs) == shape(this%latvecs)))

    this%volume = abs(determinant33(latVecs))
    this%latVecs(:,:) = latVecs

    this%tCoordsUpdated = .false.
    this%tChargesUpdated = .false.

  end subroutine updateLatVecs


  !> Print Openmmpol settings to external file
  subroutine writeOpenmmpolInfo(unit, solvation)

    !> Unit for IO
    integer, intent(in) :: unit

    !> Instance.
    type(TOpenmmpol), intent(in) :: solvation

  #:if WITH_OPENMMPOL
    !> TODO: throw in some info for good measure
    ! call ommp_print_summary_to_file(solvation%pSystem, unit)
  #:else
    call notImplementedError
  #:endif

  end subroutine writeOpenmmpolInfo
  

  !> Is the electrostic field modified by this solvent model?
  pure function isEFieldModified(this) result(isChanged)

    !> Data instance
    class(TOpenmmpol), intent(in) :: this

    !> Has the solvent model changed the electrostatic environment
    logical :: isChanged

    isChanged = .true.

  end function isEFieldModified


  !> Returns solvent region relative dielectric constant
  pure function getEpsilon_r(this) result(e_r)

    !> Data structure
    class(TOpenmmpol), intent(in) :: this

    !> epsilon_r
    real(dp) :: e_r

    e_r = 1.0_dp

  end function getEpsilon_r


  !> Distance cut off for QM/MM interaction
  function getRCutoff(this) result(cutoff)

    !> Instance.
    class(TOpenmmpol), intent(inout) :: this

    !> Resulting cutoff
    real(dp) :: cutoff

    cutoff = this%rCutoff

  end function getRCutoff


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


! #:if not WITH_OPENMMPOL
  subroutine notImplementedError
    call error("DFTB+ compiled without support for the openmmpol library")
  end subroutine notImplementedError
! #:endif
    
end module dftbp_extlibs_openmmpol