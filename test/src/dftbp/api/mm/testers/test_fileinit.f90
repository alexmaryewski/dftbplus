!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Reads in a dftb_in.hsd file then overrides some of its supplied setting, before calculating
!> properties
program test_fileinit
  use dftbplus, only : getDftbPlusApi, getDftbPlusBuild, TDftbPlus, TDftbPlus_init, TDftbPlusInput
  use dftbp_common_constants, only : AA__Bohr
  ! Only needed for the internal test system
  use testhelpers, only : writeAutotestTag
  implicit none

  ! How often should the test repeated
  integer, parameter :: nIter = 2

  ! double precision real
  integer, parameter :: dp = kind(1.0d0)

  ! reference coordinates, (xyz,:nAtom) in atomic units
  real(dp), parameter :: initialCoords(3, 2) = reshape([&
      & 0.0000000000000000_dp, 0.0000000000000000_dp, 0.0000000000000000_dp,&
      & 2.5639291987021915_dp, 2.5639291987021915_dp, 2.5639291987021915_dp], [3, 2])

  ! lattice vectors in atomic units
  real(dp), parameter :: initialLatVecs(3, 3) = reshape([&
      & 5.1278583974043830_dp, 5.1278583974043830_dp, 0.0000000000000000_dp,&
      & 0.0000000000000000_dp, 5.1278583974043830_dp, 5.1278583974043830_dp,&
      & 5.1278583974043830_dp, 0.0000000000000000_dp, 5.1278583974043830_dp], [3, 3])

  call main_()

contains

  !! Main test routine
  !!
  !! All non-constant variables must be defined here to ensure that they are all explicitely
  !! deallocated before the program finishes  (avoiding residual memory that tools like valgrind notice).
  !!
  subroutine main_()

    type(TDftbPlus), allocatable :: dftbp
    type(TDftbPlusInput), allocatable :: input

    real(dp) :: merminEnergy
    real(dp) :: coords(3, 2), latVecs(3, 3), gradients(3, 2), grossCharges(2), stressTensor(3, 3)

    integer :: devNull

    character(:), allocatable :: DftbVersion
    integer :: major, minor, patch
    integer :: iIter

    do iIter = 1, nIter
      allocate(dftbp)
      allocate(input)

      call getDftbPlusBuild(DftbVersion)
      write(*,*)'DFTB+ build: ' // "'" // trim(DftbVersion) // "'"
      call getDftbPlusApi(major, minor, patch)
      write(*,"(1X,A,1X,I0,'.',I0,'.',I0)")'API version:', major, minor, patch

      ! Note: setting the global standard output to /dev/null to suppress output and run-time error
      ! messages
      open(newunit=devNull, file="/dev/null", action="write")

      !call TDftbPlus_init(dftbp, outputUnit=devNull)
      call TDftbPlus_init(dftbp)

      ! You should provide the dftb_in.hsd and skfile found in the test/app/dftb+/non-scc/Si_2/ folder
      call dftbp%getInputFromFile("dftb_in.hsd", input)
      call dftbp%setupCalculator(input)

      ! a simple check for data values in the file
      if (dftbp%nrOfAtoms() /= size(coords, dim=2)) then
        print *, "Mismatch between expected number of atoms and dftb_in.hsd file"
        return
      end if

      ! replace the lattice vectors and coordinates in the document tree
      latVecs(:,:) = initialLatVecs
      coords(:,:) = initialCoords
      call dftbp%setGeometry(coords, latVecs)

      ! evaluate energy and forces
      call dftbp%getEnergy(merminEnergy)
      call dftbp%getGradients(gradients)
      print "(A,F15.10)", 'Obtained Mermin Energy:', merminEnergy
      print "(A,3F15.10)", 'Obtained gradient of atom 1:', gradients(:,1)

      ! make a small displacement in the lattice vectors and coordinates
      latVecs(1, 1) = latVecs(1, 1) + 0.1_dp * AA__Bohr
      coords(1, 1) = coords(1, 1) + 0.2_dp * AA__Bohr
      call dftbp%setGeometry(coords, latVecs)

      ! re-calculate energy and forces
      call dftbp%getEnergy(merminEnergy)
      call dftbp%getGradients(gradients)
      call dftbp%getGrossCharges(grossCharges)
      call dftbp%getStressTensor(stressTensor)
      print "(A,F15.10)", 'Obtained Mermin Energy:', merminEnergy
      print "(A,3F15.10)", 'Obtained gradient of atom 1:', gradients(:,1)
      print "(A,2F15.10)", 'Obtained charges:', grossCharges

      ! Deallocating instances to trigger finalizers
      deallocate(dftbp)
      deallocate(input)

    end do

    ! Write file for internal test system
    call writeAutotestTag(merminEnergy=merminEnergy, gradients=gradients, grossCharges=grossCharges,&
        & stressTensor=stressTensor)

  end subroutine main_

end program test_fileinit
