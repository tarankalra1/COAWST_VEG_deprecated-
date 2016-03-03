      SUBROUTINE ana_vegetation (ng, tile, model)
!                                                                      ! 
!! svn $Id: ana_vegetation.h 429 2015-18-05 17:00:25 Z arango $        !
!!=====================================================================!
!! Copyright (c) 2002-2014 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!================================================== John C. Warner ====!
!==================================================== Neil K. Ganju  ==!
!==================================================== Alexis Beudin  ==!
!==================================================Tarandeep S. Kalra==!
!                                                                      !
!  Vegetation Model Kernel Variables:                                  !
!  NVEG          Number of vegetation types                            !
!  NVEGP         Varying vegetation properties                         !
!  plant         Vegetation variable properties:                       !
!                   plant(:,:,:,pdiam) => diameter                     !
!                   plant(:,:,:,phght) => height                       !
!                   plant(:,:,:,pdens) => density                      !
!                   plant(:,:,:,pthck) => thickness                    !
!                   plant(:,:,:,pupbm) => above ground biomass         !
!                   plant(:,:,:,pdwbm) => below ground biomass         !
!                                                                      !
!  This routine sets initial conditions for vegetation fields          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ncparam
      USE mod_vegarr
      
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_vegetation_tile (ng, tile, model,                        &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        VEG(ng) % plant)
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(48)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_vegetation
!
!***********************************************************************
      SUBROUTINE ana_vegetation_tile (ng, tile, model,                  &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              plant)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_vegetation
      USE mod_vegarr
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: plant(LBi:,LBj:,:,:)
#else
      real(r8), intent(inout) :: plant(LBi:UBi,LBj:UBj,NVEG,NVEGP)
#endif
!
!  Local variable declarations.
!
#ifdef DISTRIBUTE
      integer :: Tstr, Tend
#endif
      integer :: i, j, k, iveg

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set initial properties for each plant 
!  To have variable properties in array->plant(x,y,iveg,iprop)
!-----------------------------------------------------------------------
!
!     
! the vegetation patch has been hardwired in a 100x100 domain to be
! spread over '11' points in the middle of the domain
! NVEG=1 (input from vegetation.in) file in /ROMS/External folder
! 
      DO iveg=1,NVEG
        plant(45:55,45:55,iveg,pdiam)=0.01_r8   ! Diameter
        plant(45:55,45:55,iveg,phght)=1.0_r8    ! Height
        plant(45:55,45:55,iveg,pdens)=2500.0_r8 ! Density
        plant(45:55,45:55,iveg,pthck)=0.0005_r8 ! Thickness
#ifdef VEGETATION_BIOMASS 
        plant(45:55,45:55,iveg,pagbm)=0.0_r8    !Above ground Biomass
        plant(45:55,45:55,iveg,pbgbm)=0.0_r8    !Below ground Biomass
#endif            
      END DO 
                                        
      RETURN

      END SUBROUTINE ana_vegetation_tile
