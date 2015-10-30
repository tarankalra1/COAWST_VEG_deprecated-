#include "cppdefs.h"
 
       MODULE vegetation_mask_hmixing_mod
#if defined VEGETATION && defined VEG_MASK_HMIXING
!
!svn $Id: vegetation_hmixing.F 429 2015-04-20 17:30:26Z arango $
!=======================================================================
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                   Alexander F. Shchepetkin   !
!================================================John C. Warner=========
!================================================Neil K. Ganju  ========
!=============================================Alexis  Beudin============
!==============================================Tarandeep S. Kalra=======
!                                                                      !
!  Calculate viscosity change at vegetation iterface and add in        !
!  hmixing.F                                                           ! 
!                                                                      !
!  References:                                                         !   
!                                                                      !
!=======================================================================
!                                                                      !
!                                                                      !
!=======================================================================

      implicit none

      PRIVATE
      PUBLIC  :: vegetation_mask_hmixing_cal

      CONTAINS
!
!***********************************************************************
      SUBROUTINE vegetation_mask_hmixing_cal (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_stepping
      USE mod_grid 
      USE mod_ocean
      USE mod_vegarr
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!

# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 16)
# endif
      CALL vegetation_mask_hmixing_tile  (ng, tile,                     &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        VEG(ng) % plant,                          &
     &                        VEG(ng) % mask_veg_bound)
             
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 16)
# endif
      RETURN
      END SUBROUTINE vegetation_mask_hmixing_cal

!***********************************************************************
      SUBROUTINE vegetation_mask_hmixing_tile (ng, tile,                &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              plant,                              &
     &                              mask_veg_bound)          
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_vegetation
      USE mod_vegarr
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: plant(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: mask_veg_bound(LBi:,LBj:)
# else
      real(r8), intent(in) :: plant(LBi:UBi,LBj:UBj,NVEG,NVEGP)
      real(r8), intent(inout) :: mask_veg_bound(LBi:UBi,LBj:UBj)
# endif

!  Local variable declarations.
!
      integer :: i, j, k, iveg
      integer, parameter :: zero = 0
      real(r8), parameter :: Inival=0.0_r8
      real(r8) :: mask_veg_bound_loc
!
# include "set_bounds.h"
!----------------------------------------------------------------------
!----------Executing the code------------------------------------------
!----------------------------------------------------------------------
! 
!   This is currently assuming that vegetation types do not overlap
! 
      mask_veg_bound=Inival
      mask_veg_bound_loc=Inival
!
      DO iveg=1,NVEG 
        DO j=Jstr,Jend
          DO i=Istr,Iend
            mask_veg_bound_loc=MIN(plant(i,j,iveg,pdens),1.0_r8)
            mask_veg_bound_loc=MIN(plant(i,j,iveg,pdens),1.0_r8)
            mask_veg_bound_loc=MIN(plant(i,j,iveg,pdens),1.0_r8)
            mask_veg_bound_loc=MIN(plant(i,j,iveg,pdens),1.0_r8)
            mask_veg_bound_loc=0.0_r8
!
!   Adding from multiple veg types
! 
            mask_veg_bound(i,j)=mask_veg_bound_loc+mask_veg_bound(i,j)
!
!   For overlapping veg types, limit mask to 1 
!
            mask_veg_bound(i,j)=MIN(1.0_r8,mask_veg_bound(i,j))
          END DO
        END DO
      END DO 
      END SUBROUTINE vegetation_mask_hmixing_tile
#endif
      END MODULE vegetation_mask_hmixing_mod
