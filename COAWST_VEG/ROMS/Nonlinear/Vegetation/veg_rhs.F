#include "cppdefs.h"

      MODULE vegetation_rhs_mod

#if defined VEGETATION && defined VEG_RHS 
!
!svn $Id: veg_rhs.F 429 2015-05-26 10:10:26Z arango $
!=======================================================================
!  Copyright (c) 2002-2014 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                   Alexander F. Shchepetkin   !
!==================================================== John C. Warner ===
!==================================================== Neil K. Ganju  ===
!==================================================== Alexis Beudin  ===
!==================================================Tarandeep S. Kalra===
!                                                                      !
!  This routine computes the effect of vegetation on the momentum      !
!  equation as the vegetation acts as a sink                           !
!                                                                      !  
!  References:                                                         !
!                                                                      !
!   Uittenbogaard, R. (2003), Modelling turbulence in vegetated        !
!   aquatic  flows, International workshop on RIParian FORest          !
!   vegetated channels: hydraulic, morphological and ecological        !
!   aspects, 20-22 February 2003, Trento, Italy.                       !
!                                                                      ! 
!   Luhar M., and H. M. Nepf (2011), Flow-induced reconfiguration of   !
!   buoyant and flexible aquatic vegetation, Limnology and             !
!   Oceanography, 56(6): 2003-2017.                                    !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: vegetation_rhs_cal

      CONTAINS
!
!***********************************************************************
      SUBROUTINE vegetation_rhs_cal (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_stepping 
      USE mod_grid
      USE mod_ocean
      USE mod_vegarr
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
      CALL vegetation_rhs_tile   (ng, tile,                             &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nrhs(ng),                                 & 
     &                        GRID(ng) % Hz,                            &
     &                        OCEAN(ng) % u,                            &
     &                        OCEAN(ng) % v,                            &
     &                        VEG(ng) % plant,                          &
# ifdef VEG_FLEX 
     &                        VEG(ng) % bend,                           &
# endif 
     &                        VEG(ng) % ru_loc_veg,                     &
     &                        VEG(ng) % rv_loc_veg,                     &
     &                        VEG(ng) % ru_veg,                         &
     &                        VEG(ng) % rv_veg,                         & 
     &                        VEG(ng) % Lveg)                      

# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 16)
# endif
      RETURN
      END SUBROUTINE vegetation_rhs_cal
!
!***********************************************************************
      SUBROUTINE vegetation_rhs_tile   (ng, tile,                       &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              nrhs,                               & 
     &                              Hz,                                 &
     &                              u, v,                               &
     &                              plant,                              &
# ifdef VEG_FLEX 
     &                              bend,                               &
# endif 
     &                              ru_loc_veg, rv_loc_veg,             &
     &                              ru_veg, rv_veg,                     &
     &                              Lveg)
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
      integer, intent(in) :: nrhs
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: plant(LBi:,LBj:,:,:)
#  ifdef VEG_FLEX 
      real(r8), intent(inout) :: bend(LBi:,LBj:,:)
#  endif 
      real(r8), intent(inout) :: ru_loc_veg(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: rv_loc_veg(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ru_veg(LBi:,LBj:,:)
      real(r8), intent(inout) :: rv_veg(LBi:,LBj:,:)
      real(r8), intent(inout) :: Lveg(LBi:,LBj:,:)
# else
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: plant(LBi:UBi,LBj:UBj,NVEG,NVEGP)
#  ifdef VEG_FLEX 
      real(r8), intent(inout) :: bend(LBi:UBi,LBj:UBj,N(ng),NVEG)
#  endif 
      real(r8), intent(inout) ::                                        &
    &               ru_loc_veg(LBi:UBi,LBj:UBj,N(ng),NVEG),         &
    &               rv_loc_veg(LBi:UBi,LBj:UBj,N(ng),NVEG)
      real(r8), intent(inout) :: ru_veg(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: rv_veg(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: Lveg(LBi:UBi,LBj:UBj,N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k, iveg
! 
      real(r8), parameter :: one_third  = 1.0_r8/3.0_r8
      real(r8) :: sma, buoy, cff1 , cff2, cff3 
      real(r8) :: theta, Uds, Ca, cflex, fveg, fac 
      real(r8), dimension(0:N(ng)) :: dab  
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vegu
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vegv
!
# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Resistance imposed on the flow by the submerged (seagrass: veg) 
!  or emergent (marsh) vegetation.
!-----------------------------------------------------------------------
!
!  Compute distance above bed - initialising dab 

      dab(0)=0.0_r8 
      VEG_LOOP: DO iveg=1,NVEG
        K_LOOP: DO k=1,N(ng)
          DO j=JstrV-1,Jend
            DO i=IstrU-1,Iend
!
!  Select the grid cell (full or part) within the canopy layer
!
              dab(k)=dab(k-1)+Hz(i,j,k)
              cff1=MIN((dab(k)-plant(i,j,iveg,phght))/Hz(i,j,k),1.0_r8)
              Lveg(i,j,k)=MIN(1.0_r8-cff1,1.0_r8)
               
# ifdef VEG_FLEX 
              sma=(plant(i,j,iveg,pdiam)*                               &
     &             plant(i,j,iveg,pthck)**3.0_r8)/12.0_r8
!
! Buoyancy parameter    
!          
              buoy=(rhow-veg_massdens(iveg,ng))*g*plant(i,j,iveg,pdiam)*&
     &              plant(i,j,iveg,pthck)*plant(i,j,iveg,phght)**3.0_r8/&
     &                                    (E_VEG(iveg,ng)*sma)
!
! Rotate the horizontal components of flow (u and v) at rho points
! into a downstream (Uds) component
!
              cff2=0.5_r8*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
              cff3=0.5_r8*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
              theta=ATAN2(cff3,cff2)
              Uds=cff2*COS(theta)+cff3*SIN(theta)
!
! Cauchy number
!
              Ca=0.5_r8*rhow*Cd_veg(iveg,ng)*plant(i,j,iveg,pdiam)*     &
     &                   Uds**2.0_r8*plant(i,j,iveg,phght)**3.0_r8/     &
     &                       (E_veg(iveg,ng)*sma)
!
              cflex=1.0_r8-((1.0_r8-0.9_r8*Ca**(-one_third))/           &
     &             (1.0_r8+(Ca**(-1.5_r8)*(8.0_r8+buoy**(1.5_r8)))))
!         
! To avoid NaN value when Ca is nul 
!
              cflex=MIN(cflex,1.0_r8)
!
! Effective blade length takes into account: 
!
              Lveg(i,j,k)=cflex*Lveg(i,j,k)

!   Compute friction force *at ask Alexis:*    
              fveg=0.5_r8*cd_veg(iveg,ng)*plant(i,j,iveg,pdiam)*        &
     &             plant(i,j,iveg,pdens)*Hz(i,j,k)*Lveg(i,j,k)
              cff2=Uds*COS(theta)
              cff3=Uds*SIN(theta)
              vegu(i,j)=fveg*cff2*ABS(cff2)
              vegv(i,j)=fveg*cff3*ABS(cff3)
! 
!  Compute blade bending angle
!
              bend(i,j,iveg)=ACOS(cflex**(-3.0_r8))
            END DO 
          END DO
! 
!   Else statement to compute for non-flexible vegetation
!       
# else 
!
!   Compute friction force at rho points          
!       
              fveg=0.5_r8*cd_veg(iveg,ng)*plant(i,j,iveg,pdiam)*        &
     &             plant(i,j,iveg,pdens)*Hz(i,j,k)*Lveg(i,j,k)
              cff2=0.5_r8*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
              cff3=0.5_r8*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
              vegu(i,j)=fveg*cff2*ABS(cff2)
              vegv(i,j)=fveg*cff3*ABS(cff3)
            END DO 
          END DO
# endif 
!
!   Interpolate at cell face   
!
          DO j=Jstr,Jend
            DO i=IstrU,Iend
              ru_loc_veg(i,j,k,iveg)=0.5_r8*(vegu(i-1,j)+vegu(i,j))
!
!   Add the ru_iveg from this veg type to another veg type
!   which can be there at the same point (i,j,k) 
!
              ru_veg(i,j,k)=ru_loc_veg(i,j,k,iveg)+ru_veg(i,j,k) 
            END DO
          END DO
          DO j=JstrV,Jend
            DO i=Istr,Iend
              rv_loc_veg(i,j,k,iveg)=0.5_r8*(vegv(i,j-1)+vegv(i,j))
!
!   Add the rv_iveg from this veg type to another veg type
!   which can be there at the same point (i,j,k) 
!
              rv_veg(i,j,k)=rv_loc_veg(i,j,k,iveg)+rv_veg(i,j,k)
            END DO
          END DO
        END DO K_LOOP
      END DO VEG_LOOP
!             
      RETURN
      END SUBROUTINE vegetation_rhs_tile
#endif
      END MODULE vegetation_rhs_mod
