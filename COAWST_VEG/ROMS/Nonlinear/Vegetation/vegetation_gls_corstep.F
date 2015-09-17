#include "cppdefs.h"

      MODULE vegetation_gls_corstep

#if defined NONLINEAR && defined VEGETATION && defined VEG_TURB
!
!svn $Id: vegetation_gls_corstep.F 429 2015-06-10 12:30:26Z arango $
!=======================================================================
!  Copyright (c) 2002-2014 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                                              !
!==================================================== John C. Warner ===
!==================================================== Neil K. Ganju  ===
!==================================================== Alexis Beudin  ===
!==================================================Tarandeep S. Kalra===
!                                                                      !
!  This routine computes the turbulent kinetic energy and length scale ! 
!  modifications due to the vegetation. This is fed back to the        !
!  gls_corstep.F eventually.                                           !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: vegetation_turb_cal

      CONTAINS
!
!***********************************************************************
      SUBROUTINE vegetation_turb_cal (ng, tile)
!***********************************************************************
!
      USE mod_stepping 
      USE mod_grid
      USE mod_ocean
      USE mod_mixing 
      USE mod_vegarr
      USE vegetation_drag_mod
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
      CALL vegetation_turb_tile  ( ng, tile,                            &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nstp(ng), nnew(ng),                       &
     &                        OCEAN(ng) % u,                            &
     &                        OCEAN(ng) % v,                            &
     &                        VEG(ng) % ru_loc_veg,                     &
     &                        VEG(ng) % rv_loc_veg,                     &
     &                        VEG(ng) % plant,                          &
# ifdef VEG_FLEX 
     &                        VEG(ng) % bend,                           &
# endif 
     &                        MIXING(ng) % gls,                         &
     &                        MIXING(ng) % tke,                         &
     &                        VEG(ng) % gls_veg,                        &
     &                        VEG(ng) % tke_veg )
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 16)
# endif
      RETURN
      END SUBROUTINE vegetation_turb_cal 
!
!***********************************************************************
      SUBROUTINE vegetation_turb_tile ( ng, tile,                       &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              nstp, nnew,                         &
     &                              u, v,                               &
     &                              ru_loc_veg, rv_loc_veg,             &
     &                              plant,                              &
# ifdef VEG_FLEX 
     &                              bend,                               &
# endif 
     &                              gls, tke,                           &
     &                              gls_veg, tke_veg )
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_vegetation 
      USE mod_vegarr
      USE vegetation_drag_mod
!
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew 
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: ru_loc_veg(LBi:,LBj:,:,:)
      real(r8), intent(in) :: rv_loc_veg(LBi:,LBj:,:,:)
      real(r8), intent(in) :: plant(LBi:,LBj:,:,:)
# ifdef VEG_FLEX 
      real(r8), intent(in) :: bend(LBi:,LBj:,:)
# endif 
      real(r8), intent(in) :: gls(LBi:,LBj:,0:,:)
      real(r8), intent(in) :: tke(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: gls_veg(LBi:,LBj:,0:)
      real(r8), intent(inout) :: tke_veg(LBi:,LBj:,0:)
# else
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),nstp)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),nstp)
      real(r8), intent(in) :: ru_loc_veg(LBi:UBi,LBj:UBj,N(ng),NVEG)
      real(r8), intent(in) :: rv_loc_veg(LBi:UBi,LBj:UBj,N(ng),NVEG)
      real(r8), intent(in) :: plant(LBi:UBi,LBj:UBj,NVEG,NVEGP)
# ifdef VEG_FLEX 
      real(r8), intent(in) :: bend(LBi:UBi,LBj:UBj,NVEG)
# endif 
      real(r8), intent(in) :: gls(LBi:UBi,LBj:UBj,0:N(ng),nnew)
      real(r8), intent(in) :: tke(LBi:UBi,LBj:UBj,0:N(ng),nnew)
      real(r8), intent(inout) :: gls_veg(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: tke_veg(LBi:UBi,LBj:UBj,0:N(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k, iveg
!
      real(r8), parameter :: one_half=1.0_r8/2.0_r8 
      real(r8), parameter :: one_third=1.0_r8/3.0_r8 
      real(r8), parameter :: Inival=0.0_r8
      real(r8), parameter :: cl_veg=1.0_r8 
      real(r8) :: wku1, wku2, wku3, wku4, wku
      real(r8) :: wkv1, wkv2, wkv3, wkv4, wkv
      real(r8) :: cmu, L
      real(r8) :: kveg
      real(r8) :: solid, eqvegT
      real(r8) :: taufree, tauveg, taueff
      real(r8) :: tke_loc_veg, gls_loc_veg
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vegu
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vegv 

# include "set_bounds.h"
!
      DO k=1,N(ng)-1
        DO j=Jstr,Jend
          DO i=Istr,Iend
            tke_veg(i,j,k)=Inival
            gls_veg(i,j,k)=Inival
          END DO 
        END DO 
      END DO 

      VEG_LOOP: DO iveg=1,NVEG
        DO k=1,N(ng)-1
          DO j=Jstr,Jend
            DO i=Istr,Iend
!
!-----------------------------------------------------------------------
!  Additional turbulence generated by the vegetation = 
!  work spent by the fluid against the plants (in m3/s3)
!-----------------------------------------------------------------------
!
              wku1=ru_loc_veg(i,j,k,iveg)*u(i,j,k,nstp)
              wku2=ru_loc_veg(i,j,k+1,iveg)*u(i,j,k+1,nstp)
              wku3=ru_loc_veg(i+1,j,k,iveg)*u(i+1,j,k,nstp)
              wku4=ru_loc_veg(i+1,j,k+1,iveg)*u(i+1,j,k+1,nstp)
              wku=0.25_r8*(wku1+wku2+wku3+wku4)
              wkv1=rv_loc_veg(i,j,k,iveg)*v(i,j,k,nstp)
              wkv2=rv_loc_veg(i,j,k+1,iveg)*v(i,j,k+1,nstp)
              wkv3=rv_loc_veg(i,j+1,k,iveg)*v(i,j+1,k,nstp)
              wkv4=rv_loc_veg(i,j+1,k+1,iveg)*v(i,j+1,k+1,nstp)
              wkv=0.25_r8*(wkv1+wkv2+wkv3+wkv4)
              tke_loc_veg=wku+wkv
!
!-----------------------------------------------------------------------
!  Dissipation due to vegetation
!  Dissipation time-scale for free turbulence
!-----------------------------------------------------------------------
!
              taufree=tke(i,j,k,nstp)/gls(i,j,k,nstp)  
!
# ifdef VEG_FLEX
!
!-----------------------------------------------------------------------
!  Equivalent thickness:horizontal projection of the bending plant
!-----------------------------------------------------------------------
! 
              eqvegT=plant(i,j,iveg,pthck)+sin(bend(i,j,iveg))*         &
     &                                       plant(i,j,iveg,phght)
# else
              eqvegT=plant(i,j,iveg,pthck)
# endif                 
!
!-----------------------------------------------------------------------
!  Solidity:cross-sectional area of a plant the number of plants per m2
!-----------------------------------------------------------------------
!
              solid=plant(i,j,iveg,pdiam)*eqvegT*plant(i,j,iveg,pdens)
!
!-----------------------------------------------------------------------
!  Eddies typical size constrained by distance in between the plants
!-----------------------------------------------------------------------
!
              L=cl_veg*((1.0_r8-MAX(solid,1.0_r8))/                     &
     &                 plant(i,j,iveg,pdens))**one_half
!
!-----------------------------------------------------------------------
!  Dissipation time-scale of eddies in between the plants
!-----------------------------------------------------------------------
!
              tauveg=(L**2.0_r8/(cmu**2.0_r8*tke_loc_veg))**one_third
!           
!-----------------------------------------------------------------------
!  Effective dissipation time-scale
!-----------------------------------------------------------------------
! 
              taueff=MIN(taufree,tauveg)
              gls_loc_veg=gls_c2(ng)*tke_loc_veg/taueff
!
!-----------------------------------------------------------------------
!  Add the tke and gls changes from all vegetation types
!-----------------------------------------------------------------------
! 
              tke_veg(i,j,k)=tke_loc_veg + tke_veg(i,j,k)
              gls_veg(i,j,k)=gls_loc_veg + gls_veg(i,j,k)
            END DO 
          END DO
        END DO 
      END DO VEG_LOOP
!
      RETURN
      END SUBROUTINE vegetation_turb_tile
#endif
      END MODULE vegetation_gls_corstep
