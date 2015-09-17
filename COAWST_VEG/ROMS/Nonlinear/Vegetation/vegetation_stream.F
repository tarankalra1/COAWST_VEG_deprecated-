#include "cppdefs.h"
 
       MODULE vegetation_stream_mod
#if defined VEGETATION && defined VEG_STREAMING 
!
!svn $Id: vegetation_stream.F 429 2015-04-20 17:30:26Z arango $
!=======================================================================
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                   Alexander F. Shchepetkin   !
!================================================John C. Warner=========
!================================================Neil K. Ganju  ========
!================================================Alexis Beudin  ========
!==============================================Tarandeep S. Kalra=======
!                                                                      ! 
!  References:                                                         !   
!                                                                      !
!=======================================================================
!                                                                      !
!=======================================================================

      implicit none

      PRIVATE
      PUBLIC  :: vegetation_stream_cal

      CONTAINS
!
!***********************************************************************
      SUBROUTINE vegetation_stream_cal (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
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
      CALL vegetation_stream_tile  (ng, tile,                           &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                         GRID(ng) % angler,                       &
# ifdef SOLVE3D
     &                         GRID(ng) % z_w,                          &
# endif
!# ifdef WAVES_HEIGHT
     &                         FORCES(ng) % Hwave,                      &
!# endif
!# ifdef WAVES_DIR
     &                         FORCES(ng) % Dwave,                      &
!# endif
!# ifdef WAVES_LENGTH
     &                         FORCES(ng) % Lwave,                      &
!# endif
     &                         VEG(ng) % plant,                         &
     &                         VEG(ng) % dissip_veg,                    &
     &                         VEG(ng) % Lveg,                          &
     &                         VEG(ng) % BWDXL_veg,                     &
     &                         VEG(ng) % BWDYL_veg)
              
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 16)
# endif
      RETURN
      END SUBROUTINE vegetation_stream_cal 

!***********************************************************************
      SUBROUTINE vegetation_stream_tile  (ng, tile,                     &
     &                              LBi, UBi, LBj, UBj,                 &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                              angler,                             &
# ifdef SOLVE3D
     &                              z_w,                                &
# endif
!# ifdef WAVES_HEIGHT
     &                              Hwave,                              &
!# endif
!# ifdef WAVES_DIR
     &                              Dwave,                              &
!# endif
!# ifdef WAVES_LENGTH
     &                              Lwave,                              &
!# endif
     &                              plant,                              &
     &                              dissip_veg, Lveg,                   &
     &                              BWDXL_veg, BWDYL_veg)  
!***********************************************************************
!
      USE mod_param
      USE mod_grid
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
      real(r8), intent(in) :: angler(LBi:,LBj:)
#  ifdef SOLVE3D
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
#  endif
!#  ifdef WAVES_HEIGHT
      real(r8), intent(in) :: Hwave(LBi:,LBj:)
!#  endif
!#  ifdef WAVES_LENGTH
      real(r8), intent(in) :: Lwave(LBi:,LBj:)
!#  endif
!#  ifdef WAVES_DIR
      real(r8), intent(in) :: Dwave(LBi:,LBj:)
!#  endif
      real(r8), intent(in) :: plant(LBi:,LBj:,:,:)
      real(r8), intent(in) :: dissip_veg(LBi:,LBj:)
      real(r8), intent(in) :: Lveg(LBi:,LBj:,:)
      real(r8), intent(inout) :: BWDXL_veg(LBi:,LBj:)
      real(r8), intent(inout) :: BWDYL_veg(LBi:,LBj:)
# else
      real(r8), intent(in)  :: angler(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,UBk)
#  endif
!#  ifdef WAVES_HEIGHT
      real(r8), intent(in)  :: Hwave(LBi:UBi,LBj:UBj)
!#  endif
!#  ifdef WAVES_LENGTH
      real(r8), intent(in)  :: Lwave(LBi:UBi,LBj:UBj)
!#  endif
!#  ifdef WAVES_DIR
      real(r8), intent(in)  :: Dwave(LBi:UBi,LBj:UBj)
!#  endif
      real(r8), intent(in) :: plant(LBi:UBi,LBj:UBj,NVEG,NVEGP)
      real(r8), intent(in) :: dissip_veg(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Lveg(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: BWDXL_veg(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: BWDYL_veg(LBi:UBi,LBj:UBj)
# endif

!  Local variable declarations.
!
      integer :: i, j, k, iveg
      real(r8) :: cff1, cff2, cff3, cff4, cff5
      real(r8) :: EWD_veg
      real(r8), parameter :: sqrt2=SQRT(2.0_r8)
      real(r8), parameter :: ks=0.03_r8
      real(r8), parameter :: awd=1.0_r8
      real(r8), parameter :: KWDmax=200.0_r8
      real(r8), parameter :: eps = 1.0E-14_r8
      real(r8), parameter :: kDmax = 5.0_r8
      real(r8), parameter :: Lwave_min = 1.0_r8

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Dstp
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: kD
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waven
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wavenx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waveny
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: sigma
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: osigma

# include "set_bounds.h"
!----------------------------------------------------------------------
!----------Executing the code------------------------------------------
!----------------------------------------------------------------------
!     
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
          Dstp(i,j)=z_w(i,j,N(ng))-z_w(i,j,0)
!----------------------------------------------------------------------
!  Compute wave amplitude (0.5*Hrms), wave number, intrinsic frequency.
!----------------------------------------------------------------------
          waven(i,j)=2.0_r8*pi/MAX(Lwave(i,j),Lwave_min)
          cff2=1.5_r8*pi-Dwave(i,j)-angler(i,j)
          wavenx(i,j)=waven(i,j)*COS(cff2)
          waveny(i,j)=waven(i,j)*SIN(cff2)
          sigma(i,j)=MIN(SQRT(g*waven(i,j)*TANH(waven(i,j)*             &
     &                                  Dstp(i,j))),2.0_r8)
          osigma(i,j)=1.0_r8/sigma(i,j)
!----------------------------------------------------------------------
!  Compute wave celerity and nonlinear water depth 
!----------------------------------------------------------------------
          kD(i,j)=MIN(waven(i,j)*Dstp(i,j)+eps,kDmax)
        END DO 
      END DO 
!
      DO k=1,N(ng)
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
!----------------------------------------------------------------------
!  Metric for vertical bottom streaming distribution     
!----------------------------------------------------------------------
            cff3=0.25_r8*sqrt2*sigma(i,j)*Hwave(i,j)/                   &
     &           (SINH(kD(i,j))+eps)

!----------------------------------------------------------------------
!  Wave friction factor (based on Soulsby, 1997).
!----------------------------------------------------------------------
            cff4=MIN(1.39_r8*(sigma(i,j)*(ks/30.0_r8)/                  &
     &               cff3)**0.52_r8,0.2_r8)
! 
#  ifdef WAVES_OCEAN
!----------------------------------------------------------------------
!   Note: Alexis - check if we need a local dissip_veg here 
!   Also Lveg is for 1 veg type only 
!----------------------------------------------------------------------
            EWD_veg=dissip_veg(i,j)
#  else
            EWD_veg=(0.5_r8/sqrt(pi))*cff4*(cff3**3.0_r8)
#  endif
            cff5=EWD_veg*osigma(i,j)*Lveg(i,j,k)
            BWDXL_veg(i,j)=cff5*wavenx(i,j)
            BWDYL_veg(i,j)=cff5*waveny(i,j)
!
          END DO
        END DO
      END DO    
!
      END SUBROUTINE vegetation_stream_tile
#endif
      END MODULE vegetation_stream_mod
