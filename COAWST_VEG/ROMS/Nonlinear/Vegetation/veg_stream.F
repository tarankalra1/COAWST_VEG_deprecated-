#include "cppdefs.h"
 
       MODULE veg_stream_mod
#if defined VEGETATION && defined VEG_STREAM 
!
!svn $Id: wave_thrust_marsh.F 429 2015-04-20 17:30:26Z arango $
!=======================================================================
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                   Alexander F. Shchepetkin   !
!================================================John C. Warner=========
!================================================Neil K. Ganju  ========
!==============================================Tarandeep S. Kalra=======
!                                                                      ! 
!  References:                                                         !   
!                                                                      !
!=======================================================================
!                                                                      !
!=======================================================================

      implicit none

      PRIVATE
      PUBLIC  :: veg_stream

      CONTAINS
!
!***********************************************************************
      SUBROUTINE veg_stream_cal (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_vegarr
      USE wec_streaming_mod

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
      CALL veg_stream_tile  (ng, tile,                                  &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                         GRID(ng) % angler,                       &
# ifdef SOLVE3D
     &                         GRID(ng) % z_w,                          &
# endif
     &                         FORCES(ng) % Hwave,                      &
     &                         FORCES(ng) % Dwave,                      &
     &                         FORCES(ng) % Lwave,                      &
     &                         VEG(ng) % dissip_veg,                    &
     &                         VEG(ng) % EWD_veg) 
              
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 16)
# endif
      RETURN
      END SUBROUTINE veg_stream_cal 

!***********************************************************************
      SUBROUTINE veg_stream_tile  (ng, tile,                            &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              angler,                             &
# ifdef SOLVE3D
     &                              z_w,                                &
# endif
     &                              Hwave, Dwave, Lwave,                &
     &                              dissip_veg, EWD_veg)           
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars
      USE mod_vegetation
      USE mod_vegarr

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
      real(r8), intent(in) :: Hwave(LBi:,LBj:)
      real(r8), intent(in) :: Lwave(LBi:,LBj:)
      real(r8), intent(in) :: Dwave(LBi:,LBj:)
      real(r8), intent(inout) :: EWD_veg(LBi:,LBj:,NVEG)
# else
      real(r8), intent(in)  :: angler(LBi:UBi,LBj:UBj)
#  ifdef SOLVE3D
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,UBk)
#  endif
      real(r8), intent(in)  :: Hwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: Lwave(LBi:UBi,LBj:UBj)
      real(r8), intent(in)  :: Dwave(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: EWD_veg(LBi:UBi,LBj:UBj,NVEG)
# endif

!  Local variable declarations.
!
      integer :: i,j
      real(r8) :: cff, cff3
      real(r8), parameter :: sqrt2=SQRT(2.0_r8)
      real(r8), parameter :: ks=0.03_r8
      real(r8), parameter :: awd=1.0_r8
      real(r8), parameter :: KWDmax=200.0_r8
      real(r8), parameter :: eps = 1.0E-14_r8
      real(r8), parameter :: kDmax = 5.0_r8
      real(r8), parameter :: Lwave_min = 1.0_r8

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: kD
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waven
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wavenx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: waveny
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: KWD
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: owd
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: EWD
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: sigma

# include "set_bounds.h"
!----------------------------------------------------------------------
!----------Executing the code------------------------------------------
!----------------------------------------------------------------------
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
          Dstp(i,j)=z_w(i,j,N(ng))-z_w(i,j,0)
!----------------------------------------------------------------------
!  Compute wave amplitude (0.5*Hrms), wave number, intrinsic frequency.
!----------------------------------------------------------------------
          waven(i,j)=2.0_r8*pi/MAX(Lwave(i,j),Lwave_min)
          cff=1.5_r8*pi-Dwave(i,j)-angler(i,j)
          wavenx(i,j)=waven(i,j)*COS(cff)
          waveny(i,j)=waven(i,j)*SIN(cff)
          sigma(i,j)=MIN(SQRT(g*waven(i,j)*TANH(waven(i,j)*Dstp(i,j))),2.0_r8)
!----------------------------------------------------------------------
!  Compute wave celerity and nonlinear water depth 
!----------------------------------------------------------------------
          kD(i,j)=MIN(waven(i,j)*Dstp(i,j)+eps,kDmax)
        END DO 
      END DO 

      VEG_LOOP: DO iveg=1,NVEG 
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
!----------------------------------------------------------------------
!  Metric for vertical bottom streaming distribution     
!----------------------------------------------------------------------
            cff=0.25_r8*sqrt2*sigma(i,j)*Hwave(i,j)/                    &
     &           (SINH(kD(i,j))+eps)

!----------------------------------------------------------------------
!  Wave friction factor (based on Soulsby, 1997).
!----------------------------------------------------------------------
            cff3=MIN(1.39_r8*(sigma(i,j)*(ks/30.0_r8)/                  &
     &               cff)**0.52_r8,0.2_r8)
 
#  ifdef WAVES_OCEAN
!----------------------------------------------------------------------
!   Note: Alexis - check if we need a local dissip_veg here 
!----------------------------------------------------------------------
            EWD_loc_veg=Dissip_veg(i,j)
#  else
            EWD_loc_veg=(0.5_r8/sqrt(pi))*cff3*(cff**3.0_r8)
#  endif
          END DO
        END DO

        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
            EWD_veg(i,j,iveg)=EWD_loc_veg+EWD_veg(i,j,iveg)
          END DO 
        END DO 

      END DO     

      END SUBROUTINE wave_thrust_marsh_tile
#endif
      END MODULE wave_thrust_marsh_mod
