!
!svn $Id: vegetation_mod.h 429 2015-06-10 17:30:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!================================================= John C. Warner =====!   
!================================================= Neil K. Ganju ======!   
!================================================= Alexis Beudin ======!   
!================================================= Tarandeep S. Kalra =!
!=======================================================================
!                                                                      !
!  Vegetation Model Kernel Variables:                                  !
!  =================================                                   !
!  NVEG          Number of vegetation types                            !
!  NVEGP         Number of vegetation array properties                 !
!  CD_VEG        Drag coefficient for each veg type                    ! 
!  E_VEG         Youngs modulus for each veg type                      !
!  VEG_MASSDEN   Mass density for each veg type                        !
!                                                                      ! 
!  Plant Property indices:                                             !
!  ======================                                              !
!  pdens         Density                                               !
!  phght         Height                                                !
!  pdiam         Diameter                                              !
!  pthck         Thickness                                             !
!  pupbm         Above ground biomass                                  !
!  pdwbm         Below ground biomass                                  !
!                                                                      !
!  Plant Property indices:                                             !
!  ======================                                              !
!  idvprp        Indices for storing plant properties                  ! 
!                                                                      !
!  Plant Property Output IDs:                                          !
!  ==========================                                          !
!  ipdens         Id to output plant density                           !
!  iphght         Id to output plant height                            !
!  ipdiam         Id to output plant diameter                          !
!  ipthck         Id to output plant thickness                         !
!  ipupbm         Id to output above ground biomass                    !
!  ipdwbm         Id to output below ground biomass                    !
!                                                                      !
!  Wave Thrust on Marsh Output:                                        !
!  ==========================                                          !
!  idTmsk        Masking for getting thrust due to waves at rho pts.   ! 
!  idTmax        Maximum thrust due to waves                           !
!  idTton        Tonelli masking based thrust due to waves             !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
      integer :: NVEG, NVEGP
      integer :: counter
      integer :: phght, pdens, pdiam, pthck
      integer :: ipdens,iphght,ipdiam,ipthck

#ifdef VEG_BIOMASS 
      integer :: pupbm, pdwbm   
      integer :: ipupbm, ipdwbm   
#endif 
      integer, allocatable :: idvprp(:)    
#ifdef WAVE_THRUST_MARSH 
      integer ::  idTmsk, idTmax, idTton 
#endif 
!
      real(r8), allocatable :: E_VEG(:,:)
      real(r8), allocatable :: CD_VEG(:,:)
      real(r8), allocatable :: VEG_MASSDENS(:,:)

      CONTAINS 
! 
      SUBROUTINE initialize_vegetation
!
      USE mod_param
!
      implicit none 
!
!     Setup property indices 
! 
       counter = 1
       pdens   = counter 
       counter = counter+1 
       phght   = counter
       counter = counter+1 
       pdiam   = counter
       counter = counter+1 
       pthck   = counter
#ifdef VEG_BIOMASS 
       counter = counter+1 
       pupbm   = counter
       counter = counter+1 
       pdwbm   = counter 
#endif 
       NVEGP = counter
       IF (.not.allocated(idvprp)) THEN
         allocate ( idvprp(NVEGP) )
       END IF
      RETURN
      END SUBROUTINE initialize_vegetation
