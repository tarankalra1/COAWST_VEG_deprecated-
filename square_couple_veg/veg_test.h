/*
** svn $Id: inlet_test.h 838 2008-11-17 04:22:18Z jcwarner $
*******************************************************************************
** Copyright (c) 2002-2008 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Vegetation Test Case, waves-ocean (SWAN/ROMS) two-way coupling.
**
** Application flag:   VEG_TEST
** Input script:       coupling_veg_test.in
**                     swan_veg_test.in
**                     ocean_veg_test.in
**                     sediment_veg_test.in
*/

#define ROMS_MODEL
#define SWAN_MODEL
#define MCT_LIB

#define UV_VIS2
#define MIX_S_UV
#define MASKING
#define UV_ADV
#undef  UV_COR
#define TS_MPDATA
#define DJ_GRADPS

#define SOLVE3D
#define SPLINES

#undef  WEC_MELLOR
#define WEC_VF
#ifdef WEC_VF
# define BOTTOM_STREAMING
# define WDISS_WAVEMOD
#endif
#define UV_KIRBY

#undef ANA_INITIAL
#define ANA_FSOBC
#define ANA_M2OBC

/* define only one of the following */
#undef UV_LOGDRAG
#define SSW_BBL
#ifdef SSW_BBL
# define SSW_CALC_ZNOT
#endif

#define VEGETATION 
# ifdef VEGETATION 
#  undef ANA_VEGETATION 
#  define VEG_RHS
#  define VEG_SWAN_COUPLING
#  define VEG_STREAMING
# endif

#ifdef SOLVE3D
# define GLS_MIXING
# ifdef GLS_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
# endif
# define SEDIMENT
# ifdef SEDIMENT
#  define SUSPLOAD
#  undef  BEDLOAD_SOULSBY
#  undef  BEDLOAD_MPM
#  undef SED_MORPH
# endif
# if defined SEDIMENT || defined SG_BBL || defined MB_BBL || defined SSW_BBL
#  undef ANA_SEDIMENT
# endif
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_BPFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define ANA_SPFLUX
# define ANA_SRFLUX
#endif

#undef AVERAGES
#undef STATIONS
#undef DIAGNOSTICS_UV
