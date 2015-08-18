/*
** svn $Id: inlet_test.h 838 2008-11-17 04:22:18Z jcwarner $
*******************************************************************************
** Copyright (c) 2002-2008 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Inlet Test Case, waves-ocean (SWAN/ROMS) two-way coupling.
**
** Application flag:   INLET_TEST
** Input script:       ocean_inlet_test.in
**                     coupling_inlet_test.in
**                     sediment_inlet_test.in
*/

#define ROMS_MODEL
#undef SWAN_MODEL
#undef MCT_LIB

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
#undef WEC_VF
#undef WDISS_WAVEMOD
#undef UV_KIRBY
#undef ANA_INITIAL
#define ANA_SMFLUX
#define ANA_FSOBC
#define ANA_M2OBC
#define WAVES_HEIGHT
#define WAVES_LENGTH
#define WAVES_DIR

/* define only one of the following 5 */
#undef  UV_LOGDRAG
#define  UV_QDRAG
#undef  MB_BBL
#undef  SG_BBL
#undef SSW_BBL
#ifdef SSW_BBL
# define SSW_CALC_ZNOT
#endif

#define VEGETATION 
# ifdef VEGETATION 
#  undef ANA_VEGETATION 
#  undef VEG_RHS
#  undef VEG_SWAN_COUPLING 
#  undef VEG_STREAMING 
#  define WAVE_THRUST_MARSH 
#   define ANA_MARSH_MASK
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
# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_BPFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define ANA_SPFLUX
# define ANA_SRFLUX
#endif

#define AVERAGES
#define STATIONS
#undef DIAGNOSTICS_UV
