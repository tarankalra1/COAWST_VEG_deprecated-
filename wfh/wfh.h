/*
** svn $Id: wfh.h 429 2009-12-20 17:30:26Z arango $
*******************************************************************************
** Copyright (c) 2002-2010 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Shore Face Planar Beach Test Case.
**
** Application flag:   WFH
** Input scripts:      ocean_wfh.in
**
*/
#define ROMS_MODEL
#define BIO_ESTUARY
#undef  ANA_BIOLOGY

#define BIO_SEDIMENT

#define ANA_SPECIR
!#define SPECTRAL_LIGHT
#define MOD_SWR_SPECT
#define CDOM_DEFAULT
!#define ALGAL_RESP
!#define SAV_MODEL
!#define DENITRIFICATION
!#define OXYGEN
!#define CARBON

#define NO_LBC_ATT

#define WET_DRY
#define OUT_DOUBLE
#define UV_ADV

#undef  TS_MPDATA_LIMIT
#define TS_U3HADVECTION
#undef  TS_SVADVECTION

#define DJ_GRADPS
#define SALINITY
#define SOLVE3D
#define SPLINES
#define CURVGRID

#undef  ANA_TOBC


#define MASKING

#define UV_QDRAG


#define BULK_FLUXES
#define SOLAR_SOURCE
#define LONGWAVE
#define ANA_RAIN
#define ANA_CLOUD

# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_BPFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define ANA_SPFLUX


# define GLS_MIXING
# if defined GLS_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
#  undef CRAIG_BANNER
#  undef CHARNOK
#  undef ZOS_HSIG
#  undef TKE_WAVEDISS
# endif

