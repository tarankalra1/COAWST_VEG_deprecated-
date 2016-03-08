# COAWST_VEG  :seedling:
This repository contains the version of COAWST accounting for vegetation:
* Wave-flow-vegetation interaction model in ROMS/Nonlinear/Vegetation 
* Wave thrust on marsh face calculation routine in ROMS/Nonlinear/Vegetation
* SAV growth model in ROMS/Nonlinear/Biology

------------------------------------------------
C-preprocessing options in application header file (.h extension)
-------------------------------------------------

/* Turning on vegetation */
* #define VEGETATION 
* #ifdef VEGETATION
* # undef ANA_VEGETATION
* # define VEG_DRAG
* # ifdef VEG_DRAG
* #  define VEG_FLEX
* #  define VEG_TURB
* # endif
* # undef VEG_SWAN_COUPLING
* # ifdef VEG_SWAN_COUPLING
* #  define VEG_STREAMING
* # endif
* #endif

/* Turning on wave thrust */
* #define WAVE_THRUST_MARSH
* #ifdef WAVE_THRUST_MARSH
* # undef ANA_MARSH_MASK
* #endif

/* Turning on SAV model */
* #define ESTUARYBGC
* #ifdef ESTUARYBGC
* # define SPECTRAL_LIGHT
* # define CARBON
* # define OXYGEN
* # define DENITRIFICATION
* # define BIO_SEDIMENT
* # define SAV_BIOMASS
* # define ALGAL_RESP
* # define DIAGNOSTICS_BIO
* #endif

--------------------------------------------------
ANAlytical or NetCDF file for initial conditions
--------------------------------------------------
* plant properties and distribution
* marsh boundary
