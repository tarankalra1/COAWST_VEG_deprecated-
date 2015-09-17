# COAWST_VEG
This repository contains the COAWST version to compute the effects of
* Vegetation in ROMS/Nonlinear/Vegetation 
* Wave thrust calculation on marsh routine in ROMS/Nonlinear/Vegetation
* Submerged Aquatic Vegetation (SAV model) in ROMS/Nonlinear/Biology

--------------------------------------------------
* 1) Analytical inputs to plant properties 
* 2) NETCDF based inputs to plant properties
* 3) NETCDF and Analytical input to initial marsh boundary

-------------------------------------------------
 Turning on vegetation or wave thrust (Executable file ./coawstM)  
* #define vegetation 
* #ifdef vegetation
*  #ifdef ana_vegetation 
*  #ifdef veg_drag
*  #ifdef veg_turb 
 *  #ifdef veg_flex
*  #ifdef veg_stream
*  #ifdef veg_swan_coupling 
*  #ifdef wave_thrust_marsh
*  #ifdef ana_marsh_mask 
* #endif 

------------------------------------------------
Turning on SAV model (Executable file ./coawstS) 
* #ifdef BIO_FENNEL
* # define SPECTRAL_LIGHT 
* # define CARBON
* # define OXYGEN
* # define DENITRIFICATION
* # define BIO_SEDIMENT
* # define SAV_BIOMASS(new flags)
* # define ALGAL_RESP (new flags) 
* # define DIAGNOSTICS_BIO
* #endif

