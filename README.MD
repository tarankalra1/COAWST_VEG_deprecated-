# COAWST_VEG
This repository contains the COAWST version to compute the effects of vegetation (along with two test cases)
and also wave thrust calculation on marsh routine
* 1) Analytical inputs to plant properties 
* 2) NETCDF based inputs to plant properties
* 3) NETCDF and Analytical input to initial marsh boundary

-------------------------------------------------
 TURNING ON VEGETATION  
* #define vegetation 
* #ifdef vegetation
*  #ifdef ana_vegetation 
*  #ifdef veg_rhs 
*  #ifdef veg_turb 
*  #ifdef veg_flex
*  #ifdef veg_stream
*  #ifdef veg_swan_coupling 
*  #ifdef wave_thrust_marsh
*   #ifdef ana_marsh_mask 
* #endif 

