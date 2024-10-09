This module is for developing a 3D interface (e.g., Top of the crust , Moho etc)

make_interface_txt2xlx.m : Script for making interface for imported interface as txt/xlx file
			  -- Preferred Format of all txt/xlx file - Longitude||Latitude||Depth/TWTT
			  -- Also can input MCS data to constrain any part
			  -- This script is using griddata, so it will take time to interpolate
			  -- Script is prepared for Cascadia subduction zone 
					where lon is decreasing from east to west(min value)
                                             and lat is decreasing from north to south(min value)

calc_complex_thickness.m : Calculate a crust of seemingly same thickness from the moho information of interface structure 


modify_interfaces.m : Script to modify interfaces and make four different ones
			-- modified interfaces will have diffrent depths
			-- Input is an existing int_3Dmat structure
			-- 'pct1' and 'pct2' are two percentage inputs that control the depth modification
			-- 'pct1' is the starting position, left/west 
				of which will be unmodified (typically constrained with some data)
			-- 'pct2' is the ending position of transition and starting postion of depth
			-- apply a low filter number (input is 'flt') to smooth corners of transition
 