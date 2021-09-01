This package includes four programs that can be used to do:

1. frequency-domain acoustic forward modeling
2. frequency-domain viscoacoustic forward modeling
3. frequency-domain acoustic reverse time migration (FA-RTM)
4. frequency-domain Q-compensated reverse time migration (FQ-RTM).

Both the forward modeling and RTM are based on frequency-domain finite-difference algorithm.
All the codes are written with MATLAB language.
The codes were tested on Windows 10, with Matlab R2014a.

------------------------------------------------------------------
Contents
------------------------------------------------------------------
The package comprises these functions
*) frequency_modeling_record.m       : frequency domain acoustic and viscoacoustic forward modeling for generating the frequency-domain records
*) frequency_modeling_wavefield.m   : frequency domain acoustic and viscoacoustic forward modeling for generating the frequency-domain wavefields
*) frequency_RTM.m        : frequency domain acoustic reverse-time migration (FA-RTM)
*) frequency_QRTM.m     : frequency-domain Q-compensation reverse-time migration (FQ-RTM) 
*) geometry.m	      : define acquisition geometry
*) impedance_matrix.m   : generate a impedance matrix use an optimal 9-point, finite-differnce method
*) PML.m		      : calculate the reciprocal of PML coefficients 
*) model_pml.m	      : extend model parameters with PML boundary 
*) fricker.m                     : creates a causal Ricker wavelet in the frequncy domain 
*) avg_amp_spectra.m    : average amplitude spectrum of a 2D plofile

------------------------------------------------------------------
Quick start
------------------------------------------------------------------
To reproduce the Figures shown in this manuscript£¬
1. Download this packages.
2. Unzip the packages.
3. Open Matlab and change to the current folder.
4. Execute the Matlab codes, e.g. figure2,  figure5_10, figure11_12, test_gas_chimney

------------------------------------------------------------------
Copyright (c) Xiong Ma 
------------------------------------------------------------------
All rights reserved.

This program is free software: you can redistribute it and/or modify it under the terms of 
the GNU General Public License as published by the Free Software Foundation, 
eitherversion 3 of the License, or (at your option) any later version. 
This program is distributed in the hope that it will be useful, 
but WITHOUT ANYWARRANTY; without even the implied warranty of MERCHANTABILITY 
or FITNESSFOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with thisprogram.

------------------------------------------------------------------
Feedback
------------------------------------------------------------------

If you have any comment, suggestion, or question, please do
contact Xiong Ma at: mx_geophysics@126.com
