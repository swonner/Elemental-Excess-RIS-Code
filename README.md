# Elemental-Excess-RIS-Code
This repository compiles code that can be used to calculate an elemental excess (atoms/nm<sup>2</sup>) from one dimensional composition profiles displaying elemental segregation. This code was primarily written in order to quantify radiation-induced segregation in proton-irradiated 316L stainless steel.

The key result from this code is the text file: "ELEMENTALEXCESSMinMaxData.txt" which reports the elemental excess for the segregation of a specific element in a line scan and the width of the integration profile utilized to calculate it. 

In this repository, you will find the core python code, which compiles all of the necessary functions utilized in the elemental excess algorithm. You will also find scripts for 4 different use cases: 
(1) A single line scan of a single element
(2) A single line scan with multiple elements
(3) Multiple line scans with a single element
(4) Multiple line scans with multiple elements

# Licence information

Elemental Excess RIS Code. Copyright (C) 2025, Sara K Wonner and Pascal Bellon

This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version. This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

# References
If you use this code, you are required to cite the following paper: 
- S. K. Wonner and P. Bellon, "Investigation of radiation-induced segregation at fully characterized coherent twin boundaries in proton-irradiated 316L stainless steel," Journal of Nuclear Materials, vol. 604, p. 155470, 2025/01/01/ 2025, doi: https://doi.org/10.1016/j.jnucmat.2024.155470.
  
This code and algorithm were developed to analyze the data from this paper. Furthermore, the data provided in the example folder of this repository was collected and reported in this manuscript. 

# Acknowledgments
This work is supported by the U.S. Department of Energy, Office of Nuclear Energy’s Nuclear Energy University Programs, under Award Number DE-NE-0008973.
