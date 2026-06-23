---
title: "Sun Shape"
---

The sun shape is an extrinsic property that depends on local atmospheric conditions, 
including atmospheric depth, precipitable water, aerosol content, and incidence angle
of the sunlight through the atmosphere.

The apparent shape of the sun is more technically described as the profile of flux 
intensity as a function of angular displacement from the center of the sun's disc. In 
other words, a measurement device with a vanishingly narrow acceptance window that 
views the sun at the center of the solar disc will see maximum relative intensity. 
As the device sweeps away from the center point, the measured intensity will decrease 
until only ambient scattered light is visible. The profile of intensity measurement 
as a function of displacement angle is called "sun shape". 

Within SolTrace, sun shape refers to the angular (normalized) intensity distribution 
of light across the sun's disk. The sun shape distribution is used to determine the 
probability of a ray being emitted at the angle.

SolTrace supports five options for defining the sun shape:
1. Gaussian
2. Pillbox
3. Buie circumsolar ratio (CSR)
4. Limb-darkened
5. User-defined

The radius of the solar disc is typically in the range of 4.60-4.75 mrad 
(0.263°-0.273°), with 4.65 mrad offered as a typical value.