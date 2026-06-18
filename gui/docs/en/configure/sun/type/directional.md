---
title: "Directional Sun"
---

The directional sun method creates a plane where rays are randomly generated and emitted towards the scene geometry. The sun plane bounds are determined by projecting scene geometry bounding boxes corner points to a plane with a normal defined by sun's azimuth and elevation angles. For the CPU runners, all rays are emitted parallel to sun plane normal and sunshape ray directional deviations are applied at the ray's first surface interaction. For the GPU runner, sunshape ray directional deviations are applied at when the ray is emitted; thereby requiring the sun plane bounds to include a "buffer" to account for rays that can be generated near the sun plane boundary and intersect scene geometry.

A directional sun position is defined by azimuth and elevation angles. Azimuth angle is bounded between zero and 360 degrees where zero degrees is due north and increases clockwise (i.e., 90 degrees due east, 180 degrees due south, and 270 degrees due west). Elevation angle is bounded between -90 and 90 degrees where 90 degrees is directly overhead (i.e., parallel to the z-axis) and zero degrees is on the "horizon" (i.e., parallel to the x-y plane).

Sun's azimuth and elevation angles can be directly input or they can be calculated using a number of built-in solar position calculators (i.e., "Set by Calculator").

SolTrace supports five options for defining the sun's position:
1. Duffie & Beckman
2. SOLPOS 2.0
3. Solar Position Algorithm (SPA)
4. Legacy (from orginal SolTrace)
5. Custom (define azimuth and elevation angle)
