---
title: "Solar Position Algorithm (SPA)"
---

Solar Position Algorithm (SPA) is a high-precision solar position model developed by NLR to calculate the sun’s apparent position from date, time, and observer location. The algorithm uses detailed astronomical models that account for Earth orbital variations, nutation, aberration, and atmospheric refraction. SPA is widely used in CSP, photovoltaic, solar resource, and tracking system applications where high pointing accuracy is required ($\pm 0.0003$ degrees). Its primary advantages are sub-arcminute accuracy and broad industry acceptance, while disadvantages include greater computational complexity and higher execution cost compared to simpler engineering formulations such as Duffie & Beckman or SOLPOS.

References:
    https://midcdmz.nlr.gov/spa/
    
    I. Reda and A. Andreas, “Solar Position Algorithm for Solar Radiation Applications (Revised),” NREL/TP-560-34302, 15003974, Jan. 2008. doi: 10.2172/15003974.