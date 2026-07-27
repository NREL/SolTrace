---
title: "Overview"
---

SolTrace is an optical ray tracing application for modeling concentrating solar power (CSP) systems and evaluating their optical performance. It was developed at the National Laboratory of the Rockies (NLR) to support systems that are too complex for simpler analytical tools, while still being useful for many general optical modeling tasks.

Use SolTrace when you need to describe a scene in terms of ray sources, optical materials, surface geometry, and staged elements, then trace how rays move through that scene. The GUI is designed around that workflow: load or create a scene, configure the ray source and optical model, run a trace, and inspect the resulting rays, intersections, flux maps, and exports.

SolTrace is especially useful when geometry, blocking, shading, reflection, transmission, refraction, slope error, specularity error, or receiver intercept matter to the result. Instead of relying only on idealized assumptions, you can build the optical system directly and evaluate how rays interact with the configured elements.

The current application adds a modern Qt interface around the SolTrace engine. It includes interactive 3D scene inspection, reusable material and geometry editors, multiple tracing backends, result management, flux analysis, export tools, and a scriptable interface for automation.

Typical reasons to use SolTrace include:

- Studying heliostats, dishes, troughs, receivers, secondary optics, or other CSP components.
- Comparing optical configurations before committing to a design.
- Investigating ray paths and intersections when a result does not match expectations.
- Generating flux maps and traced ray data for downstream engineering analysis.
- Automating repeated studies with scripts or external workflows.

For users who prefer code-driven studies, the SolTrace project also provides Python API support through `pysoltrace`. The GUI is intended for interactive model setup, inspection, debugging, and analysis.
