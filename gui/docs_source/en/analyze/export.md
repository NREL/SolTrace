---
title: "Export"
---

Export writes the selected result data to a file for external analysis. Choose the result and output format that matches the workflow you plan to use outside SolTrace.

Flux map mesh export writes a SolTrace-extended OBJ file with normal vertices, texture coordinates, normals, and faces. Each face is followed by `fa` for face area and `fv` for the face-bin ray count. Some strict OBJ readers may warn about those custom records.

Large traces can produce large exports. Reduce the ray count or filter the selected result when you only need a focused subset of the data.
