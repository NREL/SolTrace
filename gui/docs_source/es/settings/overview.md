---
title: "Resumen"
---

SolTrace es una aplicacion de trazado optico de rayos para modelar sistemas de energia solar de concentracion (CSP) y evaluar su rendimiento optico. Fue desarrollada en el National Laboratory of the Rockies (NLR) para estudiar sistemas demasiado complejos para herramientas analiticas mas simples, y tambien puede usarse en muchas tareas generales de modelado optico.

Use SolTrace cuando necesite describir una escena mediante fuentes de rayos, materiales opticos, geometria de superficies y elementos organizados en etapas, y luego trazar como se mueven los rayos por esa escena. La GUI esta organizada alrededor de ese flujo de trabajo: cargar o crear una escena, configurar la fuente de rayos y el modelo optico, ejecutar un trazado e inspeccionar los rayos, intersecciones, mapas de flujo y exportaciones resultantes.

SolTrace es especialmente util cuando la geometria, bloqueo, sombreado, reflexion, transmision, refraccion, error de pendiente, error de especularidad o intercepcion del receptor afectan el resultado. En lugar de depender solo de supuestos idealizados, puede construir directamente el sistema optico y evaluar como interactuan los rayos con los elementos configurados.

La aplicacion actual agrega una interfaz moderna en Qt alrededor del motor de SolTrace. Incluye inspeccion interactiva de escenas 3D, editores reutilizables de materiales y geometria, multiples motores de trazado, gestion de resultados, analisis de flujo, herramientas de exportacion y una interfaz de scripts para automatizacion.

Razones comunes para usar SolTrace incluyen:

- Estudiar heliostatos, platos, canaletas, receptores, optica secundaria u otros componentes CSP.
- Comparar configuraciones opticas antes de confirmar un diseno.
- Investigar trayectorias e intersecciones de rayos cuando un resultado no coincide con lo esperado.
- Generar mapas de flujo y datos de rayos trazados para analisis de ingenieria posterior.
- Automatizar estudios repetidos con scripts o flujos de trabajo externos.

Para usuarios que prefieren estudios basados en codigo, el proyecto SolTrace tambien ofrece soporte de API de Python mediante `pysoltrace`. La GUI esta pensada para configuracion interactiva de modelos, inspeccion, depuracion y analisis.
