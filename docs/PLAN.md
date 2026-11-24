# Plan de trabajo · ADORA2A Retrospective Audit

Estudio de caso reproducible sobre el receptor de adenosina $A_{2A}$ humano (`ADORA2A`, UniProt `P29274`, ChEMBL `CHEMBL251`).
Objetivo: demostrar en un repositorio limpio, rápido y verificable el criterio científico necesario para diseñar, auditar y evaluar una campaña retrospectiva de cribado virtual en GPCRs — y comunicarlo con claridad a un equipo de ingeniería de IA/datos.

---

## Principios y criterios explícitos de éxito / fallo (*go / no-go*)

- **Ningún dato sin procedencia:** toda actividad química conserva su `assay_chembl_id`, `standard_type`, `standard_relation` y `confidence_score`.
- **No mezclar peras con manzanas:** separar estrictamente ensayos de unión ($K_i$) de $IC_{50}$ / ensayos funcionales, y distinguir agonistas de antagonistas.
- **Control de sesgo del benchmark:** un modelo *dummy* basado solo en propiedades fisicoquímicas (MW, cLogP, TPSA, HBD, HBA, carga) debe dar un ROC-AUC cercano a 0,50 sobre el set de activos + decoys emparejados. Si supera ~0,60, el set de decoys está sesgado y se rehace.
- **Métricas de cribado honestas:** evaluar siempre **EF1%**, **EF5%**, **BEDROC ($\alpha = 20$)**, **ROC-AUC** y **PR-AUC** — nunca *accuracy* como métrica principal.
- **Criterio de fallo claro:** documentar explícitamente qué no funciona (p. ej. falta de correlación lineal entre score y $pK_i$, o pérdida de rendimiento al hacer split por andamio de Murcko).

---

## Fases de implementación

### Fase 1 · Auditoría de estructuras GPCR (`structures.py`)
- [ ] Consultar la API REST de **GPCRdb** (`/services/structure/protein/aa2ar_human/`) y metadatos de **RCSB PDB**.
- [ ] Clasificar estructuras por estado conformacional (`Active`, `Intermediate`, `Inactive`), resolución, método (X-ray / cryo-EM), ligando co-cristalizado (agonista vs. antagonista), proteínas de fusión (BRIL / T4L) y mutaciones termoestabilizantes.
- [ ] Generar tabla comparativa en `reports/structure_audit.md` + figura resumen (resolución vs. estado conformacional).
- [ ] **Decisión explícita:** elegir 1 estructura de referencia para antagonistas (inactiva) y 1 para agonistas (activa), justificando qué estructuras se descartan y por qué.

### Fase 2 · Curación de bioactividad y set retrospectivo (`curation.py` + `decoys.py`)
- [ ] Descargar bioactividades de `CHEMBL251` desde la API de **ChEMBL** (con caché en disco en `data/cache/` para reproducibilidad offline).
- [ ] Filtrar por `confidence_score >= 8`, `standard_relation == '='`, unidad `nM`, y separar cohortes ($K_i$ vs. $IC_{50}$).
- [ ] Normalización química con RDKit: fragmento mayor (eliminación de sales), neutralización, deduplicación por `InChIKey` (calculando mediana de $pK_i$ y señalando discrepancias > 1 log unit entre ensayos).
- [ ] Anotar alertas estructurales (PAINS/BRENK) **con contexto**, verificando que fármacos/ligandos de referencia conocidos (p. ej. cafeína, istradefilina, ZM241385) no se descartan ciegamente.
- [ ] Generar **decoys emparejados en propiedades fisicoquímicas** (MW, cLogP, HBD, HBA, enlaces rotables, carga) y con baja similitud topológica (Tanimoto ECFP4 < 0,35 respecto a cualquier activo).
- [ ] Verificar con un clasificador fisicoquímico *dummy* (LogisticRegression / RandomForest solo con descriptores 1D) que el benchmark no tiene fuga trivial.

### Fase 3 · Benchmark retrospectivo y métricas honestas (`metrics.py` + `benchmark.py`)
- [ ] Implementar y testear (`tests/test_metrics.py`) las métricas: **EF1%**, **EF5%**, **BEDROC ($\alpha=20$)**, **ROC-AUC**, **PR-AUC** y *accuracy* (para contrastar por qué engaña a tasa de hits del 1–2 %).
- [ ] Ejecutar comparativa reproducible:
  1. Baseline fisicoquímico (control negativo).
  2. Cribado 2D por similitud ECFP4 (usando ligandos co-cristalizados como consulta).
  3. Evaluación con **Murcko scaffold split** (para medir generalización a series químicas nuevas).
  4. Correlación score vs. $pK_i$ experimental en activos (demostrando que *score $\neq$ afinidad*).
- [ ] Generar figuras limpias en `reports/figures/`:
  - Curvas de enriquecimiento / ROC vs. PR.
  - Comparativa *accuracy* vs. EF1% / BEDROC.
  - Scatter score vs. $pK_i$ experimental con casos fallidos anotados.

### Fase 4 · Informe ejecutivo y cierre (`reports/executive_summary.md` + `README.md`)
- [ ] Redactar resumen ejecutivo de 1 página dirigido a ingenieros de IA/datos: pregunta terapéutica, criterios *go/no-go*, decisiones de curación, resultados honestos (incluyendo qué falló) y recomendaciones para un pipeline productivo.
- [ ] Completar `README.md` con las figuras embebidas y comandos exactos para reproducir todo el análisis.
- [ ] Verificar `pytest` y `ruff check .` en verde + commit limpio.
