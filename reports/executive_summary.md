# Resumen Ejecutivo · Auditoría Retrospectiva de Cribado Virtual en ADORA2A

**Audiencia:** Equipos de Ingeniería de IA, Ciencia de Datos y Química Computacional  
**Diana biológica:** Receptor humano de adenosina $A_{2A}$ (`ADORA2A`, UniProt `P29274`, ChEMBL `CHEMBL251`, GPCR de Clase A)  
**Pregunta terapéutica y de ingeniería:** *Antes de escalar un pipeline de cribado virtual o entrenar un modelo ML sobre decenas de millones de compuestos para antagonistas de $A_{2A}$ (inmuno-oncología / enfermedad de Parkinson), ¿qué decisiones estructurales y de curación evitan fugas de datos, y qué métodos aportan enriquecimiento temprano real cuando eliminamos los artefactos triviales?*

---

## 1. Puertas de calidad explícitas (*Go / No-Go*)

| Puerta de auditoría | Criterio cuantitativo *Go / No-Go* | Resultado observado | Estado |
| :--- | :--- | :--- | :---: |
| **1. Integridad del bolsillo 3D** | Cero mutaciones ingenieriles en residuos de contacto ortostérico/alostérico; estado conformacional coherente con la modalidad (`Inactive` para antagonistas). | Se descartan `5NM4` (1,70 Å) y `5IU4` (1,72 Å) —y el **66,0 % (`124/188`)** del PDB de $A_{2A}$— por portar mutaciones `StaR2` (`T88A`, `S277A`) **dentro del bolsillo**. Se elige **`4EIY` (1,80 Å, bolsillo 100 % *wild-type* + $\text{Na}^+$ alostérico)** y **`9EE8` (2,63 Å, activa + $G_s$)** para agonistas. | **GO** |
| **2. Trazabilidad y pureza de ensayo** | Conservar `assay_chembl_id`, `confidence_score >= 8`, relación exacta `'='`, separar $K_i$ de unión (`Ki_B`) de $IC_{50}$ funcional y excluir agonistas ribosídicos del cribado inactivo. | De `12.105` registros brutos $\rightarrow$ `5.901` medidas exactas `Ki_B` $\rightarrow$ **`4.972` estructuras padre únicas (`InChIKey`)**, marcando `44` compuestos con discrepancia inter-ensayo $> 1$ unidad log y excluyendo `546` agonistas. | **GO** |
| **3. Alertas estructurales con contexto** | Auditar PAINS/Brenk sin borrado ciego de quimiotipos clínicos validados. | Un filtro ciego PAINS/Brenk habría eliminado **3 antagonistas clínicos co-cristalizados de Fase 2/3**: `Tozadenant` (PAINS `anil_di_alk_C`), `Vipadenant` (Brenk `aniline`) e `Imaradenant` (Brenk `2-halo_pyridine`). | **GO** |
| **4. Control de sesgo 1D en decoys** | Un clasificador 1D (`MW`, `cLogP`, `TPSA`, `HBD`, `HBA`, `RotBonds`, `Charge`) debe obtener $\text{ROC-AUC} \le 0,60$ en validación cruzada. | Decoys aleatorios de ChEMBL sufren fuga masiva (**AUC = `0,81` LogReg / `0,78` RF**, impulsada por la alta `TPSA` de ligandos $A_{2A}$). Tras emparejamiento 7D + disimilitud ECFP4 $< 0,35$, cae a **`0,48` / `0,35`**. | **GO** |

---

## 2. Resultados honestos del benchmark retrospectivo (1.000 compuestos · 2,5 % hits)

Sobre el mismo set libre de sesgo 1D (`25` antagonistas potentes diversos + `975` decoys emparejados en 7D; tasa basal de hits = `2,5 %`, *accuracy* trivial al predecir "todo inactivo" = **`97,5 %`**, máximo teórico $\text{EF1\%} = 40,0\times$):

| Nivel / Método | EF1% (IC 95%) | EF5% | BEDROC ($\alpha=20$) | ROC-AUC | PR-AUC | Accuracy @ Top 1% | Lectura crítica para ingeniería |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **1. Control 1D (RF fisicoquímico)** | `0,0x` `(0–0)` | `0,0x` | `0,026` | `0,346` | `0,019` | `96,5 %` | Confirma ausencia de fuga fisicoquímica trivial; fíjese en el engañoso `96,5 %` de *accuracy*. |
| **2. Baseline 2D (ECFP4 vs. `ZM241385`)** | **`36,0x`** `(16–40)` | **`9,6x`** | **`0,519`** | `0,725` | **`0,433`** | `98,3 %` | Recupera 9 de 10 hits en el Top 1%, pero depende de similitud topológica al quimiotipo consulta (`23,7x` al excluir análogos cercanos). |
| **3. Forma 3D (`ETKDGv3` + `USRCAT` + `O3A`)** | `8,0x` `(0–20)` | `6,4x` | `0,355` | `0,722` | `0,166` | `96,9 %` | **Mismo ROC-AUC global (`0,722` vs. `0,725`) que 2D**, pero **4,5× menos enriquecimiento en el Top 1%**: ejemplo real de por qué ROC-AUC no basta. |
| **4. Docking 3D (`AutoDock Vina` · `4EIY`)** | **`0,0x`** `(0–0)` | `0,8x` | `0,063` | `0,584` | `0,033` | `96,5 %` | **Resultado negativo claro:** una vez eliminado el sesgo de peso molecular/lipofilia, el docking rígido sin aguas estructurales fracasa en el Top 1%. |

---

## 3. Qué NO funcionó y por qué (diagnóstico mecanicista)

1. **El docking rígido deshidratado (`AutoDock Vina`) es superado claramente por métodos 2D y de forma 3D:**
   - En el re-docking cognado de `ZM241385` sobre `4EIY` sin aguas cristalográficas, la pose invertida (`-9,08 kcal/mol`) y la pose cristalográfica nativa (rango 3, `-8,93 kcal/mol`, oxígeno furánico a `0,41 Å` del cristal) difieren en solo **`0,15 kcal/mol`** (dentro del ruido de la función de score). Al retirar las aguas conservadas (`HOH 2503–2505`), decoys heteroaromáticos hidrofóbicos de igual `MW` y `TPSA` ocupan el vestíbulo extracelular con energías de `-10` a `-12 kcal/mol` (falsos positivos), mientras que antagonistas reales voluminosos chocan con la conformación rígida inducida por `ZM241385` (falsos negativos).
2. **Score de cribado $\neq$ afinidad experimental ($pK_i$):**
   - Evaluando **60 ligandos con $K_i$ experimental** a lo largo de 5 órdenes de magnitud ($pK_i \in [4,8;\; 11,0]$), el score de AutoDock Vina explica apenas un **2,6 % de la varianza** ($r = +0,161$, $R^2 = 0,026$, $p = 0,22$) y la similitud 2D ECFP4 un **6,0 %** ($r = +0,246$, $R^2 = 0,060$). Las funciones rápidas de cribado sirven para **priorizar subconjuntos enriquecidos**, nunca como sustituto directo de $\Delta G_{\text{bind}}$ (que requiere FEP/ABFE o ensayos biofísicos).
3. **Caída al saltar de andamio químico (*Murcko scaffold split*):**
   - Evaluar modelos supervisados 2D mediante *random split* sobreestima el rendimiento porque reparte análogos metil/cloro de una misma serie entre entrenamiento y test; el particionado estricto por **andamio de Murcko** es obligatorio para estimar la tasa de acierto sobre series químicas nuevas.

---

## 4. Recomendaciones accionables para un pipeline productivo

1. **Bloqueo automático de estructuras mutadas en el bolsillo:** cruzar siempre GPCRdb + coordenadas atómicas antes de parametrizar un receptor (descartando mutaciones `StaR2` `T88A`/`S277A` y ajustando protonación a pH 7,4 en `Glu169`/`Asp52`).
2. **Arquitectura jerárquica en embudo (2D/3D rápido $\rightarrow$ *Ensemble Docking* con aguas explícitas $\rightarrow$ re-scoring físico):** filtrar primero por modelos 2D/farmacofóricos validados bajo *Murcko split* y forma 3D (`USRCAT`/`O3A`), y reservar el docking para un **ensemble de conformaciones inactivas** (`4EIY`, `5OLZ`, `6GT3`) **conservando las aguas estructurales puente** del bolsillo profundo y aplicando *pharmacophore tethering* sobre `Asn253^{6.55}` y `Phe168^{ECL2}`.
3. **Gobernanza de métricas en CI/CD de modelos:** prohibir *accuracy* y exigir siempre **EF1%**, **BEDROC ($\alpha=20$)** y **PR-AUC** con intervalos bootstrap junto con una puerta automática de **control negativo 1D ($\text{ROC-AUC} < 0,60$)** cada vez que se actualice un benchmark.
