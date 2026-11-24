# ADORA2A Retrospective Virtual Screening & Curation Audit

A reproducible case study demonstrating structure selection, bioactivity curation, decoy-bias control, and honest retrospective evaluation on a Class A GPCR: the human **Adenosine $A_{2A}$ receptor** (`ADORA2A`, UniProt [`P29274`](https://www.uniprot.org/uniprotkb/P29274)).

## Why this project exists

Virtual screening pipelines and ML models frequently fail prospectively because of avoidable upstream artefacts:

1. **Wrong conformational state** — docking or evaluating agonists against inactive GPCR structures (or vice versa), or ignoring engineered fusion proteins / missing pocket loops.
2. **Uncurated bioactivity data** — mixing incompatible assay types ($K_i$ vs. $IC_{50}$, binding vs. functional), keeping duplicate salts, or dropping provenance (`assay_chembl_id`).
3. **Trivial benchmark leakage** — evaluating against unmatched decoys where a simple physicochemical baseline (MW / LogP) already separates actives from decoys.
4. **Misleading metrics** — reporting accuracy on imbalanced screens or treating virtual screening scores as binding affinities.

This repository turns those failure modes into an explicit, reproducible audit with clear go/no-go criteria.

## Scope & Structure

- `src/adora2a_audit/` — reusable modules for GPCRdb/PDB structure auditing, ChEMBL curation with provenance, property-matched decoy generation, and enrichment metrics (EF1%, EF5%, BEDROC, ROC-AUC, PR-AUC).
- `notebooks/` — walkthrough analysis notebook.
- `reports/` — executive summary for ML/data engineering collaborators and generated figures.
- `tests/` — unit tests for curation, deduplication, and enrichment metric calculations.
- `docs/PLAN.md` — implementation plan and explicit success/failure criteria.

## Status

Work in progress — see [docs/PLAN.md](docs/PLAN.md).
