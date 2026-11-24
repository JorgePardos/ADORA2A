"""Unit tests for chemical standardization, deduplication, and contextual alert auditing."""

from __future__ import annotations

from adora2a_audit.curation import (
    audit_reference_ligand_alerts,
    standardize_smiles,
)


def test_salt_stripping_and_neutralization_converge_to_same_inchikey():
    # ZM241385 free base vs dihydrochloride salt vs deprotonated phenolate
    zm_free = "Oc1ccc(cc1)CCNc1nc(N)n2c(n1)nc(n2)c1ccco1"
    zm_hcl = "Oc1ccc(cc1)CCNc1nc(N)n2c(n1)nc(n2)c1ccco1.Cl.Cl"
    zm_phenolate = "[O-]c1ccc(cc1)CCNc1nc(N)n2c(n1)nc(n2)c1ccco1.[Na+]"

    s1 = standardize_smiles(zm_free)
    s2 = standardize_smiles(zm_hcl)
    s3 = standardize_smiles(zm_phenolate)

    assert s1 is not None and s2 is not None and s3 is not None
    assert s1.inchikey == s2.inchikey == s3.inchikey
    assert "." not in s2.canonical_smiles
    assert "." not in s3.canonical_smiles


def test_invalid_or_inorganic_inputs_rejected_cleanly():
    assert standardize_smiles("") is None
    assert standardize_smiles("not_a_smiles") is None
    # Sodium chloride / tiny fragment (< 6 heavy atoms) rejected
    assert standardize_smiles("[Na+].[Cl-]") is None
    assert standardize_smiles("CCO") is None


def test_contextual_alert_audit_preserves_reference_drugs():
    """Verify that reference A2A ligands (Caffeine, Istradefylline, ZM241385, Tozadenant,
    Vipadenant, Imaradenant) are audited rather than blindly dropped: Phase 3 Tozadenant
    triggers PAINS (anil_di_alk_C) and Phase 2 Vipadenant/Imaradenant trigger Brenk alerts!
    """
    df = audit_reference_ligand_alerts()
    names = set(df["name"])
    assert {"Caffeine", "Istradefylline", "ZM241385", "Tozadenant", "Vipadenant", "Imaradenant"}.issubset(names)

    toza = df[df["name"] == "Tozadenant"].iloc[0]
    assert bool(toza["has_pains"]) is True

    imara = df[df["name"] == "Imaradenant"].iloc[0]
    assert bool(imara["has_brenk"]) is True

    caff = df[df["name"] == "Caffeine"].iloc[0]
    assert bool(caff["has_pains"]) is False
