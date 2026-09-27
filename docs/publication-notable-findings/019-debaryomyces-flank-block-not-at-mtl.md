# 019. In Debaryomyces the PAP1-OBP1-PIK1 block is not at the MTL locus

- **Category:** biology; assembly-or-annotation artefact (a spurious call it caused)
- **Status:** verified (CBS767 chromosome coordinates)
- **Lineage:** Ascomycota, Serinales (Debaryomycetaceae)

## Summary
In most Serinales the MTL flanks PAP1, OBP1 and PIK1 sit at the MTL. In D.
hansenii CBS767 the block is about 0.7 Mb away, so a weak MTLA2 fragment next to
it produced spurious second calls.

## Evidence
- CBS767 `NC_006047.2`: PAP1-OBP1-PIK1 at 0.88 Mb; a1/a2/alpha1 locus at
  1.59 Mb. The spurious call rested on a 29% MTLA2 fragment 8 kb from the
  block and turned 42 `homothallic` genotypes into `a+homothallic`.
- Flank-carried calls in the Serinales-wide scan: 131 of 2,647; 52 are
  Debaryomyces (46 D. hansenii).
- Source: `docs/notes/2026-09-26_polish-cap-measured-and-serinales-scan.md`;
  annotation report D9.

## Method that found it
Flank-carried audit of the Serinales-wide scan (`flank_carried_audit.py`).

## Verification done / still open
Done: the flank-carried rule withholds all 52 (commits `3684c60`, `b0898d9`).

## Limits
Checked in CBS767 only; other genera not mapped.

## Related
Entries 004, 007.
