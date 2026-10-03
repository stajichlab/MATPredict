# Confidence rules
Status: decided

## Results
- Allele-absent rule replay (`results/2026-09-27_tier_rule_replay/NOTE.md`):
  unguarded variant B unsafe (108 strong-gene promotions).
- Implemented with whole-call guard 1e69622 (`results/2026-09-27_tier_rule_implemented/NOTE.md`):
  86 risers vs 89 predicted; homothallic calls excluded (dab2eea).
- Fallback identity floor (`results/2026-09-27_fallback_confidence_floor/NOTE.md`):
  identity does not separate good from bad fallback calls; no floor.
