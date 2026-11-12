# Decisions and recommendations on the open questions of the junction-panel design (#321)

Date 2026-11-12. The curator answered all six questions in section 5 of `docs/notes/2026-11-09_junction-panel-design.md` on 2026-11-12, accepting each of the recommendations below (Q2 wording as drafted; Q3 relabel the three now, check the other two later; Q5 controls approved, C. albicans source to be verified first, #343). The recommendations were also posted as a comment on PR #321. Follow-ups: #195 (hand-off for the reads pages), #343 (truth sets).

1. **Active mating type for switching (HO+) species.** No. For wild *Saccharomyces* the honest output is "silent cassettes: not typable" (or "a/alpha or switching"), not a frequency. That also means steps 4-7 of the design are not needed for v0.9.0.
2. **Site wording.** Approve the wording table in section 2 once it has been read. It is the curator's to edit; there is no evidence here that changes it.
3. **Relabel *K. marxianus*, *K. phaffii* and *S. kudriavzevii* to `silent_cassettes`.** Yes: the note's own measurements describe them as cassette yeasts. Re-check *H. anomala* and *M. pulcherrima* before relabelling them; their reads have not been looked at.
4. **Ratio test A first, junction panel B only if A fails the controls.** Yes. A needs no new reads where a panel already types. Ratio and junction typing are already moved to v0.10.0 (#322, #323).
5. **Truth sets.** Yes for *S. cerevisiae* BY4741/2/3 as controls (they stay off the site) and yes for *S. pombe* 972, 975 and 968. For the published *C. albicans* MTL zygosity table no source was found: two PubMed searches returned nothing, so none is named here. It needs a literature search by someone with access.
6. **What closes v0.9.0.** The 2026-11-09 scope split already answers it: v0.9.0 closes on the gate reasons and wording, the cause TSV and *C. fimbriata*; ratio and junction typing move to v0.10.0.
