# MATPredict logo concepts (fable)

Three concepts, each as a square mark (`<name>-mark.svg`, viewBox 0 0 256 256, no text)
and a horizontal lockup (`<name>-lockup.svg`, mark + wordmark). `preview.html` shows every
file at 192/96 px, 48 px, 32 px and 16 px on white and on #17171a.

All files are plain SVG: no external references, scripts, filters, gradients or raster.
Every mark is a handful of paths; the lockups are ~4 KB because the wordmark is outlines.

## Wordmark

"MATPredict" is traced to paths (HarfBuzz shaping, so kerning is baked in) from Avenir Next,
the geometric sans that ships with macOS: **MAT** in Bold, **Predict** in Medium, so the
gene-name part carries the weight. Because the letters are outlines, the lockups need no
font at all, anywhere. The wordmark is filled with `currentColor` and the root `<svg>`
carries `color="#1d1d1f"`, so:

- used as `<img>` or a file, the wordmark is #1d1d1f (light theme);
- inlined in the site's HTML, a CSS rule like `.dark svg { color: #f2f2f4 }` flips it;
- for `<img>` use on a dark background, use `lens-lockup-dark.svg` (dark tokens + #f2f2f4).

If a libre typeface is preferred for redistribution, the lockups can be regenerated from any
font file (the tracing is ~30 lines of fontTools + uharfbuzz; it was kept out of the repo);
Avenir Next is only the outline source, nothing links to it.

## fit — compatible idiomorphs

**Idea.** Two gene arrows from a gene map, one per mating type: the blue arrow docks into
the complementary notch of the magenta one, so "these two fit" reads without a word of
biology. Corners are rounded by stroking each shape in its own colour.

**Colours.** Blue #0a5da8, magenta #c2418c. Mono: fill the arrow solid and draw the socket
as an outline of the same colour (the site's own solid-vs-outline language).

**Weaknesses.** Out of context it can read as a "skip forward" or pennant glyph rather than
genes; it is the least obviously fungal of the three; the two mid-value colours sit close in
lightness, so on #17171a the blue is a little dim and the seam between the pieces vanishes
below 32 px (the shape still reads, the "docking" does not).

## cap — a mushroom whose spots are genes

**Idea.** The fly-agaric silhouette everyone recognises, with the spots replaced by gene
arrows: two green flank genes on the sides and the yellow called core gene on top, so the
cap carries the locus. Playful, and unmistakably fungal at every size.

**Colours.** Cap blue #0a5da8, stem sky #4fb3e8, flank spots green #2fc596 (the dark-theme
green, chosen because #0f8a5f disappears on the blue cap), core spot yellow #e8d44d. Mono:
solid silhouette with the three arrows knocked out (transparent).

**Weaknesses.** It says "fungus" much louder than "mating type", and a mushroom is a
basidiomycete fruiting body while most of the users work on yeasts and moulds; the gene
arrows only register as arrows above ~48 px (at 16 px it is a mushroom with a yellow dot);
flat-bottomed cap and no gills, so it is more icon than illustration.

## lens — find the locus (recommended)

**Idea.** A lens with the found locus inside, drawn in the site's own gene-map vocabulary:
green flank genes, the winning core MAT gene in blue, on the pale-yellow called span. It
shows the thing the tool does (search a genome, call the locus) and the mark on the site
matches the gene maps below it.

**Colours.** Lens and core gene blue #0a5da8, flanks green #0f8a5f, called span yellow
#e8d44d. `lens-mark-dark.svg` / `lens-lockup-dark.svg` use the dark tokens (#4c9be8,
#2fc596). `lens-mono.svg` is one colour via `currentColor`: lens and core gene solid,
flanks at 50 %, span at 18 %; it works as black on white or white on dark.

**Weaknesses.** A magnifier is common iconography, so the silhouette alone is not unique;
the interior is only readable from ~32 px up, and at 16 px it is a ring with a yellow-and-blue
dot (still distinct from a plain search icon, but only just); the "set aside" competing gene
had to be dropped from the interior because it turned the small sizes to mush.

## Recommendation

**lens.** It is the most legible of the three at every size and on both backgrounds, it has
a real one-colour form, it carries the product's actual drawing language rather than a
metaphor, and the "find" verb is the product. `cap` is the runner-up if the site wants to
lead with warmth rather than precision; `fit` is the strongest idea but the weakest read.
