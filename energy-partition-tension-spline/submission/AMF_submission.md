# Applied Mathematical Finance (Taylor & Francis) — submission & help

**Manuscript file:** `swap_curve_tension_AMF.tex` (coded for the T&F `interact` class)
**Publisher:** Taylor & Francis · **Submission system:** ScholarOne Manuscripts
**Why this journal:** Hagan & West's "Interpolation methods for curve construction" — which
this paper benchmarks against and cites — was published in AMF, so a curve-construction method
is squarely in scope and well-precedented. AMF has no tight word limit, so the full-length paper
fits.

---

## 1. What you must download (I cannot ship these)

The Taylor & Francis class and bibliography style files are hosted on the journal's pages and are
copyrighted by T&F, so they are **not** included here. Get them from the AMF
"Instructions for authors → Formatting and templates / LaTeX" page (search:
*Applied Mathematical Finance instructions for authors*), or from the T&F Overleaf gallery:

1. **`interact.cls`** — the Interact layout class (article-based, single-column for review).
2. The AMF **reference style** `.bst` (T&F call these "Reference Style …"). AMF's instructions
   page states which one to use — confirm it there and download the matching `.bst`
   (e.g. a Chicago author-date style such as `tfcad.bst`, or the numbered style if AMF specifies
   one). The provided manuscript uses a **manual** `thebibliography`, so it compiles **without**
   a `.bst`; you only need the `.bst` if you prefer to switch to BibTeX (see §4).

Place `interact.cls` (and any `.bst`) in the **same directory** as `swap_curve_tension_AMF.tex`
— T&F explicitly say you do not need to install it into your TeX distribution.

## 2. Compile

```bash
pdflatex swap_curve_tension_AMF.tex
pdflatex swap_curve_tension_AMF.tex      # twice for cross-references
```

The manuscript already uses the Interact front-matter commands:
- `\title{...}`, `\author{\name{Ashutosh Sinha\thanks{Email: ...}} \affil{Independent Researcher}}`,
  `\maketitle`
- `\begin{abstract} ... \end{abstract}`
- `\begin{keywords} ... \end{keywords}`
- `\begin{jelcode} C63; G12; E43; C61; G17 \end{jelcode}` (Interact provides this environment)
- `\articletype{Research Article}`

If a package clash appears on first compile (rare), it will be between `interact` and one of
`amsthm`/`hyperref`; the fix is to load `hyperref` last (already done) and, if needed, drop the
`\theoremstyle{plain}` line.

## 3. Figures

All figures referenced by the manuscript (`fig_curves.png`, `fig_fbd.png`, `fig_phi.png`,
`fig_locality.png`, `fig_forwards.png`, `fig_fwdspace.png`, `fig_nu_phi.png`, `fig_nu_sigma.png`,
`fig_nu_curve.png`, `fig_fwd_nonunif.png`, `fig_multicurve.png`, `fig_jacobian.png`,
`fig_benchmark.png`) must sit in the directory. T&F prefer **vector** figures (EPS/PDF) for
production; the PNGs are fine for review, but regenerate as PDF/EPS for the final version
(`benchmark.py`, `multicurve.py`, etc. can save with `.savefig('name.pdf')`).

## 4. Optional: switch to BibTeX

If you prefer BibTeX over the manual list, create a `.bib` from the entries, set
`\bibliographystyle{<AMF style>}` and `\bibliography{<file>}`, and remove the manual
`thebibliography`. Keep the manual list otherwise — it already matches the cited works.

## 5. Submission steps (ScholarOne)

1. Anonymise if AMF uses double-blind review (check the instructions). To anonymise in Interact,
   move `\maketitle` to sit *between* `\title` and `\author` so the name/affiliation are hidden,
   and remove the email `\thanks` and any self-identifying acknowledgements.
2. Prepare: the PDF, the `.tex` source, `interact.cls`, any `.bst`, and all figure files.
3. Submit at the AMF ScholarOne site (linked from the instructions page). Provide title,
   abstract, keywords, and a cover letter (see template below).
4. Suggested/excluded reviewers if requested.
5. On acceptance, supply editable source + figures (vector) as T&F request.

## 6. Cover letter (paste/edit)

> Dear Editors,
>
> Please consider the enclosed manuscript, "Swap Curve Construction and Optimization of the
> Tension Parameter: An Energy-Partition Formulation," for publication in *Applied Mathematical
> Finance*.
>
> The paper addresses a long-standing gap in term-structure interpolation: the tension spline is
> well established, but the choice of its tension parameter has remained heuristic. We place the
> tension spline on a variational footing as the equilibrium of a thin elastic beam under axial
> tension, and we show that the curvature and first-derivative penalties are the two halves of a
> single stored elastic energy. From this we derive a dimensionless, scale-invariant selection
> criterion — the bending-energy fraction φ(σ) — with a transparent equipartition default, and we
> connect the resulting tension to hedge locality and risk. All empirical results use the real
> U.S. Treasury par-yield curve of 22 May 2026; a fully reproducible reference implementation
> accompanies the paper. The method sits naturally alongside Hagan and West (2006), published in
> this journal, against which we benchmark.
>
> The manuscript is original, not under consideration elsewhere, and has no competing interests.
> A preprint is posted on SSRN. I am the sole author.
>
> Thank you for your consideration.
> Ashutosh Sinha — ajsinha@gmail.com

## 7. Checklist
- [ ] `interact.cls` (+ optional AMF `.bst`) in the manuscript directory.
- [ ] Compiles cleanly twice; cross-references resolved.
- [ ] Figures present (vector for final).
- [ ] Keywords + JEL present (already in the file).
- [ ] Cover letter ready; reviewers list if requested.
- [ ] Anonymised version if double-blind.
