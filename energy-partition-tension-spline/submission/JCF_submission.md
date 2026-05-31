# Journal of Computational Finance (Risk Journals / Infopro Digital) — submission & help

**Publisher:** Infopro Digital (Risk.net) · **Editor-in-Chief:** Christoph Reisinger (Oxford)
**Submission:** by email to **technical@infopro-digital.com** · **Preferred format: PDF**
**Why this journal:** its scope is exactly numerical/computational techniques in pricing, hedging
and risk management — a tension-spline construction method fits well.

---

## 1. Important: there is no JCF LaTeX class, and there is a word limit

Two things differ from the Taylor & Francis route:

- **No submission class.** Unlike AMF (which uses `interact.cls`), Risk Journals do **not** provide
  or require a LaTeX class at submission. They want a **PDF** for review; only **on publication**
  do they ask for the LaTeX source **including the `.bbl` file** and charts in **EPS** (or XLS for
  data charts). Word files are also accepted. So you submit the compiled PDF — any clean class is
  fine for that PDF.
- **Length.** Risk Journals state a **maximum recommended length of ~5,500 words**, plus some
  allowance for charts and formulas. The current full manuscript is roughly **double** that
  (≈10,000+ words of body text across 22 pages). **The paper must be condensed for JCF** — this is
  the binding constraint, more than any formatting. See §3.

## 2. House reference style (use `jcf_references.tex`)

JCF/Risk uses **author–date** citations, alphabetical by author, last names + initials, with DOIs
where available. Example of the required form:

> Hagan, P. S., and West, G. (2006). Interpolation methods for curve construction.
> *Applied Mathematical Finance* 13(2), 89–129. https://doi.org/10.1080/13504860500396032.

I have reformatted the paper's full reference list into this house style in
**`jcf_references.tex`** — drop it in to replace the numbered list. You also need to switch
**in-text** citations from the numbered `[n]` form to author–date, e.g. `\cite{Hagan}` →
"Hagan and West (2006)". The cleanest way is to load `natbib` with an author–date style and use
`\citet`/`\citep`; the reference keys in `jcf_references.tex` match the existing `\cite` keys, so
the substitution is mechanical.

## 3. Condensing plan to reach ~5,500 words

The paper was deliberately written at full, teaching-level granularity. For JCF, keep the
contribution and empirics; move the expository derivation detail out of the main text:

- **Keep in main text:** the variational set-up and the energy-partition criterion (Sec. IV–V),
  the φ properties result, the validation/locality results (Sec. VI), the multi-curve example
  (Sec. X), the benchmark comparison and φ* characterisation (Sec. XII). These are the contribution.
- **Move to a technical appendix or online supplement:** the step-by-step algebra now in
  Sec. III (cubic), the long-form construction algebra, the closed-form energy integrals
  (App. B), the stable-evaluation forms (App. C), and the non-uniform tridiagonal derivation
  (App. A). JCF allows supplementary material, and the reproducible code already documents these.
- **Compress** the fixed-income basics (Sec. II) to a short paragraph with references, and tighten
  the discussion. Target the 5,500-word body, with the derivations and code as supplements.

This is an editorial pass rather than new work; I can produce the condensed JCF manuscript on
request.

## 4. Figures

- Keep the number of figures/tables to a minimum.
- Figures must be in the main PDF **and** supplied as **separate editable vector files** — EPS
  (Illustrator) or MATLAB-style SVG. Regenerate the matplotlib figures as EPS/SVG, e.g.
  `plt.savefig('fig_locality.eps')`, for the final files.
- Code and data may be submitted as **supplementary material** — point to the reproducible
  package (and the Zenodo DOI once minted).

## 5. Submission steps

1. Condense to ~5,500 words (§3); compile a clean PDF (any standard class is fine — the `interact`
   AMF file or the original both work for producing the PDF; JCF only needs the PDF for review).
2. Apply the author–date references (`jcf_references.tex`) and author–date in-text citations.
3. Email the PDF to **technical@infopro-digital.com** with a brief cover note (template below).
4. Expect **8–12 weeks** for single-blind peer review by two or more referees.
5. On acceptance: supply LaTeX source **with the `.bbl`**, figures as EPS/XLS, plus a short bio
   and headshot for the journal website.

## 6. Cover note (email body)

> Dear Editor,
>
> Please find attached, for consideration in the *Journal of Computational Finance*, the manuscript
> "Swap Curve Construction and Optimization of the Tension Parameter: An Energy-Partition
> Formulation."
>
> The paper introduces a dimensionless, scale-invariant criterion for selecting the tension
> parameter in tension-spline curve construction, derived from the elastic-beam energy the spline
> minimises, and links the choice to hedge locality and risk. It includes a multi-curve SOFR/OIS
> construction, analytic risk Jacobians, and an even-handed benchmark against Hagan–West
> monotone-convex, Smith–Wilson, log-discount and piecewise-flat methods. All results use the real
> U.S. Treasury curve of 22 May 2026; a fully reproducible implementation accompanies the paper and
> can be provided as supplementary material.
>
> The manuscript is original and not under consideration elsewhere; I am the sole author and
> declare no competing interests. A preprint is available on SSRN.
>
> Kind regards,
> Ashutosh Sinha — ajsinha@gmail.com

## 7. Checklist
- [ ] Condensed to ~5,500 words (derivations/appendices moved to supplement).
- [ ] Author–date references (`jcf_references.tex`) + author–date in-text citations.
- [ ] Clean compiled PDF.
- [ ] Figures as separate EPS/SVG files; kept to a minimum.
- [ ] Code/data prepared as supplementary material (+ Zenodo DOI).
- [ ] Cover email to technical@infopro-digital.com.
