# How to publish: SSRN preprint + Zenodo code DOI

A short, practical walk-through. Do the code DOI first if you want to reference it in the SSRN
abstract; otherwise the order does not matter and you can add the DOI to SSRN later.

## A. Post the paper to SSRN

1. Create a free account at https://www.ssrn.com and choose **Submit a paper**.
2. Upload `swap_curve_tension_ssrn.pdf` (the SSRN copy, which carries the Keywords and JEL block
   on the title page). The original `swap_curve_tension.pdf` is preserved separately and unchanged.
3. Fill the form from `SSRN_submission.md`:
   - **Title** (section 1), **Abstract** (section 2, plain text), **Keywords** (section 3),
     **JEL Classification** (section 4).
   - **Author**: Ashutosh Sinha; affiliation "Independent Researcher" (or as you prefer);
     email ajsinha@gmail.com.
   - **Classification / networks**: Financial Engineering; Econometrics: Mathematical Methods;
     Capital Markets: Asset Pricing & Valuation.
4. Choose the license you want on SSRN and confirm you hold the rights.
5. Submit. SSRN review typically takes a few business days; you then receive a permanent
   **SSRN abstract ID** and URL (`https://ssrn.com/abstract=XXXXXXX`).
6. (Optional) Once you have the Zenodo DOI (below), edit the SSRN abstract to add a line:
   "Code and reproducibility: https://doi.org/10.5281/zenodo.XXXXXXX".

Note on journals: SSRN posting is compatible with later journal submission to Applied
Mathematical Finance, the Journal of Computational Finance, or Quantitative Finance. If you target
a Risk.net journal, check its preprint policy first, as some are stricter than the academic
journals.

## B. Mint a Zenodo DOI for the code

**Easiest route (GitHub + Zenodo):**

1. Create a public GitHub repo and push the contents of this `code/` folder (it already contains
   `README.md`, `LICENSE`, `requirements.txt`, `CITATION.cff`, `.zenodo.json`, and the modules).
2. Sign in to https://zenodo.org with GitHub, open **Settings -> GitHub**, and toggle the repo
   **On**.
3. On GitHub, create a **Release** (e.g. tag `v1.0.0`). Zenodo automatically archives that release
   and mints a DOI. It reads `.zenodo.json` / `CITATION.cff` for the metadata, so the title,
   authors, and description are filled in for you.
4. Zenodo issues two DOIs: a **version DOI** (this exact release) and a **concept DOI** (always
   points to the latest). Cite the version DOI in the paper for exact reproducibility.

**No-GitHub route:** zip this `code/` folder, upload it directly at https://zenodo.org ->
**New upload**, set Upload type = Software and License = MIT, paste the description from
`.zenodo.json`, and **Publish** to get the DOI.

5. Add the DOI to `CITATION.cff` (replace the placeholder) and, optionally, to the SSRN abstract
   and the paper's footer/acknowledgements on the next compile.

## C. Optional finishing touches

- Add a one-line "Code availability" sentence to the paper pointing at the Zenodo DOI, then
  recompile, before or after SSRN posting.
- If you later submit to a journal, keep the SSRN version and the journal version in sync, and add
  the journal reference to `CITATION.cff` once accepted.
