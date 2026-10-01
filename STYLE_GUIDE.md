# Style guide: writing in the voice of Ivanov, Adkins, Kennedy et al. (JPP)

Reference text: P. G. Ivanov, T. Adkins, D. Kennedy, M. Giacomin, M. Barnes and A. A. Schekochihin,
"Suppression of temperature-gradient-driven turbulence by sheared flows in fusion plasmas"
(arXiv 2405.00854; JPP). Source: `reference_style/ivanov2025_source/ivanov25shear_arxiv.tex`,
with `maths.tex`, `format.tex`, `ref.tex`, `tikz.tex`, `jpp.cls`, `jpp.bst`.

All quotations below are verbatim from the source (LaTeX macros shown as they render, or left as
macros where the macro is the point). Counts are from a grep of the main `.tex` file. Equation and
section numbers inside quotes, e.g. "(3.12)" or "§4.1", are reconstructed from the `\label` order and
have not been checked against a compiled PDF; the `\label` names in the source are authoritative.

---

## 0. The voice in one paragraph

A theorist explaining a scaling argument to a capable colleague. The text is in the first-person
plural, patient and careful, and it develops its argument step by step. It states assumptions
explicitly and numbers them (i)–(iii). Each new result is tied back, by equation number, to the
assumptions it used. Confidence comes from the scaling argument and the agreement with data, and
the paper does not oversell. Every claim is hedged exactly as far as the evidence warrants:
"consistent with", "likely due to", "appears to be". Mild dry wit is allowed ("an unsurprising
outcome", "even simplistic"), but the paper never hypes. British spelling throughout.

---

## 1. Voice and person

- **"We" throughout** (93 lower-case "we", 20 "We", 21 "our", 7 "Our"). "We" is the authors, and
  it is also the authors and reader working together. The paper never uses "I" and never writes
  "the authors" in the body.
  - "In this article, we consider the effects of an imposed perpendicular flow shear on saturated electrostatic gyrokinetic (GK) turbulence."
  - "We do not know how to determine λ theoretically but we can confirm numerically that the above arguments are sound and that λ > 0."
- **"Let us" moves the argument forward** (5 times). It starts a new step and invites the reader along.
  - "Let us analyse the physics of both regimes, starting with the weak-shear one."
  - "Let us revisit the current understanding of how such an energy cascade determines the properties of the saturated turbulence."
  - "Nevertheless, the theory derived in §3.2 makes numerically falsifiable predictions for the momentum transport of turbulence with imposed flow shear. Let us investigate them here."
- **The reader is addressed in the third person**, usually to head off an objection before a
  referee raises it:
  - "A cautious reader may have spotted a potential clash between having A°(0) ≪ 1 at the outer scale and the theory of the energy cascade laid out in §3.1"
  - "A careful reader may spot another issue, which is also resolved by using k_x,tilt instead of k_x."
  - "The reader might also wonder why we consider k_x = 0 given that the outer-scale radial wavenumber is not zero."
  - "Here, leaving the reader cognisant of these recent developments, we shall nevertheless focus on the case A ≲ 1."
  - "...remind the reader of the standard results..."
- **"We shall" / "we will"** for what comes next: "We shall distinguish two different regimes of sheared turbulence".
- **"One" as the impersonal subject** when a procedure is general: "One can then perform a coordinate transformation...", "one cannot make a similar statement about the radial wavenumber".
- **Active voice is the default.** The passive is kept for established facts and for numerical
  observations: "has been shown numerically to give rise to", "is found to be approximately".

---

## 2. Sentence length and rhythm

- Prose sentences are long. The median is about 25–30 words (rough count with equations stripped),
  and the top decile runs past 55 words. The long ones are built with semicolons, "viz.," clauses,
  parenthetical asides and inline equations, and they stay readable because the logic runs
  strictly left to right.
- **Short sentences are punctuation.** They close a paragraph or turn the argument:
  - "Thus, we require additional assumptions. There are multiple ways to proceed."
  - "Therefore, this approach does not work."
  - "This argument is, however, incorrect."
  - "The agreement with figure 3 is evident."
  - "Its precise form will not be needed here."
- **Semicolons** join a claim to its gloss or consequence: "...pass the energy injected at large scales down to dissipative ones (...); the energy of the fluctuations is then thermalised at these small scales, heating the plasma." / "Crucially, the nonlinear interactions (3.3) and the linear drive (3.4) have the same form in both the laboratory frame and the shearing frame; therefore, (3.7) captures completely the effects of flow shear in the shearing frame."
- **Colons** introduce an equation or a list of consequences: "The electrostatic GK equation is closed by the quasineutrality condition:".
- **Spaced em-dashes `---`** (9 times) set off an aside or a verdict at the end of a sentence:
  - "...will not reduce the turbulent transport by more than an order-unity amount --- an unsurprising outcome."
  - "However, this occurs in a surprising and nontrivial way --- turbulence becomes radially localised into disjoint turbulent patches"
  - "...(...constitute a small sample) --- it is thus a natural testbed for any theory aspiring to tokamak relevance."
- **"In other words,"** restates a result in plain physical language after the algebra: "In other words, if the unsheared turbulence has A°(0) ∼ 1 at the outer scale, shearing it with any γ_E < γ°(0) will not reduce the turbulent transport by more than an order-unity amount".
- **Physical intuition follows the derivation**, often with "This is intuitively clear:" or "This reflects the intuitive notion that...":
  - "This is intuitively clear: stronger flow shear pushes turbulence towards smaller (and thus faster) scales since the larger (and slower) eddies are more strongly affected by the shear."
  - "This reflects the intuitive notion that radially elongated fluctuations should be more susceptible to sheared poloidal flows."

---

## 3. How sections open and close

**Openings** state the scope or purpose of the section in one sentence and often point back to an
earlier section:
- "To test the validity of the theory presented in §3.2, we consider two different models of turbulence."
- "In this section, we report numerical simulations in a triply periodic domain..."
- "We now explore the validity of our theory for plasma turbulence in axisymmetric toroidal geometry."
- "We have thus far focused on the influence of imposed flow shear on the turbulent heat transport. In sheared systems, another important quantity of interest is..."
- "Magnetised-plasma turbulence exhibits a broad range of different saturation mechanisms, and so a universal theory of turbulent saturation under the influence of flow shear is not feasible. Instead, here we focus on one particular type of saturated state, viz., ..."
- "For the remainder of this article, we assume an equilibrium shear flow..."
- Appendix: "Here we develop the theory of the transition range ..., as promised in §3.2.1."
- Appendix: "We want to compute the radial flux of poloidal momentum, [eq] in the fluid ETG model used in §4.1 and §5."

**Closings** do one of three things:
1. A short bridging sentence to the next subsection: "Let us analyse the physics of both regimes, starting with the weak-shear one."
2. A summary pointer to a figure: "Figure 3 summarises the expected dependence of the heat flux on γ_E in both regimes."
3. A statement of scope or an open question: "We find no such structures and it remains an open question exactly what the determining factor for their appearance is."

They often begin with "Finally, let us mention that..." or "Finally, we note that...":
- "Finally, let us mention that, while here we shall consider only cases where A ≲ 1, this is not necessarily satisfied in all instances of fusion-relevant turbulence."
- "Finally, we note that the appearance of coherent structures ... can naturally lead to nonunique saturation."

Each section also states the condition under which its results hold and says what stays general:
"Thus, our arguments will hold regardless of whether the zero-shear outer scale is chosen à la Barnes et al. (2011), through a grand critical balance, or otherwise."

---

## 4. Stating and hedging results

**Predictions** are stated as consequences of the equations: "X implies", "we conclude that", "we find", "we expect".
- "Then, (3.12) implies that the radial turbulent heat flux satisfies [eq]"
- "Assuming a local, constant-flux cascade, ..., we conclude that [eq]"
- "Thus, we propose that, for γ_E ≫ γ°(0), the outer scale will be governed by the balance of nonlinear, injection, and shearing rates:"

**Agreement with numerics** is stated plainly. The paper calibrates the strength of each claim
carefully:
- Strong: "The agreement with figure 3 is evident." / "holds very well up to γ̂_E ≈ 100" / "in excellent agreement with numerical simulations" (abstract only) / "Given the remarkably good fit for Q as a function of γ_E, ..."
- Medium: "is followed reasonably well for both values of the ion-temperature gradient" / "There is good agreement with the weak-shear scaling" / "approximately linear, as expected from (3.17)"
- Weak: "roughly consistent with λ = 2" / "Nevertheless, our numerical results are consistent with it."
- "consistent with" is the main hedge (10 times). Use it for anything short of a quantitative fit.

**Discrepancies** are reported openly, together with a likely cause and a check of that cause:
- "However, the prediction that Π should be constant in the strong-shear regime is not observed, likely due to finite-hyperviscosity effects and the finite extent of the inertial range. Our theory does not take into account the finite value of hyperviscosity..."
- "...consistent with the hypothesis that the hyperviscous cutoff is responsible for the discrepancy with the theoretical prediction (5.8)."
- "...the small discrepancy is likely due to our estimates of k_x° and k_y° ... being an imperfect measure of the outer scale."

**Hedging vocabulary.** "likely", "appears to be", "suggests", "may", "might", "(approximately)",
"(at least approximately)", "at least in some regimes", "roughly", "of the order of". Parenthetical
qualifiers sit inside the sentence instead of in a separate caveat sentence:
- "Π is (approximately) proportional to γ_E"
- "...the system can saturate at (at least) two different levels of heat transport."
- "we are (at this stage) unable to draw any conclusions about the more general applicability..."
- "...a scale-independent, even if numerically small, number..."

**Assumptions are stated before they are used** and marked "Assuming that ...":
- "Assuming that the heat flux Q_s is dominated by contributions from the outer scale, we can estimate it..."
- "Assuming that the rate of energy injection is determined by the linear-instability growth rate γ_k and that the latter satisfies γ_k ∝ k_y, ..."

**"Note that"** (17 times) adds a necessary qualification or a check on consistency:
- "Note that at no step leading to (3.13) did we use any formulae from §3.1 that relied on isotropy..."
- "Note that the precise definition of φ̄ is not important because the phenomenological theory that is to follow predicts only scalings"

**"Crucially,"** marks the key step or the key result (used sparingly, twice): "Crucially, in both regimes, [Π/Q ∝ γ_E]."

**"Recall that"** reminds the reader of an earlier result before a number is compared with it: "Recall that the theory of §3.2.1 predicts that the transition between the two regimes should occur at γ_E ∼ γ°(0)..."

---

## 5. Equations: introduction, punctuation, reference

- **Equations are part of the sentence.** Each display ends with a comma or full stop according to
  the grammar. The sentence usually carries on after it with "where", "which", "and so" or "i.e.,".
  - "...satisfies [eq], where φ̄ is a measure of the characteristic amplitude..."
  - "...the shearing frame (Newton et al. 2010; ...): [eq]." (colon before, full stop after)
  - "[Π(γ_E) ∝ γ_E^0,] i.e., the momentum flux is independent of the imposed flow shear"
- **Lead-in phrases** for a display, all used in the paper: "satisfies", "is given by", "can be written as", "becomes", "implies", "we find", "we conclude that", "leading to", "as follows:", "viz.,", "of the form", "defined respectively as".
- **"where" clauses** define every new symbol right after the equation. When there are several,
  they run as one long comma-separated list (see the paragraph after the GK equation (2.1)):
  - "where k_⊥ and k_∥ are the typical perpendicular and parallel (to the mean magnetic field) wavenumbers, ρ_s and Ω_s are the Larmor radius and frequency of the charged particles of species s, ..."
  - "where we have defined the fluctuation aspect ratio at scale k_y as A ≡ k_x/k_y."
  - "where we have introduced the *critical shearing rate*"
- **≡ for definitions**, ∼ for order-of-magnitude estimates, ∝ for scalings, ≈ for numerical
  values. These are kept strictly apart, and ∼ is never used where ∝ is meant.
- **Referring to equations.** Use `\cref`, which renders as a bare "(3.5)" with no word
  "equation". A sentence may start with "Equation (3.7) tells us that..." or "Expressions (3.13) and (3.14) predict that..." because a sentence cannot start with "(". Within a sentence, write "(3.12) implies", "by (3.13)", "via (3.27)", "as per (2.4)", "in view of (2.2) and (2.3)", "together with (3.7) and (3.9)". Ranges render as "(4.1)–(4.3)".
- **A named equation is the noun, with the number in apposition**: "the quasineutrality condition (2.5)", "the wavenumber drift (3.7)", "the GK equation (2.1)", "the 'grand critical balance' (3.6)".
- **Short inline maths is wrapped in `\mbox{\(...\)}`** to prevent line breaks, e.g. `\mbox{\(\gammaE < \gammaoO\)}`.
- **Multi-line derivations** use `align` and `\nonumber`. They are short, and each line follows from
  the last.
- **Labels are descriptive**: `eq:heatflux_weakly_sheared`, `eq:kxo_weak_shear`, `fig:tripleplot`, `sec:weak_shear`, `appendix:aniso_to_iso`, `footnote:pvg`.

---

## 6. Footnotes

- **Many footnotes** (13 in the main text and appendices). They hold technical caveats,
  definitions of jargon, and replies to "what about...?" objections, which keeps the main line of
  argument clean.
- They are often long (up to about 150 words), in full sentences, with citations and equation
  references.
- They are labelled and cross-referenced when used again: `\footnote{\label{footnote:pvg}...}`, then "(see also footnote 4)" or "Given the discussion in footnote 3, we must mention that...".
- Typical openings: "Strictly speaking, ...", "Note that ...", "Throughout this article, we use 'GK' to refer to ...", "In general, ...", "Indeed, ...", "This assumption can be made weaker: ...", "A careful reader may spot another issue...", "Notice that, ...".
  - "Strictly speaking, (2.7) contains another injection term that is associated with the radial gradient of u. Here, we assume that this can be neglected (see also footnote 4)."
  - "This assumption can be made weaker: we will only need γ_k to be approximately independent of k_x for k_x < k_y."
- A footnote marker sits after punctuation and is attached directly to it: `as follows,\footnote{...}` or `...satisfies\footnote{...}` just before a display.

---

## 7. Figures, tables and captions

- **Captions are long, self-contained paragraphs.** They define every line style, colour,
  normalisation and panel, so that the figure can be read without the text. Fig. 4's caption is
  about 270 words.
- **Panel labels go first in parentheses**: "(a) Time-averaged, saturated radial turbulent heat flux, ... (b) The outer-scale wavenumbers ..." A panel can also be referred to as "panel (a)".
- **The first sentence is a noun phrase** that says what is plotted. It often has no main verb:
  - "An illustration of the relationship between the nonlinear mixing rate τ_nl⁻¹, the energy-injection rate γ_k, and the location of the outer scale, where τ_nl⁻¹ ∼ γ_k."
  - "A qualitative diagram of the heat flux Q as a function of the flow shear γ_E in the case of (a) ... and (b) ..."
  - "Snapshots of φ (top row) and δT_e/T_e (bottom row) in the (x, y) plane for Sim1 simulations with four different values of γ_E, as specified above each column."
  - "Radial localisation of turbulent perturbations at very large values of flow shear."
- **Line keys are given in words**: "The black dashed and dash-dotted lines show the theoretical predictions (3.13) and (3.18), respectively", "The vertical black dotted line marks...", "shown as a black dashed line", "hollow triangles".
- **Captions give interpretation and cross-references** to equations and sections: "In the former, k_y° is (approximately) pinned to k_y°(0) but k_x° increases linearly with γ_E."
- **Fitted values are stated in the caption**: "plotted using γ̂_c ≈ 39, found by fitting to the data presented here."
- **In the text, cite figures with `\cref`**, which renders "figure 4(a)" in lower case mid-sentence and `\Cref` gives "Figure 4(a)" at the start of a sentence. Write "figures 2–3" for ranges. Figure sentences start "Figure 4(a) shows ...". A pointer can also go in parentheses: "(see figure 1)", or in brackets inside a parenthetical maths context: "[see figure 3(b)]".
- Schematic figures (tikz) come before data figures. The theory gets a qualitative diagram first,
  and the simulation figure is later compared against that diagram ("The agreement with figure 3 is evident.").
- Tables carry a caption that defines all symbols: "A summary of the simulation parameters used in §4.1."

---

## 8. Citing the literature in-line

- **natbib, author–year** (`jpp.bst`). Use `\citep` (56 times) for parenthetical support and `\citet` (25 times) when the authors are the subject of the sentence.
  - Many citations are grouped in one `\citep`, appended to a noun phrase: "the impact of sheared flows on the turbulence (Artun & Tang 1992; Synakowski et al. 1997; ...)".
  - `\citet` as subject: "In the absence of flow shear, Barnes et al. (2011) posit (i) that ...", "Recent numerical and analytical work by Nies et al. (2024) suggests that ...".
  - `\citealt` inside an existing parenthesis: "(same for all species, see \citealt{abel13})", "(\citealt{lin99,...} constitute a small sample)".
  - "e.g., \citet{abel13} or \citet{catto2019}" points to derivations.
- **A citation goes right after the concept it supports**, not at the end of the sentence:
  "the so-called 'quench' rule \citep{...}, according to which ...", "a free-energy conservation law \citep{abel13} of the form".
- **The paper is generous to earlier work and explicit about the differences**: "There are some parallels that can be drawn between these theories and our approach: ... However, the derivation of our theory ... is not related to these decorrelation theories and produces different scalings..."
- Similar results elsewhere are acknowledged directly: "Similar bistability in gyrokinetic turbulence with mean flow shear has been reported by \citet{christen22}."
- Credit is given to unpublished sources: "the same as that used by \citet{lithwick07_shear} (and attributed by him to Gordon Ogilvie)".
- **Codes are set in `\texttt{}` with a citation**: "the GK code \texttt{GENE} \citep{jenko00, jenko00GENE}", "\texttt{stella}", "\texttt{GS2}", "\texttt{GKW}".

---

## 9. Notation and macros (maths.tex conventions)

Load `format.tex`, `ref.tex` and `maths.tex` in the same way: `\def\closesymbol{\!}` before `\input{maths.tex}`.

- **Bold vectors**: `\vec{x}` is redefined to `\boldsymbol`. Unit vectors use `\uvec{b}`, giving b̂. Write `\vr, \vk, \vv, \vR, \vE, \vB`, `\vkperp` = **k**_⊥. Use `\bcdot` for a bold dot product and `\grad` for a bold ∇.
- **Upright subscript labels** via `\text`: `\vth` = v_th, `\vthe`, `\vthi`, `\mfp`, `\ope`, `\lDe`. Superscript labels are upright too: τ_nl^o is `\tau_\text{nl}^\text{o}`.
- **Upright differential d**: `\rmd` (`\frac{\rmd W}{\rmd t}`). Partial derivatives are `\pt, \px, \py, \pz` = ∂_t etc., and there is also `\partd[n]{f}{x}`.
- **Species-generic subscripts**: `\s` = s. `\fs, \Fs, \dfs` give f_s, F_s, δf_s, with `\closesymbol` tightening the space in "δ f". There are electron and ion versions: `\fe, \Foe, \dfe, \fion` (not `\fi`, which is TeX's `\fi`), `\dfi`.
- **Perpendicular and parallel**: `\kperp, \kpar, \vperp, \vpar, \vvperp, \hvpar`; `\kperprhoisq` etc.
- **Integrals**: `\intr, \intv, \intw, \intR` give ∫ d³**r** with a trailing thin space.
- **Averages**: `\avgR{}` gives ⟨·⟩_R (gyroaverage) and `\avgr{}` gives ⟨·⟩_r. Use the `...inline` versions to avoid `\left/\right` sizing in running text.
- **Poisson bracket**: `\pbra{f}{g}` gives {f, g} (`\pbrainline` in text).
- **Order symbol**: `\order{}` and `\orderinline{}` give O(·).
- **ExB in prose**: `\exb{}` gives "**E**×**B**" (it includes its own math delimiters, so follow it with `{}`).
- **Paper-specific macros are defined in the preamble of the main file**, not in maths.tex: `\kxo, \kyo, \gammao, \taunlo, \AoO, \gammaE, \gammaani`. Outer-scale quantities carry an upright superscript "o", and the zero-shear value is written as an argument "(0)", as in γ°(0). The paper says so explicitly: "Here and in what follows, the superscript '^o' denotes quantities associated with the outer scale."
- **Normalised quantities wear a hat**: Q̂, γ̂ (with defining equations given). Calligraphic symbols are used for derived dimensionless ratios: A ≡ k_x/k_y.
- **Inline maths is `\( ... \)`, not `$...$`**. The source mostly uses `\(...\)`, with a few `$`.
- **cleveref formats** (format.tex): sections render as §3.2, figures as "figure 4", equations as "(3.5)", equation ranges as "(4.1)–(4.3)" and figure ranges as "figures 2–3". Always use `\cref`/`\Cref` and never hand-type "Eq." or "Fig.".
- **Defined terms are italicised on first use** with `\textit` (10 times): *outer scale*, *dissipation scale*, *inertial range*, *fluctuation aspect ratio*, *shearing frame*, *weak-shear regime*, *critical shearing rate*. Emphasis in an argument also uses `\textit`: "at *any* scale".
- **Single quotes for coined or loose terms**: `` `slab' ``, `` `streamers' ``, `` `quench' rule ``, `` `grand critical balance' ``, `` `ferdinons' ``, `` `natural' ``, `` `isotropise' ``. LaTeX open quote ` and close '. Double quotes are never used.
- Write "(i)~...; and (ii)~..." for enumerated assumptions, with a tie `~` after the label, and then refer back as "assumption (ii)" and "(i)–(iii) imply".

---

## 10. Abstract and introduction structure

**Abstract** (one paragraph, about 190 words, no citations, no equations, no section numbers,
`\noindent` first):
1. **One long opening sentence** that starts from the premise ("Starting from the assumption that ..., we formulate a detailed phenomenological theory for ...").
2. **What the theory introduces**: "Our theory introduces two distinct regimes, called the weak-shear and strong-shear regimes, each with its own set of scaling laws for ..."
3. **The key discovery**, in first person: "We discover that the ratio of ... (i.e., their aspect ratio) at the outer scale plays a central role in determining ..."
4. **Validation**, naming the models with an inline (i)/(ii) list: "Our theoretical predictions are found to be in excellent agreement with numerical simulations of two paradigmatic models ...: (i) an electrostatic fluid model of slab electron-scale turbulence, and (ii) Cyclone-base-case gyrokinetic ion-scale turbulence."
5. **Wider implication**, hedged: "Additionally, our theory envisions a potential mechanism for ..."

**Introduction** (three paragraphs, about 650 words):
1. **Context, from broad to narrow.** It starts with fusion ("The quest for controlled fusion as a viable and sustainable energy source has been a long-standing scientific and engineering challenge."), then turbulence limits performance, then "therefore crucial", then the specific topic, which is given a long `\citep` list. It ends on why the topic matters physically: "Sheared flows can modify the size and shape of the fluctuations, and thus have a direct impact on the transport properties of the plasma."
2. **Intellectual framing.** It says what kind of theory is possible and gives its lineage: "Despite the absence of a rigorous theory ..., it is still possible to develop phenomenological models that, at least in some regimes, capture its essential features and allow us to make falsifiable, qualitative, and sometimes even quantitative, predictions ..." This leads to K41, the cascade picture, and how it sets transport. It closes with a one-line hook: "An imposed or self-generated sheared flow plays a nontrivial role in all of this."
3. **Roadmap.** "In this article, we consider ..." It then walks through every section with `\cref` in narrative order: "We first give ... in §2, and then, in §3.1, remind the reader of ... In §3.2, we proceed to develop ... To verify our theoretical predictions, in §4, we present ... Then, in §5, we discuss ..., before finally summarising and discussing our results in §6." The roadmap also previews the main result in words: "The effect of this shear is to suppress the turbulent fluctuations and, in turn, the turbulent heat flux according to a certain scaling with the size of the shear."

There is no bulleted "contributions" list and no "the main results of this paper are".

**Summary and discussion** (the final section):
- It opens by restating the premise and what was done, in the present perfect tense: "Starting from the standard picture of turbulent saturation via a local energy cascade (§3.1), we have developed a theory for ..."
- It then has one paragraph per regime or result, each restating the scaling in words with an equation reference in apposition, followed by validation ("are confirmed to hold over a range of four orders of magnitude for the flow shear in idealised fluid ETG simulations (§4.1)").
- Next it places the work against competing theories, in a measured way.
- It ends on an application or speculation, framed as future work: "Whether this mechanism ... is indeed realised in multiscale plasma turbulence appears to be a promising subject for future work."

**Back matter order**: Acknowledgements ("We thank ... for inspiring discussions and invaluable feedback."), then Funding (with grant numbers and CSD3/DiRAC wording), then Declaration of interests ("The authors report no conflict of interest."), then appendices.

**Appendices** hold long derivations, numerical-method details and secondary theory. Each opens by
saying why it exists and what it will show. The main text forward-references them: "In Appendix A, we develop a simple theory for the transition region."

---

## 11. Limitations and open questions

- **Scope is set at the start and repeated where it matters**:
  - "a universal theory of turbulent saturation under the influence of flow shear is not feasible. Instead, here we focus on one particular type of saturated state, viz., ..."
  - "Our analysis depends only on the injection rate being a function solely of the poloidal wavenumber, while the precise relationship between the injection rate and the linear growth rate is outside of the scope of the current work."
- **Stock phrases** (use these and avoid "beyond the scope"):
  - "...falls outside the scope of this paper." / "is outside of the scope of this paper." / "falls outside of the range of validity of the theory presented in §3"
  - "...whose detailed investigation is the subject of our ongoing work"
  - "it remains an open question exactly what the determining factor for their appearance is."
  - "We do not know how to determine λ theoretically but ..."
  - "...we are (at this stage) unable to draw any conclusions about ..."
  - "...appears to be a promising subject for future work."
- **Admit model limits plainly, then say why the result is still useful**:
  - "While this model is extremely simple, even simplistic, the benefit of using it is that ..."
  - "Of course, our oversimplified model of electron-scale turbulence cannot be applied directly to any experimental studies. Nevertheless, our theory suggests that ..."
  - "Unfortunately, because the aspect ratio A°(0) cannot be varied in the fluid model, we are unable to vary the size of the transition range ... Nevertheless, our numerical results are consistent with it."
- **Unrealistic idealisations are flagged in a footnote** with the physics they omit: "In general, a pure perpendicular linear shear is not realistic: e.g., u is purely toroidal in axisymmetric devices ... Here we assume that there is no PVG instability (or at least that it is irrelevant for the saturated state, which is reasonable if the shear is not too large; ...)".
- **Unexpected observations are reported as observations** and not over-interpreted: "We find no such structures and it remains an open question..."

---

## 12. Spelling and typography

- **British spelling, -ise forms**: normalised (17), localised, magnetised, initialised, thermalised, realised, maximise, characterised, summarising, idealised, parameterised, linearised, utilise, isotropise. Also behaviour, colour, centred, gyrocentre and analyse. No -ize forms appear in prose. The only -or words are standard ones (factor, vector, operator, major, minor).
- **Closed prefixes, no hyphen**: nonlinear, nonzero, nontrivial, nonmonotonic, nonunique, nonexponentially, quasineutrality, quasistatic, reemerges, multiscale. The exception is "non-interacting" / "non-negligible", which keep the hyphen.
- **"Time scale" and "length scale" as two words.** "Wavenumber" is one word.
- **Compound adjectives are hyphenated**: "temperature-gradient-driven turbulence", "free-energy injection rate", "weak-shear regime", "outer-scale eddies", "local-energy-cascade phenomenology", "order-unity amount", "Cyclone-base-case gyrokinetic ion-scale turbulence".
- **Latin abbreviations** always take a following comma: "i.e.," (23), "e.g.," (15), "viz.," (14). "etc" is written without a stop in one place. "à la" and "na\"ive" appear with diacritics.
- **"viz."** introduces an exact specification ("namely"). It is a signature of this voice: "the gyroradii of the main ion species ρ_i and of the electrons ρ_e", "viz., L_ns ∼ L_Ts ∼ L".
- **Section signs**: §3.2 (via cleveref). Write "Appendix A" for appendices.
- **"vs."** in captions ("Π vs. flow shear"). Use "versus" in a caption's opening phrase: "Radial turbulent heat flux versus time".
- **"data" is singular**: "The data from all four sets overlays", "whose data is presented".
- **Ordinal connectives: use "Firstly," / "Secondly,"** (author preference). The source is
  inconsistent here: Appendix C has "First, the ion gyroradius is large, ..." followed by
  "Secondly, the electron distribution function is expressed as ...", and Appendix A has "First, let us show that ...". The new paper should use "Firstly, ..." / "Secondly, ..." / "Finally, ..." and never "First, ..." / "Second, ...".

---

## 13. Connective words (with approximate frequency in source)

| Function | Words used |
|---|---|
| Consequence | Therefore (11), Thus (7), Consequently (3), and so, hence (3), follows immediately |
| Contrast | However (12), In contrast (4), Nevertheless (6), Despite, While, Instead, unlike |
| Emphasis/confirmation | Indeed (4), Crucially (2), In particular, Specifically, Of course |
| Specification | viz. (14; the paper never writes "namely"), i.e. (23), e.g. (15), In other words, To be more specific, Specifically |
| Pointing back/forward | Recall that, As discussed in §, as promised in §, as per, in view of, see also, Here and in what follows, For the remainder of this article, thus far |
| Qualification | Note that (17), Strictly speaking, at least approximately, at least in some regimes, provided, as long as |
| Sequence | Firstly / Secondly (preferred over First / Second), Then, Finally, before finally |
| Addition | Additionally, Furthermore, Also, Similarly, Just like we did for |
| Evaluative (rare) | Fortunately (1), Unfortunately (1), surprising, unsurprising, remarkably good |

Other phrases characteristic of this voice: "It is instructive to consider", "it is meaningful to
distinguish", "is readily generalisable to", "the so-called", "a prime candidate", "a natural
testbed", "paradigmatic", "borne out by", "manifestly", "na\"ively", "at the cost of".

---

## 14. Things the paper never does

- It never uses "I", "the present authors", "this author" or "the authors" in the body.
- It never uses hype words: "novel", "groundbreaking", "for the first time", "unprecedented" and
  "Interestingly," have zero occurrences. "significant(ly)" (5 times) and "clearly" (once) are used
  only literally ("too weak to influence the saturated state significantly"), never as rhetoric.
- It never uses contractions (no "don't", "it's", "we'll").
- It never uses bullet points or numbered lists in the body. Lists are inline: (i)~..., (ii)~....
- It never writes "Eq.", "Fig.", "Eqs. (3)-(5)" or "equation (3)" mid-sentence. It always uses
  cleveref's bare "(3.5)" and lower-case "figure 4".
- It never uses double quotation marks.
- It never uses American spelling.
- It never leaves an equation unpunctuated, and never leaves a new symbol undefined after its first
  display.
- It never claims agreement stronger than the data shows, and it never hides a discrepancy.
- It never uses "beyond the scope". It writes "outside (of) the scope" instead.
- It never uses "First," / "Second," as ordinals in the new paper (author preference:
  "Firstly," / "Secondly,").
- It never puts citations or equations in the abstract.
- It never has a "Contributions" list or a "Conclusions" section separate from "Summary and
  discussion".
- It never reports numerical parameters in the main-text theory sections. Simulation details are
  kept in the numerical-results section, tables and appendices.

---

## 15. Quick checklist before submitting a section

1. Does the section open by stating what it does, with a `\cref` back to what it builds on?
2. Is every assumption stated as "Assuming that ..." before it is used?
3. Does each display end in "," or "." and get followed by a "where" clause defining new symbols?
4. Are relations ≡ / ∼ / ∝ / ≈ used correctly?
5. Is each numerical comparison hedged to match the evidence (evident / reasonably well / consistent with)?
6. Are caveats in footnotes and not in the main line?
7. Spelling -ise, "nonlinear", "time scale", "i.e.,", "viz.,", single quotes, "Firstly/Secondly"?
8. Do the captions stand alone?
9. Does the section end with a bridge, a figure summary or an explicit open question?
