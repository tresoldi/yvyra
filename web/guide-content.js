/*
 * Guide content: 8-step interactive walkthrough for yvyra.
 * Each step has: title, HTML content, and a YAML example.
 */

const GUIDE_STEPS = [
{
title: "What makes yvyra different",
content: `
<p>Most phylogenetic software was written for molecular sequence data: thousands
of sites, symmetric substitution models, automated alignment. Scholars working
with linguistic, stemmatological, or cultural data face a different situation:
30 to 100 characters where every coding reflects a comparative judgment.
Directionality matters (lenition is common, fortition is rare; morphological
loss is easy, spontaneous complexification is not), and different characters
carry different evidential weight.</p>

<p><strong>Asymmetric models.</strong> yvyra treats custom rate matrices as an
edge case buried in the NEXUS syntax. In yvyra, directional cost matrices are
the main modeling tool. Built-in templates cover common patterns (irreversible
loss, lenition, grammaticalization, ordered change), or you can write arbitrary
cost matrices. A shared innovation under an asymmetric model contributes more to the
likelihood than a shared retention, which is how comparativists have
always treated the evidence informally.</p>

<p><strong>Named states.</strong> Characters use meaningful labels
(<code>conservative</code>/<code>lenited</code>, not <code>0</code>/<code>1</code>)
from input through to ancestral state reconstruction output, where you read
<code>p(lenited)</code> instead of <code>p(1)</code>.</p>

<p><strong>Confidence-weighted tips.</strong> Fragmentary attestations are not
binary (known or missing). A damaged Tocharian tablet, a contested manuscript
reading, a late recording of a dying language: these carry graded uncertainty.
Partial confidence (<code>{lost: 0.8}</code>) and polymorphic codings
(<code>{state_a: 0.6, state_b: 0.4}</code>) pass through the likelihood
computation rather than being collapsed to missing data.</p>

<p><strong>Estimated character weights.</strong> Setting <code>weight: estimated</code>
treats a character's weight as a free parameter sampled during MCMC under a
Gamma prior. Fixed weights are also available when the analyst has prior
knowledge about relative informativeness.</p>

<p><strong>Per-model-group ascertainment correction.</strong> When models with
different numbers of states coexist in one analysis (binary loss characters
alongside three-state allomorphy systems), the Lewis (2001) Mkv correction must
be computed separately for each model group. yvyra applies a single global
correction using the symmetric Mk model for all characters. yvyra computes the
correction per model group, using each group's actual Q matrix.</p>

<p><strong>Per-character diagnostics.</strong> A posterior probability tells you
what the tree recovered. Diagnostics tell you which characters are responsible:</p>
<ul>
<li><strong>Clade-character support</strong>: signed per-character contribution
to a clade's posterior (positive = supports, negative = conflicts)</li>
<li><strong>Sensitivity</strong>: clade posterior when each character is removed,
computed by importance sampling without rerunning MCMC</li>
<li><strong>Consistency index</strong>: Fitch parsimony on the MAP tree, flagging
homoplastic characters</li>
<li><strong>ASR entropy</strong>: ancestral state probabilities at internal nodes,
with Shannon entropy and credible intervals</li>
</ul>

<p><strong>Beyond linguistics.</strong> Stemmatology and cultural phylogenetics
involve the same kind of data: discrete characters, directional change,
differential weight. The examples in this guide include manuscript
transmission and pottery traditions alongside the linguistic datasets.</p>

<p><strong>YAML input.</strong> Analyses are specified in YAML rather than NEXUS
command blocks. A single file holds models, weights, confidence values, and
diagnostics, and can be versioned alongside the data.</p>

<p><strong>Browser execution.</strong> The MCMC engine compiles to WebAssembly
and runs in this browser tab. No installation, no server.</p>
`,
example: null
},

{
title: "Your first analysis",
content: `
<p>A fictional language family inspired by Ancient Greek dialects: six languages,
14 phonological and morphological characters, all binary, under the default
symmetric Mk model.</p>

<p>The data has clear phylogenetic signal: five characters are shared by Arcadian,
Ionic, and Attic; three more link Ionic and Attic specifically; and three group
Dorian with Aeolic. The tree should recover these groupings if the signal
is strong enough.</p>

<p>Click <strong>Load and Run this example</strong> below.
The analysis runs 25,000 iterations and should complete in a few seconds.</p>

<p><strong>What to look for in the output:</strong></p>
<ul>
<li>The <em>consensus tree</em>: does it group (Ionic, Attic) and (Dorian, Aeolic)?</li>
<li>The <em>clade posteriors</em>: which groupings have strong support (&gt;80%)?</li>
<li>The <em>best log-likelihood</em>: a measure of how well the model fits the data</li>
<li>The <em>convergence assessment</em>: is the ESS above 200?</li>
</ul>
`,
example: `# Step 2: Your first analysis
# A fictional family inspired by Ancient Greek dialects
# 6 languages, 14 characters (phonological and morphological)

name: first_analysis

taxa:
  - id: koine
    name: Koine
  - id: dorian
    name: Dorian
  - id: aeolic
    name: Aeolic
  - id: arcadian
    name: Arcadian
  - id: ionic
    name: Ionic
  - id: attic
    name: Attic

characters:
  # Shared innovations in Arcadian, Ionic, Attic
  - id: h_loss
    data:
      koine: preserved
      dorian: preserved
      aeolic: preserved
      arcadian: lost
      ionic: lost
      attic: lost
  - id: w_loss
    data:
      koine: preserved
      dorian: preserved
      aeolic: preserved
      arcadian: lost
      ionic: lost
      attic: lost
  - id: long_vowel_shift
    data:
      koine: absent
      dorian: absent
      aeolic: absent
      arcadian: present
      ionic: present
      attic: present
  - id: dat_pl_merger
    data:
      koine: distinct
      dorian: distinct
      aeolic: distinct
      arcadian: merged
      ionic: merged
      attic: merged
  - id: gen_s_ending
    data:
      koine: absent
      dorian: absent
      aeolic: absent
      arcadian: present
      ionic: present
      attic: present
  # Shared innovations in Ionic, Attic only
  - id: vowel_contraction
    data:
      koine: absent
      dorian: absent
      aeolic: absent
      arcadian: absent
      ionic: present
      attic: present
  - id: quantitative_metathesis
    data:
      koine: absent
      dorian: absent
      aeolic: absent
      arcadian: absent
      ionic: present
      attic: present
  - id: psilosis
    data:
      koine: absent
      dorian: absent
      aeolic: absent
      arcadian: absent
      ionic: present
      attic: present
  # Shared innovations in Dorian, Aeolic
  - id: first_aor_part
    data:
      koine: absent
      dorian: present
      aeolic: present
      arcadian: absent
      ionic: absent
      attic: absent
  - id: dat_en_ending
    data:
      koine: absent
      dorian: present
      aeolic: present
      arcadian: absent
      ionic: absent
      attic: absent
  - id: patronymic_suffix
    data:
      koine: absent
      dorian: present
      aeolic: present
      arcadian: absent
      ionic: absent
      attic: absent
  # Innovations separating Koine as outgroup
  - id: dual_loss
    data:
      koine: preserved
      dorian: lost
      aeolic: lost
      arcadian: lost
      ionic: lost
      attic: lost
  - id: optative_reduction
    data:
      koine: preserved
      dorian: lost
      aeolic: lost
      arcadian: lost
      ionic: lost
      attic: lost
  # Noisy character
  - id: accent_shift
    data:
      koine: conservative
      dorian: innovative
      aeolic: conservative
      arcadian: innovative
      ionic: conservative
      attic: innovative
analysis:
  coding: variable
  seed: 1305
  iterations: 25000
  sample_frequency: 50
`
},

{
title: "Named states",
content: `
<p>Opaque state labels, whether numeric ("0" and "1") or arbitrary ("x" and "y"),
obscure what the data represents. Named states make the YAML readable and carry
through to the output: ancestral state reconstruction shows
<code>p(lenited)</code> rather than <code>p(1)</code>.</p>

<p>This example uses Romance languages, with characters encoding phonological
and morphological changes from Latin. States are labeled
<code>preserved</code>/<code>lenited</code>/<code>lost</code> for intervocalic
consonants, <code>preposed</code>/<code>postposed</code> for the definite
article. One three-state character (<code>intervocalic_p</code>) shows
multi-state encoding.</p>

<p><strong>What to look for:</strong> The state names appear in the ASR entropy
table columns. Compare the readability with the numeric output from the
previous step.</p>
`,
example: `# Step 3: Named states
# Romance languages with meaningful state labels
# 6 languages, 14 characters including a 3-state character

name: named_states

taxa:
  - id: lat
    name: Latin
  - id: fra
    name: French
  - id: spa
    name: Spanish
  - id: ita
    name: Italian
  - id: por
    name: Portuguese
  - id: ron
    name: Romanian

characters:
  # All daughters share (Latin as outgroup)
  - id: future_tense
    states: [synthetic, periphrastic]
    data:
      lat: synthetic
      fra: periphrastic
      spa: periphrastic
      ita: periphrastic
      por: periphrastic
      ron: periphrastic
  - id: articles
    states: [absent, present]
    data:
      lat: absent
      fra: present
      spa: present
      ita: present
      por: present
      ron: present
  # Italo-Western (all except Romanian)
  - id: cl_palatalization
    states: [preserved, palatalized]
    data:
      lat: preserved
      fra: palatalized
      spa: palatalized
      ita: palatalized
      por: palatalized
      ron: preserved
  - id: case_loss
    states: [retained, lost]
    data:
      lat: retained
      fra: lost
      spa: lost
      ita: lost
      por: lost
      ron: retained
  - id: neuter_loss
    states: [retained, lost]
    data:
      lat: retained
      fra: lost
      spa: lost
      ita: lost
      por: lost
      ron: retained
  - id: preposed_article
    states: [absent, present]
    data:
      lat: absent
      fra: present
      spa: present
      ita: present
      por: present
      ron: absent
  # Western Romance (French, Spanish, Portuguese)
  - id: final_s
    states: [lost, preserved]
    data:
      lat: preserved
      fra: preserved
      spa: preserved
      ita: lost
      por: preserved
      ron: lost
  - id: voicing
    states: [preserved, voiced]
    data:
      lat: preserved
      fra: voiced
      spa: voiced
      ita: preserved
      por: voiced
      ron: preserved
  - id: intervocalic_lenition
    states: [preserved, weakened]
    data:
      lat: preserved
      fra: weakened
      spa: weakened
      ita: preserved
      por: weakened
      ron: preserved
  # Ibero-Romance (Spanish, Portuguese)
  - id: initial_f
    states: [preserved, lost]
    data:
      lat: preserved
      fra: preserved
      spa: lost
      ita: preserved
      por: lost
      ron: preserved
  - id: pl_palatalization
    states: [preserved, palatalized]
    data:
      lat: preserved
      fra: preserved
      spa: palatalized
      ita: preserved
      por: palatalized
      ron: preserved
  # Gallo-Romance (French)
  - id: front_rounding
    states: [absent, present]
    data:
      lat: absent
      fra: present
      spa: absent
      ita: absent
      por: absent
      ron: absent
  - id: liaison
    states: [absent, present]
    data:
      lat: absent
      fra: present
      spa: absent
      ita: absent
      por: absent
      ron: absent
  # 3-state character: demonstrates multi-state named encoding
  - id: intervocalic_p
    states: [preserved, lenited, lost]
    data:
      lat: preserved
      fra: lost
      spa: lenited
      ita: preserved
      por: lenited
      ron: preserved

analysis:
  coding: variable
  seed: 1305
  iterations: 25000
  sample_frequency: 50
  sitelikes: yes

diagnostics:
  - asr_entropy
`
},

{
title: "Asymmetric models",
content: `
<p>The symmetric Mk model assigns equal probability to forward and reverse
transitions. But lenition is far more frequent than fortition, and
morphological categories are lost more readily than they are created.
An asymmetric rate matrix gives shared innovations more weight in the
likelihood than shared retentions.</p>

<p>Celtic is a natural test case: lenition pervades the family's phonological
history. This example combines two built-in templates, <code>loss</code>
(cost 1 forward, 100 reverse) and <code>lenition</code> (cost 1 forward,
2 reverse), with a custom <code>one_way</code> model for innovations that
do not reverse. Old Irish serves as outgroup.</p>

<p><strong>What to look for:</strong> Goidelic (Scottish Gaelic + Manx) and
Brythonic (Welsh + Breton + Cornish) should both receive posterior support.
The character diagnostics show which features contribute to the Welsh+Breton
subgrouping and which conflict with it.</p>
`,
example: `# Step 4: Asymmetric models
# Celtic languages with directional evolution models
# 6 languages, 12 characters, 3 model types

name: asymmetric_celtic

taxa:
  - id: oir
    name: Old_Irish
  - id: wel
    name: Welsh
  - id: bre
    name: Breton
  - id: sga
    name: Scottish_Gaelic
  - id: cor
    name: Cornish
  - id: mnx
    name: Manx

models:
  # Custom model: innovation easy, reversal nearly impossible
  one_way:
    costs:
      simple:    0      1.0
      extended:  100.0  0

characters:
  # All share (Old Irish as archaic outgroup)
  - id: initial_mutation
    states: [absent, present]
    model: loss
    data:
      oir: absent
      wel: present
      bre: present
      sga: present
      cor: present
      mnx: present
  - id: verb_initial
    states: [free, fixed]
    model: loss
    data:
      oir: free
      wel: fixed
      bre: fixed
      sga: fixed
      cor: fixed
      mnx: fixed
  # Brythonic (Welsh, Breton, Cornish)
  - id: p_lenition
    states: [preserved, lenited]
    model: lenition
    data:
      oir: preserved
      wel: lenited
      bre: lenited
      sga: preserved
      cor: lenited
      mnx: preserved
  - id: kw_preservation
    states: [preserved, lost]
    model: loss
    data:
      oir: lost
      wel: preserved
      bre: preserved
      sga: lost
      cor: preserved
      mnx: lost
  - id: internal_inflection
    states: [preserved, lost]
    model: loss
    data:
      oir: preserved
      wel: lost
      bre: lost
      sga: preserved
      cor: lost
      mnx: preserved
  - id: preposition_conjugation
    states: [full, reduced]
    model: lenition
    data:
      oir: full
      wel: reduced
      bre: reduced
      sga: full
      cor: reduced
      mnx: full
  # Goidelic (Scottish Gaelic, Manx)
  - id: lenition_system
    states: [simple, extended]
    model: one_way
    data:
      oir: simple
      wel: simple
      bre: simple
      sga: extended
      cor: simple
      mnx: extended
  - id: palatalization
    states: [preserved, simplified]
    model: lenition
    data:
      oir: preserved
      wel: preserved
      bre: preserved
      sga: simplified
      cor: preserved
      mnx: simplified
  - id: dual_number
    states: [preserved, lost]
    model: loss
    data:
      oir: preserved
      wel: lost
      bre: lost
      sga: lost
      cor: lost
      mnx: lost
  # Welsh + Breton closer
  - id: aspirate_mutation
    states: [absent, present]
    model: one_way
    data:
      oir: absent
      wel: present
      bre: present
      sga: absent
      cor: absent
      mnx: absent
  # Mixed signal
  - id: nasal_mutation
    states: [productive, residual]
    model: lenition
    data:
      oir: productive
      wel: productive
      bre: residual
      sga: residual
      cor: residual
      mnx: residual
  - id: case_system
    states: [retained, lost]
    model: loss
    data:
      oir: retained
      wel: lost
      bre: lost
      sga: lost
      cor: lost
      mnx: lost

analysis:
  coding: variable
  seed: 1305
  iterations: 50000
  sample_frequency: 50
  sitelikes: yes

diagnostics:
  - chardiag: [Welsh, Breton]
`
},

{
title: "Character weights",
content: `
<p>Characters in a phylogenetic dataset rarely carry equal amounts of signal.
A shared morphological restructuring is more informative than a lexical
retention; a complex technological innovation outweighs an ornamental
choice. The question is how much more, and whether the analyst or the
data should decide.</p>

<p>The example here is from stemmatology, where the distinction between
significant and insignificant variants has been debated since Lachmann.
Six manuscript witnesses, 14 variant readings: conjunctive errors (shared
mistakes that could only arise through copying) get weight 3.0, separative
errors 2.0, and orthographic variants get <code>weight: estimated</code>,
so the MCMC sampler determines their contribution under a Gamma(2,2) prior.</p>

<p><strong>What to look for:</strong> The stemma should recover clear
transmission groups. Estimated weights are in the <code>.w</code> file;
orthographic variants should receive lower posterior weights than the
conjunctive errors. The convergence diagnostics report weight ESS alongside
tree length ESS.</p>
`,
example: `# Step 5: Manuscript tradition with weighted variants
# 6 witnesses, 14 readings: conjunctive errors, separative errors,
# and orthographic variants with estimated weights

name: manuscript_weights

taxa:
  - id: A
    name: Codex_A
  - id: B
    name: Codex_B
  - id: C
    name: Codex_C
  - id: D
    name: Codex_D
  - id: E
    name: Codex_E
  - id: F
    name: Codex_F

characters:
  # Conjunctive errors: strong evidence of copying (weight 3.0)
  - id: lacuna_3v
    states: [present, absent]
    weight: 3.0
    data:
      A: present
      B: present
      C: absent
      D: absent
      E: absent
      F: absent
  - id: transposition_ch7
    states: [original, transposed]
    weight: 3.0
    data:
      A: original
      B: original
      C: transposed
      D: transposed
      E: transposed
      F: transposed
  - id: interpolation_12r
    states: [absent, present]
    weight: 3.0
    data:
      A: absent
      B: absent
      C: present
      D: present
      E: present
      F: present
  # Separative errors: distinguish sub-branches (weight 2.0)
  - id: omission_5v
    states: [present, absent]
    weight: 2.0
    data:
      A: present
      B: present
      C: present
      D: absent
      E: absent
      F: absent
  - id: substitution_8r
    states: [original, variant]
    weight: 2.0
    data:
      A: original
      B: original
      C: original
      D: variant
      E: variant
      F: variant
  - id: gloss_marginal
    states: [absent, incorporated]
    weight: 2.0
    data:
      A: absent
      B: absent
      C: absent
      D: incorporated
      E: incorporated
      F: incorporated
  - id: rewriting_15v
    states: [original, rewritten]
    weight: 2.0
    data:
      A: original
      B: original
      C: original
      D: original
      E: rewritten
      F: rewritten
  - id: rubric_change
    states: [original, modified]
    weight: 2.0
    data:
      A: original
      B: original
      C: original
      D: original
      E: modified
      F: modified
  - id: abbreviation_style
    states: [expanded, contracted]
    weight: 1.5
    data:
      A: expanded
      B: contracted
      C: contracted
      D: expanded
      E: expanded
      F: expanded
  # Orthographic variants: let the data decide (estimated weight)
  - id: spelling_ae
    states: [ae, e]
    weight: estimated
    data:
      A: ae
      B: ae
      C: e
      D: e
      E: ae
      F: e
  - id: spelling_ci
    states: [ci, ti]
    weight: estimated
    data:
      A: ti
      B: ci
      C: ci
      D: ti
      E: ci
      F: ci
  - id: word_order_9r
    states: [order_a, order_b]
    weight: estimated
    data:
      A: order_a
      B: order_a
      C: order_b
      D: order_a
      E: order_b
      F: order_b
  # Shared by all descendants (A as archaic witness)
  - id: colophon
    states: [present, absent]
    data:
      A: present
      B: absent
      C: absent
      D: absent
      E: absent
      F: absent
  - id: chapter_division
    states: [original, modified]
    data:
      A: original
      B: modified
      C: modified
      D: modified
      E: modified
      F: modified

analysis:
  coding: variable
  seed: 1305
  iterations: 100000
  sample_frequency: 50
  sitelikes: yes
`
},

{
title: "Uncertain data",
content: `
<p>Linguistic attestations are rarely complete or unambiguous. Tocharian, known
only from fragmentary Buddhist manuscripts in Central Asia, is the extreme
case: damaged tablets, partial paradigms, disputed readings. But the same
problem arises whenever a datum is not certain but not entirely absent either.</p>

<p>Three kinds of uncertainty can be encoded:</p>
<ul>
<li><code>~</code> — missing data (marginalized over all states)</li>
<li><code>{state: 0.8}</code> — partial confidence (80% probability for this state)</li>
<li><code>{state_a: 0.6, state_b: 0.4}</code> — polymorphic or contested coding</li>
</ul>

<p>This example uses six Indo-European branches with Hittite as the Anatolian
outgroup. Several Tocharian characters are coded as missing or with partial
confidence, reflecting the fragmentary record. Greek and Latin carry some
uncertain codings for disputed features as well.</p>

<p><strong>What to look for:</strong> The Tocharian A+B clade should receive
support despite the incomplete data. Partial confidence passes more information
through the likelihood than <code>~</code> (missing) alone.</p>
`,
example: `# Step 6: Uncertain data — fragmentary Tocharian
# 6 IE branches, 14 characters with missing, partial, and polymorphic codings

name: uncertain_ie

taxa:
  - id: hit
    name: Hittite
  - id: grk
    name: Greek
  - id: lat
    name: Latin
  - id: san
    name: Sanskrit
  - id: toA
    name: Tocharian_A
  - id: toB
    name: Tocharian_B

characters:
  # Core IE vs Anatolian (Hittite as outgroup)
  - id: augment
    states: [absent, present]
    data:
      hit: absent
      grk: present
      lat: present
      san: present
      toA: present
      toB: present
  - id: thematic_vowel
    states: [absent, present]
    data:
      hit: absent
      grk: present
      lat: present
      san: present
      toA: present
      toB: present
  - id: feminine_gender
    states: [absent, present]
    data:
      hit: absent
      grk: present
      lat: present
      san: present
      toA: present
      toB: present
  # Tocharian A+B clade
  - id: vowel_rotation
    states: [absent, present]
    data:
      hit: absent
      grk: absent
      lat: absent
      san: absent
      toA: present
      toB: present
  - id: case_agglutination
    states: [absent, present]
    data:
      hit: absent
      grk: absent
      lat: absent
      san: absent
      toA: present
      toB: present
  - id: verb_class_merger
    states: [distinct, merged]
    data:
      hit: distinct
      grk: distinct
      lat: distinct
      san: distinct
      toA: merged
      toB: merged
  # Graeco-Aryan signal
  - id: prohibitive_ma
    states: [absent, present]
    data:
      hit: absent
      grk: present
      lat: absent
      san: present
      toA: ~                                # damaged tablet
      toB: ~                                # unattested
  - id: augment_preterite
    states: [absent, present]
    data:
      hit: absent
      grk: present
      lat: absent
      san: present
      toA: ~                                # missing
      toB: {present: 0.6, absent: 0.4}     # partially attested
  # Latin innovation
  - id: gen_sg_i
    states: [absent, present]
    data:
      hit: absent
      grk: absent
      lat: present
      san: absent
      toA: absent
      toB: absent
  # Fragmentary Tocharian: partial confidence
  - id: perfect_tense
    states: [retained, lost]
    data:
      hit: lost
      grk: retained
      lat: retained
      san: retained
      toA: {lost: 0.8}                     # probably lost (damaged text)
      toB: ~                                # unreadable
  - id: dual_number
    states: [retained, lost]
    data:
      hit: lost
      grk: retained
      lat: lost
      san: retained
      toA: ~                                # not attested
      toB: {retained: 0.6, lost: 0.4}      # contested
  - id: middle_voice
    states: [retained, lost]
    data:
      hit: retained
      grk: retained
      lat: {retained: 0.7}                 # partially retained
      san: retained
      toA: lost
      toB: lost
  # Noisy / uncertain
  - id: laryngeal_reflex
    states: [colored, lost]
    data:
      hit: colored
      grk: {colored: 0.3, lost: 0.7}       # disputed reflex
      lat: lost
      san: lost
      toA: lost
      toB: lost
  - id: subjunctive
    states: [retained, lost]
    data:
      hit: lost
      grk: retained
      lat: retained
      san: retained
      toA: {retained: 0.5, lost: 0.5}      # ambiguous forms
      toB: retained

analysis:
  coding: variable
  seed: 1305
  iterations: 100000
  sample_frequency: 50
  sitelikes: yes
`
},

{
title: "Diagnostics: who drives what",
content: `
<p>A posterior probability tells you what the tree recovered. The diagnostics
tell you which characters are responsible, and which ones conflict:</p>

<ul>
<li><strong>chardiag</strong>: per-character likelihood contribution and signed
clade support (positive = supports, negative = conflicts)</li>
<li><strong>sensitivity</strong>: clade posterior when each character is removed
(importance sampling, no rerun)</li>
<li><strong>consistency index</strong>: Fitch parsimony on the MAP tree, flagging
characters with homoplasy</li>
</ul>

<p>The example here moves outside linguistics to Bronze Age Mediterranean pottery.
Manufacturing techniques (slip, burnishing, kiln firing) tend to be
transmitted within traditions, while decorative features (painted bands,
pilgrim flasks) spread through trade. The diagnostics target the
Egyptian+Levantine clade, where trade-diffused characters pull against
the tree.</p>

<p><strong>What to look for:</strong> <code>lotus_decoration</code> supports the
clade; <code>base_ring_form</code> (shared by Cyprus and the Levant through
trade) conflicts. The sensitivity table shows the clade dropping from ~73% to
~46% if the lotus decoration character is removed.</p>
`,
example: `# Step 7: Bronze Age Mediterranean pottery traditions
# 6 traditions, 13 characters, diagnostics for Egyptian+Levantine clade

name: pottery_diagnostics

taxa:
  - id: min
    name: Minoan
  - id: myc
    name: Mycenaean
  - id: cyp
    name: Cypriot
  - id: hit
    name: Hittite
  - id: egy
    name: Egyptian
  - id: lev
    name: Levantine

characters:
  # Aegean shared technology
  - id: slip_technique
    states: [unslipped, slipped]
    data:
      min: slipped
      myc: slipped
      cyp: slipped
      hit: unslipped
      egy: unslipped
      lev: unslipped
  - id: burnishing
    states: [absent, present]
    data:
      min: present
      myc: present
      cyp: present
      hit: absent
      egy: absent
      lev: absent
  - id: stirrup_jar
    states: [absent, present]
    data:
      min: present
      myc: present
      cyp: present
      hit: absent
      egy: absent
      lev: absent
  # Minoan + Mycenaean exclusive
  - id: marine_motifs
    states: [absent, present]
    data:
      min: present
      myc: present
      cyp: absent
      hit: absent
      egy: absent
      lev: absent
  - id: palace_workshop
    states: [absent, present]
    data:
      min: present
      myc: present
      cyp: absent
      hit: absent
      egy: absent
      lev: absent
  - id: fresco_influence
    states: [absent, present]
    data:
      min: present
      myc: present
      cyp: absent
      hit: absent
      egy: absent
      lev: absent
  # Near Eastern shared
  - id: coil_building
    states: [absent, dominant]
    data:
      min: absent
      myc: absent
      cyp: absent
      hit: dominant
      egy: dominant
      lev: dominant
  - id: geometric_stamp
    states: [absent, present]
    data:
      min: absent
      myc: absent
      cyp: absent
      hit: present
      egy: present
      lev: present
  # Egyptian + Levantine closer
  - id: lotus_decoration
    states: [absent, present]
    data:
      min: absent
      myc: absent
      cyp: absent
      hit: absent
      egy: present
      lev: present
  - id: faience_glaze
    states: [absent, present]
    data:
      min: absent
      myc: absent
      cyp: absent
      hit: absent
      egy: present
      lev: present
  # Trade-diffused features (create homoplasy)
  - id: painted_bands
    states: [absent, present]
    data:
      min: present
      myc: present
      cyp: present
      hit: absent
      egy: present
      lev: present
  - id: pilgrim_flask
    states: [absent, present]
    data:
      min: absent
      myc: present
      cyp: present
      hit: absent
      egy: present
      lev: present
  # Cypriot trade contact with Levant (conflicts with tree)
  - id: base_ring_form
    states: [absent, present]
    data:
      min: absent
      myc: absent
      cyp: present
      hit: absent
      egy: absent
      lev: present

analysis:
  coding: variable
  seed: 1305
  iterations: 200000
  sample_frequency: 50
  sitelikes: yes

diagnostics:
  - chardiag: [Egyptian, Levantine]
  - sensitivity: [Egyptian, Levantine]
`
},

{
title: "The full picture: Anatolian languages",
content: `
<p>This example brings together everything from the previous steps in a single
analysis. The dataset is the full Anatolian character matrix from Billing and
Elgh (2023): 5 languages (Hittite, Palaic, Lydian, Luwian, Lycian), 27
characters (10 phonological, 17 morphological), 6 asymmetric cost matrices,
differential weights, missing data, and one polymorphic coding.</p>

<p>The six models encode different types of linguistic change: irreversible
innovation (<code>one_way</code>), lenition bias (<code>fort_len</code>),
a three-state genitive allomorphy system (<code>third_pl</code>), and an
unusual case where forward change is costlier than reversal
(<code>four_w</code>). Twelve cells are coded as missing (<code>~</code>),
mostly in Lydian and Palaic where the attestation is fragmentary. Palaic
<code>gen_ers</code> is polymorphic: 60% <code>no_ers</code>, 40%
<code>ers_primary</code>.</p>

<p>The diagnostics target the broader Luwic clade (Lydian, Luwian, Lycian)
rather than the Luwian+Lycian pair, because the latter receives near-100%
support and cannot be meaningfully decomposed.</p>

<p><strong>What to look for:</strong></p>
<ul>
<li>The Luwian+Lycian clade should receive strong support; the broader
(Lydian, Luwian, Lycian) grouping should be lower</li>
<li>Which characters drive the broader clade? (clade-character support)</li>
<li>Is the grouping robust to removing any single character? (sensitivity)</li>
<li>Which characters show homoplasy? (consistency index)</li>
<li>The ASR entropy table for three-state characters like <code>gen_ers</code>
and <code>encl_du</code></li>
</ul>
`,
example: `# Step 8: Full Anatolian analysis — Billing & Elgh (2023)
# 5 taxa, 27 characters, 6 asymmetric models, differential weights
# Adapted from anat2_210823.nex

name: anatolian_full

taxa:
  - id: hit
    name: Hittite
  - id: pal
    name: Palaic
  - id: lyd
    name: Lydian
  - id: luw
    name: Luwian
  - id: lyc
    name: Lycian

models:
  one_way:
    costs:
      conservative:  0      1.0    1.0
      intermediate:  100.0  0      1.0
      innovative:    100.0  100.0  0
  fort_len:
    costs:
      conservative:  0    1.0
      lenited:       2.0  0
  nti:
    costs:
      conservative:  0    1.0
      innovative:    2.0  0
  third_pl:
    costs:
      no_ers:         0       1.0    3.0
      ers_secondary:  100.0   0      3.0
      ers_primary:    100.0   100.0  0
  four_w:
    costs:
      conservative:  0    4.0
      innovative:    1.0  0
  tu_pronoun:
    costs:
      type_a:    0      2.0    4.0
      type_b:    100.0  0      4.0
      type_c:    100.0  2.0    0

characters:
  # Phonological innovations
  - id: len_fric
    model: fort_len
    data:
      hit: conservative
      pal: conservative
      lyd: lenited
      luw: ~
      lyc: lenited
  - id: eh1_a
    model: one_way
    states: [conservative, innovative]
    data:
      hit: conservative
      pal: conservative
      lyd: innovative
      luw: ~
      lyc: innovative
  - id: ye_yi
    model: one_way
    states: [conservative, innovative]
    weight: 2.0
    data:
      hit: conservative
      pal: conservative
      lyd: innovative
      luw: innovative
      lyc: innovative
  - id: centum
    model: one_way
    states: [conservative, innovative]
    data:
      hit: innovative
      pal: innovative
      lyd: innovative
      luw: conservative
      lyc: conservative
  - id: irreg_s_t
    model: one_way
    states: [conservative, innovative]
    weight: 4.0
    data:
      hit: conservative
      pal: conservative
      lyd: conservative
      luw: innovative
      lyc: innovative
  - id: g_y_front
    model: one_way
    states: [conservative, innovative]
    weight: 2.0
    data:
      hit: conservative
      pal: conservative
      lyd: ~
      luw: innovative
      lyc: innovative
  - id: init_gw_w
    model: one_way
    states: [conservative, innovative]
    data:
      hit: conservative
      pal: conservative
      lyd: ~
      luw: innovative
      lyc: innovative
  - id: med_gw_w
    model: one_way
    states: [conservative, innovative]
    data:
      hit: conservative
      pal: conservative
      lyd: innovative
      luw: innovative
      lyc: innovative
  - id: tn_nn
    model: one_way
    states: [conservative, innovative]
    data:
      hit: innovative
      pal: innovative
      lyd: ~
      luw: conservative
      lyc: ~
  - id: cop
    model: one_way
    states: [conservative, innovative]
    weight: 3.0
    data:
      hit: conservative
      pal: conservative
      lyd: ~
      luw: innovative
      lyc: innovative
  - id: lar_phon
    model: fort_len
    data:
      hit: conservative
      pal: conservative
      lyd: lenited
      luw: lenited
      lyc: lenited
  - id: sh_s
    model: one_way
    states: [conservative, innovative]
    data:
      hit: conservative
      pal: conservative
      lyd: innovative
      luw: conservative
      lyc: innovative
  # Morphological innovations
  - id: innov_nsi
    model: one_way
    states: [conservative, innovative]
    data:
      hit: conservative
      pal: conservative
      lyd: innovative
      luw: innovative
      lyc: innovative
  - id: gen_ha
    model: one_way
    states: [conservative, innovative]
    data:
      hit: conservative
      pal: innovative
      lyd: conservative
      luw: innovative
      lyc: innovative
  - id: gen_u
    model: one_way
    states: [conservative, innovative]
    data:
      hit: conservative
      pal: conservative
      lyd: innovative
      luw: innovative
      lyc: innovative
  - id: gen_ers
    model: third_pl
    weight: 3.0
    data:
      hit: no_ers
      pal: {no_ers: 0.6, ers_primary: 0.4}
      lyd: ers_secondary
      luw: ers_primary
      lyc: ers_primary
  - id: gen_to
    model: one_way
    states: [conservative, innovative]
    data:
      hit: conservative
      pal: conservative
      lyd: ~
      luw: innovative
      lyc: innovative
  - id: ablins_Vdi
    model: one_way
    states: [conservative, innovative]
    data:
      hit: conservative
      pal: conservative
      lyd: conservative
      luw: innovative
      lyc: innovative
  - id: nza_abl
    model: nti
    data:
      hit: innovative
      pal: innovative
      lyd: conservative
      luw: conservative
      lyc: conservative
  - id: lose_ont
    model: one_way
    states: [conservative, innovative]
    data:
      hit: conservative
      pal: conservative
      lyd: innovative
      luw: innovative
      lyc: innovative
  - id: encl_du
    model: tu_pronoun
    weight: 4.0
    data:
      hit: type_a
      pal: type_c
      lyd: type_b
      luw: type_c
      lyc: type_c
  - id: i_mut
    model: one_way
    states: [conservative, innovative]
    data:
      hit: conservative
      pal: conservative
      lyd: innovative
      luw: innovative
      lyc: innovative
  - id: a_o_stem_merge
    model: one_way
    states: [conservative, innovative]
    data:
      hit: innovative
      pal: innovative
      lyd: conservative
      luw: conservative
      lyc: conservative
  - id: len_e_stem
    model: four_w
    weight: 2.0
    data:
      hit: conservative
      pal: conservative
      lyd: innovative
      luw: conservative
      lyc: innovative
  - id: nu_stem_hi
    model: one_way
    states: [conservative, innovative]
    data:
      hit: conservative
      pal: ~
      lyd: conservative
      luw: innovative
      lyc: innovative
  - id: gen_genadj
    model: one_way
    states: [conservative, innovative]
    weight: 3.0
    data:
      hit: conservative
      pal: conservative
      lyd: conservative
      luw: innovative
      lyc: innovative
  - id: lex_tuw
    model: one_way
    states: [conservative, innovative]
    data:
      hit: conservative
      pal: ~
      lyd: innovative
      luw: innovative
      lyc: innovative

analysis:
  coding: variable
  seed: 1305
  iterations: 200000
  sample_frequency: 100
  sitelikes: yes

diagnostics:
  - chardiag: [Lydian, Luwian, Lycian]
  - sensitivity: [Lydian, Luwian, Lycian]
  - asr_entropy
`
}
];
