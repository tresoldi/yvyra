# yvyra guide

## What makes yvyra different

Most phylogenetic software was written for molecular sequence data: thousands
of sites, symmetric substitution models, automated alignment. Scholars working
with linguistic, stemmatological, or cultural data face a different situation:
30 to 100 characters where every coding reflects a comparative judgment.
Directionality matters (lenition is common, fortition is rare; morphological
loss is easy, spontaneous complexification is not), and different characters
carry different evidential weight.

**Asymmetric models.** Traditional phylogenetic tools treat custom rate matrices
as an edge case buried in the NEXUS syntax. In yvyra, directional cost matrices
are the main modeling tool. Built-in templates cover common patterns (irreversible
loss, lenition, grammaticalization, ordered change), or you can write arbitrary
cost matrices. A shared innovation under an asymmetric model contributes more to
the likelihood than a shared retention, which is how comparativists have always
treated the evidence informally.

**Named states.** Characters use meaningful labels (`conservative`/`lenited`,
not `0`/`1`) from input through to ancestral state reconstruction output, where
you read `p(lenited)` instead of `p(1)`.

**Confidence-weighted tips.** Fragmentary attestations are not binary (known or
missing). A damaged Tocharian tablet, a contested manuscript reading, a late
recording of a dying language: these carry graded uncertainty. Partial
confidence (`{lost: 0.8}`) and polymorphic codings (`{state_a: 0.6, state_b:
0.4}`) pass through the likelihood computation rather than being collapsed to
missing data.

**Estimated character weights.** Setting `weight: estimated` treats a
character's weight as a free parameter sampled during MCMC under a Gamma prior.
Fixed weights are also available when the analyst has prior knowledge about
relative informativeness.

**Per-model-group ascertainment correction.** When models with different numbers
of states coexist in one analysis (binary loss characters alongside three-state
allomorphy systems), the Lewis (2001) Mkv correction must be computed separately
for each model group. yvyra computes the correction per model group, using each
group's actual Q matrix.

**Per-character diagnostics.** A posterior probability tells you what the tree
recovered. Diagnostics tell you which characters are responsible:

- *Clade-character support*: signed per-character contribution to a clade's
  posterior (positive = supports, negative = conflicts)
- *Sensitivity*: clade posterior when each character is removed, computed by
  importance sampling without rerunning MCMC
- *Consistency index*: Fitch parsimony on the MAP tree, flagging homoplastic
  characters
- *ASR entropy*: ancestral state probabilities at internal nodes, with Shannon
  entropy and credible intervals

**Beyond linguistics.** Stemmatology and cultural phylogenetics involve the same
kind of data: discrete characters, directional change, differential weight. The
examples in this guide include manuscript transmission and pottery traditions
alongside the linguistic datasets.

**YAML input.** Analyses are specified in YAML rather than NEXUS command blocks.
A single file holds models, weights, confidence values, and diagnostics, and
can be versioned alongside the data.

**Browser execution.** The MCMC engine compiles to WebAssembly and runs in the
browser. No installation, no server.


## Your first analysis

A fictional language family inspired by Ancient Greek dialects: six languages,
14 phonological and morphological characters, all binary, under the default
symmetric Mk model.

The data has clear phylogenetic signal: five characters are shared by Arcadian,
Ionic, and Attic; three more link Ionic and Attic specifically; and three group
Dorian with Aeolic. The tree should recover these groupings if the signal is
strong enough.

What to look for in the output:

- The consensus tree: does it group (Ionic, Attic) and (Dorian, Aeolic)?
- The clade posteriors: which groupings have strong support (>80%)?
- The best log-likelihood: a measure of how well the model fits the data
- The convergence assessment: is the ESS above 200?

```yaml
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
```


## Named states

Opaque state labels, whether numeric ("0" and "1") or arbitrary ("x" and "y"),
obscure what the data represents. Named states make the YAML readable and carry
through to the output: ancestral state reconstruction shows `p(lenited)` rather
than `p(1)`.

This example uses Romance languages, with characters encoding phonological and
morphological changes from Latin. States are labeled
`preserved`/`lenited`/`lost` for intervocalic consonants,
`preposed`/`postposed` for the definite article. One three-state character
(`intervocalic_p`) shows multi-state encoding.

What to look for: the state names appear in the ASR entropy table columns.
Compare the readability with numeric output.

```yaml
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
  # ... 9 more characters (see examples/romance_named_states.yaml)
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
```


## Asymmetric models

The symmetric Mk model assigns equal probability to forward and reverse
transitions. But lenition is far more frequent than fortition, and morphological
categories are lost more readily than they are created. An asymmetric rate
matrix gives shared innovations more weight in the likelihood than shared
retentions.

Celtic is a natural test case: lenition pervades the family's phonological
history. This example combines two built-in templates, `loss` (cost 1 forward,
100 reverse) and `lenition` (cost 1 forward, 2 reverse), with a custom
`one_way` model for innovations that do not reverse. Old Irish serves as
outgroup.

What to look for: Goidelic (Scottish Gaelic + Manx) and Brythonic (Welsh +
Breton + Cornish) should both receive posterior support. The character
diagnostics show which features contribute to the Welsh+Breton subgrouping and
which conflict with it.

### Defining a custom model

A model is a cost matrix where rows and columns are states. The diagonal is
always zero. Higher costs mean less probable transitions:

```yaml
models:
  one_way:
    costs:
      simple:    0      1.0     # simple -> extended costs 1.0
      extended:  100.0  0       # extended -> simple costs 100 (nearly impossible)
```

### Built-in templates

Three templates are available without defining a cost matrix:

- `loss`: forward 1, reverse 100 (irreversible loss)
- `lenition`: forward 1, reverse 2 (weakening bias)
- `analogy`: 3-state asymmetric change

Use them by name:

```yaml
characters:
  - id: initial_mutation
    model: loss
    data:
      oir: absent
      wel: present
      bre: present
      sga: present
```

The full Celtic example is in `examples/celtic_asymmetric.yaml`.


## Character weights

Characters in a phylogenetic dataset rarely carry equal amounts of signal. A
shared morphological restructuring is more informative than a lexical retention;
a complex technological innovation outweighs an ornamental choice. The question
is how much more, and whether the analyst or the data should decide.

The example here is from stemmatology, where the distinction between significant
and insignificant variants has been debated since Lachmann. Six manuscript
witnesses, 14 variant readings: conjunctive errors (shared mistakes that could
only arise through copying) get weight 3.0, separative errors 2.0, and
orthographic variants get `weight: estimated`, so the MCMC sampler determines
their contribution under a Gamma(2,2) prior.

### Fixed weights

```yaml
characters:
  - id: lacuna_3v
    weight: 3.0       # conjunctive error: strong signal
    data: { ... }
  - id: omission_5v
    weight: 2.0       # separative error: moderate signal
    data: { ... }
```

### Estimated weights

```yaml
characters:
  - id: spelling_ae
    weight: estimated  # let the data decide
    data: { ... }
```

What to look for: estimated weights are in the `.w` file; orthographic variants
should receive lower posterior weights than the conjunctive errors. The
convergence diagnostics report weight ESS alongside tree length ESS.

The full manuscript example is in `examples/manuscript_weights.yaml`.


## Uncertain data

Linguistic attestations are rarely complete or unambiguous. Tocharian, known
only from fragmentary Buddhist manuscripts in Central Asia, is the extreme case:
damaged tablets, partial paradigms, disputed readings. But the same problem
arises whenever a datum is not certain but not entirely absent either.

Three kinds of uncertainty can be encoded:

- `~` -- missing data (marginalized over all states)
- `{state: 0.8}` -- partial confidence (80% probability for this state)
- `{state_a: 0.6, state_b: 0.4}` -- polymorphic or contested coding

This example uses six Indo-European branches with Hittite as the Anatolian
outgroup. Several Tocharian characters are coded as missing or with partial
confidence, reflecting the fragmentary record. Greek and Latin carry some
uncertain codings for disputed features as well.

```yaml
characters:
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
  - id: perfect_tense
    states: [retained, lost]
    data:
      hit: lost
      grk: retained
      lat: retained
      san: retained
      toA: {lost: 0.8}                     # probably lost (damaged text)
      toB: ~                                # unreadable
```

What to look for: the Tocharian A+B clade should receive support despite the
incomplete data. Partial confidence passes more information through the
likelihood than `~` (missing) alone.

The full Indo-European example is in `examples/uncertain_tocharian.yaml`.


## Diagnostics

A posterior probability tells you what the tree recovered. The diagnostics tell
you which characters are responsible, and which ones conflict:

- **chardiag**: per-character likelihood contribution and signed clade support
  (positive = supports, negative = conflicts)
- **sensitivity**: clade posterior when each character is removed (importance
  sampling, no rerun)
- **consistency index**: Fitch parsimony on the MAP tree, flagging characters
  with homoplasy

The example here moves outside linguistics to Bronze Age Mediterranean pottery.
Manufacturing techniques (slip, burnishing, kiln firing) tend to be transmitted
within traditions, while decorative features (painted bands, pilgrim flasks)
spread through trade. The diagnostics target the Egyptian+Levantine clade,
where trade-diffused characters pull against the tree.

### Requesting diagnostics

Add a `diagnostics` section at the end of the YAML file. The clade is specified
by listing the taxa it contains:

```yaml
diagnostics:
  - chardiag: [Egyptian, Levantine]
  - sensitivity: [Egyptian, Levantine]
  - asr_entropy
```

What to look for: `lotus_decoration` supports the clade; `base_ring_form`
(shared by Cyprus and the Levant through trade) conflicts. The sensitivity
table shows the clade dropping from ~73% to ~46% if the lotus decoration
character is removed.

The full pottery example is in `examples/pottery_diagnostics.yaml`.


## Worked example: Anatolian languages

This example brings together everything from the previous sections in a single
analysis. The dataset is the full Anatolian character matrix from Billing and
Elgh (2023): 5 languages (Hittite, Palaic, Lydian, Luwian, Lycian), 27
characters (10 phonological, 17 morphological), 6 asymmetric cost matrices,
differential weights, missing data, and one polymorphic coding.

The six models encode different types of linguistic change: irreversible
innovation (`one_way`), lenition bias (`fort_len`), a three-state genitive
allomorphy system (`third_pl`), and an unusual case where forward change is
costlier than reversal (`four_w`). Twelve cells are coded as missing (`~`),
mostly in Lydian and Palaic where the attestation is fragmentary. Palaic
`gen_ers` is polymorphic: 60% `no_ers`, 40% `ers_primary`.

The diagnostics target the broader Luwic clade (Lydian, Luwian, Lycian) rather
than the Luwian+Lycian pair, because the latter receives near-100% support and
cannot be meaningfully decomposed.

What to look for:

- The Luwian+Lycian clade should receive strong support; the broader (Lydian,
  Luwian, Lycian) grouping should be lower
- Which characters drive the broader clade? (clade-character support)
- Is the grouping robust to removing any single character? (sensitivity)
- Which characters show homoplasy? (consistency index)
- The ASR entropy table for three-state characters like `gen_ers` and `encl_du`

The full example is in `examples/anatolian_full.yaml`.


## YAML reference

### Top-level structure

```yaml
name: <string>            # required, analysis identifier
taxa: <list>              # required, at least 4
characters: <list>        # required, at least 1
models: <mapping>         # optional, custom cost matrices
analysis: <mapping>       # optional, MCMC settings
diagnostics: <list>       # optional, post-run analyses
```

### Taxa

```yaml
taxa:
  - id: hit               # required, short identifier used in data
    name: Hittite          # optional, display name
```

### Characters

```yaml
characters:
  - id: feature_name       # required
    states: [a, b, c]      # optional, explicit state names
    model: model_name       # optional, reference to models section or built-in
    weight: 2.0             # optional, fixed weight (default 1.0)
    weight: estimated       # optional, sampled under Gamma(2,2) prior
    data:                   # required, taxon -> state mapping
      taxon_id: state_name
      taxon_id: ~                          # missing
      taxon_id: {state: 0.8}              # partial confidence
      taxon_id: {a: 0.6, b: 0.4}         # polymorphic
```

### Models

```yaml
models:
  model_name:
    costs:
      state_a: [0,   1.0]     # row: costs from state_a to each state
      state_b: [2.0, 0  ]     # diagonal must be 0
```

Built-in templates (no definition needed): `loss`, `lenition`, `analogy`.

### Analysis

```yaml
analysis:
  coding: variable            # required for yvyra
  seed: 1305                  # optional, PRNG seed
  iterations: 10000           # optional, MCMC generations
  sample_frequency: 50        # optional, record every Nth
  runs: 1                     # optional, independent runs
  chains: 1                   # optional, chains per run
  print_frequency: 1000       # optional, screen output interval
  diagnostic_frequency: 5000  # optional, diagnostics interval
  sitelikes: yes              # optional, per-site log-likelihoods
```

### Diagnostics

```yaml
diagnostics:
  - chardiag: [Taxon1, Taxon2]      # clade-character support
  - sensitivity: [Taxon1, Taxon2]   # leave-one-out clade posterior
  - asr_entropy                      # ancestral state probabilities
```

Clade members are specified by taxon *name* (not id).


## Running yvyra

### Command line

```sh
# Build
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)

# Run
./build/src/yvyra analysis.yaml
```

### Browser

Open the web interface, paste or load a YAML file, and click Run. The MCMC
engine runs as WebAssembly in the browser tab.

### Output files

For an analysis named `mydata`, yvyra produces:

- `mydata.nex.p` -- parameter samples (log-likelihood, tree length, alpha)
- `mydata.nex.t` -- tree samples (Newick format)
- `mydata.nex.w` -- character weight samples (if estimated weights used)
- `mydata.nex.slk` -- per-site log-likelihoods (if `sitelikes: yes`)
- `mydata.nex.con.tre` -- majority-rule consensus tree
