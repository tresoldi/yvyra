# Changelog

## 0.2.0 (2026-04-01)

### Bug fixes

- **Fix rooting artifact in YAML analyses.** In unrooted trees, the first
  taxon in the matrix was always pinned as the calculation root (a MrBayes
  implementation detail), making it appear as sister to all other taxa in the
  consensus tree regardless of the data. For linguistic datasets without a
  natural outgroup, this produced spurious topologies -- e.g., three languages
  with identical feature values would never all cluster together because one
  was silently locked at the root.

  For unrooted analyses, a dummy outgroup taxon (`_outgroup_`) with all-missing
  data is now automatically injected and pruned from output. However, the
  default tree model is now `clock` (see below), which avoids the issue
  entirely.

### New features

- **Rooted clock trees (new default).** YAML analyses now use a rooted
  birth-death tree with an IGR (independent gamma rates) relaxed clock by
  default. The root is inferred by the model -- no outgroup needed. This
  produces more intuitive results for linguistic data: identical taxa form
  clades, and branch lengths represent relative time. Set
  `tree_model: unrooted` in the `analysis:` section for the old behavior.

- **Clock model selection.** Set `clock: strict|igr|iln|tk02|cpp|wn` in the
  `analysis:` section to choose a relaxed clock model. Default is `igr`.

- **Explicit outgroup support (unrooted only).** Set `outgroup: <taxon_id>`
  in the `analysis:` section to designate a specific outgroup for unrooted
  trees.

- **Auto-fallback for tipweights.** Datasets with confidence-weighted or
  polymorphic codings automatically use unrooted trees, because clock tree
  likelihoods are not yet compatible with tipweights. A note is printed when
  this happens.

Thanks to Deepthi Gopal for reporting the rooting bug.

## 0.1.0 (2026-03-14)

Initial public release. Working fork of MrBayes 3.2.7a with molecular models
removed and linguistic phylogenetics extensions (YAML input, cost matrices,
tipweights, character diagnostics, WASM interface).
