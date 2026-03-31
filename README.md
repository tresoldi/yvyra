# yvyra

Bayesian inference of phylogeny for linguistics and culture.

yvyra (Guarani: "tree") adapts the Mk model framework (Lewis 2001) and
the MrBayes MCMC sampler (Ronquist et al. 2012) for historical linguistics,
stemmatology, and cultural phylogenetics. It adds asymmetric evolution models,
character weights, confidence-weighted tips, and per-character diagnostics.

Forked from [MrBayes 3.2.7a](https://github.com/NBISweden/MrBayes),
stripped of all molecular models.

## Features

- **YAML configuration** with named states, cost matrices, confidence values
- **Asymmetric rate matrices**: directional cost matrices per character, with
  22 built-in templates (loss, lenition, grammaticalization, etc.)
- **Character weights**: fixed or estimated via MCMC under a Gamma prior
- **Confidence-weighted tips**: partial confidence and polymorphic codings
- **Per-character diagnostics**: clade-character support, leave-one-out
  sensitivity (importance sampling), consistency index, ASR entropy
- **Per-model-group ascertainment correction** for heterogeneous USERTYPE models
- Lewis (2001) Mk model with `coding=variable`
- Gamma rate variation, ordered/unordered characters
- Metropolis-coupled MCMC (MC³) with convergence diagnostics
- **Browser interface**: full engine compiled to WebAssembly, no installation

## Quick start: command line

```bash
# Build
cmake -B build
cmake --build build -j$(nproc)

# Run with YAML
./build/src/yvyra analysis.yaml

# Run with NEXUS
./build/src/yvyra analysis.nex

# Test
python3 testing/run_tests.py --binary ./build/src/yvyra
```

## Quick start: web interface

The web interface runs at [yvyra.tresoldi.org](https://yvyra.tresoldi.org)
with no installation required. To serve it locally:

```bash
# Build WASM (requires Emscripten SDK)
source ~/emsdk/emsdk_env.sh
./web/build-wasm.sh

# Serve
cd web && python3 -m http.server 8080
```

Then open `http://localhost:8080`.

## YAML input format

yvyra accepts YAML configuration files with named states, cost matrices,
confidence values, and diagnostic commands:

```yaml
name: example
taxa:
  - id: hit
    name: Hittite
  - id: luw
    name: Luwian
  - id: lyc
    name: Lycian

models:
  one_way:
    costs:
      conservative:  0      1.0
      innovative:    100.0  0

characters:
  - id: lenition
    model: one_way
    data:
      hit: conservative
      luw: innovative
      lyc: innovative

analysis:
  coding: variable
  seed: 1305
  iterations: 50000

diagnostics:
  - chardiag: [Luwian, Lycian]
```

YAML files are converted to NEXUS internally. Standard NEXUS input is
also supported.

## Example: Anatolian languages

The full Billing & Elgh (2023) dataset: 27 characters across 5 Anatolian
languages with 6 asymmetric evolution models, differential weights,
missing data, and polymorphic codings.

```
   /--- Hittite
   |
   |--- Palaic
   +
   |          /--- Lydian
   \---83-----+
              |          /--- Luwian
              \---100----+
                         \--- Lycian
```

Luwian-Lycian clade at 100%, Lydian grouping with Luwic at 83%.

## Cost matrices

Linguists reason in terms of costs (higher = harder change). The `usercost`
command (or `costs:` in YAML) takes a cost matrix and converts to rates
internally (`rate = 1/cost`):

| Pattern | Cost matrix | Meaning |
|---------|-------------|---------|
| One-way innovation | `0 1 / 100 0` | Forward cost 1, reversal cost 100 |
| Fortition/lenition | `0 1 / 2 0` | Lenition cost 1, fortition cost 2 |
| Irreversible 3-state | `0 1 1 / 100 0 1 / 100 100 0` | Innovation can't reverse |

Built-in templates: `loss`, `lenition`, `grammaticalization`, `word_order`,
`sound_merger`, and others. See the web interface guide for the full list.

## Requirements

- C99 compiler (gcc, clang)
- CMake 3.16+
- Optional: readline/libedit (interactive mode)
- Optional: Emscripten SDK (WASM build)
- Python 3 (test suite, stdlib only)

## License

GNU General Public License v3 or later.

Copyright (C) 2026 Tiago Tresoldi.
Based on MrBayes 3.2.7a, Copyright (C) 2002-2023 John P. Huelsenbeck,
Fredrik Ronquist et al.

## Citation

If you use yvyra, please cite both the tool and the underlying methodology:

- Tresoldi, T. (2026). yvyra: Bayesian phylogenetic inference for linguistic
  and cultural data. Software available at https://yvyra.tresoldi.org.
- Ronquist, F., et al. (2012). MrBayes 3.2: Efficient Bayesian phylogenetic
  inference and model choice across a large model space. *Systematic Biology*,
  61(3), 539-542.
- Lewis, P. O. (2001). A likelihood approach to estimating phylogeny from
  discrete morphological character data. *Systematic Biology*, 50(6), 913-925.
