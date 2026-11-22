# AGENTS.md

Instructions for coding agents working in this repository.

## Project

PBASEX-Julia is a Julia implementation of the PBASEX (Polar Basis Set Expansion)
algorithm for inverting 2D photoelectron/photoion images from velocity map
imaging experiments. It is a port of CPBASEX (MATLAB/Python), after Garcia,
Nahon & Powis, *Rev. Sci. Instrum.* 75, 4989 (2004).

Developed against Julia 1.11; 1.10 is also supported. The module is `pBASEX`,
defined in `src/pbasex_core.jl`.

This repository is a Julia project **environment**, not a registered package:
`Project.toml` intentionally carries no `name`/`uuid`. See the comment at the top
of `Project.toml` for the reason. The public API is therefore loaded with:

```julia
include("src/pbasex_core.jl")
using .pBASEX
```

## Setup

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

`Manifest.toml` is generated locally and is gitignored; do not commit it.

## Tests

```bash
julia --project=. test/runtests.jl
```

`Pkg.test()` does not apply, because there is no package name/UUID.

The unit tests deliberately avoid `findG`, the line-of-sight quadrature over
every pixel and basis function that dominates runtime (minutes to hours at
512x512). Keep it that way: tests must stay fast. If you add coverage that needs
a real basis, use a tiny grid and a cached basis file under a temporary
directory.

## Git workflow

Branch from main using feature/<name> or fix/<name>. The owner has
authorized agents to commit, push the branch, and open a PR. Do not merge the
PR or fast-forward local main; after opening a PR, wait for the owner to
review and merge. All changes reach main through PRs. Run CI only when
requested. Preserve unrelated edits and dataset settings.

If direct GitHub access fails, use `HTTPS_PROXY=http://127.0.0.1:7897`.

## Conventions

- Keep changes scoped to what the task asks for. Do not reformat, rename, or
  "clean up" source files that the task does not concern.
- `data/`, `*.h5`, `*.hdf5`, `Manifest.toml`, and `*.txt` (except
  `examples/image.txt`) are gitignored. Cached basis files and analysis output
  must never be committed.
- `examples/image.txt` is a real 1023x1023 experimental image, about 25 MB of
  text. Treat it as read-only input data: do not regenerate, resample, trim, or
  replace it, and do not add other dataset files to the repository.
- Each `src/*.jl` file is its own module and several are `include`d from more
  than one place, so load order matters. Adding a new source file means adding
  it to the `include` list in `src/pbasex_core.jl`.
- Docstrings use Julia's `"""..."""` convention with a signature line, an
  `# Arguments` section, and a `# Returns` section. Match the surrounding style.
