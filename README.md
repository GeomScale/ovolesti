# ovolesti

[![CI](https://github.com/GeomScale/ovolesti/actions/workflows/ci.yml/badge.svg)](https://github.com/GeomScale/ovolesti/actions/workflows/ci.yml)
[![Octave](https://img.shields.io/badge/Octave-%3E%3D9.1-blue)](https://www.octave.org)
[![License](https://img.shields.io/badge/License-GPLv3%2B-green)](./LICENSE)

Octave binding for the [volesti](https://github.com/GeomScale/volesti) C++ library —
high-dimensional sampling and volume computation for convex polytopes.

---

## 1. About

`ovolesti` wraps volesti's templated C++ algorithms behind Octave classdef
objects (`HPolytope`, *VPolytope*, *Zonotope*, *Spectahedron* — the latter three
are planned).  Each class exposes `sample` and `volume` methods that dispatch
through compiled oct-file bindings to the underlying C++ library.

The package is in **early development**.  The current release provides:

| Class | Methods | Status |
|-------|---------|--------|
| `HPolytope` | `sample`, `volume` | ✓ implemented |

---

## 2. Requirements

- **GNU Octave ≥ 9.1**
- **cmake**, **git**, **wget** (or **curl**) — needed only to run the
  dependency-fetching script; not needed at `pkg install` time if installing
  from a release tarball.

All C++ dependencies (Eigen, Boost, Spectra, and the volesti library itself)
are vendored by the fetch script at pinned versions — no system packages
beyond Octave are required at build time.

---

## 3. Install

### 3.1. Latest dev version

Install the latest dev version from the Octave command prompt by typing
```
pkg install "https://github.com/GeomScale/ovolesti/archive/refs/heads/main.zip"
```
Load the package by typing

```
pkg load ovolesti
```

---

## 4. Usage

```octave
pkg load ovolesti

## A unit cube in R²
A = [1 0; -1 0; 0 1; 0 -1];
b = [1; 1; 1; 1];
P = HPolytope (A, b);

## Volume estimate
vol = volume (P)

## Draw 100 uniform samples
pts = sample (P, 10, 100, [0; 0]);
plot (pts(1,:), pts(2,:), '.');
```

Run the test suite:

```octave
pkg test ovolesti
```


---

## 5. Documentation

Every public function and class method carries a texinfo help block.  In Octave:

```octave
help HPolytope
help HPolytope.sample
```

---

## 6. Provide feedback

Please report bugs, missing features, or volesti capabilities not yet exposed
through the Octave binding on the [GitHub Issues
tracker](https://github.com/GeomScale/ovolesti/issues).  Include enough detail
to reproduce the problem: Octave version, platform, and a minimal code snippet.

---

## 7. Contribute

`ovolesti` is open source under GPLv3+.  Contributions are welcome —
oct-file bindings for new polytope types, additional volesti algorithms,
better tests, or documentation improvements. 