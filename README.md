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
| `HPolytope` | `sample`, `volume`, `is_in` | ✓ implemented |

---

## 2. Requirements

`pkg install` compiles the oct-file bindings from source, so you need:

- **GNU Octave ≥ 9.1**, including its development files (`mkoctfile` and
  headers).  On Debian/Ubuntu these are in the `octave-dev` package, on
  Fedora in `octave-devel`.
- **A C++17 compiler** (e.g. `g++`) and **make**.

All C++ dependencies (the volesti library, Eigen, Boost, and Spectra) are
vendored in the repository at pinned versions, so installing needs no network
access beyond the download itself and no other system libraries.

CI tests the package on Linux with Octave 9.1 and 11.3.

*Maintainers only:* re-vendoring the dependencies with `tools/fetch_deps.sh`
additionally requires **cmake**, **git**, and **wget** (or **curl**).

---

## 3. Install

### 3.1. Latest release

Install the latest release (v0.0.1) from the Octave command prompt by typing
```
pkg install "https://github.com/GeomScale/ovolesti/archive/refs/tags/v0.0.1.tar.gz"
```

### 3.2. Latest dev version

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

## The square [-1, 1]² in R²
A = [1 0; -1 0; 0 1; 0 -1];
b = [1; 1; 1; 1];
P = HPolytope (A, b);

## Volume estimate
vol = volume (P)

## Test point membership
is_in (P, [0; 0])    # true (interior)
is_in (P, [2; 0])    # false (outside)

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
to reproduce the problem: Octave version, compiler/platform, and a minimal code snippet.

---

## 7. Contribute

`ovolesti` is open source under GPLv3+.  Contributions are welcome: 
oct-file bindings for new polytope types, additional volesti algorithms,
better tests, or documentation improvements. 

---

## 8. Funding

This project's page at NLnet: <https://nlnet.nl/project/Ovolesti/>.

<a href="https://nlnet.nl"><img src="https://nlnet.nl/logo/banner.svg" alt="Logo NLnet: abstract logo of four people seen from above" height="60"></a>
&nbsp;
<a href="https://nlnet.nl/core"><img src="https://nlnet.nl/image/logos/NGI0CommonsFund_tag.svg" alt="Logo NGI Zero Commons Fund: letterlogo shaped like a tag" height="60"></a>

This project is funded through the [NGI0 Commons Fund](https://nlnet.nl/commonsfund),
a fund established by [NLnet](https://nlnet.nl) with financial support from the
European Commission's [Next Generation Internet](https://ngi.eu) programme, under
the aegis of [DG Communications Networks, Content and Technology](https://commission.europa.eu/about-european-commission/departments-and-executive-agencies/communications-networks-content-and-technology_en)
under grant agreement N<sup>o</sup> [101135429](https://cordis.europa.eu/project/id/101135429).
Additional funding is made available by the
[Swiss State Secretariat for Education, Research and Innovation](https://www.sbfi.admin.ch/sbfi/en/home.html) (SERI).
