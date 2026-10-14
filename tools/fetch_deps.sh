#!/usr/bin/env bash
#
# fetch_deps.sh — download volesti source + pinned C++ dependencies
#
# Fetches everything needed to build the ovolesti oct-files so that the
# resulting tree is self-contained and installs offline via "pkg install".
#
# Sources downloaded:
#   src/volesti/                        volesti C++ library (tarball from GitHub)
#   src/volesti/external/_deps/eigen-src    Eigen 3.4.0
#   src/volesti/external/_deps/boost-src    Boost 1.84.0 (stripped to ~15 libs)
#   src/volesti/external/_deps/spectra-src  Spectra v1.1.0
#
# The volesti commit is pinned below — bump it to update the binding's view
# of the C++ library.  Eigen / Boost / Spectra versions are pinned inside
# volesti's own cmake-files/ modules.
#
# Prerequisites: bash, cmake, git, wget (or curl), tar
#
# Usage:
#   tools/fetch_deps.sh           # first fetch
#   tools/fetch_deps.sh           # subsequent runs are no-ops
#   tools/fetch_deps.sh --clean   # wipe everything, then fetch fresh

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC_DIR="$ROOT/src"
VOLESTI_DIR="$SRC_DIR/volesti"
DEPS_DIR="$VOLESTI_DIR/external/_deps"
CMAKE_DIR="$VOLESTI_DIR/external/cmake-files"

# ---------------------------------------------------------------------------
#  Pin the volesti version here
# ---------------------------------------------------------------------------
VOLESTI_COMMIT="75b76ccd33fe27aaf6a8ed6216fac521889908d2"
VOLESTI_URL="https://github.com/GeomScale/volesti/archive/${VOLESTI_COMMIT}.tar.gz"

# ---------------------------------------------------------------------------
#  --clean: wipe everything so the next run starts from scratch
# ---------------------------------------------------------------------------
if [ "${1:-}" = "--clean" ]; then
  echo "Removing $VOLESTI_DIR ..."
  rm -rf "$VOLESTI_DIR"
fi

# ---------------------------------------------------------------------------
#  Step 1 — download volesti source (if not already present)
# ---------------------------------------------------------------------------
if [ ! -d "$VOLESTI_DIR" ]; then
  echo "Downloading volesti @ ${VOLESTI_COMMIT:0:7} ..."

  TMPDIR=$(mktemp -d)
  trap "rm -rf $TMPDIR" EXIT

  ARCHIVE="$TMPDIR/volesti.tar.gz"

  if command -v wget >/dev/null 2>&1; then
    wget -q --show-progress -O "$ARCHIVE" "$VOLESTI_URL"
  else
    curl -L -o "$ARCHIVE" "$VOLESTI_URL"
  fi

  mkdir -p "$SRC_DIR"
  tar xzf "$ARCHIVE" -C "$TMPDIR"
  mv "$TMPDIR/volesti-${VOLESTI_COMMIT}" "$VOLESTI_DIR"

  trap - EXIT
  rm -rf "$TMPDIR"

  # --- Strip volesti — keep only include/ and external/ ----------------------
  echo "Stripping volesti source ..."
  for item in "$VOLESTI_DIR"/* "$VOLESTI_DIR"/.*; do
    name=$(basename "$item")
    case "$name" in
      .|..|include|external) ;;
      *) rm -rf "$item" ;;
    esac
  done

  echo "volesti source installed at $VOLESTI_DIR"
else
  echo "volesti source already present at $VOLESTI_DIR"
fi

# ---------------------------------------------------------------------------
#  Step 2 — fetch pinned Eigen, Boost, Spectra via volesti's CMake modules
# ---------------------------------------------------------------------------
mkdir -p "$DEPS_DIR"

TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

cat > "$TMPDIR/CMakeLists.txt" << CMAKE_EOF
cmake_minimum_required(VERSION 3.14)
project(ovolesti-fetch-deps NONE)

# Suppress system-path searches so find_path() only looks in _deps/.
# This forces FetchContent to download on the first run regardless of
# whether the host has libeigen3-dev installed.
set(CMAKE_FIND_USE_SYSTEM_ENVIRONMENT_PATH OFF)
set(CMAKE_FIND_USE_CMAKE_SYSTEM_PATH        OFF)

# Volesti's own dependency-fetching modules (version pins live here).
include("${CMAKE_DIR}/Eigen.cmake")
include("${CMAKE_DIR}/Boost.cmake")
include("${CMAKE_DIR}/Spectra.cmake")

GetEigen()
GetBoost()
GetSpectra()
CMAKE_EOF

echo "Fetching Eigen, Boost, Spectra at pinned versions ..."
cmake -S "$TMPDIR" -B "$TMPDIR/build"

# --- Verify _deps/ populated --------------------------------------------------
if [ ! -d "$DEPS_DIR/eigen-src" ] \
    || [ ! -d "$DEPS_DIR/boost-src" ] \
    || [ ! -d "$DEPS_DIR/spectra-src" ]; then
  echo ""
  echo "ERROR: _deps/ not fully populated after cmake run." >&2
  echo "Check the cmake output above for errors." >&2
  exit 1
fi

# --- Clean up cmake build leftovers -------------------------------------------
rm -rf "$DEPS_DIR"/*-build "$DEPS_DIR"/*-subbuild

# --- Delete VCS and CI metadata from cloned dependencies ----------------------
find "$DEPS_DIR" \( -name ".git" -o -name ".github" -o -name ".gitlab" \) \
  -type d -exec rm -rf {} + 2>/dev/null || true
find "$DEPS_DIR" \( -name ".gitignore" -o -name ".gitattributes" \
  -o -name ".hgeol" -o -name ".travis.yml" \) \
  -type f -delete 2>/dev/null || true

# ---------------------------------------------------------------------------
#  Step 3 — strip each dependency to the minimum needed at compile time
# ---------------------------------------------------------------------------

# Eigen — keep only the Eigen/ header tree and the find_path marker file
EIGEN_SRC="$DEPS_DIR/eigen-src"
if [ -d "$EIGEN_SRC" ]; then
  for item in "$EIGEN_SRC"/*; do
    name=$(basename "$item")
    case "$name" in
      Eigen|signature_of_eigen3_matrix_library) ;;
      *) rm -rf "$item" ;;
    esac
  done
fi

# Boost — keep only the ~15 libraries our oct-file actually includes.
# Regenerate this list after changing the oct-file's volesti include set:
#
#   mkoctfile -E -H src/__volesti_HPolytope__.cc 2>&1 1>/dev/null \
#     | grep -oP 'boost/\K[^/]+' | sort -u
BOOST_SRC="$DEPS_DIR/boost-src"
rm -rf "$BOOST_SRC"/{libs,doc,tools,more,status}

BOOST_KEEP_DIRS="
  array        assert       concept      config       core
  detail       exception    integer      io           iterator
  mpl          preprocessor random       range        type_traits
  utility
"

if [ -d "$BOOST_SRC/boost" ]; then
  for d in "$BOOST_SRC/boost"/*/; do
    name=$(basename "$d")
    keep=0
    for k in $BOOST_KEEP_DIRS; do
      [ "$k" = "$name" ] && keep=1
    done
    if [ "$keep" -eq 0 ]; then
      rm -rf "$d"
    fi
  done
fi

# Spectra — keep only the include/ directory
SPECTRA_SRC="$DEPS_DIR/spectra-src"
if [ -d "$SPECTRA_SRC" ]; then
  for item in "$SPECTRA_SRC"/*; do
    name=$(basename "$item")
    case "$name" in
      include) ;;
      *) rm -rf "$item" ;;
    esac
  done
fi

# ---------------------------------------------------------------------------
#  Report
# ---------------------------------------------------------------------------
VOLESTI_SIZE=$(du -sh "$VOLESTI_DIR" 2>/dev/null | cut -f1)
EIGEN_SIZE=$(du -sh "$EIGEN_SRC" 2>/dev/null | cut -f1)
BOOST_SIZE=$(du -sh "$BOOST_SRC" 2>/dev/null | cut -f1)
SPECTRA_SIZE=$(du -sh "$SPECTRA_SRC" 2>/dev/null | cut -f1)
echo ""
echo "All dependencies fetched and stripped:"
echo "  $VOLESTI_DIR  ($VOLESTI_SIZE)"
echo "  └── external/_deps/eigen-src    ($EIGEN_SIZE)"
echo "  └── external/_deps/boost-src    ($BOOST_SIZE)"
echo "  └── external/_deps/spectra-src  ($SPECTRA_SIZE)"
echo ""
echo "The package is ready to build (make -C src)."
