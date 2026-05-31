# WENOExt Patch Notes — GCC 15 / OpenFOAM-v2512 Compatibility

## Environment

| Item | Version |
|------|---------|
| OS | Linux 7.0.0-22-generic |
| Compiler | GCC 15.2.0 |
| OpenFOAM | v2512 (ESI) |
| CMake | 4.2 |
| Catch2 | bundled in `tests/Catch2/` |

---

## Background

Building WENOExt against OpenFOAM-v2512 with GCC 15.2.0 requires fixes in four areas:

1. **C++ standard upgrade** — OpenFOAM-v2512 uses C++17 features throughout its headers (`std::is_same_v`, `if constexpr`, `std::string_view`, structured bindings, `std::align_val_t`, etc.). WENOExt was previously configured for C++14.
2. **Missing `<cstdint>` in Catch2** — GCC 15 enforces stricter header self-sufficiency. The bundled Catch2 relied on `uint8_t`/`uint32_t`/`uint64_t` being transitively available, which GCC 15 no longer guarantees.
3. **Deprecated OpenFOAM API** — Several OpenFOAM APIs used in WENOExt were deprecated in recent versions (`Pstream::scatterList`, `autoPtr::set`, `argList::optionReadIfPresent`).
4. **Code quality fixes** — Dangling reference bug, uninitialized variables, signed/unsigned comparison, array bounds bug in tests.

---

## Fix 1 — Raise C++ standard to 17

**File:** `CMakeLists.txt`

| Line | Before | After |
|------|--------|-------|
| 43 | `# Enforce C++ 14 standard` | `# Enforce C++ 17 standard` |
| 44 | `# C++ 14 is required for blaze` | `# C++ 17 is required for OpenFOAM-v2312+` |
| 45 | `set(CMAKE_CXX_STANDARD 14)` | `set(CMAKE_CXX_STANDARD 17)` |

Blaze 3.8 is fully compatible with C++17.

---

## Fix 2 — Suppress third-party warnings in compile flags

**File:** `CMakeLists.txt` — `add_compile_options(...)` line

Added two flags to suppress warnings from OpenFOAM and Blaze headers that cannot be fixed in WENOExt:

| Flag | Reason |
|------|--------|
| `-Wno-overloaded-virtual` | GCC 15 extended `-Woverloaded-virtual` to `operator=`; triggers on OpenFOAM's `fvPatchField` hierarchy (intentional design, not a bug) |
| `-Wno-dangling-reference` | Blaze 3.8 internal `DMatDVecMultExpr.h` triggers this; cannot be fixed without modifying Blaze |

---

## Fix 3 — Add `#include <cstdint>` to Catch2 sources

GCC 15 no longer provides `uint8_t`/`uint32_t`/`uint64_t` via transitive includes. Added `#include <cstdint>` after the last existing `#include` in each of the following 24 files:

| File |
|------|
| `tests/Catch2/src/catch2/internal/catch_string_manip.hpp` |
| `tests/Catch2/src/catch2/catch_test_case_info.hpp` |
| `tests/Catch2/src/catch2/internal/catch_xmlwriter.cpp` |
| `tests/Catch2/src/catch2/catch_totals.cpp` |
| `tests/Catch2/src/catch2/catch_config.hpp` |
| `tests/Catch2/src/catch2/catch_timer.cpp` |
| `tests/Catch2/src/catch2/catch_test_case_info.cpp` |
| `tests/Catch2/src/catch2/catch_config.cpp` |
| `tests/Catch2/src/catch2/internal/catch_random_seed_generation.cpp` |
| `tests/Catch2/src/catch2/reporters/catch_reporter_streaming_base.hpp` |
| `tests/Catch2/src/catch2/matchers/catch_matchers_floating_point.hpp` |
| `tests/Catch2/src/catch2/internal/catch_floating_point_helpers.cpp` |
| `tests/Catch2/src/catch2/internal/catch_test_case_registry_impl.cpp` |
| `tests/Catch2/src/catch2/internal/catch_run_context.cpp` |
| `tests/Catch2/src/catch2/reporters/catch_reporter_combined_tu.cpp` |
| `tests/Catch2/src/catch2/internal/catch_random_number_generator.cpp` |
| `tests/Catch2/src/catch2/interfaces/catch_interfaces_config.hpp` |
| `tests/Catch2/src/catch2/reporters/catch_reporter_multi.hpp` |
| `tests/Catch2/src/catch2/reporters/catch_reporter_compact.cpp` |
| `tests/Catch2/src/catch2/interfaces/catch_interfaces_reporter.hpp` |
| `tests/Catch2/src/catch2/reporters/catch_reporter_console.cpp` |
| `tests/Catch2/src/catch2/reporters/catch_reporter_multi.cpp` |
| `tests/Catch2/src/catch2/reporters/catch_reporter_cumulative_base.hpp` |
| `tests/Catch2/src/catch2/reporters/catch_reporter_event_listener.hpp` |

---

## Fix 4 — Deprecated OpenFOAM API in WENOExt library code

### 4a — `Pstream::scatterList` → `Pstream::broadcastList`

Deprecated since OpenFOAM-v2512 (2025-03). Replaced in three locations:

| File | Line |
|------|------|
| `libWENOEXT/WENOBase/WENOBase.C` | 346 |
| `libWENOEXT/WENOBase/globalfvMesh.C` | 124 |
| `libWENOEXT/WENOBase/globalfvMesh.C` | 193 |

### 4b — `autoPtr::set()` → `autoPtr::reset()`

Deprecated since OpenFOAM 2022-01. Replaced in two locations:

| File | Lines |
|------|-------|
| `libWENOEXT/WENOBase/WENOBase.C` | 527, 537 |

### 4c — `argList::optionReadIfPresent()` → `argList::readIfPresent()`

Deprecated since OpenFOAM 2018-01. Replaced in three locations:

| File | Line |
|------|------|
| `utilities/writeWENOStats/writeWENOStats.C` | 21 |
| `utilities/writeStencilCells/writeStencilCells.C` | 30 |
| `utilities/writeStencilCells/writeStencilCells.C` | 31 |

---

## Fix 5 — Dangling reference bug in WENOUpwindFit

**File:** `libWENOEXT/WENOUpwindFit/WENOUpwindFit.C` — line 130

`patchNeighbourField()` returns `tmp<Field<Type>>` (a temporary object). Binding a `const&` directly to `()` on the temporary caused a dangling reference after the expression ended.

```cpp
// Before (bug: tmp destroyed at end of expression, vfN dangles)
const Field<Type>& vfN = (vf.boundaryField()[patchI].patchNeighbourField())();

// After (tmp lifetime extended by named variable)
const tmp<Field<Type>> tvfN = vf.boundaryField()[patchI].patchNeighbourField();
const Field<Type>& vfN = tvfN();
```

---

## Fix 6 — Uninitialized variables in realEigenValues.H

**File:** `libWENOEXT/WENOBase/geometryWENO/realEigenValues.H` — line 222

GCC 15's `-Wmaybe-uninitialized` flagged `r` and `q` in the QR iteration loop. The variables are logically always initialized before use (set in the preceding `m`-loop), but the compiler cannot trace the control flow. Explicit initialization eliminates the warning and is safe.

```cpp
// Before
double z,y,x,w,v,u,t,s,r,q,p,anorm=0.0;

// After
double z=0,y=0,x=0,w=0,v=0,u=0,t=0,s=0,r=0,q=0,p=0,anorm=0.0;
```

---

## Fix 7 — Signed/unsigned comparison in matrixDB.C

**File:** `libWENOEXT/WENOBase/matrixDB.C` — line 120

`blaze::size()` returns `std::size_t` (unsigned); `A.size()` returns `Foam::label` (signed `int`). Added explicit cast.

```cpp
// Before
if (blaze::size(cmpA) == A.size())

// After
if (blaze::size(cmpA) == static_cast<std::size_t>(A.size()))
```

---

## Fix 8 — Array out-of-bounds in test main

**File:** `tests/src/main.C` — line 79

`malloc(sizeof(char*))` allocated space for only 1 pointer, but the code wrote to both `argvOF[0]` and `argvOF[1]`.

```cpp
// Before (undefined behaviour: only 1 slot allocated, 2 written)
char **argvOF = static_cast<char**>(malloc(sizeof(char*)));

// After
char **argvOF = static_cast<char**>(malloc(2 * sizeof(char*)));
```

---

## Fix 9 — Signed/unsigned comparison in List3D-Test.C

**File:** `tests/src/List3D-Test.C` — lines 61, 64

Loop variables `int i`/`int j` compared against `std::vector::size_type` (`size_t`).

```cpp
// Before
for (int i = 0; i< vecMatrix.size(); ++i)
    for (int j=0; j < vecMatrix[i].size(); j++)

// After
for (size_t i = 0; i < vecMatrix.size(); ++i)
    for (size_t j = 0; j < vecMatrix[i].size(); ++j)
```

---

## Fix 10 — Deprecated API in test code

**File:** `tests/src/globalFvMesh-Test.C` — line 104

```cpp
// Before
Pstream::scatterList(allCellCenters);

// After
Pstream::broadcastList(allCellCenters);
```

---

## Build result

Final build: **100% success, zero warnings, zero errors.**

All targets built and installed:

| Target | Location |
|--------|---------|
| `libWENOEXT.so` | `$FOAM_USER_LIBBIN/` |
| `writeWENOStats` | `$FOAM_USER_APPBIN/` |
| `writeStencilCells` | `$FOAM_USER_APPBIN/` |
| `WENO_TEST` | `tests/src/` |
| `performanceRun.exe` | `tests/Cases/performanceTest/src/` |

---

## Fix 11 — CI support for OpenFOAM v2506 and v2512

The existing CI only covered v1912, v2006, v2012 (ESI) and OF5/7/8 (Org). Added CI support for v2506 and v2512.

### New files

#### `CI/Dockerfile.OFv2506` and `CI/Dockerfile.OFv2512`

Both use the same structure. Key differences from older Dockerfiles:

| Item | Old (v2012 and earlier) | New (v2506/v2512) |
|------|------------------------|-------------------|
| Base image | `ubuntu:18.04` (GCC 7) | `ubuntu:22.04` (GCC 11) |
| OpenFOAM install | Pre-built Docker image or source build | Official ESI apt repository |
| C++ standard support | C++14 only | C++17 (required by OpenFOAM v2312+) |
| `foamDotFile` path | `/home/gitlab/OpenFOAM/OpenFOAM-vXXXX/etc/bashrc` | `/usr/lib/openfoam/openfoamXXXX/etc/bashrc` |

Installation uses the official ESI repository:
```dockerfile
RUN curl -s https://dl.openfoam.com/add-debian-repo.sh | bash \
 && apt-get install -y openfoam2512-dev
```

### Modified files

#### `Makefile`

Added two new tag variables and targets, updated `all` and `clean`:

```makefile
DOCKER_TAG_OF_v2506=wenotest:v2506
DOCKER_TAG_OF_v2512=wenotest:v2512

runTestsOFv2506:
    docker build --rm -t ${DOCKER_TAG_OF_v2506} -f CI/Dockerfile.OFv2506 .

runTestsOFv2512:
    docker build --rm -t ${DOCKER_TAG_OF_v2512} -f CI/Dockerfile.OFv2512 .
```

Also fixed a pre-existing typo in the `clean` target (`$ ` spurious character removed before `${DOCKER_TAG_OF_8}`).

#### `.github/workflows/c-ofESI.yml`

Added two new jobs after `build-OF2012`:

```yaml
build-OF2506:
  name: OpenFOAM v2506
  runs-on: ubuntu-latest
  steps:
  - uses: actions/checkout@v2
  - name: OpenFOAMv2506
    run: make runTestsOFv2506

build-OF2512:
  name: OpenFOAM v2512
  runs-on: ubuntu-latest
  steps:
  - uses: actions/checkout@v2
  - name: OpenFOAMv2512
    run: make runTestsOFv2512
```
