# Project-local number systems

This directory holds number systems that extend the
[Universal](https://github.com/stillwater-sc/universal) number library and are
maintained in this repository rather than in the upstream submodule. The
`include/` tree mirrors Universal's layout (`include/sw/universal/number/...`),
so the headers are reached exactly like any built-in Universal type:

```cpp
#include <universal/number/tau/tau.hpp>
```

The include root `src/number-systems/include/sw` is added to the build by the
top-level `CMakeLists.txt` (`MPADAO_NUMBER_SYSTEMS_INCLUDE_DIR`).

## `tau` — the balanced tau (phi²) number system

`tau<nbits, rbits, bt>` implements the *ternary phi-representation* of the
Generalized Phi Number System (GPNS) introduced in

> S. Agaian and Y. Zhou, *"Generalized Phi Number System and its Applications
> for Image Decomposition and Enhancement"*, Proc. SPIE-IS&T Electronic
> Imaging, Vol. 7881, 78810M, 2011.

It is a positional number system whose base is the square of the golden ratio

```
tau = phi^2 = (3 + sqrt(5)) / 2  ≈ 2.6180339887...
```

(`tau` is the positive root of `x² − 3x + 1 = 0`, hence the defining identity
`tau + 1/tau = 3`) and whose digit set is the balanced-ternary set
`{ −1, 0, +1 }`. A value is

```
X = Σ_p  d_p · tau^(p − rbits),   d_p ∈ {−1, 0, +1}
```

where `p ∈ [0, nbits)` indexes the trits, the `rbits` least-significant trits
form the fractional part and the remaining `nbits − rbits` trits form the
integer part.

### Template parameters

| parameter | meaning                                            |
|-----------|----------------------------------------------------|
| `nbits`   | total number of balanced-ternary trits             |
| `rbits`   | number of fractional trits (positions −1 … −rbits) |
| `bt`      | block storage type (default `std::uint32_t`)       |

Each trit is packed into two bits of a `blockbinary` (magnitude bit + sign bit).

### Properties

- **Symmetric range.** Because the digit set is symmetric there is no separate
  sign field and `maxneg() == -maxpos()`.
- **`constexpr`-capable.** Construction, conversion, arithmetic and comparison
  can all be evaluated in constant expressions (the C++23 requirement for an
  arithmetic type), using a `constexpr` integer power of the irrational base.
- **Rounded arithmetic.** `tau` is an irrational base, so the representable set
  is not closed under the field operations. Arithmetic is evaluated in extended
  precision and rounded to the nearest representable value, saturating at the
  extremes — the trit encoding remains the canonical source of truth.

### Following the Universal pattern

The header set mirrors the structure of an existing Universal type such as
`cfloat`/`lns`:

```
tau/
├── tau.hpp            umbrella header
├── tau_fwd.hpp        forward declarations
├── tau_impl.hpp       class template + operators
├── tau_traits.hpp     is_tau / enable_if_tau
├── exceptions.hpp     arithmetic / internal exceptions
├── numeric_limits.hpp std::numeric_limits specialization
├── manipulators.hpp   type_tag, to_triplets, to_binary, range, color_print, ...
├── attributes.hpp     sign/scale/range helpers
└── mathlib.hpp        elementary math functions
```

A runnable demonstration lives in `src/apps/mixed-precision/tau/` and regression
tests in `test/src/tau_tests.cpp`.
