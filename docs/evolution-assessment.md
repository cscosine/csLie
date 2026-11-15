# csLie — Evolution Assessment

Status: draft for discussion · Date: 2026-11-15 · Scope: library at `dev` @ `b271322`

This document checks csLie against its goals, measures the current numerical
behaviour, compares it with existing libraries, and proposes an architecture
and a phased roadmap. It is meant to be kept in the repository and updated as
decisions are made.

---

## 1. Goals

| # | Goal | Short name |
|---|------|------------|
| G1 | Best available **computation speed** for SO3, SE3, S2, SO2, SE2 | *Speed* |
| G2 | Best available **numerical accuracy and correct edge cases**, in **float32 and float64** | *Accuracy* |
| G3 | Scope: **SO3, SE3, S2** for 3D and **SO2, SE2** for 2D, all first-class (Sim3 etc. out of scope) | *Scope* |
| G4 | **CPU vectorization** (SIMD), assessed and favoured | *SIMD* |
| G5 | Usable **on GPU** as well | *GPU* |
| G6 | Header-only vs static library: choose whatever gives the best speed and optimization | *Packaging* |

---

## 2. Where the code is today

### 2.1 Content

~1.2 kLOC of headers plus Catch2 tests (float32 and float64 for every test).

| Group | Available | Missing for a complete Lie API |
|-------|-----------|-------------------------------|
| SO3 | `so3Exp`, `so3Log`, left Jacobian `so3J`, `so3Jinv`, `skewSymm`, RPY conversions | right Jacobians, compose/inverse/act with Jacobians, plus/minus (⊕/⊖), adjoint, interpolation, quaternion conversions, projection back onto SO3 (re-orthonormalization) |
| SE3 | `se3Exp`, `se3Log`, `se3Adjoint`, point-action Jacobian at zero | left/right Jacobians (the Q block), Jacobians of compose/inverse/act/log/exp, ⊕/⊖, interpolation |
| S2 | `s2Exp`, `s2oplus`, tangent basis (`s2toZAlignedSO3`), Jacobian of ⊕ at δ=0 | **log / ⊖** (needed by any estimator), Jacobians of ⊕/⊖ at δ≠0, distance, basis Jacobians |
| SO2 / SE2 | exp/log (SO2), exp + xyYaw (SE2) | `se2Log` (TODO in code) |

API style: **free functions on Eigen matrices** (`SO3<S>` = `Matrix3`, `SE3<S>` =
`Eigen::Transform<S,3,Isometry>`, algebras = `Vector3`/`Vector6`). Tangent
ordering for se3 is `[translation; rotation]`, the same as Sophus and manif.

### 2.2 Packaging fact: it is already header-only

`csLie/CMakeLists.txt` declares `cs_add_library(csLie INTERFACE)`. The library
is **already header-only**: there is no compiled object, and everything is
inlined into the caller. So the question in G6 is really "should it stay
header-only?" (answered in §3.6).

### 2.3 Measured numerical behaviour

A probe compared csLie against a long double reference: 200 random axes per
angle, maximum absolute error in units of machine epsilon of the respective type.
The probe was built with GCC `-O2` (RelWithDebInfo asserts) and is not committed.
§6 proposes turning it into a permanent precision suite.

| θ [rad] | so3Exp f32 | so3Log f32 | se3 log∘exp f32 | so3Exp f64 | so3Log f64 | se3 log∘exp f64 |
|--------:|----:|----:|----:|----:|----:|----:|
| 1e-10 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | **193 554** |
| 1e-7  | 0.0 | 0.0 | 0.2 | 0.3 | 0.0 | **154 007** |
| 1e-5  | 0.0 | 0.0 | 0.2 | 0.3 | 0.0 | **1 588** |
| 1e-3  | 0.2 | 0.0 | 0.5 | 0.3 | 0.0 | 30 |
| 0.017 | 0.3 | 7.3 | 7.3 | 0.4 | 38.6 | 38.8 |
| 0.05  | 0.3 | 11.4 | 11.4 | 0.4 | 15.9 | 3.3 |
| 1     | 0.8 | 1.3 | 1.5 | 1.0 | 1.3 | 1.0 |
| 3     | 2.4 | 4.9 | 12.0 | 2.5 | 6.1 | 16.0 |
| 3.14159 | 2.8 | 3.3 | 4.0 | 2.7 | 2.4 | 3.0 |
| **exactly π** (`diag(1,−1,−1)`) | | **NaN** | | | **NaN** | |

Findings:

1. **`so3Log` returns NaN for an exact 180° rotation** (both precisions). In the
   near-π branch the sign is taken from `signum` of the skew part, which is
   exactly 0 at π, so the code divides by zero. Exact 180° rotations are common
   in practice (e.g. flipped frames, `diag(1,−1,−1)`).
2. **`se3Exp` loses ~1e-11 absolute accuracy (float64) for small angles.** The
   closed form `(1−cos θ)/θ²` is used down to θ = ε, but `1−cos θ` cancels
   completely below θ ≈ 1e-8 (it evaluates to 0 instead of θ²/2). The rotation
   block hides this, because the coefficient is multiplied by θ², but the
   translation (the V matrix) only multiplies it by θ. Isolated test at
   θ = 1e-10, |t| ≈ 1: translation error 6.5e-11, and `se3Log` alone is exact.
   Small angles are the **most common** case in estimation (small increments),
   so this matters.
3. **The float32 path partly computes in double.** Some calls are unqualified
   (`sin`, `cos`) and some are `std::`-qualified. For `float` arguments the
   unqualified calls resolve to `::cos(double)`. The float instantiations link
   both `sincos` (double) and `sincosf`. This is slower, and float accuracy then
   depends on accidental promotion, which is why float32 does not show finding 2.
   It also breaks custom scalar types (Ceres `Jet`, SIMD batches) in the
   `std::`-qualified places.
4. **`so3Log` uses `acos` in the mid range.** The error peaks right above the
   switch point (`cos θ > 0.999856`, θ ≈ 1°): ~39 ε in f64. `acos` amplifies
   input rounding by 1/sin θ. Not catastrophic, but avoidable: `atan2(sin θ, cos θ)`
   is uniformly accurate.
5. **Thresholds are ad-hoc** ("magic numbers" in `so3Log`, `se3Log`, `s2Exp`),
   not derived per precision from truncation error vs. cancellation error.
6. **No input projection**: logs assume a valid rotation (checked only in
   Paranoid). Real inputs drift (products of many rotations); results on
   slightly non-orthonormal input are not specified.

### 2.4 Performance-relevant observations

- **SO3 stored as a 3×3 matrix** (9 scalars), **SE3 as `Eigen::Transform`**
  (16 scalars including the constant bottom row). Composition costs 27 multiply-adds
  for SO3 (quaternion: 16), and normalization needs SVD or Gram–Schmidt
  (quaternion: one rsqrt).
- **The S2 basis uses `atan2`, `cos`, `sin`** where `c = y/r`, `s = z/r`
  gives the same result without transcendental functions.
- Every function branches between Taylor and closed form, which is fine for
  scalar CPU code but costly for SIMD lanes and GPU warps (§3.4, §3.5).
- **No benchmarks** exist, so there is no baseline to claim "best" against.

### 2.5 Strengths worth keeping

- float32 **and** float64 tested everywhere (template tests); Paranoid config
  with UB sanitizer and runtime precondition checks.
- Eigen-native, generic over `MatrixBase` (works with blocks/maps).
- Explicit, documented Jacobian conventions (left Jacobian, ⊕ at δ = 0).
- A tested S2 implementation (few libraries offer S2 at all).

---

## 3. Considerations, goal by goal

### 3.1 G2 — Accuracy and edge cases (first, because speed without it is useless)

**Principle: write every formula so it has no cancellation anywhere, then pick
Taylor thresholds per precision at the crossover point of truncation error vs.
rounding error, and prove it with a precision suite (§6).**

Concrete recommendations:

| Quantity | Robust formulation |
|----------|-------------------|
| Rotation angle from R | `θ = atan2(‖vee(R−Rᵀ)‖/2, (tr R − 1)/2)`, uniformly accurate on [0, π]; or via quaternion `θ = 2·atan2(‖v‖, w)` |
| `(1−cos θ)/θ²` | `½·sinc²(θ/2)` (no cancellation at any θ) |
| `sin θ/θ` | `sinc` with Taylor below a precision-dependent threshold |
| `(θ − sin θ)/θ³`, `(1/θ²)(1 − θ sin θ / (2(1−cos θ)))` (Jl⁻¹, SE3 V⁻¹, Q) | Taylor below the threshold, or a rewritten form using half-angle identities; **never** `1 − sinc` divided by θ² near 0 |
| Log near π | axis from the largest diagonal element of `(R+Rᵀ)/2 − cos θ·I`, or matrix → quaternion with **Shepperd's method** (pick the largest of w,x,y,z), then quaternion log; at **exactly π** pick a deterministic sign (no division by a zero skew part) |
| S2 log / distance | `atan2(‖a×b‖, a·b)`, never `acos(a·b)` |
| S2 tangent basis | algebraic (`y/r`, `z/r`) or GTSAM-style cross product with the least-aligned axis; document the unavoidable discontinuity (hairy-ball theorem) and keep it far from the working point |

Threshold selection: for a Taylor series truncated after degree k, use it while
`θ^(k+1)·|c_{k+1}| < ε/2` (truncation below rounding), and verify that the
closed form above the threshold has error ≲ a few ε. Thresholds differ for
float32 and float64 and must be `constexpr` per scalar type, not shared.

Other rules:

- **Input projection**: provide `normalize()` (quaternion: rescale; matrix:
  polar/SVD or one Newton step `R ← ½(R + R⁻ᵀ)`), and document that log/act
  assume valid input. Optionally `checked*` variants returning `std::optional`.
- **No `-ffast-math`** for this code: it removes NaN/Inf handling and allows
  reassociation that breaks the carefully ordered formulas. `-fno-math-errno`
  is safe and lets `sqrt` inline.
- **Reproducibility**: CPU vs GPU results differ (FMA contraction, different
  libm). Specify accuracy in ε/ulps per function, never bitwise equality.
- **Scalar genericity**: call math through one internal layer
  (`using std::sin; sin(x)` or a traits struct), so float stays float and
  custom scalars (Ceres `Jet`, SIMD batch, CUDA) resolve correctly. This also fixes finding 3.

### 3.2 G3 — Scope and representation

**Scope.** SO3, SE3, S2 for 3D and SO2, SE2 for 2D, all first-class: same API,
same Jacobian conventions, same accuracy and speed requirements. No
Sim3/SE_2(3)/RxSO3.

**SO2 and SE2 (2D).** They are needed for planar robots and odometry, 2D SLAM
and scan matching, particle filters, and image-plane motion (csImgTracker,
csVisOpenGL 2D). They are cheap, but they are often evaluated in very large
numbers (particles, scan points, image patches), so the batch/SIMD path (§3.4)
matters for them as much as for 3D.

| Item | Recommendation |
|------|----------------|
| SO2 storage | unit complex number `(c, s)`: compose = 4 FMA, act = 4 FMA, normalize = 1 rsqrt; the angle is not stored (no wrap-around) |
| SO2 log | `atan2(s, c)`, result in (−π, π] |
| SE2 storage | unit complex + translation (4 scalars) |
| se2 tangent ordering | `[x, y, θ]` (translation first, consistent with se3 `[t; ω]`) |
| SE2 exp | `V = [[A, −B], [B, A]]` with `A = sinθ/θ` and `B = (1−cosθ)/θ = sin(θ/2)·sinc(θ/2)` (cancellation-free); Taylor near 0 with per-precision thresholds |
| SE2 log (missing today) | `V⁻¹ = [[a, b], [−b, a]]` with `a = (θ/2)·cot(θ/2)` (Taylor `1 − θ²/12` near 0) and `b = θ/2`; well-defined up to ±π |
| Jacobians | left/right Jacobians and Jacobians of compose/inverse/act/⊕/⊖, all closed form and small (3×3) |
| 2D/3D embedding | helpers SE2 ↔ SE3 (planar pose at z = 0, yaw about z), used by planar robots in 3D maps |

**Representation, the main speed/accuracy design decision:**

| Storage | Compose | Act on point | Exp/Log | Normalize | Size |
|---------|---------|--------------|---------|-----------|------|
| SO3 as unit quaternion | 16 FMA | ~15–18 FMA (`v + 2w(q×v) + 2q×(q×v)`) | simpler, one `sincos` / one `atan2` | 1 rsqrt | 4 |
| SO3 as 3×3 matrix | 27 FMA | 9 FMA | Rodrigues | SVD / Newton | 9 |

Recommendation: **quaternion as canonical SO3 storage** (and SE3 = quaternion +
translation, 7 scalars), with cheap `matrix()` conversion. To transform many
points, convert once to a matrix and apply it to the whole batch. This is what
Sophus and manif do. Matrix-based free functions stay available for
interoperability with existing code that holds `Matrix3`.

**Strong types.** Today `SO3<S>`, `SO3J<S>` and `SO3Jinv<S>` are the same type
(`Matrix3`), and `S2<S>` and `SO3Alg<S>` are both `Vector3`. Overloading is
impossible and mistakes compile silently. Introduce value classes
(`SO3<S>`, `SE3<S>`, `S2<S>`, tangent types) on top of the free-function kernels.

**Full API** (manif-style naming, "A micro Lie theory" conventions; Jacobians
optional so they cost nothing when not requested):

```cpp
SO3f R = SO3f::exp(w, &J_R_w);       // optional Jacobian out-params
Vec3f w = R.log(&J_w_R);
SO3f C = A.compose(B, &J_C_A, &J_C_B);
Vec3f p2 = R.act(p, &J_p2_R, &J_p2_p);
SO3f Rp = R.plus(dw);  Vec3f d = A.minus(B);   // right ⊕/⊖ (left variants explicit)
Mat3f Jr = SO3f::Jr(w), JrInv = SO3f::JrInv(w), Jl = ..., JlInv = ...;
Mat6f Ad = T.adj();
SO3f I = slerp(A, B, t);
```

Conventions to fix once and document: tangent ordering `[ρ; θ]` (current),
**right** perturbation as default (with left also available), Hamilton
quaternions, `w ≥ 0` canonical form (optional).

### 3.3 G1 — Speed of single operations

For one element, the cost is dominated by **transcendental functions** (each
`sin`/`cos`/`atan2` is ~20–100 cycles) and by branch mispredictions, not by the
~30–150 multiply-adds. Gains come from:

1. **Fewer transcendentals**: one `sincos` per exp; one `atan2` per log
   (no `acos` + `sqrt` + extra `sin`); a trig-free S2 basis.
2. **Shared subexpressions**: exp + Jacobian computed together share `sincos`,
   `sinc` and `θ²`; an API that returns both avoids recomputation.
3. **Quaternion storage** for composition and normalization (§3.2).
4. **Inlining** (header-only, §3.6) so Eigen's fixed-size code is fully unrolled.
5. Horizontal SIMD inside one 3×3 / 4×4 operation gives **little**: 3 is not a
   multiple of the vector width, and Eigen's fixed-size 3×3 is mostly scalar.
   Do not chase it.

### 3.4 G4 — CPU vectorization (the big lever)

The large speedups come from **vertical SIMD**: processing N elements at once
(AVX2: 8 floats / 4 doubles; AVX-512: 16 / 8; NEON: 4 / 2), with **SoA layout**
(all x components contiguous, then all y, …) and **vectorized math functions**.
For transcendental-bound kernels, typical gains are roughly proportional to the
lane count (≈4–8× for AVX2 float). This must be measured (§6), not assumed.

**Key design: write each kernel once, scalar-generic and branch-light**

```cpp
template <class T>                  // T = float, double, simd batch, Jet, CUDA scalar
CSLIE_HD void so3ExpKernel(T wx, T wy, T wz, T (&R)[9]) {
  T th2 = wx*wx + wy*wy + wz*wz;
  auto small = th2 < Thr<T>::exp2;  // mask, not branch
  T th = sqrt(th2);
  auto [s, c] = sincos(th);
  T a = select(small, taylorSinc(th2), s / th);
  ...
}
```

- `select(mask, a, b)` instead of `if`. With SIMD both sides are evaluated and
  blended, which is cheap here because both sides are short. For scalar code it compiles to a
  conditional move or a well-predicted branch.
- Division by zero in the unused lane must be harmless: guard denominators
  (`th = select(small, 1, th)`). The symforce approach of adding a tiny ε
  (`sqrt(θ² + ε)`) is an alternative that removes branches entirely, at the
  price of a bounded O(ε) bias near 0. Use it only where it is provably below rounding.
- **Batch API** over spans / SoA containers:
  `exp(span<const Vec3f> in, span<SO3f> out)`, plus SoA containers for
  zero-transposition hot paths. AoS → SoA transposition in registers is cheap for
  4- or 8-lane blocks.
- **Reusing Eigen code with a batch scalar**: `Eigen::Matrix<Batch,3,3>` with a
  SIMD batch as the scalar type turns Eigen's unrolled fixed-size code into
  vertical SIMD for free (needs a `NumTraits` specialization). The Rust port
  sophus-rs follows the same idea with batch scalar/vector/matrix types
  (behind a `simd` feature that currently needs nightly Rust).

**SIMD backend options** (C++17 project):

| Option | License | Header-only | Runtime CPU dispatch | Vector math (sin/cos/atan2) | Notes |
|--------|---------|:---:|:---:|:---:|-------|
| **xsimd** | BSD-3 | ✅ | manual | ✅ (few-ulp) | used by Arrow/Pythran; good fit for header-only |
| **Google Highway** | Apache-2.0 | ❌ (compiled for dispatch) | ✅ built-in | ✅ (`hwy/contrib/math`) | best portability across x86/ARM/RISC-V/SVE |
| `std::experimental::simd` / C++26 `std::simd` | std | ✅ | ❌ | partial | libstdc++ only today; the long-term standard |
| EVE | BSL-1.0 | ✅ | ❌ | ✅ | requires C++20 |
| SLEEF | BSL-1.0 | ❌ | ✅ | ✅ (reference quality) | math only; pair with another layer |
| Eigen packets | MPL-2.0 | ✅ | ❌ | partial | internal API, not intended for this |
| OpenMP `simd` + branch-free code | — | ✅ | ❌ | depends on compiler/libmvec | zero dependency; fragile, compiler-dependent |

Recommendation: **an internal math/traits layer** (`sqrt`, `rsqrt`, `sincos`,
`atan2`, `select`, `fma`) so the backend is swappable. Start with **xsimd**
(header-only, keeps G6) for the batch API. Consider **Highway** later if
runtime dispatch (one binary for AVX2/AVX-512/NEON) becomes necessary. Its
dispatched kernels then live in an optional compiled component (§3.6).

### 3.5 G5 — GPU

- **Same headers on device**: annotate kernels with `CSLIE_HD`
  (`__host__ __device__` under NVCC/HIP, empty otherwise). Eigen fixed-size types
  work in CUDA (`EIGEN_DEVICE_FUNC`, used by Sophus as `SOPHUS_FUNC`). One
  element per GPU thread is the natural SIMT mapping, and the same branch-light
  kernels avoid warp divergence.
- **No host-only dependencies in kernels**: csCore asserts (libassert,
  exceptions) cannot run on device. Device builds compile asserts out or use a
  trap. The kernel layer must not depend on csCore at runtime.
- **float32 is the GPU precision**: consumer GPUs run float64 at 1/32–1/64 of
  float32 throughput. That makes float32 accuracy (§3.1) a hard requirement,
  not a nice-to-have.
- **Batch GPU kernels** (exp/log/compose/act over device arrays) belong in an
  optional `csLie::cuda` component. Data layout is the same SoA as the CPU batch API.
- Reference: **lietorch** (CUDA kernels for SO3/SE3 under PyTorch) and
  **symforce** (generated CUDA code) show the approach works. NVIDIA Warp and
  PyPose are Python-side references.
- HIP/SYCL portability follows from the same macro approach. Not a goal now.

### 3.6 G6 — Header-only vs static library

It is **already header-only** (§2.2), and that is the right choice for the
core:

- Small functions on fixed-size types only reach full speed when inlined and
  unrolled at the call site. A static library would block that (unless LTO is used
  everywhere).
- Templates over the scalar (float, double, Jet, SIMD batch, CUDA) need headers anyway.

Where compiled parts make sense (optional components, not the core):

| Component | Why compiled |
|-----------|-------------|
| `csLie::batch` dispatched kernels | runtime CPU dispatch (Highway or per-target TUs): one binary for AVX2/AVX-512/NEON |
| `csLie::cuda` | `.cu` batch kernels |
| explicit instantiations (`extern template` for float/double) | compile-time reduction only. Opt-in, and only for heavy non-hot functions |

### 3.7 Dependencies

- **Eigen**: keep as the public interface (interoperability), not as the kernel
  engine. Kernels operate on scalars so they work with batches and on device.
- **csCore**: today used for asserts and shape checks. Keep it for the host-side
  API (Paranoid checks are valuable), but the kernel layer stays dependency-free
  (header-only, no exceptions, no allocation).
- **Relation to csCore: separate repository, depending only on csCore's base
  layer.** csLie is not merged into csCore: it has its own goals (precision
  suite, benchmarks, SIMD/GPU, possibly a public release), its own optional
  dependencies (xsimd, CUDA) that must not land in csCore, and csCamera depends
  on it. csLie links only csCore's base components (assert, types, math traits),
  so when csCore grows utilities (properties, serialization, logging) csLie does
  not pull them. See the csCore assessment, §4.3.

### 3.8 Python bindings (optional)

**Worth adding, as a thin optional layer, once the C++ API is stable (after
phase 2).** The aim is not to compete with jaxlie/PyPose, but to use *this*
implementation from Python:

- **Validation**: the precision suite (§6) can compare against mpmath/SymPy
  references and plot error-vs-angle curves with matplotlib far more easily than in C++.
- **Tools**: notebooks, dataset conversion, quick experiments with the same
  numerics as the C++ code.

How:

- **nanobind** (smaller and faster bindings than pybind11, same author).
- A **batch API on NumPy arrays** (`(N,4)` quaternions, `(N,3)` tangents, …):
  calling one function per element from Python is slow, while batches map
  directly onto the SoA/SIMD path (§3.4).
- Later, **DLPack** to pass GPU arrays to and from PyTorch/JAX without copies (§3.5).
- Packaging with scikit-build-core and cibuildwheel, in a separate optional
  target or package. The C++ build never depends on Python.

Cost: wheel builds, a Python matrix in CI, and keeping the binding API stable.

---

## 4. Existing alternatives

| Library | Groups (of interest) | Language / form | Storage | Jacobians | float32 focus | SIMD batch | GPU | Edge-case care | License |
|---------|---------------------|-----------------|---------|-----------|:---:|:---:|:---:|----------------|---------|
| **csLie (today)** | SO3, SE3, S2, SO2, SE2 | C++17 header-only, Eigen | matrix | partial (left J, S2 at 0) | tested | ❌ | ❌ | issues (§2.3) | own |
| **Sophus** | SO2/3, SE2/3 (+Sim) | C++ header-only, Eigen | quaternion | exp/log Jacobians, left J | ⚠️ | ❌ | ⚠️ `SOPHUS_FUNC` device-annotated | good; Taylor thresholds | MIT |
| **Sophus 2 (C++, in farm-ng-core)** | + unit vectors (S2), splines, camera models | C++ rewrite | quaternion | yes | ⚠️ | not confirmed | ⚠️ | good | verify (farm-ng-core LICENSE); Sophus 1 itself is now in maintenance mode |
| **sophus-rs** (Rust port) | SO2/3, SE2/3 (+ opt., autodiff) | Rust | quaternion | via dual numbers | ⚠️ | ✅ batch types (`simd` feature, nightly) | ❌ | good | MIT / Apache-2.0; API declared unstable |
| **manif** | SO2/3, SE2/3, SE_2(3), bundles | C++11 header-only (+Python) | quaternion | **complete** (right/left for every op) | ⚠️ | ❌ | ❌ | good | MIT |
| **smooth** | SO2/3, SE2/3, bundles | C++20 header-only | quaternion | complete, autodiff | ⚠️ | ❌ | ❌ | good | MIT |
| **Lie++** (aau-cns) | SO3, SE3, SE_2(3), group products | C++ header-only, Eigen | matrix/quat | for filter design (equivariant filters, MSCEqF) | ⚠️ | ❌ | ❌ | not assessed | check repo |
| **GTSAM** (Rot3/Pose3/Unit3) | SO3, SE3, **S2 (Unit3)** | C++ compiled library | matrix (quat opt.) | complete (`OptionalJacobian`) | ❌ double | ❌ | ❌ | **very careful**: `ExpmapFunctor`/`DexpFunctor` with Taylor terms, configurable near-0 and near-π thresholds | BSD |
| **symforce** | Rot3, Pose3, Unit3 | Python → generated C++ (and CUDA via Caspar, 2026) | quaternion | generated, exact | ✅ | ⚠️ (flattened, compiler-vectorizable) | ✅ generated | **ε-regularized, branch-free** | Apache-2.0 |
| **wrenfold** | quaternion / rotation types as symbolic building blocks | Python → generated C++17 / Rust | user-defined | generated | ✅ | ⚠️ (flattened code) | ❌ (no CUDA target found) | depends on expressions | check repo; JOSS 2025 |
| **Pinocchio** (spatial) | SE3, SO3 | C++ (large library) | matrix | Jexp/Jlog | ⚠️ | ❌ | ❌ | careful precision-aware Taylor | BSD-2 |
| **GSLAM transforms** | SO3, SE3, Sim3 | C++ | quaternion | partial | ⚠️ | ❌ | ❌ | not assessed | check repo; reports scalar exp/log faster than Sophus |
| **lietorch** | SO3, SE3 (+Sim3) | PyTorch CUDA kernels | quaternion | via autograd (tangent-space backprop) | ✅ | n/a | ✅ batched | good | BSD-3 |
| **PyPose** (LieTensor) | SO3, SE3, Sim3 | PyTorch | quaternion | autograd, batched Jacobians | ✅ | n/a | ✅ batched | good | BSD-style (verify) |
| **Theseus / torchlie** (Meta) | SO3, SE3 | PyTorch | matrix | autograd | ✅ | n/a | ✅ batched | good | MIT |
| **jaxlie** | SO2/3, SE2/3 | JAX (Sophus-inspired) | quaternion (+t) | AD (forward/reverse) | ✅ | ✅ via `vmap`/XLA | ✅ via XLA | good | MIT |
| **liealg** (Rust) | SO3, SE3 | Rust | matrix | ❌ | ⚠️ | ❌ | ❌ | basic | check crate |
| **Klein** | rigid motions as PGA motors | C++ SSE | motor (8 floats) | ❌ | ✅ only | ✅ SSE-packed | ❌ | n/a | MIT; low activity |
| **Ceres `rotation.h`** | conversions only | C++ header | various | via Jet | ✅ | ❌ | ❌ | careful (Jet-friendly) | BSD |

### 4.1 Verdict on alternatives

- **No library combines all six goals.** Each is strong on one axis:
  manif/smooth on API completeness, GTSAM on careful numerics (but double-only
  and heavy), symforce/wrenfold on branch-free generated code, lietorch / PyPose /
  Theseus / jaxlie on GPU batches (but Python-side), sophus-rs on SIMD batch types (Rust).
- **The C++ gap is real.** A dedicated search found **no C++ library with SIMD
  batches of SO3/SE3 exp/log in SoA layout**. The batched implementations are all in
  Python frameworks (PyTorch/JAX) or Rust. Layout studies in other fields
  agree that choosing SoA matters more than hand-written intrinsics.
- **S2 is rare**: only GTSAM (Unit3), symforce (Unit3) and Sophus 2
  (unit vectors) offer it. That is a reason in itself to keep csLie.
- **Codegen as a tool, not a dependency**: symforce, wrenfold and Caspar show
  that generated, flattened, branch-free code is very fast. csLie can use
  SymPy/symforce/wrenfold *offline* to derive and verify formulas, Taylor
  coefficients and Jacobians, while shipping plain hand-checked headers.
- **Continuing csLie is justified** if it sets out to be the
  *"accurate in float32, batch-SIMD and GPU-ready"* library, which none of the
  above fully is. Borrow, don't reinvent:
  - API and Jacobian conventions from **manif** (Solà et al., "A micro Lie theory").
  - Numerically robust formulas and Unit3 basis from **GTSAM**.
  - Branch-free / ε-regularized formulation and codegen-based verification from **symforce** / **wrenfold**.
  - Batch-scalar design from **sophus-rs**; batched API ergonomics from **jaxlie** / **PyPose**.
- **Use them as benchmarks and oracles**: Sophus, manif, GTSAM and GSLAM in the
  CPU benchmark suite; lietorch / PyPose as GPU throughput references (§6).

Facts about third-party projects come from my knowledge plus web searches made
on 2026-11-15 (see Sources). Licenses and project activity marked "verify" or
"check repo" were not confirmed.

### 4.2 Sources

- manif: <https://github.com/artivis/manif>
- smooth: <https://pettni.github.io/smooth/>, <https://github.com/pettni/smooth>
- Lie++: <https://github.com/aau-cns/Lie-plusplus>
- Sophus: <https://github.com/strasdat/Sophus> (README: Sophus 1 maintenance mode, Sophus 2 in farm-ng-core)
- sophus-rs: <https://github.com/farm-ng/sophus-rs>, <https://lib.rs/crates/sophus_core/features>
- GTSAM SO3 numerics: <https://borglab.github.io/gtsam/so3>, DexpFunctor PR <https://github.com/borglab/gtsam/pull/2099>
- symforce: <https://github.com/symforce-org/symforce>, paper <https://arxiv.org/abs/2204.07889>
- Caspar (CUDA from SymForce): <https://arxiv.org/pdf/2605.30583>
- wrenfold: <https://wrenfold.org/>, <https://github.com/wrenfold/wrenfold>
- lietorch: <https://github.com/princeton-vl/lietorch>, paper <https://arxiv.org/pdf/2103.12032>
- PyPose: <https://pypi.org/project/pypose/>
- torchlie / Theseus: <https://pypi.org/project/torchlie>
- jaxlie: <https://brentyi.github.io/jaxlie>
- liealg: <https://docs.rs/crate/liealg/latest>
- GSLAM transforms: <https://zdzhaoyong.github.io/GSLAM/transform.html>
- Rotation representations survey: <https://arxiv.org/html/2605.08086v1>
- SoA vs AoS vs intrinsics study: <https://arxiv.org/pdf/2603.18052>

---

## 5. Target architecture

```
┌ csLie (header-only, public) ───────────────────────────────────────────────┐
│ types:   SO3<S>, SE3<S>, S2<S>, SO2<S>, SE2<S>, tangent types              │
│          exp/log/compose/inverse/act/plus/minus/adj/Jr/Jl/… + opt. Jacobians│
│ compat:  current free functions on Eigen matrices (thin wrappers)           │
│ checks:  Paranoid preconditions via csCore asserts (host only)              │
└──────────────────────────────┬─────────────────────────────────────────────┘
┌ csLie kernels (header-only, dependency-free, CSLIE_HD) ─────────────────────┐
│ scalar-generic, branch-light formulas; per-precision constexpr thresholds   │
│ math traits: sqrt/rsqrt/sincos/atan2/select/fma for float, double, Jet,     │
│              SIMD batch, CUDA                                               │
└──────────┬───────────────────────────┬─────────────────────────────────────┘
┌ csLie::batch (optional) ┐   ┌ csLie::cuda (optional) ┐
│ span/SoA API, xsimd     │   │ device-array kernels    │
│ (Highway if dispatch)   │   │ same kernels, SoA       │
└─────────────────────────┘   └─────────────────────────┘
```

---

## 6. Proof: precision and performance suites

"Best available" needs evidence, not claims.

**Precision suite** (extends the probe of §2.3):

- Reference in long double / `__float128` / MPFR (boost::multiprecision).
- Grids: θ logarithmically from 1e-12 to π, **π − 10⁻ᵏ**, **exactly 0 and π**,
  random axes, axis-aligned axes, subnormals, inputs with injected
  non-orthonormality (1e-7, 1e-4).
- Metrics: max and p99 error in ε (or ulps) per function × precision;
  NaN/Inf must never appear for valid input.
- Jacobians: compared with central differences in extended precision, or with
  autodiff (`Jet`).
- Properties: `exp(log(X)) = X`, `log(exp(ξ)) = ξ` (for |θ| < π), `J·J⁻¹ = I`,
  adjoint identities, group axioms.
- Tolerances checked in CI, so regressions fail the build.

**Performance suite** (Google Benchmark or nanobench):

- ns/op for single ops (exp, log, compose, act, exp + Jacobian) and throughput
  for batches (1k, 1M elements), float32 and float64.
- Baselines: current csLie, Sophus, manif, GTSAM, Eigen `AngleAxis`/`Quaternion`,
  symforce-generated code.
- Flags documented: `-O3`, `-march=x86-64-v3` (AVX2) and `-march=native`, no `-ffast-math`.
- GPU: elements/s on representative cards for the CUDA component.

---

## 7. Risks and open questions

| Topic | Risk | Stance |
|-------|------|--------|
| API break | Moving from matrix free-functions to types breaks callers | Keep free functions as wrappers for a transition release |
| Branch-free formulas | Evaluating both branches costs scalar performance | Measure; keep a scalar path with real branches where it wins |
| ε-regularization | Bias near 0 | Only where provably below rounding; covered by the precision suite |
| SIMD dependency | xsimd/Highway adds a dependency | Only in the optional `batch` component; core stays dependency-free |
| GPU + Eigen | Compile time, warnings, missing features on device | Kernels on scalars; Eigen only at the API layer |
| Determinism | CPU vs GPU vs compilers differ | Accuracy specified in ε/ulps, never bitwise |

Open questions:

1. Is **Ceres `Jet`** (or another autodiff) compatibility required?
2. **Right or left** perturbation as default? (Right is the most common in VIO/SLAM.)
3. Quaternion storage: **Hamilton** convention and canonical `w ≥ 0`?
4. Target CPUs: x86-64 only (AVX2 baseline?) or also ARM (NEON/SVE, e.g. Jetson)?
   This decides xsimd vs Highway.
5. Minimum GPU: NVIDIA only (CUDA), or HIP as well?
6. Raise the C++ standard (C++20 for concepts / `std::span`) or stay on C++17?

---

## 8. Proposed roadmap

**Phase 0: Fix and measure.**
Fix the NaN at exactly π, the small-angle `se3Exp` cancellation, the
unqualified/qualified math calls (float promotion), and the `acos` usage. Add the
precision suite (§6) and a benchmark baseline against Sophus and manif.

**Phase 1: Kernel layer.**
Scalar-generic, branch-light kernels with math traits and per-precision
thresholds; `CSLIE_HD`; no csCore dependency inside kernels. Existing free functions
re-implemented on top, with no API change.

**Phase 2: Types and complete API for SO3, SE3, S2, SO2, SE2.**
Quaternion-based `SO3`/`SE3`, `S2` with log/⊖ and full Jacobians; complex-based
`SO2`/`SE2` with `se2Log` completed and full Jacobians; left and right
Jacobians; compose/inverse/act/plus/minus with optional Jacobians; slerp;
normalization; SE2 ↔ SE3 helpers. Optional Ceres `Jet` support.

**Phase 3: CPU batch SIMD.**
`csLie::batch` with span/SoA API on xsimd; measure vs scalar and alternatives;
Highway dispatch if multi-target binaries are needed.

**Phase 4: GPU.**
Device-compilable headers verified under NVCC; `csLie::cuda` batch kernels;
float32 accuracy validated on device.

**Phase 5: Hardening.**
Extended-precision references in CI, fuzzing of edge inputs, documentation of
conventions and accuracy per function.

**Optional, any time after phase 2: Python bindings** (§3.8), nanobind with a
batch NumPy API; DLPack after phase 4.

---

## 9. Summary

- csLie is a **small, tested, header-only** seed with a rare **S2** implementation,
  but it is incomplete for SO3/SE3/S2 (and SE2 lacks `log`) and has **measurable numerical defects**:
  NaN at exactly π, ~1e-11 small-angle error in float64 `se3Exp`, silent
  double promotion in float32, and `acos`-based logs.
- **Accuracy first**: cancellation-free formulas, `atan2`-based angles,
  Shepperd-style near-π handling, per-precision thresholds, proven by a precision suite.
- **Speed**: quaternion storage, fewer transcendentals, shared subexpressions.
  The real multiplier is **vertical SIMD on SoA batches**, written once as
  scalar-generic, branch-light kernels.
- **The same kernels serve GPU** (`__host__ __device__`), where float32 accuracy becomes essential.
- **Stay header-only** for the core. Only optional batch-dispatch and CUDA
  components are compiled.
- **No existing library meets all goals**, and no C++ library offers SIMD SoA
  batches of SO3/SE3 exp/log. Borrow API from manif, numerics from GTSAM,
  branch-free and codegen ideas from symforce/wrenfold, batch types from
  sophus-rs, and benchmark against them.
- **Scope**: SO3, SE3, S2 for 3D and SO2, SE2 for 2D, all first-class.
- **Stays a separate repository** on csCore's base layer. **Python bindings**
  are an optional nanobind layer once the API is stable.
