// Quad-precision (__float128, ~1e-34) reference of the csLie functions, used by generateReferenceData.cpp only
// (__float128 and libquadmath: GCC on x86-64 Linux): the tests read the values it generates (data/), so
// they run on every platform. Every formula is accurate to far below double precision for any angle: Taylor series (to
// theta^8) below 1e-3, closed forms above; no formula subtracts nearly equal quad values by more than 1e-6.
#pragma once
#include <quadmath.h>

#include <array>
#include <cmath>
#include <cstddef>

namespace ref {

using Q = __float128;

template <int R, int C> struct M {
  std::array<Q, static_cast<std::size_t>(R *C)> a{};
  Q &operator()(int r, int c) { return a[static_cast<std::size_t>(r * C + c)]; }
  Q operator()(int r, int c) const { return a[static_cast<std::size_t>(r * C + c)]; }
};
using V3 = M<3, 1>;
using M3 = M<3, 3>;
using V2 = M<2, 1>;
using M2 = M<2, 2>;
using V6 = M<6, 1>;

// pi from atanq, not M_PIq: its Q literal suffix is a GNU extension (rejected by -std=c++17)
inline const Q PI = 4 * atanq(1);

template <int R, int K, int C> M<R, C> mul(M<R, K> const &a, M<K, C> const &b) {
  M<R, C> o;
  for (int r = 0; r < R; ++r)
    for (int c = 0; c < C; ++c) {
      Q s = 0;
      for (int k = 0; k < K; ++k)
        s += a(r, k) * b(k, c);
      o(r, c) = s;
    }
  return o;
}
template <int R, int C> M<R, C> add(M<R, C> a, M<R, C> const &b, Q sb = 1) {
  for (std::size_t i = 0; i < a.a.size(); ++i)
    a.a[i] += sb * b.a[i];
  return a;
}
template <int R, int C> M<C, R> tr(M<R, C> const &a) {
  M<C, R> o;
  for (int r = 0; r < R; ++r)
    for (int c = 0; c < C; ++c)
      o(c, r) = a(r, c);
  return o;
}
template <int R, int C> Q norm(M<R, C> const &a) {
  Q s = 0;
  for (auto v : a.a)
    s += v * v;
  return sqrtq(s);
}
inline M3 eye3() {
  M3 o;
  o(0, 0) = o(1, 1) = o(2, 2) = 1;
  return o;
}
inline M3 hat(V3 const &w) {
  M3 o;
  o(0, 1) = -w(2, 0);
  o(0, 2) = w(1, 0);
  o(1, 0) = w(2, 0);
  o(1, 2) = -w(0, 0);
  o(2, 0) = -w(1, 0);
  o(2, 1) = w(0, 0);
  return o;
}
inline V3 cross(V3 const &a, V3 const &b) { return mul(hat(a), b); }
inline Q dot3(V3 const &a, V3 const &b) { return a(0, 0) * b(0, 0) + a(1, 0) * b(1, 0) + a(2, 0) * b(2, 0); }

// A = sin t / t, B = (1 - cos t) / t^2, C = (t - sin t) / t^3, D = (1 - A / (2B)) / t^2
struct Coeffs {
  Q A, B, C, D;
};
inline Coeffs coeffs(Q t) {
  Q const t2 = t * t;
  Coeffs c;
  if (t < Q(1e-3)) {
    c.A = 1 - t2 / 6 * (1 - t2 / 20 * (1 - t2 / 42 * (1 - t2 / 72)));
    c.B = Q(0.5) * (1 - t2 / 12 * (1 - t2 / 30 * (1 - t2 / 56 * (1 - t2 / 90))));
    c.C = (1 - t2 / 20 * (1 - t2 / 42 * (1 - t2 / 72 * (1 - t2 / 110)))) / 6;
    c.D = Q(1) / 12 + t2 * (Q(1) / 720 + t2 * (Q(1) / 30240 + t2 * (Q(1) / 1209600 + t2 / Q(47900160))));
  } else {
    Q const s = sinq(t), sh = sinq(t / 2);
    c.A = s / t;
    c.B = 2 * sh * sh / t2;
    c.C = (t - s) / (t2 * t);
    c.D = (1 - c.A / (2 * c.B)) / t2;
  }
  return c;
}

// ---- SO3 ----
inline M3 so3Exp(V3 const &w) {
  Coeffs const c = coeffs(norm(w));
  M3 const W = hat(w);
  return add(add(eye3(), W, c.A), mul(W, W), c.B);
}
inline M3 so3Jl(V3 const &w) {
  Coeffs const c = coeffs(norm(w));
  M3 const W = hat(w);
  return add(add(eye3(), W, c.B), mul(W, W), c.C);
}
inline M3 so3JlInv(V3 const &w) {
  Coeffs const c = coeffs(norm(w));
  M3 const W = hat(w);
  return add(add(eye3(), W, Q(-0.5)), mul(W, W), c.D);
}
// Shepperd's method to an unnormalized quaternion, then theta = 2 atan2(|v|, w).
// Valid for any (also slightly non-orthonormal) input; at exactly pi the sign is the pivot's.
inline V3 so3Log(M3 const &R) {
  Q const trc = R(0, 0) + R(1, 1) + R(2, 2);
  Q w;
  V3 v;
  if (trc > 0) {
    w = 1 + trc;
    v(0, 0) = R(2, 1) - R(1, 2);
    v(1, 0) = R(0, 2) - R(2, 0);
    v(2, 0) = R(1, 0) - R(0, 1);
  } else {
    int i = 0;
    if (R(1, 1) > R(i, i))
      i = 1;
    if (R(2, 2) > R(i, i))
      i = 2;
    int const j = (i + 1) % 3, k = (i + 2) % 3;
    Q const r = sqrtq(1 + 2 * R(i, i) - trc);
    v(i, 0) = r;
    v(j, 0) = (R(j, i) + R(i, j)) / r;
    v(k, 0) = (R(k, i) + R(i, k)) / r;
    w = (R(k, j) - R(j, k)) / r;
    if (w < 0) {
      w = -w;
      for (auto &x : v.a)
        x = -x;
    }
  }
  Q const n = norm(v);
  V3 o;
  if (!(n > 0))
    return o;
  Q const th = 2 * atan2q(n, w);
  for (int i = 0; i < 3; ++i)
    o(i, 0) = th / n * v(i, 0);
  return o;
}

// ---- SE3: tangent [rho; omega], T = [R, Jl(omega) rho] ----
struct SE3q {
  M3 R;
  V3 t;
};
inline SE3q se3Exp(V6 const &xi) {
  V3 rho, w;
  for (int i = 0; i < 3; ++i)
    rho(i, 0) = xi(i, 0), w(i, 0) = xi(i + 3, 0);
  return {so3Exp(w), mul(so3Jl(w), rho)};
}
inline V6 se3Log(SE3q const &T) {
  V3 const w = so3Log(T.R);
  V3 const rho = mul(so3JlInv(w), T.t);
  V6 o;
  for (int i = 0; i < 3; ++i)
    o(i, 0) = rho(i, 0), o(i + 3, 0) = w(i, 0);
  return o;
}

// ---- SO2 / SE2: tangent [x, y, theta] ----
inline M2 so2Exp(Q t) {
  M2 o;
  o(0, 0) = o(1, 1) = cosq(t);
  o(1, 0) = sinq(t);
  o(0, 1) = -o(1, 0);
  return o;
}
inline Q so2Log(M2 const &R) { return atan2q(R(1, 0), R(0, 0)); }
struct SE2q {
  M2 R;
  V2 t;
};
// V = [[A, -B t], [B t, A]]
inline SE2q se2Exp(V3 const &xi) {
  Q const t = xi(2, 0);
  Coeffs const c = coeffs(fabsq(t));
  Q const a = c.A, b = c.B * t;
  V2 tr2;
  tr2(0, 0) = a * xi(0, 0) - b * xi(1, 0);
  tr2(1, 0) = b * xi(0, 0) + a * xi(1, 0);
  return {so2Exp(t), tr2};
}
inline V3 se2Log(SE2q const &T) {
  Q const t = so2Log(T.R);
  Coeffs const c = coeffs(fabsq(t));
  Q const a = c.A, b = c.B * t, det = a * a + b * b;
  V3 o;
  o(0, 0) = (a * T.t(0, 0) + b * T.t(1, 0)) / det;
  o(1, 0) = (-b * T.t(0, 0) + a * T.t(1, 0)) / det;
  o(2, 0) = t;
  return o;
}

// ---- S2: geodesic exp / log at p (|p| = 1), tangent u in R^3 with u.p = 0 ----
inline V3 s2Exp(V3 const &p, V3 const &u) {
  Q const n = norm(u);
  Coeffs const c = coeffs(n);
  V3 o = p;
  for (auto &x : o.a)
    x *= cosq(n);
  return add(o, u, c.A);
}
inline V3 s2Log(V3 const &p, V3 const &q) {
  V3 const pq = cross(p, q);
  Q const s = norm(pq), c = dot3(p, q);
  Q const th = atan2q(s, c);
  // direction: q - (p.q) p, normalized; its norm equals s for unit p, q
  V3 d = add(q, p, -c);
  Q const dn = norm(d);
  V3 o;
  if (!(dn > 0))
    return o;
  for (int i = 0; i < 3; ++i)
    o(i, 0) = th / dn * d(i, 0);
  return o;
}

} // namespace ref
