// Generates the reference data of the csLie accuracy tests, one header per group of functions (data/<group>.hpp, the
// types in ReferenceTypes.hpp): inputs and their exact results, computed in quad precision (QuadReference.hpp).
//   generateReferenceData <output folder>
// Built where __float128 and libquadmath exist (GCC on x86-64 Linux), not installed. The test "csLie reference data
// up to date" fails when its output differs from the committed data/ folder; to regenerate it after changing this
// file or QuadReference.hpp:
//   cmake --build <build dir> --config <config> --target csLieRegenerateReferenceData
// Each input is rounded to the tested precision first; its reference is computed from that rounded value, so a
// test measures the error of the function alone. Values are written as hexadecimal floating literals (exact).
#include "QuadReference.hpp"

#include <cstddef>
#include <cstdio>
#include <string>
#include <vector>

namespace {

using ref::Q;

// the file being written
std::FILE *out = nullptr;

struct Angle {
  char const *label;
  Q theta;
};

std::vector<Angle> const &angles() {
  static std::vector<Angle> const a{
      {"0", 0},
      {"1e-12", Q(1e-12)},
      {"1e-10", Q(1e-10)},
      {"1e-8", Q(1e-8)},
      {"1e-6", Q(1e-6)},
      {"1e-4", Q(1e-4)},
      {"1e-3", Q(1e-3)},
      {"1deg", ref::PI / 180},
      {"0.1", Q(0.1)},
      {"0.5", Q(0.5)},
      {"1", Q(1)},
      {"2", Q(2)},
      {"2pi/3", 2 * ref::PI / 3},
      {"3", Q(3)},
      {"pi-1e-2", ref::PI - Q(1e-2)},
      {"pi-1e-4", ref::PI - Q(1e-4)},
      {"pi-1e-6", ref::PI - Q(1e-6)},
      {"pi-1e-8", ref::PI - Q(1e-8)},
  };
  return a;
}

struct Axis {
  char const *label;
  ref::V3 v;
};

std::vector<Axis> const &axes() {
  static std::vector<Axis> const a = [] {
    std::vector<Axis> o;
    auto add = [&](char const *label, Q x, Q y, Q z) {
      ref::V3 v;
      v(0, 0) = x, v(1, 0) = y, v(2, 0) = z;
      Q const n = ref::norm(v);
      for (auto &e : v.a)
        e /= n;
      o.push_back({label, v});
    };
    add("x", 1, 0, 0);
    add("y", 0, 1, 0);
    add("z", 0, 0, 1);
    add("-z", 0, 0, -1);
    add("xyz", 1, 1, 1);
    add("x-y", 1, -1, 0);
    add("r1", Q(0.3), Q(-0.8), Q(0.52));
    add("r2", Q(-0.91), Q(0.12), Q(-0.4));
    add("r3", Q(0.05), Q(0.6), Q(-0.79));
    return o;
  }();
  return a;
}

// the angles of the exp tests: those of the log, exactly pi and beyond pi (the exp takes any angle)
std::vector<Angle> const &expAngles() {
  static std::vector<Angle> const a = [] {
    std::vector<Angle> o = angles();
    o.push_back({"pi", ref::PI});
    o.push_back({"4", Q(4)});
    o.push_back({"2pi-1e-3", 2 * ref::PI - Q(1e-3)});
    o.push_back({"10", Q(10)});
    return o;
  }();
  return a;
}

// the angles of the coefficient tests: those of the exp, plus tiny ones and the decades between
std::vector<Angle> const &coeffAngles() {
  static std::vector<Angle> const a = [] {
    std::vector<Angle> o = expAngles();
    o.push_back({"1e-30", Q(1e-30)});
    o.push_back({"1e-20", Q(1e-20)});
    o.push_back({"1e-5", Q(1e-5)});
    o.push_back({"3e-4", Q(3e-4)});
    o.push_back({"1e-2", Q(1e-2)});
    o.push_back({"0.3", Q(0.3)});
    return o;
  }();
  return a;
}

// the translations of the SE2 tests
std::vector<std::pair<char const *, std::array<Q, 2>>> const &translations2() {
  static std::vector<std::pair<char const *, std::array<Q, 2>>> const t{
      {"x", {1, 0}}, {"y", {0, 1}}, {"r1", {Q(0.3), Q(-0.7)}}, {"r2", {Q(-2.5), Q(1.25)}}, {"far", {Q(100), Q(-40)}},
  };
  return t;
}

// rotations by exactly pi, exactly representable
std::vector<std::pair<char const *, std::array<int, 9>>> const &exactPi() {
  static std::vector<std::pair<char const *, std::array<int, 9>>> const m{
      {"pi x", {1, 0, 0, 0, -1, 0, 0, 0, -1}},  {"pi y", {-1, 0, 0, 0, 1, 0, 0, 0, -1}},
      {"pi z", {-1, 0, 0, 0, -1, 0, 0, 0, 1}},  {"pi x+y", {0, 1, 0, 1, 0, 0, 0, 0, -1}},
      {"pi x+z", {0, 0, 1, 0, -1, 0, 1, 0, 0}},
  };
  return m;
}

template <class T> Q roundTo(Q v) { return Q(static_cast<T>(v)); }

std::string hex(Q v) {
  char b[64];
  std::snprintf(b, sizeof b, "%a", static_cast<double>(v));
  return b;
}

template <int R, int C> std::string list(ref::M<R, C> const &m) {
  std::string s;
  for (std::size_t i = 0; i < m.a.size(); ++i)
    s += (i ? ", " : "") + hex(m.a[i]);
  return s;
}

template <class T> void so3LogCases(char const *name) {
  std::fprintf(out, "inline constexpr So3LogCase %s[] = {\n", name);
  auto emit = [](std::string const &label, ref::M3 const &R) {
    std::fprintf(out, "    {\"%s\", {%s}, {%s}},\n", label.c_str(), list(R).c_str(), list(ref::so3Log(R)).c_str());
  };
  for (auto const &an : angles())
    for (auto const &ax : axes()) {
      ref::V3 w;
      for (int i = 0; i < 3; ++i)
        w(i, 0) = an.theta * ax.v(i, 0);
      ref::M3 R = ref::so3Exp(w);
      for (auto &e : R.a)
        e = roundTo<T>(e);
      emit(std::string(an.label) + " " + ax.label, R);
      if (!(an.theta > 0))
        break; // one identity is enough
    }
  for (auto const &[label, m] : exactPi()) {
    ref::M3 R;
    for (std::size_t i = 0; i < m.size(); ++i)
      R.a[i] = m[i];
    emit(label, R);
  }
  std::fprintf(out, "};\n\n");
}

// a 3x3 function of the rotation vector (fn: ref::so3Exp, ref::so3Jl, ref::so3JlInv) at the angles of the exp
template <class T, class Fn> void so3MatrixCases(char const *name, Fn fn) {
  std::fprintf(out, "inline constexpr So3MatrixCase %s[] = {\n", name);
  for (auto const &an : expAngles())
    for (auto const &ax : axes()) {
      ref::V3 w;
      for (int i = 0; i < 3; ++i)
        w(i, 0) = roundTo<T>(an.theta * ax.v(i, 0));
      std::fprintf(out, "    {\"%s %s\", {%s}, {%s}},\n", an.label, ax.label, list(w).c_str(), list(fn(w)).c_str());
      if (!(an.theta > 0))
        break; // one zero is enough
    }
  std::fprintf(out, "};\n\n");
}

template <class T> void coeffCases(char const *name) {
  std::fprintf(out, "inline constexpr CoeffCase %s[] = {\n", name);
  for (auto const &an : coeffAngles()) {
    Q const t = roundTo<T>(an.theta);
    ref::Coeffs const c = ref::coeffs(t);
    std::fprintf(out, "    {\"%s\", %s, %s, %s, %s, %s, %s},\n", an.label, hex(t).c_str(), hex(c.A).c_str(),
                 hex(c.B).c_str(), hex(c.C).c_str(), hex(c.D).c_str(), hex(c.A / (2 * c.B)).c_str());
  }
  std::fprintf(out, "};\n\n");
}

// SE2 exp: every angle of the exp, both signs, every translation
template <class T> void se2ExpCases(char const *name) {
  std::fprintf(out, "inline constexpr Se2ExpCase %s[] = {\n", name);
  for (auto const &an : expAngles())
    for (Q const sign : {Q(1), Q(-1)}) {
      for (auto const &[tl, tv] : translations2()) {
        ref::V3 xi;
        xi(0, 0) = roundTo<T>(tv[0]), xi(1, 0) = roundTo<T>(tv[1]), xi(2, 0) = roundTo<T>(sign * an.theta);
        ref::SE2q const e = ref::se2Exp(xi);
        ref::M3 M;
        M(0, 0) = e.R(0, 0), M(0, 1) = e.R(0, 1), M(0, 2) = e.t(0, 0);
        M(1, 0) = e.R(1, 0), M(1, 1) = e.R(1, 1), M(1, 2) = e.t(1, 0);
        M(2, 2) = 1;
        std::fprintf(out, "    {\"%s%s %s\", {%s}, {%s}},\n", sign < 0 ? "-" : "", an.label, tl, list(xi).c_str(),
                     list(M).c_str());
      }
      if (!(an.theta > 0))
        break; // one zero is enough
    }
  std::fprintf(out, "};\n\n");
}

// SE2 log: the exp of the angles of the log (up to pi - 1e-8), both signs, every translation, rounded
template <class T> void se2LogCases(char const *name) {
  std::fprintf(out, "inline constexpr Se2LogCase %s[] = {\n", name);
  for (auto const &an : angles())
    for (Q const sign : {Q(1), Q(-1)}) {
      for (auto const &[tl, tv] : translations2()) {
        ref::V3 xi;
        xi(0, 0) = tv[0], xi(1, 0) = tv[1], xi(2, 0) = sign * an.theta;
        ref::SE2q e = ref::se2Exp(xi);
        for (auto &v : e.R.a)
          v = roundTo<T>(v);
        for (auto &v : e.t.a)
          v = roundTo<T>(v);
        ref::M3 M;
        M(0, 0) = e.R(0, 0), M(0, 1) = e.R(0, 1), M(0, 2) = e.t(0, 0);
        M(1, 0) = e.R(1, 0), M(1, 1) = e.R(1, 1), M(1, 2) = e.t(1, 0);
        M(2, 2) = 1;
        std::fprintf(out, "    {\"%s%s %s\", {%s}, {%s}},\n", sign < 0 ? "-" : "", an.label, tl, list(M).c_str(),
                     list(ref::se2Log(e)).c_str());
      }
      if (!(an.theta > 0))
        break; // one zero is enough
    }
  std::fprintf(out, "};\n\n");
}

// RPY: pitches from 0 to the gimbal lock (exactly +-pi / 2, rounded) and pairs of roll, yaw
template <class T> void rpyCases(char const *name) {
  std::fprintf(out, "inline constexpr RpyCase %s[] = {\n", name);
  std::vector<Angle> const pitches{
      {"0", 0},
      {"0.3", Q(0.3)},
      {"1.2", Q(1.2)},
      {"pi/2-1e-2", ref::PI / 2 - Q(1e-2)},
      {"pi/2-1e-4", ref::PI / 2 - Q(1e-4)},
      {"pi/2-1e-6", ref::PI / 2 - Q(1e-6)},
      {"pi/2-1e-8", ref::PI / 2 - Q(1e-8)},
      {"pi/2", ref::PI / 2},
  };
  std::vector<std::pair<char const *, std::array<Q, 2>>> const rollYaw{
      {"r0.2 y-0.7", {Q(0.2), Q(-0.7)}},
      {"r-2.9 y3", {Q(-2.9), Q(3)}},
      {"r1.1 y0", {Q(1.1), 0}},
      {"r0 y-1.9", {0, Q(-1.9)}},
  };
  for (auto const &p : pitches)
    for (Q const sign : {Q(1), Q(-1)}) {
      for (auto const &[label, ry] : rollYaw) {
        Q const roll = roundTo<T>(ry[0]), pitch = roundTo<T>(sign * p.theta), yaw = roundTo<T>(ry[1]);
        ref::V3 rpy;
        rpy(0, 0) = roll, rpy(1, 0) = pitch, rpy(2, 0) = yaw;
        std::fprintf(out, "    {\"p%s%s %s\", {%s}, {%s}},\n", sign < 0 ? "-" : "", p.label, label, list(rpy).c_str(),
                     list(ref::rollPitchYawToSo3(roll, pitch, yaw)).c_str());
      }
      if (!(p.theta > 0))
        break; // one zero is enough
    }
  std::fprintf(out, "};\n\n");
}

// writes <folder>/<group>.hpp: the header, the tables written by body, the footer
template <class Body> bool writeGroup(std::string const &folder, char const *group, Body body) {
  std::string const path = folder + "/" + group + ".hpp";
  out = std::fopen(path.c_str(), "w");
  if (out == nullptr) {
    std::fprintf(stderr, "generateReferenceData: cannot write %s\n", path.c_str());
    return false;
  }
  std::fprintf(out, "// Generated by generateReferenceData.cpp: do not edit, regenerate.\n"
                    "// clang-format off\n"
                    "#pragma once\n\n"
                    "#include \"../ReferenceTypes.hpp\"\n\n"
                    "namespace csLie::testReference {\n\n");
  body();
  std::fprintf(out, "} // namespace csLie::testReference\n"
                    "// clang-format on\n");
  bool const ok = std::fclose(out) == 0;
  out = nullptr;
  return ok;
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 2) {
    std::fprintf(stderr, "usage: generateReferenceData <output folder>\n");
    return 2;
  }
  std::string const folder = argv[1];
  bool const ok = writeGroup(folder, "SO3", [] {
    so3LogCases<float>("so3LogCasesF32");
    so3LogCases<double>("so3LogCasesF64");
    so3MatrixCases<float>("so3ExpCasesF32", ref::so3Exp);
    so3MatrixCases<double>("so3ExpCasesF64", ref::so3Exp);
    so3MatrixCases<float>("so3JlCasesF32", ref::so3Jl);
    so3MatrixCases<double>("so3JlCasesF64", ref::so3Jl);
    so3MatrixCases<float>("so3JlInvCasesF32", ref::so3JlInv);
    so3MatrixCases<double>("so3JlInvCasesF64", ref::so3JlInv);
    rpyCases<float>("rpyCasesF32");
    rpyCases<double>("rpyCasesF64");
    coeffCases<float>("coeffCasesF32");
    coeffCases<double>("coeffCasesF64");
  });
  bool const okSe2 = writeGroup(folder, "SE2", [] {
    se2ExpCases<float>("se2ExpCasesF32");
    se2ExpCases<double>("se2ExpCasesF64");
    se2LogCases<float>("se2LogCasesF32");
    se2LogCases<double>("se2LogCasesF64");
  });
  return ok && okSe2 ? 0 : 1;
}
