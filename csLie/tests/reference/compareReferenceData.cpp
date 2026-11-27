// Compares a committed reference data file with the generator's output: the same text, except that two hexadecimal
// floating literals may differ by one unit in the last place. The generator's quad-precision results are exact to far
// below double precision, but a value that is the cancellation of larger terms can round to either of two doubles
// depending on libquadmath's last bit, which differs between GCC versions.
//   compareReferenceData <committed file> <generated file>
// Exit 0 when they match, 1 when they do not (the first difference on stderr), 2 on a usage or read error.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>

namespace {

bool read(char const *path, std::string &text) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    return false;
  }
  std::ostringstream s;
  s << in.rdbuf();
  text = s.str();
  return true;
}

// the length of the hexadecimal floating literal at text[i] (an optional '-', "0x", digits and '.', 'p', exponent),
// 0 when there is none
std::size_t hexFloatLength(std::string const &text, std::size_t const i) {
  std::size_t j = i;
  if (j < text.size() && text[j] == '-') {
    ++j;
  }
  if (j + 1 >= text.size() || text[j] != '0' || text[j + 1] != 'x') {
    return 0;
  }
  j += 2;
  while (j < text.size() && (std::isxdigit(static_cast<unsigned char>(text[j])) != 0 || text[j] == '.')) {
    ++j;
  }
  if (j >= text.size() || text[j] != 'p') {
    return 0;
  }
  ++j;
  if (j < text.size() && (text[j] == '+' || text[j] == '-')) {
    ++j;
  }
  while (j < text.size() && std::isdigit(static_cast<unsigned char>(text[j])) != 0) {
    ++j;
  }
  return j - i;
}

std::size_t lineOf(std::string const &text, std::size_t const i) {
  std::size_t line = 1;
  for (std::size_t k = 0; k < i && k < text.size(); ++k) {
    line += text[k] == '\n' ? std::size_t{1} : std::size_t{0};
  }
  return line;
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 3) {
    std::fprintf(stderr, "usage: compareReferenceData <committed file> <generated file>\n");
    return 2;
  }
  std::string a;
  std::string b;
  if (!read(argv[1], a) || !read(argv[2], b)) {
    std::fprintf(stderr, "compareReferenceData: cannot read %s or %s\n", argv[1], argv[2]);
    return 2;
  }
  std::size_t i = 0;
  std::size_t j = 0;
  while (i < a.size() && j < b.size()) {
    std::size_t const la = hexFloatLength(a, i);
    std::size_t const lb = hexFloatLength(b, j);
    if (la > 0 && lb > 0) {
      double const x = std::strtod(a.substr(i, la).c_str(), nullptr);
      double const y = std::strtod(b.substr(j, lb).c_str(), nullptr);
      // equal (+0 and -0 too), or neighbours
      if (!(x <= y && x >= y) && !(std::nextafter(x, y) <= y && std::nextafter(x, y) >= y)) {
        std::fprintf(stderr, "line %zu: %s against %s\n", lineOf(a, i), a.substr(i, la).c_str(),
                     b.substr(j, lb).c_str());
        return 1;
      }
      i += la;
      j += lb;
    } else if (a[i] == b[j]) {
      ++i;
      ++j;
    } else {
      std::fprintf(stderr, "line %zu: the text differs\n", lineOf(a, i));
      return 1;
    }
  }
  if (i != a.size() || j != b.size()) {
    std::fprintf(stderr, "line %zu: one file is longer\n", lineOf(a, i));
    return 1;
  }
  return 0;
}
