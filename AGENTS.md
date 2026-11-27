# AGENTS.md

<!-- cs-template:begin (this block is owned by cs-project-template: do not edit; the rest of the file is the project's) -->
Entry point for coding agents (and people) working in csLie. The rules of every cscosine repository are
in [.claude/cscosine-rules.md](.claude/cscosine-rules.md): read them first. What is specific to this repository
follows this block.

Commands:

- Checks: `pre-commit run --all-files` (`pre-commit install` once).

Who owns which file (template: [cs-project-template](https://github.com/cscosine/cs-project-template); change it
there, then `copier update --defaults`; never edit its parts here):

- Template: `open-code.sh`, `open-code.ps1`, `.clang-format`,
  `.vscode/`, `.github/` (the workflows: `checks.yml` only), `CLAUDE.md`, `.claude/cscosine-rules.md`, `.claude/settings.json`,
  `LICENSE`.
- Template and project: `.pre-commit-config.yaml` (project hooks
  between the `project: begin` / `project: end` markers), `.gitignore`, `.gitattributes` (project lines after
  `--- project ---`), `README.md` and this file (the
  `cs-template` block).
- The project: everything else.
<!-- cs-template:end -->

## csLie

Lie groups and algebras (SE(3), SE(2)) in C++ on Eigen

Rules of this library:

- Every target builds with `cscmake_strict_warnings` (`DEFAULT_COMPILE_PROFILES` in `CMakeLists.txt`): GCC and Clang
  both clean.
- Every public function has a test over `float32_t, float64_t` (`TEMPLATE_TEST_CASE`): a header is a template, its body
  is compiled (and warned about) only where a test instantiates it.
- Generic scalar code: `using std::sin; using std::cos; using std::sqrt; ...` in the function and unqualified calls (a
  float stays float, a scalar type with its own overloads finds them); literals as `S(...)`, never a bare double.
