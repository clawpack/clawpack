# AGENTS.md

Guidance for AI coding tools working in this repository. Humans: see
[AI_POLICY.md](AI_POLICY.md). This file relaxes nothing there; it makes the
same rules legible to tooling.

- A human owns and must understand every change. Do not open a PR the human
  has not read and cannot explain.
- Never add `Signed-off-by:` or list yourself as an author or co-author. When
  you meaningfully assist a commit, add `Assisted-by: <tool> <model>`; the
  human adds the sign-off.
- Do not write the human summary at the top of a PR description, or draft
  review replies for the human to paste back. You may draft a detailed
  description below the summary, marked as AI-generated.
- Numerical core (Riemann solvers, fluxes, limiters, source terms,
  well-balancing, time stepping, AMR flagging and interpolation): you may
  implement a decision the human has made, write tests, build verification
  drivers, or refactor with the math held fixed. STOP and flag before making
  a numerical decision. Methods come from cited literature, never from
  recalled code.
- Results-affecting changes need verification evidence (manufactured/exact
  solution, convergence rate, conservation/well-balancing checks). A run that
  looks right is not evidence.
- Clawpack is BSD-licensed. Do not introduce code that may derive from
  copyleft or license-incompatible sources.
- Keep diffs small and focused.
