# Clawpack AI-Assisted Contribution Policy

> **TL;DR:** We are responsible for any contributions and interactions, especially those that touch core numerical code.
> - AI-assisted contributions are welcome, and you are fully responsible for their content.
> - You are responsible for the numerical decisions and verification that they are correct.
> - You disclose AI help with `git commit -s` plus an `Assisted-by:` trailer, and your PR opens with your own summary, followed by AI content clearly marked.

> _Drafted with AI assistance (Claude), directed and edited by Kyle Mandli; under review by the Clawpack maintainers._

Clawpack welcomes AI-assisted contributions. Using an AI tool never absolves a contributor of their responsibility for their contributions, however they are produced. This applies to code, tests, docs, examples, and build/CI configuration in any Clawpack repository.

## The four commitments

**1. Vet: you own it and can explain it.** Submit only what you understand and could defend yourself. A human who understands the change is always in the loop; we don't accept unattended agent PRs or "please review what the agent wrote." Keep PRs small and focused, and for anything non-trivial start from an issue that says what problem is being solved. Answer review questions in your own words, not by relaying them to a model.

**2. You do the numerics.** You make every numerical decision: any mathematical choice that changes what the code computes (scheme, derivation, formulation, flux, limiter, source-term or well-balancing treatment, AMR flagging and interpolation). AI can help you implement a decision you've made, learn a method, check your understanding, write tests, or refactor with the math held fixed, but it doesn't make the decision. Going from Cartesian to spherical brings new metric and source terms, so it's yours to derive; an affine relabeling isn't. Results-affecting changes ship with verification evidence: a manufactured or exact solution, a convergence rate at design order, or conservation/well-balancing/positivity checks. A plausible-looking plot is not verification.

**3. Attribute: be transparent.** When AI meaningfully helped with a commit, disclose it with an `Assisted-by:` trailer and sign off (`-s`) to certify that you own that help:

```
git commit -s --trailer "Assisted-by: Claude Code (claude-opus-4-8)"
```

Autocomplete in otherwise hand-written code doesn't need disclosing; generated or agentic work does. Tools never sign off or appear as co-authors. Your PR description opens with a short summary you wrote, in Clawpack's terminology, of what problem is solved and how; any longer AI-written detail goes below it, clearly marked. Disclosure informs review and is never held against you.

**4. Protect: license and IP.** Clawpack is BSD. Methods come from cited papers, not code an AI reproduced from memory. Don't submit output that may derive from copyleft or license-incompatible sources, check that your tool's terms don't conflict with BSD, and never put proprietary, embargoed, or NDA'd material into a prompt.

## For reviewers

AI use doesn't by itself change how a PR is reviewed. Scrutiny scales, as it always has, with how well we know the contributor, and the sign-off puts responsibility on the submitter. Ask for verification evidence on results-affecting changes before merging.

---

Sign-off is currently required only on AI-assisted commits and is on the honor system; a project-wide [DCO](https://developercertificate.org/) or a check bot can be revisited later. This policy will be revisited as tools and norms change.
