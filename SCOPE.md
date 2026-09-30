# What this repository covers, and what it does not

A short note for the corresponding author and for reviewers, so that the
boundary of this deposit is explicit rather than inferred.

## The generalized estimating equation models are not here

The paper reports generalized estimating equation models. **No GEE is fitted
anywhere in this repository, and none was fitted by the codebase these scripts
were extracted from.** This was checked rather than assumed: the full set of
statistical routines invoked across every analysis contributing to Figure 5 and
Extended Data Figure 7 is ordinary least squares regression, Spearman
correlation, Benjamini-Hochberg and Holm multiplicity correction, Wilcoxon
signed-rank, Friedman, exact binomial, Fisher's exact test, logistic
regression, a permutation test, and recovery of a z statistic from a p-value.

The GEE models therefore belong to a different analysis, most plausibly a
repeated-measures analysis of the longitudinal behavioural or clinical data.
Whoever fitted them holds that code and should deposit it, either alongside
this repository or as a second deposit referenced from the data availability
statement. It cannot be reconstructed here, and reconstructing a model of that
kind without its original data and specification would produce something that
looks like the published analysis without being it.

## Also not here

- The mouse experiments: surgery, behaviour, tissue collection, telomere
  measurement, FACS, qPCR, immunohistochemistry and its quantification.
- The proximity-labelling proteomics and the cell-line work.
- The Extended Data Fig. 7F osteoporosis replication meta-analysis, which uses
  summary statistics taken directly from the GWAS Catalog and from BioBank
  Japan rather than through code in this repository. The fixed-effect
  inverse-variance meta-analysis and Cochran's Q are three lines of arithmetic
  on published betas and standard errors, reported in the figure legend.
- The dorsal root ganglion analysis of Figure 5E and the longitudinal whole
  blood analysis of Figure 5F. These use published cohort data and are not
  included, because in both cases the analysis depends on participant-level
  phenotype assignments that are not fully recoverable from the public
  deposits. This is a limitation of those two panels rather than an omission
  from this repository. The second dorsal root ganglion cohort's pain
  phenotype is not in the public record, and the deposited Figure 5E metadata
  carries no donor identifier, so repeat samples from the same donor cannot be
  collapsed from public data alone.

## What the assertions are for

Each script ends by checking its own results against the values in the
published figures, and fails loudly if they do not match. This is deliberate.
Three of the pinned resources will drift: the CZ CELLxGENE Census `"stable"`
alias advances, GTEx gene identifiers are release-specific, and FinnGen
endpoint definitions change between releases. A future re-run that quietly
produces different numbers is worse than one that stops, so the scripts stop.

If an assertion fails after a release change, that is information rather than a
bug: the pinned release in `src/common.py` no longer matches what the service
returns, and the fix is to pin the earlier release or to report the new
numbers, not to relax the check.
