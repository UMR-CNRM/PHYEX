<!--
The coding norms ('norms' workflow, or phyex-check_norms locally) and the compilation and
execution of the tests ('compil/exec' jobs of the 'tests' workflow, run for the pull requests
on master, or phyex-check_testprogs/phyex-check_phyex1d locally) must pass: a pull request
cannot be integrated otherwise. The 'reproduction' jobs fail when the results differ from the
references, which is expected when the results change (see 'Impact on the results').

Title: it must start with 'bugfix: ', 'feature: ' or 'bugfix+feature: ', in any case (this is
checked by the 'PR title' workflow), for instance 'bugfix: wrong ice fall speed'. Separate pull
requests are preferred when a bug fix and a new feature can be split.
-->

## Description

<!-- What is modified and why (bug, new feature...), the related issues and pull requests (IAL, Meso-NH...) -->

## Impact on the results

<!-- Tick the box which applies -->
- [ ] The results are unchanged
- [ ] The results change, as expected because: <!-- explain which configurations/schemes are affected and why -->
