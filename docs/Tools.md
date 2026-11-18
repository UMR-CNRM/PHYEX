# PHYEX TOOLS DOCUMENTATION

## About this section

This document is intended for persons who want to use the tools provided with the PHYEX package: preparation of the source code for the different models and checks of the contributions (bit-reproducibility of the test cases between two commits, coding norms...).

## Installation

The tools are provided as a python package (pyphyextools, in the tools directory) which must be installed in editable mode (the tools use the files of the repository), for instance in a virtual environment:

```
python3 -m venv phyex.env
. phyex.env/bin/activate
pip install -e PHYEX/tools
```

The python dependencies (listed in tools/pyproject.toml) are installed automatically. In addition, some tools need:

  - phyex-check\_ial: the gmkpack build system must be installed and configured (the ial-to\_pack command, coming from the ial\_build python package, is used to create the packs)
  - all the phyex-check\_\* tools: FORTRAN and C compilers, CMake, git, wget, make and an internet connection (the model source codes, the test case data and some dependencies are downloaded)

The tools can be used from a repository (the version of the tools is the version of the repository) to check another commit or another local repository.

## User directories and environment variables

${PHYEXCONF} (the PHYEXCONF environment variable defaults to ${HOME}/.phyex) is the user directory containing:

  - testprogs\_data: input data of the test programs (see the [Offline documentation](./Offline.md))
  - arch: user architecture files (searched before those of the arch directory of the repository by phyex-check\_testprogs)
  - phyex1d: grids, namelists and DEPHY-SCM cases of the 1D model (see the [phyex1d documentation](./phyex1d.md))

The phyex-check\_\<model\> tools create their compilation directories in:

  - $HOMEPACK (default: $HOME/pack) for IAL; the PHYEXROOTPACK environment variable, if set, is used as root pack for the incremental packs
  - $MNHPACK (default: $HOME/MesoNH/PHYEX) for Meso-NH
  - $LMDZPACK (default: $HOME/LMDZ/PHYEX) for LMDZ
  - $TESTPROGSDIR (default: $HOME/TESTPROGS) for the test programs
  - $PHYEX1DDIR (default: $HOME/PHYEX1D) for the 1D model (phyex1d)

The PHYEXREPOuser (default: UMR-CNRM), PHYEXREPOprotocol (https or ssh) and PHYEXREPOrepo environment variables give the GitHub repository from which the commits are cloned.

## Tools

All the commands have a -h option describing their arguments.

### phyex-prep\_code

This command prepares the source code for a model: extraction of the source code (from a commit or from a local repository), transformations (with pyfortool), generation of the interface modules, and, optionally, push of the result as a new branch.
It is used to build the model specific branches (see the [Integrator documentation](./Integrator.md)) and by the phyex-check\_\<model\> tools.

### phyex-check\_ial, phyex-check\_mesonh, phyex-check\_lmdz and phyex-check\_testprogs

These commands compile a model (IAL, Meso-NH, LMDZ) or the test programs using the PHYEX source code of a given commit (or of a local repository), execute the test cases and compare the results against a reference simulation.

The steps are selected with options: -p (creation of the compilation directory), -c (compilation), -r (execution), -C (comparison), --remove (suppression of the compilation directory). The first argument is the commit (or the directory) to test, the second one is the reference: a commit, a directory or REF (the references given, test case by test case, in the 'testing' entry of the src/phyex/\<model\>\_version.json file).

The version of the model to use is described in the src/phyex/\<model\>\_version.json file (for instance, repository and commit of IAL in ial\_version.json); the model source code is downloaded automatically.

Example for the test programs, with the source code of a local repository:
```
phyex-check_testprogs -p -c -r -C --computeRefIfNeeded $HOME/PHYEX REF
```
The --computeRefIfNeeded option computes the reference simulations that do not exist yet.

### phyex-check\_phyex1d

This command checks the 1D model (see the [phyex1d documentation](./phyex1d.md)) with the same steps and arguments
as phyex-check\_testprogs: the source code is copied (or cloned) in $PHYEX1DDIR/\<name\>/PHYEX, a virtual environment
is created next to it (-p), the pyphyex package is compiled and installed in this virtual environment
(-c, the architecture is given with -a as for phyex-check\_testprogs), the 1D model is run on DEPHY-SCM cases
(-r, ARMCU/REF by default, -t to give other cases; the cases are installed in ${PHYEXCONF}/phyex1d if needed)
and the results are compared to the reference (-C). The reference commits are given, case by case, in the
'testing' entry of the src/phyex/phyex1d\_version.json file.
As for phyex-check\_testprogs, the --expand and --shumantrans options apply the corresponding source code
transformations when building the package (PHYEX\_TRANSFORMATION CMake variable); the reference simulations
are always computed without transformation.

Example, with the source code of a local repository:
```
phyex-check_phyex1d -p -c -r -C --computeRefIfNeeded $HOME/PHYEX REF
```

### Reference hashes (phyex-check\_testprogs and phyex-check\_phyex1d)

After each simulation, phyex-check\_testprogs and phyex-check\_phyex1d write the hash of the results in a
hash\_result.txt file next to them (for the test programs, the hash of the output listing starting at
the first 'IBL' line, normalised as for the comparison; for phyex1d, the hash of the datasets of the
result file, the attributes being ignored). Two results have the same hash if, and only if, the
comparison step considers them identical.

With the --refHashFile FILE option, the comparison step (-C) uses a json file containing, for each
test case, the reference commit and the hash of its results:

  - if the reference commit of a case (given on the command line or in the 'testing' entry of the
    version json file, compared as strings: the full SHAs must be used in the json files) is the
    commit stored in the file, the results are compared to the stored hash; the reference simulation
    is neither needed nor computed. In case of difference, the detailed comparison is obtained by
    running the comparison step again without the --refHashFile option.
  - otherwise, the reference simulation is used (and computed if --computeRefIfNeeded is set) and its
    commit and hash are stored in the file (created if it does not exist).

The hashes depend on the configuration (architecture, compiler, precision...): a file must be used
for a single configuration. This mechanism is used by the GitHub workflow (.github/workflows/tests.yaml),
the files being stored in the GitHub cache (see .github/actions/refhash/action.yml): the references are
compiled and run only when their hash is unknown.

### phyex-check\_norms

This command checks the coding norms (described in the [CodingNorms documentation](./CodingNorms.md)) of the source code (for instance: phyex-check\_norms -v src/phyex).

### phyex-check\_consistency

This command checks the consistency of the repository: the version (src/phyex/VERSION file) is the version of the last release tag, the minimal pyfortool version is the same everywhere.

### phyex-check\_commit

This command performs all the checks (consistency, coding norms, documentation generation, python wheel
creation and the test cases of IAL, Meso-NH, testprogs, LMDZ and phyex1d) on the repository containing the tools
(with its uncommitted modifications); the --only option restricts the checks and the --parallel option
checks the different models simultaneously.
The report is written on the standard output; it can also be sent by email in case of an error (--mail)
and attached as a status to a github commit (--github-status, the API url to use must be given).

### clone\_and\_test.sh

This script, at the root of the repository, is designed to be run periodically by cron. It clones (or updates) the PHYEX repository in ${PHYEXWORKDIR} (default: $HOME/PHYEXTESTING), checks out the last commit of a branch (or a given commit), installs the tools and runs phyex-check\_commit, the result being attached as a status to the commit on GitHub (and sent by email in case of an error).

### Others

  - phyex-compare: comparison of model outputs (used by the phyex-check\_\<model\> tools)
  - phyex-plot\_perf: plot of the performance statistics (files produced with the --perf option of the phyex-check\_\<model\> tools)
