#!/bin/bash

# Clones (or updates) the PHYEX repository, installs the tools and runs
# phyex-check_commit on a commit, the report being attached to the commit on
# GitHub. See the usage function.

set -e
set -o pipefail #abort if left command on a pipe fails

function usage {
  echo "Usage: $0 [-h] [--repo-user USER] [--repo-protocol PROTOCOL] [--repo-repo REPO]"
  echo "                         [--commit SHA | --ref REF] [--no-remove] [--hostname HOSTNAME] [PHYEX-CHECK_COMMIT OPTIONS]"
  echo "--repo-user USER"
  echo "                user hosting the PHYEX repository on github,"
  echo "                defaults to the env variable PHYEXREPOuser (or UMR-CNRM)"
  echo "--repo-protocol PROTOCOL"
  echo "                protocol (https or ssh) to reach the PHYEX repository on github,"
  echo "                defaults to the env variable PHYEXREPOprotocol (or ssh)"
  echo "--repo-repo REPO"
  echo "                repository name, defaults to the env variable PHYEXREPOrepo (or PHYEX)"
  echo "--commit SHA    commit to check"
  echo "--ref REF       check the last commit of REF (defaults to refs/heads/master)"
  echo "--no-remove     do not remove the compilation directories of the models (by default, the --remove"
  echo "                option is given to phyex-check_commit)"
  echo "--hostname HOSTNAME"
  echo "                hostname used to build the log file name and the context of the github status,"
  echo "                (continuous-integration/HOSTNAME) instead of the real hostname."
  echo "                This can be useful when running on a cluster node which can vary from run to run"
  echo ""
  echo "The other options are given to the phyex-check_commit script (e.g. --only, --mail, --force, --perf,"
  echo "--insecure, --options-<model>, options for the model checks...), see its -h option for details."
  echo ""
  echo "This script is designed to be run periodically by cron to test the last commit of a branch"
  echo "(eg '00 22 * * * bash -l -c \"PHYEXWORKDIR=~/PHYEXTESTING ~/PHYEXTESTING/PHYEX/clone_and_test.sh \\"
  echo "                             --repo-user UMR-CNRM --repo-protocol ssh --repo-repo PHYEX --mail user@domain\"')."
  echo "Steps:"
  echo "  - the commit is checked out in the \${PHYEXWORKDIR}/PHYEX clone (created if needed); if this"
  echo "    version of the clone_and_test.sh script differs from the running one, it is executed instead;"
  echo "  - the tools are installed (pip install -e) in the \${PHYEXWORKDIR}/venv virtual environment"
  echo "    (created if needed);"
  echo "  - phyex-check_commit is run, with the report attached as a status to the commit on github"
  echo "    (the github token is read in the .netrc file). Nothing is done if the commit already has"
  echo "    a status with the same context (unless --force is given)."
  echo "All the output is written in \${PHYEXWORKDIR}/logfile_HOSTNAME (previous files are kept with the"
  echo ".old and .old2 suffixes). If an error occurs before running phyex-check_commit, an error status"
  echo "is attached to the commit (and a mail is sent if --mail is given)."
  echo ""
  echo "All the work is done within the \${PHYEXWORKDIR} directory (it defaults to ~/PHYEXTESTING)."
  echo "The TESTPROGSDIR, HOMEPACK, MNHPACK and LMDZPACK environment variables, used by the model checks,"
  echo "are set to the TESTPROGS, pack, MesoNH and LMDZ subdirectories of \${PHYEXWORKDIR}."
  echo "It may be necessary to fill these directories with base source code (for arome and/or mesonh),"
  echo "with data (testprogs) and with the reference simulations."
}

origargs=("$@")
#The repository variables are exported (they are used by the model checks)
export PHYEXREPOuser=${PHYEXREPOuser:=UMR-CNRM}
export PHYEXREPOprotocol=${PHYEXREPOprotocol:=ssh}
export PHYEXREPOrepo=${PHYEXREPOrepo:=PHYEX}
WORKDIR=${PHYEXWORKDIR:=${HOME}/PHYEXTESTING}
commit=""
ref="refs/heads/master"
contextHostname=$(hostname)
mail=""
insecure=""
remove="--remove"
checkCommitOptions=()
while [ $# -gt 0 ]; do
  case "$1" in
    '-h') usage; exit;;
    '--repo-user') export PHYEXREPOuser=$2; shift;;
    '--repo-protocol') export PHYEXREPOprotocol=$2; shift;;
    '--repo-repo') export PHYEXREPOrepo=$2; shift;;
    '--commit') commit=$2; shift;;
    '--ref') ref=$2; shift;;
    '--no-remove') remove="";;
    '--hostname') contextHostname=$2; shift;;
    # Options also given to phyex-check_commit, used here to report errors
    '--mail') mail=$2; checkCommitOptions+=("$1" "$2"); shift;;
    '--insecure') insecure="--insecure"; checkCommitOptions+=("$1");;
    *) checkCommitOptions+=("$1");;
  esac
  shift
done
if [ "${PHYEXREPOprotocol}" == 'ssh' ]; then
  giturl="git@github.com:${PHYEXREPOuser}/${PHYEXREPOrepo}.git"
elif [ "${PHYEXREPOprotocol}" == 'https' ]; then
  giturl="https://github.com/${PHYEXREPOuser}/${PHYEXREPOrepo}.git"
else
  echo "--repo-protocol must be ssh or https" >&2
  exit 1
fi
context="continuous-integration/${contextHostname}"

#Output redirection (when the script has been re-executed after an update, the log file
#is continued)
mkdir -p "${WORKDIR}"
logfile="${WORKDIR}/logfile_${contextHostname}"
if [ "${CLONE_AND_TEST_UPDATED}" == 1 ]; then
  exec >> "${logfile}" 2>&1
else
  if [ -f "${logfile}" ]; then
    [ -f "${logfile}.old" ] && mv "${logfile}.old" "${logfile}.old2"
    mv "${logfile}" "${logfile}.old"
  fi
  exec > "${logfile}" 2>&1
fi
echo "$(date): $0 ${origargs[*]}"

#Error reporting (before phyex-check_commit is run)
SHA=""
function fail {
  message="${context}\n\n$(date)\n__ ABNORMAL EXIT __\n$1\n\nMore information can be found in $(hostname):${logfile}"
  echo -e "${message}"
  if [ "${mail}" != "" ]; then
    if command -v mail > /dev/null; then
      echo -e "${message}" | mail -s "${context}" "${mail}"
    else
      echo -e "Subject: ${context}\n\n${message}" | sendmail "${mail}"
    fi
  fi
  if [ "${SHA}" != "" ]; then
    if TOKEN=$(python3 -c "import netrc; print(netrc.netrc().authenticators('github.com')[2])"); then
      #The description is limited to 140 characters
      data=$(description="$1 (see $(hostname):${logfile})" context="${context}" python3 -c \
             "import json, os; print(json.dumps({'state': 'error', 'context': os.environ['context'],
                                                 'description': os.environ['description'][:140]}))")
      curl -sS -L ${insecure} -X POST \
           -H "Accept: application/vnd.github+json" \
           -H "Authorization: Bearer ${TOKEN}" \
           -H "X-GitHub-Api-Version: 2022-11-28" \
           "https://api.github.com/repos/${PHYEXREPOuser}/${PHYEXREPOrepo}/statuses/${SHA}" \
           -d "${data}" || echo "Cannot attach the error status to the commit"
    else
      echo "No github token found in the .netrc file, the error status cannot be attached to the commit"
    fi
  fi
  exit 1
}

#Commit to test
if [ "${commit}" != "" ]; then
  SHA=${commit}
else
  echo "Getting last commit hash of ${ref}"
  SHA=$(git ls-remote "${giturl}" "${ref}" | cut -f1) || fail "Cannot get the last commit of ${ref}"
  [ "${SHA}" == "" ] && fail "${ref} not found in ${giturl}"
fi
echo "Commit hash is ${SHA}"

#Clone or update the repository, and checkout the commit
PHYEXDIR="${WORKDIR}/PHYEX"
if [ ! -d "${PHYEXDIR}" ]; then
  echo "Cloning ${giturl} in ${PHYEXDIR}"
  git clone "${giturl}" "${PHYEXDIR}" || fail "Cannot clone ${giturl}"
else
  echo "Fetching ${giturl} in ${PHYEXDIR}"
  git -C "${PHYEXDIR}" fetch --tags --force "${giturl}" '+refs/heads/*:refs/remotes/origin/*' || \
    fail "Cannot fetch ${giturl}"
fi
echo "Checkout commit ${SHA}"
git -C "${PHYEXDIR}" checkout --force --detach "${SHA}" || fail "Cannot checkout ${SHA}"

#Use the version of this script of the tested commit
if [ "${CLONE_AND_TEST_UPDATED}" != 1 -a -f "${PHYEXDIR}/clone_and_test.sh" ]; then
  if ! cmp -s "${BASH_SOURCE[0]}" "${PHYEXDIR}/clone_and_test.sh"; then
    echo "Script has changed, running the new version"
    CLONE_AND_TEST_UPDATED=1 exec bash "${PHYEXDIR}/clone_and_test.sh" "${origargs[@]}"
  fi
fi

#Installation of the tools (and of their dependencies, among them pyfortool, if
#they are not already installed)
venv="${WORKDIR}/venv"
if [ ! -f "${venv}/bin/activate" ]; then
  echo "Creating the virtual environment ${venv}"
  python3 -m venv "${venv}" || fail "Cannot create the virtual environment ${venv}"
fi
echo "Installing the tools"
"${venv}/bin/pip" install --prefer-binary -e "${PHYEXDIR}/tools" || fail "Cannot install the tools"

#Directories used by the model checks
export TESTPROGSDIR="${WORKDIR}/TESTPROGS"
export HOMEPACK="${WORKDIR}/pack"
export MNHPACK="${WORKDIR}/MesoNH"
export LMDZPACK="${WORKDIR}/LMDZ"
for d in "${TESTPROGSDIR}" "${HOMEPACK}" "${MNHPACK}" "${LMDZPACK}"; do
  if [ ! -d "${d}" ]; then
    echo "Creating directory ${d}"
    mkdir -p "${d}"
  fi
done

#Checks
exec "${venv}/bin/phyex-check_commit" \
     --github-status "https://api.github.com/repos/${PHYEXREPOuser}/${PHYEXREPOrepo}/statuses/${SHA}" \
     --github-target-url "https://${PHYEXREPOuser}.github.io/${PHYEXREPOrepo}/displayparam.html?" \
     --hostname "${contextHostname}" ${remove} "${checkCommitOptions[@]}"
