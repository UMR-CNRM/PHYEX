"""
List of the changes between two commits of the PHYEX repository, sorted by
category according to the title of the pull request (PR) which introduced each
commit (see the PR rules in docs/Developer.md):
  - "Bugfixes": PR title beginning with "bugfix:";
  - "Features": PR title beginning with "feature:";
  - "Bugfixes and features": PR title beginning with "bugfix+feature:";
  - "Others": commits which do not come from a PR, or which come from a PR whose
    title does not follow the rule (PR merged before the rule was set up).
The case of the title prefix is not significant.

For each commit, a line gives the PR number (if any), the date (author date),
the abbreviated hash, the author and the first line of the commit message.

The commits are those reachable from the second commit but not from the first
one (git log FIRST..SECOND), merge commits excluded (unless --merges is given);
they must be available in the local repository. The PR of each commit is
obtained from the GitHub API (pull request which introduced the commit in the
repository): the commits must have been pushed on GitHub.

The GitHub token, optional for a public repository but recommended (without
token, the GitHub API is limited to 60 requests per hour, one request being
done per commit), is read in the GITHUB_TOKEN or GH_TOKEN environment variable
or, otherwise, in the .netrc file (machine api.github.com or github.com).

Example:
  phyex-changelog v0.9.3 v0.9.4
  phyex-changelog --format markdown v0.9.4 HEAD
"""

import argparse
import concurrent.futures
import json
import netrc
import os
import re
import ssl
import subprocess
import sys
import urllib.error
import urllib.request

from pyphyextools import REPO_DIR

# Same rule as in the .github/workflows/pr_title.yaml workflow
TITLE_RE = re.compile(r'^(bugfix|feature|bugfix\+feature):\s', re.IGNORECASE)

CATEGORIES = {'bugfix': 'Bugfixes',
              'feature': 'Features',
              'bugfix+feature': 'Bugfixes and features',
              None: 'Others'}


def get_token():
    """
    :return: the GitHub token (environment variable or .netrc file), or None
    """
    for var in ('GITHUB_TOKEN', 'GH_TOKEN'):
        if os.environ.get(var):
            return os.environ[var]
    try:
        net = netrc.netrc()
    except (FileNotFoundError, netrc.NetrcParseError):
        return None
    for machine in ('api.github.com', 'github.com'):
        authenticators = net.authenticators(machine)
        if authenticators is not None:
            return authenticators[2]
    return None


def get_commits(git_dir, first, second, merges=False):
    """
    :param git_dir: git repository
    :param first, second: commits (hash, tag, branch...)
    :param merges: True to include the merge commits
    :return: list of commits (dictionaries with the 'sha', 'short', 'date',
             'author' and 'subject' keys), the oldest first
    """
    sep = '\x1f'
    fmt = sep.join(['%H', '%h', '%as', '%an', '%s'])
    cmd = ['git', '-C', git_dir, 'log', '--reverse', f'--format={fmt}']
    if not merges:
        cmd.append('--no-merges')
    cmd.append(f'{first}..{second}')
    output = subprocess.run(cmd, check=True, capture_output=True,
                            encoding='utf-8').stdout
    commits = []
    for line in output.splitlines():
        sha, short, date, author, subject = line.split(sep, 4)
        commits.append({'sha': sha, 'short': short, 'date': date,
                        'author': author, 'subject': subject})
    return commits


def get_pr(repo, sha, token, insecure=False):
    """
    :param repo: GitHub repository (OWNER/NAME)
    :param sha: commit hash
    :param token: GitHub token (or None)
    :param insecure: True to disable the SSL certificate verification
    :return: the PR (dictionary from the GitHub API) which introduced the
             commit, or None
    """
    headers = {'Accept': 'application/vnd.github+json',
               'X-GitHub-Api-Version': '2022-11-28'}
    if token is not None:
        headers['Authorization'] = f'Bearer {token}'
    request = urllib.request.Request(
        f'https://api.github.com/repos/{repo}/commits/{sha}/pulls', headers=headers)
    ssl_context = ssl._create_unverified_context() if insecure else None  # pylint: disable=protected-access
    try:
        with urllib.request.urlopen(request, context=ssl_context) as response:
            prs = json.loads(response.read().decode('utf-8'))
    except urllib.error.HTTPError as error:
        if error.code == 422:
            # Commit unknown on GitHub
            print(f'Warning: commit {sha} not found on GitHub ({repo})', file=sys.stderr)
            return None
        if error.code in (403, 429):
            raise RuntimeError(
                f'GitHub API request refused ({error.code}), probably because of the '
                'rate limit: give a GitHub token (GITHUB_TOKEN environment variable '
                'or .netrc file)') from error
        raise
    # Merged PR which introduced the commit; otherwise (commit not in the
    # default branch), an open PR containing the commit
    merged = sorted([pr for pr in prs if pr.get('merged_at')], key=lambda pr: pr['merged_at'])
    if len(merged) > 0:
        return merged[0]
    return prs[0] if len(prs) > 0 else None


def category(pr):
    """
    :param pr: PR (dictionary from the GitHub API) or None
    :return: category key (key of CATEGORIES)
    """
    if pr is None:
        return None
    match = TITLE_RE.match(pr['title'])
    return None if match is None else match.group(1).lower()


def changelog(first, second, *, git_dir=REPO_DIR,  # pylint: disable=too-many-arguments
              repo='UMR-CNRM/PHYEX', merges=False, token=None, insecure=False, jobs=8):
    """
    :param first, second: commits (hash, tag, branch...)
    :param git_dir: git repository
    :param repo: GitHub repository (OWNER/NAME)
    :param merges: True to include the merge commits
    :param token: GitHub token (or None)
    :param insecure: True to disable the SSL certificate verification
    :param jobs: number of simultaneous GitHub API requests
    :return: dictionary {category title: list of commits}, each commit being a
             dictionary (see get_commits) with an additional 'pr' key (PR
             number or None); the categories are in the CATEGORIES order
    """
    commits = get_commits(git_dir, first, second, merges)
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
        prs = list(executor.map(lambda commit: get_pr(repo, commit['sha'], token, insecure),
                                commits))
    result = {title: [] for title in CATEGORIES.values()}
    for commit, pr in zip(commits, prs):
        commit['pr'] = None if pr is None else pr['number']
        result[CATEGORIES[category(pr)]].append(commit)
    for title in result:
        # Sorted by PR (commits not coming from a PR at the end), then by date
        result[title].sort(key=lambda commit: (commit['pr'] is None, commit['pr'] or 0,
                                               commit['date']))
    return result


def format_changelog(result, fmt='text'):
    """
    :param result: output of changelog
    :param fmt: 'text' or 'markdown'
    :return: the changelog as a string (empty categories are omitted)
    """
    lines = []
    for title, commits in result.items():
        if len(commits) == 0:
            continue
        if fmt == 'markdown':
            lines.extend([f'## {title}', ''])
            for commit in commits:
                pr = '' if commit['pr'] is None else f"#{commit['pr']} "
                lines.append(f"- {pr}({commit['date']}, {commit['short']}, "
                             f"{commit['author']}): {commit['subject']}")
        else:
            lines.extend([title, '-' * len(title)])
            prwidth = max(len(str(commit['pr'] or '')) for commit in commits) + 1
            authorwidth = max(len(commit['author']) for commit in commits)
            for commit in commits:
                pr = '' if commit['pr'] is None else f"#{commit['pr']}"
                lines.append(f"{pr:{prwidth}}  {commit['date']}  {commit['short']}  "
                             f"{commit['author']:{authorwidth}}  {commit['subject']}")
        lines.append('')
    return '\n'.join(lines)


def main():
    """
    Command line entry point
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('first', help='first commit (hash, tag...), excluded')
    parser.add_argument('second', help='second commit (hash, tag...), included')
    parser.add_argument('--format', choices=['text', 'markdown'], default='text',
                        help='output format (default: text)')
    parser.add_argument('--repo', default='UMR-CNRM/PHYEX',
                        help='GitHub repository (OWNER/NAME, default: UMR-CNRM/PHYEX)')
    parser.add_argument('--git-dir', default=REPO_DIR,
                        help='local git repository (default: the repository containing '
                             'these tools)')
    parser.add_argument('--merges', action='store_true', default=False,
                        help='include the merge commits')
    parser.add_argument('--insecure', action='store_true', default=False,
                        help='disable the SSL certificate verification')
    parser.add_argument('-j', '--jobs', type=int, default=8,
                        help='number of simultaneous GitHub API requests (default: 8)')
    args = parser.parse_args()

    result = changelog(args.first, args.second, git_dir=args.git_dir, repo=args.repo,
                       merges=args.merges, token=get_token(), insecure=args.insecure,
                       jobs=args.jobs)
    print(format_changelog(result, args.format), end='')


if __name__ == '__main__':
    main()
