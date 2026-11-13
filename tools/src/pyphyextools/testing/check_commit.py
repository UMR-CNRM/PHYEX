"""
Global check of the PHYEX repository containing these tools (with its
uncommitted modifications): consistency, coding norms, documentation
generation, python wheel creation and the test cases of the different models
(IAL, Meso-NH, testprogs and LMDZ).

The report is written on the standard output. Optionally:
  - it is sent by e-mail in case of failure (--mail);
  - it is attached, as a status, to a commit on GitHub (--github-status, the
    GitHub token is read from the .netrc file). If the commit already has a
    status with the same context, the checks are not performed (unless
    --force is given).
To check a given commit, the repository must be cloned, the commit checked
out and the tools installed from this clone before running this script, e.g.:
  git clone https://github.com/UMR-CNRM/PHYEX.git && cd PHYEX && \\
  git checkout <sha> && python3 -m venv venv && . venv/bin/activate && \\
  pip install -e tools && \\
  phyex-check_commit --mail user@domain \\
    --github-status https://api.github.com/repos/UMR-CNRM/PHYEX/statuses/<sha> \\
    --github-target-url 'https://UMR-CNRM.github.io/PHYEX/displayparam.html?'

The model checks (phyex-check_<model>) are run on the source code of the
repository, after the other checks; with --parallel, the different models are
checked simultaneously (their output lines are then prefixed by [<model>]). For each model, the test cases and their references are read in
the src/phyex/<model>_version.json file ('testing' entry); without such an
entry, the default test cases are compared against the 'REF' reference.
The model checks use their own environment variables (HOMEPACK, MNHPACK,
TESTPROGSDIR, LMDZPACK...).
"""

import argparse
import concurrent.futures
import copy
import datetime
import json
import netrc
import os
import shutil
import socket
import ssl
import subprocess
import sys
import tempfile
import threading
import traceback
import urllib.parse
import urllib.request

from pyphyextools.testing.check_consistency import REPO_DIR

MODELS = ('ial', 'mesonh', 'testprogs', 'lmdz')
CHECKS = ('consistency', 'norms', 'doc', 'wheel') + MODELS

# Options of the different steps of the model checks
MODEL_STEPS = {
    # model: (compilation, execution, comparison (None to disable))
    'ial': (['-p', '-c'], ['-r'], ['-C', '--computeRefIfNeeded']),
    'mesonh': (['-p', '-c'], ['-r'], ['-C', '--computeRefIfNeeded']),
    'testprogs': (['-p', '-c'], ['-r'], ['-C', '--computeRefIfNeeded']),
    'lmdz': (['-p', '-c', '--nofcm'], ['-r', '--nofcm'], None),
}


class Reporter:
    """
    Collects the results and reports them (standard output and, optionally,
    e-mail and GitHub status)
    """

    def __init__(self, args):
        self.args = args
        self.context = f'continuous-integration/{args.hostname}'
        self.message = []  # lines of the report
        self.log_message = []  # all the logged lines
        self.token = None
        if args.github_status is not None:
            # Token read at the beginning to fail early
            authenticators = netrc.netrc().authenticators('github.com')
            if authenticators is None:
                raise RuntimeError('No GitHub token found (github.com) in the .netrc file')
            self.token = authenticators[2]

    def log(self, message, report=True):
        """
        Logs a message
        :param report: True if the message must also be in the report
        """
        print(message)
        sys.stdout.flush()
        self.add(message, report)

    def add(self, message, report=True):
        """
        Adds a message (already printed)
        :param report: True if the message must also be in the report
        """
        if report:
            self.message.append(message)
        self.log_message.append(message)

    def github(self, url, data=None):
        """
        GitHub API request
        :param data: dictionary to post (GET request if None)
        :return: (decoded json response, url of the next page or None)
        """
        request = urllib.request.Request(
            url, data=None if data is None else json.dumps(data).encode('utf-8'),
            headers={'Accept': 'application/vnd.github+json',
                     'Authorization': f'Bearer {self.token}',
                     'X-GitHub-Api-Version': '2022-11-28'})
        ssl_context = None
        if self.args.insecure:
            ssl_context = ssl._create_unverified_context()  # pylint: disable=protected-access
        with urllib.request.urlopen(request, context=ssl_context) as response:
            next_url = None
            for link in response.headers.get('Link', '').split(','):
                if 'rel="next"' in link:
                    next_url = link.split(';')[0].strip().strip('<>')
            return json.loads(response.read().decode('utf-8')), next_url

    def already_reported(self):
        """
        True if the commit already has a status with the same context (GET on
        the --github-status url, which lists the statuses of the commit)
        """
        url = self.args.github_status + '?per_page=100'
        while url is not None:
            statuses, url = self.github(url)
            if any(status.get('context') == self.context for status in statuses):
                return True
        return False

    def add_status(self, state, message):
        """Attaches a status to the commit (POST on the --github-status url)"""
        data = {'state': state, 'context': self.context}
        if self.args.github_target_url is not None:
            data['target_url'] = (self.args.github_target_url +
                                  urllib.parse.quote('<pre>' + message + '</pre>', safe=':/='))
        self.github(self.args.github_status, data)

    def send_mail(self, message):
        """Sends the report by e-mail"""
        if shutil.which('mail'):
            cmd = ['mail', '-s', self.context, self.args.mail]
        else:
            cmd = ['sendmail', self.args.mail]
            message = f'Subject: {self.context}\n\n{message}'
        subprocess.run(cmd, input=message, text=True, check=True)

    def report(self, success, error=None):
        """
        Reports the result
        :param success: True if all the checks succeeded
        :param error: abnormal exit message (traceback), if any
        """
        lines = [self.context, '', f'Repository {REPO_DIR}']
        try:
            lines.append('Commit ' + subprocess.run(['git', 'describe', '--always', '--dirty'],
                                                    cwd=REPO_DIR, check=True, text=True,
                                                    capture_output=True).stdout.strip())
        except (OSError, subprocess.CalledProcessError):
            pass
        lines.append(self.args.start)
        if error is not None:
            lines.extend(['__ ABNORMAL EXIT __'] + self.log_message + [error])
        else:
            lines.extend(self.message)
        lines.append('')
        if not success:
            lines.extend(['XXXXXXXXXXXXXXXXXXXX ERROR XXXXXXXXXXXXXXXXXXXX', ''])
        lines.append(str(datetime.datetime.now()))
        message = '\n'.join(lines)

        print('\n' + '=' * 70 + '\n' + message)
        sys.stdout.flush()
        if self.args.mail and not success:
            self.send_mail(message)
        if self.args.github_status is not None:
            self.add_status('success' if success else ('error' if error else 'failure'),
                            message)


class CheckCommit:
    """Runs the different checks on a PHYEX directory"""

    def __init__(self, args, reporter, phyexdir):
        self.args = args
        self.reporter = reporter
        self.log = reporter.log
        self.phyexdir = phyexdir
        self.prefix = None  # prefix of the output lines (parallel execution)
        self.lock = threading.Lock()  # for the output of the parallel executions

    def run(self, cmd, cwd=None, input_string=None):
        """
        Runs a command (output not captured, except in parallel execution where
        the output lines are prefixed)
        :return: True in case of success
        """
        self.log('### ' + ' '.join(cmd), report=False)
        sys.stdout.flush()
        sys.stderr.flush()
        try:
            if self.prefix is None:
                result = subprocess.run(cmd, cwd=cwd, input=input_string,
                                        text=input_string is not None, check=False)
                returncode = result.returncode
            else:
                with subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE,
                                      stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                      text=True, errors='replace') as process:
                    for line in process.stdout:
                        with self.lock:
                            print(self.prefix + line, end='')
                            sys.stdout.flush()
                    returncode = process.wait()
        except OSError as err:
            self.log(f'Command failed: {err}', report=False)
            return False
        return returncode == 0

    def global_result(self, name, success):
        """Logs the global result of a check"""
        if success:
            self.log(f'..... global result for {name}: OK')
        else:
            self.log(f'XXXXX global result for {name}: ERROR')
        return success

    def check_consistency(self):
        """Consistency checks (version, pyfortool version)"""
        self.log('Test consistency')
        success = self.run([sys.executable, '-m', 'pyphyextools.testing.check_consistency',
                            '--repo-dir', self.phyexdir])
        self.log('  consistency: ' + ('OK' if success else 'error'))
        return self.global_result('consistency', success)

    def check_norms(self):
        """Coding norms"""
        self.log('Test coding conventions')
        success = self.run([sys.executable, '-m', 'pyphyextools.testing.coding_norms',
                            os.path.join(self.phyexdir, 'src', 'phyex')])
        self.log('  check coding conventions: ' + ('OK' if success else 'non-conforming'))
        return self.global_result('coding conventions', success)

    def check_doc(self):
        """Documentation generation (same configuration as the doxygen workflow)"""
        self.log('Test doc generation')
        doxygen = shutil.which('doxygen')
        if doxygen is None:
            self.log('  doc generation with doxygen: error (doxygen not found)')
            return self.global_result('doc generation', False)
        with open(os.path.join(self.phyexdir, 'docs', 'doxygen', 'doxygen_config'),
                  encoding='utf-8') as config_file:
            config = config_file.read()
        with open(os.path.join(self.phyexdir, 'src', 'phyex', 'VERSION'),
                  encoding='utf-8') as version_file:
            version = version_file.read().strip()
        with tempfile.TemporaryDirectory() as outdir:
            # The documentation is generated outside of the repository
            config += f'\nPROJECT_NUMBER = {version}\nOUTPUT_DIRECTORY = {outdir}\n'
            success = self.run([doxygen, '-'], cwd=os.path.join(self.phyexdir, 'docs', 'doxygen'),
                               input_string=config)
        self.log('  doc generation with doxygen: ' + ('OK' if success else 'error'))
        return self.global_result('doc generation', success)

    def check_wheel(self):
        """Python wheel creation"""
        self.log('Test python wheel creation')
        with tempfile.TemporaryDirectory() as wheeldir:
            success = self.run([sys.executable, '-m', 'pip', 'wheel', '--no-deps',
                                '-w', wheeldir, self.phyexdir])
            if success:
                wheels = os.listdir(wheeldir)
                self.log(f'  wheel creation: OK ({", ".join(wheels)})')
            else:
                self.log('  wheel creation: error')
        return self.global_result('python wheel creation', success)

    def get_cases(self, model):
        """
        Test cases and their references
        :return: list of (case, reference), case is None for the default cases
        """
        jsonfile = os.path.join(self.phyexdir, 'src', 'phyex', f'{model}_version.json')
        cases = []
        if os.path.exists(jsonfile):
            with open(jsonfile, encoding='utf-8') as json_file:
                cases = list(json.load(json_file).get('testing', {}).items())
        if len(cases) == 0:
            self.log(f'No cases found in {jsonfile}, we only test the default cases',
                     report=False)
            cases = [(None, 'REF')]
        return cases

    def check_model(self, model):
        """Compilation, execution, comparison and cleaning for a model"""
        self.log(f'Tests for model {model}')
        compilation, execution, comparison = MODEL_STEPS[model]
        if self.args.suppress:
            compilation = compilation + ['-s']
        cmd = [sys.executable, '-m', f'pyphyextools.testing.check_{model}']
        if self.args.name is not None:
            cmd += ['--name', self.args.name]
        cmd += self.args.check_options + getattr(self.args, f'options_{model}')
        if self.args.perf is not None:
            execution = execution + ['--perf', self.args.perf]

        # Compilation
        success = self.run(cmd + compilation + [self.phyexdir])
        self.log(f'  {model} compilation: ' + ('OK' if success else 'error'))

        if success:
            for case, ref in self.get_cases(model):
                if case is None:
                    casearg = []
                    casedescr = 'default case(s)'
                else:
                    casearg = ['-t', case]
                    casedescr = f'{case} (ref={ref})'
                    self.log(f'Testing case {case} against reference {ref}', report=False)

                # Execution
                result = self.run(cmd + execution + casearg + [self.phyexdir])
                if result:
                    self.log(f'  {model} {casedescr}: execution OK (but status not reliable)')
                else:
                    self.log(f'  {model} {casedescr}: execution error')
                    success = False

                # Comparison
                if result and comparison is not None:
                    result = self.run(cmd + comparison + casearg + [self.phyexdir, ref])
                    self.log(f'  {model} {casedescr}: comparison ' + ('OK' if result else 'error'))
                    success = success and result

            # Cleaning
            if self.args.remove:
                result = self.run(cmd + ['--remove', self.phyexdir])
                self.log(f'  {model}: cleaning ' + ('OK' if result else 'error'))
                success = success and result

        return self.global_result(f'model {model}', success)

    def check_all(self):
        """
        Runs all the requested checks
        :return: True if all the checks succeeded
        """
        success = True
        for check in [check for check in CHECKS
                      if check in self.args.only and check not in MODELS]:
            success = getattr(self, f'check_{check}')() and success

        models = [model for model in MODELS if model in self.args.only]
        if self.args.parallel and len(models) > 1:
            with concurrent.futures.ThreadPoolExecutor(max_workers=len(models)) as executor:
                futures = [executor.submit(self.check_model_parallel, model) for model in models]
                # Results (and report lines) in the order of the models
                for future in futures:
                    result, lines = future.result()
                    for message, report in lines:
                        self.reporter.add(message, report)
                    success = result and success
        else:
            for model in models:
                success = self.check_model(model) and success
        return success

    def check_model_parallel(self, model):
        """
        check_model executed in a thread: output lines are prefixed by the
        model name, and the logged messages are returned instead of being
        added to the report (to keep the messages of a model together)
        :return: (result of check_model, list of (message, report))
        """
        lines = []
        worker = copy.copy(self)
        worker.prefix = f'[{model}] '

        def log(message, report=True):
            with self.lock:
                print(worker.prefix + message)
                sys.stdout.flush()
            lines.append((message, report))
        worker.log = log
        return worker.check_model(model), lines


def parse_arguments():
    """Parses the command line arguments"""
    parser = argparse.ArgumentParser(
        description='Global check of the PHYEX repository containing these tools',
        epilog='Unrecognized options are given to all the model check scripts.\n\n' + __doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--only', action='append', choices=CHECKS, default=None,
                        help='performs only this check (option can be provided several times), '
                             'by default all the checks are performed')
    parser.add_argument('-s', action='store_true', dest='suppress',
                        help='suppress the existing compilation directories of the models '
                             '(-s option given to the model checks at the compilation step)')
    parser.add_argument('--remove', action='store_true',
                        help='remove the compilation directories of the models after their '
                             'checks (when the compilation succeeded)')
    parser.add_argument('--parallel', action='store_true',
                        help='check the different models simultaneously (after the other '
                             'checks)')
    parser.add_argument('--name', default=None,
                        help='name given to the model checks (--name option), used for their '
                             'compilation directories and in the performance file; it should be '
                             'the commit hash, to be able to use these directories as references '
                             '(by default, the model checks build a name from the path of the '
                             'repository)')
    parser.add_argument('--perf', metavar='FILE', default=None,
                        help='add performance statistics in file FILE')
    parser.add_argument('--mail', default=None,
                        help='comma-separated list of e-mail addresses (no spaces) to which '
                             'the report is sent in case of failure')
    parser.add_argument('--github-status', metavar='URL', default=None,
                        help='attach the report, as a status, to a GitHub commit, URL is the '
                             'GitHub API url to use (https://api.github.com/repos/<user>/<repo>/'
                             'statuses/<sha>)')
    parser.add_argument('--github-target-url', metavar='URL', default=None,
                        help='target url of the GitHub status, the report (url-quoted) is '
                             'appended to URL (e.g. https://<user>.github.io/<repo>/'
                             'displayparam.html?); without this option, the status has no '
                             'target url')
    parser.add_argument('--force', action='store_true',
                        help='with --github-status, perform the checks even if the commit '
                             'already has a status with the same context')
    parser.add_argument('--hostname', default=socket.gethostname(),
                        help='the context of the report (continuous-integration/HOSTNAME) is '
                             'built using the provided hostname instead of the real hostname. '
                             'This can be useful when running on a cluster node which can vary '
                             'from run to run')
    parser.add_argument('--insecure', action='store_true',
                        help='do not verify the SSL certificates for the GitHub API request')
    for model in MODELS:
        parser.add_argument(f'--options-{model}', metavar='OPTIONS', default='',
                            help=f'options to give to the check script for model {model}')

    # The --options-<model> values start with '-' and would be rejected by
    # argparse: they are extracted before parsing the other arguments
    argv = sys.argv[1:]
    options = {model: [] for model in MODELS}
    remaining = []
    while argv:
        arg = argv.pop(0)
        model = arg[len('--options-'):] if arg.startswith('--options-') else None
        if model in MODELS:
            if not argv:
                parser.error(f'argument {arg}: expected one argument')
            options[model].extend(argv.pop(0).split())
        else:
            remaining.append(arg)
    args, check_options = parser.parse_known_args(remaining)
    args.check_options = check_options
    for model in MODELS:
        setattr(args, f'options_{model}', options[model])
    if args.only is None:
        args.only = list(CHECKS)
    if args.github_target_url is not None and args.github_status is None:
        parser.error('--github-target-url is only allowed with --github-status')
    if args.force and args.github_status is None:
        parser.error('--force is only allowed with --github-status')
    if args.insecure and args.github_status is None:
        parser.error('--insecure is only allowed with --github-status')
    return args


def main():
    """Main function"""
    args = parse_arguments()
    args.start = str(datetime.datetime.now())
    reporter = Reporter(args)
    if args.github_status is not None and not args.force and reporter.already_reported():
        print(f'The commit already has a {reporter.context} status, nothing is done '
              '(use --force to perform the checks anyway)')
        return
    try:
        success = CheckCommit(args, reporter, REPO_DIR).check_all()
    except Exception:  # pylint: disable=broad-exception-caught
        error = traceback.format_exc()
        print(error)
        reporter.report(False, error)
        sys.exit(1)
    reporter.report(success)
    if not success:
        sys.exit(1)


if __name__ == '__main__':
    main()
