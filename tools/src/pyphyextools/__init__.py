"""
pyphyextools python package
"""

import os
import shutil
import ssl
import subprocess
import sys
import tarfile
import tempfile
import urllib.request

# Root of the PHYEX repository containing this package (tools/src/pyphyextools)
REPO_DIR = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        '..', '..', '..'))

# Check if current package is installed in editable mode
if not os.path.exists(os.path.join(REPO_DIR, 'tools', 'pyproject.toml')):
    raise RuntimeError('pyphyextools package must be installed in editable mode')

# Version of the tools: the PHYEX version, stored in src/phyex/VERSION (also
# used by tools/pyproject.toml)
with open(os.path.join(REPO_DIR, 'src', 'phyex', 'VERSION'),
          encoding='utf-8') as _version_file:
    __version__ = _version_file.read().strip()


def run_command(cmd, cwd=None, out_path=None, mode='w', env=None, input_string=None,
                encoding='utf8', check=True, display='on-error', banner=True, err_path=True):
    """Run *cmd* in *cwd*, optionally teeing its output to *out_path*.

    The called command is echoed to standard output (unless *banner* is
    False) and its output (merged with its standard error by default, see
    *err_path*) is always returned as a :class:`~subprocess.CompletedProcess`.

    Parameters
    ----------
    cmd : list of str
        Command to execute.
    cwd : str or None
        Directory in which the command is run.
    out_path : str or None
        Path to which the output is appended/overwritten. If None,
        no file is written.
    mode : str
        Open mode for *out_path* and *err_path* (default ``'w'``).
    env : dict or None
        Environment for the subprocess.
    input_string : str or None
        Text written to the process stdin.
    encoding : str
        Encoding used for text I/O (default ``'utf8'``).
    check : bool
        If True and the process exits with non-zero status, raise a
        :class:`~subprocess.CalledProcessError` carrying the captured output
        (default True).
    display : {'never', 'always', 'on-error'}
        Whether to echo the captured output to the terminal: ``'never'``
        suppresses it, ``'always'`` shows it regardless of the exit status,
        and ``'on-error'`` only shows it when the process fails.
    banner : bool
        If True (default), the command is echoed (``### <command>``) before
        being executed.
    err_path : bool or str
        Standard error of the process: merged with the standard output (True,
        default), discarded (False) or written in the file *err_path* (it is
        then neither captured nor displayed).

    Returns
    -------
    subprocess.CompletedProcess
        Result with ``returncode`` and the output in ``stdout``.
    """
    if display not in ('never', 'always', 'on-error'):
        raise ValueError(
            f"Invalid display value '{display}': must be one of "
            "'never', 'always' or 'on-error'")
    if not isinstance(err_path, (bool, str)):
        raise ValueError(f"Invalid err_path value '{err_path}': must be a boolean or a string")

    if banner:
        print('### ' + ' '.join(map(str, cmd)))
        sys.stdout.flush()

    out = None
    err = None
    if out_path is not None:
        out = open(out_path, mode, encoding=encoding)
    try:
        if err_path is True:
            stderr = subprocess.STDOUT
        elif err_path is False:
            stderr = subprocess.DEVNULL
        else:
            err = open(err_path, mode, encoding=encoding)
            stderr = err
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=stderr,
            bufsize=1, text=True, encoding=encoding, cwd=cwd, env=env,
            stdin=subprocess.PIPE if input_string is not None else None)
        if input_string is not None:
            proc.stdin.write(input_string)
            proc.stdin.close()
        captured = []
        for line in proc.stdout:
            captured.append(line)
            if out is not None:
                out.write(line)
            if display == 'always':
                print(line, end='')
        proc.wait()
        if display == 'on-error' and proc.returncode != 0:
            for line in captured:
                print(line, end='')
    finally:
        if out is not None:
            out.close()
        if err is not None:
            err.close()
    result = subprocess.CompletedProcess(
        args=cmd, returncode=proc.returncode, stdout=''.join(captured))
    if check and proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd,
                                            output=result.stdout)
    return result


def download_and_extract(url, directory):
    """Download the tar archive *url* and extract it in *directory*.

    The SSL certificate of the server is not verified (as with
    ``wget --no-check-certificate``). The archive is extracted as the ``tar``
    command would do (``'tar'`` extraction filter) and is then removed.

    Parameters
    ----------
    url : str
        URL of the archive.
    directory : str
        Directory in which the archive is extracted.
    """
    print(f'### Download {url}')
    sys.stdout.flush()
    with tempfile.NamedTemporaryFile(dir=directory, suffix='.tar', delete=False) as archive:
        try:
            with urllib.request.urlopen(url, context=ssl._create_unverified_context()) as response:  # pylint: disable=protected-access
                shutil.copyfileobj(response, archive)
            archive.close()
            with tarfile.open(archive.name) as tar:
                if hasattr(tarfile, 'tar_filter'):
                    tar.extractall(directory, filter='tar')
                else:  # python without extraction filters
                    tar.extractall(directory)
        finally:
            archive.close()
            os.remove(archive.name)
