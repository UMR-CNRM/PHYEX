"""
pyphyextools python package
"""

import os
import subprocess

__version__ = "0.9.0"

# Check if current package is installed in editable mode
pyproject = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         '..', '..', 'pyproject.toml')
if not os.path.exists(pyproject):
    raise RuntimeError('pyphyextools package must be installed in editable mode')


def run_command(cmd, cwd=None, out_path=None, mode='w', env=None, input_string=None,
                encoding='utf8', check=True, display='on-error'):
    """Run *cmd* in *cwd*, optionally teeing its output to *out_path*.

    The called command is always echoed to standard output and the merged
    standard output/error is always returned as a
    :class:`~subprocess.CompletedProcess`.

    Parameters
    ----------
    cmd : list of str
        Command to execute.
    cwd : str or None
        Directory in which the command is run.
    out_path : str or None
        Path to which the merged output is appended/overwritten. If None,
        no file is written.
    mode : str
        Open mode for *out_path* (default ``'w'``).
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

    Returns
    -------
    subprocess.CompletedProcess
        Result with ``returncode`` and the merged output in ``stdout``.
    """
    if display not in ('never', 'always', 'on-error'):
        raise ValueError(
            f"Invalid display value '{display}': must be one of "
            "'never', 'always' or 'on-error'")

    print('### ' + ' '.join(map(str, cmd)))

    out = None
    if out_path is not None:
        out = open(out_path, mode, encoding=encoding)
    try:
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
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
    result = subprocess.CompletedProcess(
        args=cmd, returncode=proc.returncode, stdout=''.join(captured))
    if check and proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd,
                                            output=result.stdout)
    return result
