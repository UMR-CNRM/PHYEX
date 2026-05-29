"""Commit-check tool for the LMDZ model.

Compiles LMDZ using a specific commit for the externalised physics,
runs RICO and ARM-CU 1D cases.
"""

import glob as glob_mod
import json
import os
import shutil
import subprocess
import time

from pyphyextools.testing.check_common import (
    CheckCommitBase, CheckCommitError, escape_commit, mvdiff, run_tool)
from pyphyextools.prep_code import prep_code


class CheckCommitLmdz(CheckCommitBase):
    """Check a commit against the LMDZ reference."""

    default_expand = False
    default_buildSystem = 'fcm'
    alternative_buildSystem = 'make'

    def __init__(self, prepCodeOpts='', buildSys=None, **kwargs):
        super().__init__(**kwargs)
        self.defaultTest = ["rico"]
        self.allowedTests = ["rico", "arm_cu"]
        self.prepCodeOpts = prepCodeOpts
        self.buildSys = buildSys or self.default_buildSystem
        self.LMDZPACK = os.environ.get('LMDZPACK', os.path.expanduser('~/LMDZ/PHYEX'))
        self.version = None
        self.rad = None
        self.install_arg = None
        self.lmdzdir = None
        self.phyexdir = None
        self.compilecmd = None
        self.link = 0
        self.L = 79

        if os.path.isdir(os.path.join(self.commit, 'src')):
            self.model_ready = False
            json_path = os.path.join(self.commit, 'src', 'lmdz', 'lmdz_version.json')
        else:
            self.model_ready = True
            json_path = os.path.join(self.commit, 'lmdz_version.json')
        if os.path.isfile(json_path):
            with open(json_path, encoding='utf-8') as f:
                self.json_content = json.load(f)
        else:
            self.json_content = {}

        testing = self.json_content.get('testing', '')
        all_tests = []
        if testing:
            for t in self.allowedTests:
                ref = testing.get(t, '')
                if ref:
                    all_tests.append(t)
        if not self.tests:
            self.tests = list(self.defaultTest)
        elif 'ALL' in self.tests:
            expanded = []
            for t in self.tests:
                expanded.extend(all_tests if t == 'ALL' else [t])
            self.tests = expanded
        if not self.name:
            self.name = escape_commit(self.commit)

        self.version = self.json_content.get('version', '')
        self.rad = self.json_content.get('rad', '')
        self.install_arg = self.json_content.get('install_arg', '')
        main_exe = "lmdz1d"
        self.lmdzdir = os.path.join(self.LMDZPACK, self.name, 'LMDZ')
        self.phyexdir = os.path.join(self.LMDZPACK, self.name, 'PHYEX')
        self.compilecmd = (f"./compile -L {self.L} -rad {self.rad} "
                           f"-cosp 0 -main {main_exe}")

    @classmethod
    def parse_arguments(cls):
        """Build the argument parser, add LMDZ-specific args, parse and return the dict."""
        parser = super().parse_arguments()
        parser.description = (
            "Compile the LMDZ model using a specific commit for the "
            "externalised physics, run RICO and ARM-CU 1D cases")
        parser.add_argument('--prep_code-opts', dest='prepCodeOpts', default='',
                            help='options added to prep_code call')
        parser.add_argument(
            f'--no{cls.default_buildSystem}',
            action='store_const', dest='buildSys',
            const=cls.alternative_buildSystem,
            help=f'Do not use {cls.default_buildSystem}')
        args = parser.parse_args()
        return vars(args), cls._build_commitcmd(parser)

    def pack_creation(self):
        """Create the LMDZ pack via *install_lmdz.sh*."""
        packdir = os.path.join(self.LMDZPACK, self.name)

        if self.suppress and os.path.isdir(packdir):
            shutil.rmtree(packdir)

        if os.path.isdir(packdir) and self.onlyIfNeeded:
            return

        if os.path.isdir(packdir):
            raise CheckCommitError(
                f"Pack already exists ({packdir}),\n"
                "suppress it to be able to compile it again "
                "(or use the -s option to automatically suppress it)", 5)

        print(f"### Compilation of commit {self.commit}")
        os.makedirs(packdir, exist_ok=True)

        with open(os.path.join(packdir, 'arch-mylocal.env'), 'w', encoding='utf-8') as f:
            f.write('')

        with open(os.path.join(packdir, 'arch-mylocal.fcm'), 'w', encoding='utf-8') as f:
            f.write('%COMPILER            gfortran\n')
            f.write('%LINK                gfortran\n')
            f.write('%FPP                 cpp\n')
            f.write('%AR                  ar\n')
            f.write('%ARFLAGS             rU\n')
            f.write('%MAKE                make\n')
            f.write('%FPP_FLAGS           -P -traditional\n')
            f.write('%FPP_DEF             NC_DOUBLE\n')
            f.write('%BASE_FFLAGS          -cpp -ffree-line-length-0 -fdefault-real-8 '
                    '-DNC_DOUBLE\n')
            f.write('%PROD_FFLAGS         -O3 -fallow-argument-mismatch\n')
            f.write('%DEV_FFLAGS          -Wall -fbounds-check  -fallow-argument-mismatch\n')
            f.write(('%DEBUG_FFLAGS        -g3 -Wall -fbounds-check -ffpe-trap=invalid,zero,'
                     'overflow -O0 -fstack-protector-all -fbacktrace -finit-real=snan  '
                     '-fallow-argument-mismatch\n'))
            f.write('%MPI_FFLAGS\n')
            f.write('%OMP_FFLAGS\n')
            f.write('%BASE_LD\n')
            f.write('%MPI_LD\n')
            f.write('%OMP_LD\n')

        with open(os.path.join(packdir, 'arch-mylocal.path'), 'w', encoding='utf-8') as f:
            f.write('NETCDF_INCDIR="-I/usr/include"\n')
            f.write('NETCDF_LIBDIR=""\n')
            f.write('NETCDF_LIB="-lnetcdff -lnetcdf"\n\n')
            f.write('NETCDF95_INCDIR=-I$LMDGCM/../../include/\n')
            f.write('NETCDF95_LIBDIR=-L$LMDGCM/../../lib\n')
            f.write('NETCDF95_LIB=-lnetcdf95\n\n')
            f.write('IOIPSL_INCDIR="-I$LMDGCM/../../lib -I$LMDGCM/../IOIPSL/inc"\n')
            f.write('IOIPSL_LIBDIR="-L$LMDGCM/../../lib -L$LMDGCM/../IOIPSL/lib"\n')
            f.write('IOIPSL_LIB="-lioipsl"\n\n')
            f.write('XIOS_INCDIR="-I$LMDGCM/../XIOS/inc"\n')
            f.write('XIOS_LIBDIR="-L$LMDGCM/../XIOS/lib"\n')
            f.write('XIOS_LIB="-lxios -lstdc++"\n\n')
            f.write('ORCH_INCDIR="-I$LMDGCM/../../lib"\n')
            f.write('ORCH_LIBDIR="-L$LMDGCM/../../lib"\n\n')
            f.write('OASIS_INCDIR="-I$LMDGCM/../../oasis3-mct/BLD/build/lib/psmile.MPI1"\n')
            f.write('OASIS_LIBDIR="-L$LMDGCM/../../oasis3-mct/BLD/lib"\n')
            f.write('OASIS_LIB="-lpsmile.MPI1 -lscrip -lmct -lmpeu"\n\n')
            f.write('INCA_INCDIR="-I$LMDGCM/../INCA/build/inc"\n')
            f.write('INCA_LIBDIR="-L$LMDGCM/../INCA/build/lib"\n')
            f.write('INCA_LIB="-lchimie"\n')

        subprocess.run([
            'wget', 'https://lmdz.lmd.jussieu.fr/pub/install_lmdz.sh',
            '-O', 'install_lmdz.sh'], cwd=packdir, check=True)
        subprocess.run([
            'bash', 'install_lmdz.sh',
            '-v', self.version,
            *(self.install_arg.split() if self.install_arg else []),
            '-bench', '0',
            '-rad', self.rad,
            '-name', 'LMDZ',
            '-arch_dir', packdir,
            '-arch', 'mylocal'], cwd=packdir, check=True)

        subprocess.run([
            'wget', 'https://lmdz.lmd.jussieu.fr/pub/1D/1D.tar.gz'], cwd=self.lmdzdir, check=True)
        subprocess.run(['tar', 'xf', '1D.tar.gz'], cwd=self.lmdzdir, check=True)

    def pack_update(self):
        """Copy and prepare source files with *prep_code*."""
        packdir = os.path.join(self.LMDZPACK, self.name)

        if self.link:
            if self.packupdate:
                raise CheckCommitError(
                    "link option not compatible with the update option", 10)
            target = os.path.join(self.lmdzdir, 'modipsl', 'modeles', 'LMDZ',
                                  'libf', 'phylmd')
            common_glob = glob_mod.glob(os.path.expanduser('~/PHYEX/src/common/*/*'))
            for src in common_glob:
                os.symlink(src, os.path.join(target, os.path.basename(src)))
            lmdz_glob = glob_mod.glob(os.path.expanduser('~/PHYEX/src/lmdz/*/*'))
            for src in lmdz_glob:
                dst = os.path.join(target, os.path.basename(src))
                if os.path.lexists(dst):
                    os.unlink(dst)
                os.symlink(src, dst)
        else:
            if not os.path.isdir(packdir):
                raise CheckCommitError(
                    f"Pack directory doesn't exist ({packdir})", 9)

            if self.packupdate:
                shutil.move(os.path.join(packdir, 'PHYEX'),
                            os.path.join(packdir, 'PHYEXori'))

            subs = ['turb', 'micro', 'aux', 'ext']

            print(f"Copy {self.commit}")
            os.makedirs(os.path.join(packdir, 'PHYEX'), exist_ok=True)
            subprocess.run(['scp', '-q', '-r', os.path.join(self.commit, 'src'),
                            os.path.join(packdir, 'PHYEX')], check=True)

            prep_kwargs = self._parse_prep_code_opts(self.prepCodeOpts)
            prep_code(
                directory=os.path.join(packdir, 'PHYEX'),
                model='lmdz',
                mnh_expand=self.useexpand,
                subs=subs,
                no_raise_on_coding_norms=True,
                pyfortool_options=['--shumanFUNCtoCALL', '--removeACC'],
                **prep_kwargs)

            if self.packupdate:
                #Update only modified files
                phyex_dir = os.path.join(packdir, 'PHYEX')
                for subdir in subs:
                    sdir = os.path.join(phyex_dir, subdir)
                    if not os.path.isdir(sdir):
                        continue
                    for fname in os.listdir(sdir):
                        fpath = os.path.join(sdir, fname)
                        if os.path.isfile(fpath):
                            mvdiff(fpath, os.path.join(packdir, 'PHYEXori', subdir, fname))
                shutil.rmtree(os.path.join(packdir, 'PHYEX'))
                shutil.move(os.path.join(packdir, 'PHYEXori'),
                            os.path.join(packdir, 'PHYEX'))
            else:
                #Put PHYEX source code in the LMDZ source tree
                target = os.path.join(self.lmdzdir, 'modipsl', 'modeles', 'LMDZ',
                                      'libf', 'phylmd')
                for subdir in os.listdir(self.phyexdir):
                    sub_full = os.path.join(self.phyexdir, subdir)
                    if os.path.isdir(sub_full):
                        for f in os.listdir(sub_full):
                            src = os.path.join(self.phyexdir, subdir, f)
                            dst = os.path.join(target, f)
                            if os.path.isfile(src):
                                if os.path.exists(dst):
                                    os.unlink(dst)
                                os.symlink(src, dst)

        target = os.path.join(self.lmdzdir, 'modipsl', 'modeles', 'LMDZ', 'libf', 'phylmd')
        if not self.packupdate:
            shutil.copytree(
                target,
                os.path.join(os.path.dirname(target), 'phylmdorig'),
                dirs_exist_ok=True, symlinks=True)

        if self.buildSys == 'make':
            if os.path.isfile(os.path.join(target, 'modd_dimphyexn.F90')):
                mvdiff(os.path.join(target, 'modd_dimphyexn.F90'),
                       os.path.join(target, 'modd_dimphyex.F90'))
            result = subprocess.run(
                ['grep', '-i', 'END MODULE'] + [f for f in os.listdir(target)
                                                 if f.endswith('n.F90')],
                capture_output=True, text=True, cwd=target, check=True)
            for line in result.stdout.split('\n'):
                parts = line.split(':')
                if len(parts) >= 1:
                    name = parts[0].replace('n.F90', '')
                    if name:
                        src = os.path.join(target, f"{name}n.F90")
                        dst = os.path.join(target, f"{name}_n.F90")
                        if os.path.isfile(src):
                            mvdiff(src, dst)
            for src, dst in [('hypgeo.F90', 'modi_hypgeo.F90'),
                             ('hypser.F90', 'modi_hypser.F90'),
                             ('momg.F90', 'modi_momg.F90'),
                             ('tools.F90', 'mode_tools.F90')]:
                src_full = os.path.join(target, src)
                dst_full = os.path.join(target, dst)
                if os.path.isfile(src_full):
                    mvdiff(src_full, dst_full)
            for src, dst in [('shuman_mf.F90', 'modi_shuman_mf.F90'),
                             ('shuman_phy.F90', 'mode_shuman_phy.F90')]:
                src_full = os.path.join(target, src)
                dst_full = os.path.join(target, dst)
                if os.path.isfile(src_full):
                    mvdiff(src_full, dst_full)

        if not self.packupdate and self.rad != 'ecrad':
            for pattern in ('ecrad/yom*', 'ecrad/abor1.F90', 'ecrad/abor1.intfb.h',
                            'ecrad/parkind1.F90',
                            'ecrad/*/yom*', 'ecrad/*/abor1.F90',
                            'ecrad/*/abor1.intfb.h', 'ecrad/*/parkind1.F90'):
                for f in glob_mod.glob(os.path.join(target, pattern)):
                    basename = os.path.basename(f)
                    bpath = os.path.join(target, basename)
                    if not os.path.isfile(bpath):
                        os.symlink(f, bpath)

    def compilation_step(self):
        """Compile LMDZ with the chosen build system."""
        packdir = os.path.join(self.LMDZPACK, self.name)
        if not self.onlyIfNeeded or not os.path.isfile(
                os.path.join(packdir, 'bin', 'MASTERODB')):
            print(f"### Compilation of commit {self.commit}")
            bin_dir = os.path.join(self.lmdzdir, '1D', 'bin')

            exec_names = [os.path.join(bin_dir, 'lmdz1d.e'),
                          os.path.join(self.lmdzdir, 'modipsl', 'modeles', 'LMDZ', 'bin',
                                       f'lmdz1d_{self.L}_phylmd_{self.rad}_seq.e')]
            for ex in exec_names:
                if os.path.isfile(ex):
                    os.remove(ex)

            if self.buildSys == 'fcm':
                compile_file = os.path.join(bin_dir, 'compile')
                with open(compile_file, encoding='utf-8') as f:
                    content = f.read()
                content = content.replace('fcm=0', 'fcm=1')
                with open(compile_file, 'w', encoding='utf-8') as f:
                    f.write(content)

            compile_cmd = self.compilecmd.split()
            with open(os.path.join(self.lmdzdir, 'compilation.log'), 'w', encoding='utf-8') as out:
                subprocess.run(compile_cmd, cwd=bin_dir, stdout=out, stderr=subprocess.STDOUT,
                               check=False)

            if self.buildSys == 'fcm':
                print("Using fcm, compilation exits with error even if everything is OK")

    def execution(self):
        """Run RICO and ARM-CU test cases."""
        print(f"### Execution of commit {self.commit}")
        input_dir = os.path.join(self.lmdzdir, '1D', 'INPUT', 'PHYS')

        with open(os.path.join(input_dir, 'physiq.def_6A'), encoding='utf-8') as f:
            phys_def_6a = f.read()
        with open(os.path.join(input_dir, 'physiq.def_PHYLMD'), 'w', encoding='utf-8') as f:
            f.write('iflag_physiq=1\n' + phys_def_6a)
        with open(os.path.join(input_dir, 'physiq.def_PHYEX'), 'w', encoding='utf-8') as f:
            f.write('iflag_physiq=2\n' + phys_def_6a)

        for t in self.tests:
            for DEF in ('PHYEX', 'PHYLMD'):
                d = os.path.join(self.lmdzdir, '1D', 'EXEC', f'{DEF}L{self.L}', t)
                os.makedirs(d, exist_ok=True)

                oldcases = os.path.join(self.lmdzdir, '1D', 'OLDCASES', t)
                for item in os.listdir(oldcases):
                    src = os.path.join(oldcases, item)
                    dst = os.path.join(d, item)
                    if os.path.isdir(src):
                        os.symlink(src, dst)
                    else:
                        os.symlink(src, dst)
                for item in os.listdir(os.path.join(self.lmdzdir, '1D', 'INPUT', 'DEF')):
                    src = os.path.join(self.lmdzdir, '1D', 'INPUT', 'DEF', item)
                    dst = os.path.join(d, item)
                    if os.path.isfile(src) and not os.path.isfile(dst):
                        os.symlink(src, dst)

                shutil.copy2(os.path.join(self.lmdzdir, '1D', 'INPUT', 'PHYS',
                                          f'physiq.def_{DEF}'), os.path.join(d, 'physiq.def'))

                if self.rad == 'oldrad':
                    subprocess.run(
                        ['sed', '-i', '-e', 's/iflag_rrtm=.*$/iflag_rrtm=0/',
                         '-e', 's/NSW=.*$/NSW=2/', 'physiq.def'],
                        cwd=d, check=False)

                if self.rad == 'ecrad':
                    shutil.copy2(os.path.join(self.lmdzdir, 'modipsl', 'modeles', 'LMDZ',
                                              'DefLists', 'namelist_ecrad'), d)
                    data_dir = os.path.join(self.lmdzdir, 'modipsl', 'modeles', 'LMDZ',
                                            'libf', 'phylmd', 'ecrad', 'data')
                    if os.path.isdir(data_dir):
                        shutil.copytree(data_dir, os.path.join(d, 'data'),
                                        dirs_exist_ok=True, symlinks=True)
                    subprocess.run(
                        ['sed', '-i', '-e', 's@iflag_rrtm=1@iflag_rrtm=2@', 'physiq.def'],
                        cwd=d, check=True)

                vert_dir = os.path.join(self.lmdzdir, '1D', 'INPUT', 'VERT', f'L{self.L}')
                for item in os.listdir(vert_dir):
                    src = os.path.join(vert_dir, item)
                    dst = os.path.join(d, item)
                    if os.path.isfile(src) and not os.path.isfile(dst):
                        os.symlink(src, dst)

                os.symlink(f'L{self.L}.def', os.path.join(d, 'vert.def'))

                for ext in ('.dat', '.data', '.def'):
                    for f in os.listdir(oldcases):
                        if f.endswith(ext):
                            dst = os.path.join(d, f)
                            if os.path.isfile(dst) or os.path.islink(dst):
                                os.unlink(dst)
                            shutil.copy2(os.path.join(oldcases, f), d)

                compile_script = os.path.join(d, 'compile.sh')
                with open(compile_script, 'w', encoding='utf-8') as fh:
                    fh.write(f"cd {self.lmdzdir}/1D/bin\n{self.compilecmd}\n")
                os.chmod(compile_script, 0o755)

                if self.buildSys == 'make':
                    os.symlink(os.path.join(self.lmdzdir, '1D', 'bin', 'lmdz1d.e'),
                               os.path.join(d, 'lmdz1d.e'))
                else:
                    exec_name = os.path.join(
                        self.lmdzdir, 'modipsl', 'modeles', 'LMDZ', 'bin',
                        f'lmdz1d_{self.L}_phylmd_{self.rad}_seq.e')
                    if os.path.isfile(exec_name):
                        os.symlink(exec_name, os.path.join(d, 'lmdz1d.e'))
                    else:
                        alt = os.path.join(
                            self.lmdzdir, 'modipsl', 'modeles', 'LMDZ',
                            f'lmdz1d_{self.rad}_L{self.L}.e')
                        if os.path.isfile(alt):
                            os.symlink(alt, os.path.join(d, 'lmdz1d.e'))

                if DEF == 'PHYEX':
                    gcm_file = os.path.join(d, 'gcm1d.def')
                    if os.path.isfile(gcm_file):
                        with open(gcm_file, encoding='utf-8') as f:
                            content = f.read()
                        content = content.replace('day_step=144$', 'day_step=1440')
                        with open(gcm_file, 'w', encoding='utf-8') as f:
                            f.write(content)

                t1_us = time.time_ns() // 1000
                with open(os.path.join(d, 'execution.log'), 'w', encoding='utf-8') as out:
                    subprocess.run(['./lmdz1d.e'], cwd=d, stdout=out, stderr=subprocess.STDOUT,
                                   check=False)
                t2_us = time.time_ns() // 1000
                elapsed_ms = (t2_us - t1_us) // 1000

                if self.perffile:
                    with open(self.perffile, 'a', encoding='utf-8') as pf:
                        pf.write(f"{self.commit} lmdz {t} {elapsed_ms}\n")

    def comparison(self):
        """Comparison is not yet implemented for LMDZ."""
        raise CheckCommitError(
            "Check functionality is not yet implemented for LMDZ because:\n"
            "  1) the PHYEX interface will evolve and bit-reproducibility is not guaranteed\n"
            "  2) there is no surface scheme plugged in LMDZ-PHYEX, a recompilation\n"
            "     must be done to change the surface fluxes", 6)

    def cleaning(self):
        """Remove the LMDZ pack directory."""
        packdir = os.path.join(self.LMDZPACK, self.name)
        print(f"### Remove model directory for commit {self.commit}")
        if os.path.isdir(packdir):
            shutil.rmtree(packdir)


def main():
    """Entry point for the ``phyex-check_lmdz`` command."""
    run_tool(CheckCommitLmdz, module_name='pyphyextools.testing.check_lmdz')


if __name__ == '__main__':
    main()
