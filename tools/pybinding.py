#!/usr/bin/env python3

"""
This script contains a uniq function that reads a FORTRAN file and tries
to create a new FORTRAN SUBROUTINE that can be used to access the
original routine from python through the ctypesForFortran utility.

This function is designed to work with PHYEX codes by recognizing
the name of certain structures. In addition, it performs special
processing for ini_phyex (adding a file name for the namelist and
flushing the listing file).
"""

import os
import re
from pyfortool import PYFT
from pyfortool.util import n2name

class PybindingError(Exception):
    """
    Exceptions for pybinding
    """

def pybinding(fortran_in, scope, fortran_out, python_out, libso,
              tpfileIsNam=False, Findexing=False):
    """
    :param fortran_in: FORTRAN input file
    :param scope: scope to deal with
    :param fortran_out: FORTRAN output file
    :param python_out: python output file
    :param libso: path to the shared lib (relative to the python_out file)
    :param tpfileIsNam: True if tpfile is a namelist
    :param Findexing: True to keep the same index order as in the FORTRAN subroutine
    """

    def findcomponents(bound):
        """
        :param bound: a string representing an array bound
        :returns: a list of components
        If bound is '5+(X+3)', result is ['5', 'X', '3']
        """
        if bound is None:
            return []
        elemlist = bound
        for op in ('+', '-', '*', '/', '(', ')'):
            elemlist = elemlist.replace(op, ' ')
        return elemlist.split(' ')

    #From the original subroutine, two FORTRAN subroutines will be created:
    # - the main one which receives the dummy arguments from python
    #   transforms them and call the original subroutine
    # - an helper subroutine to help dealing with kinds and some
    #   dimensions (read in module)
    #A python code to call the main FORTRAN subroutine is also created.

    #Open input FORTRAN file and get the right scope
    pftin = PYFT(fortran_in)
    scopeNode = pftin.getScopeNode(scope, excludeContains=True)
    kind, name = scope.split('/')[-1].split(':')
    name = name.upper()
    kind = {'sub': 'SUBROUTINE',
            'func': 'FUNCTION'}[kind]

    #Loop on all the dummy arguments and prepare their use
    argList1 = [] #argument list of the main public interface
    argList2 = [] #argument list of the existing FORTRAN routine
    moduleList = [] #declarative module to add to the main subroutine
    moduleListHelper = [] #declarative module to add to the helper subroutine
    declList = [] #variable declaration list (dummy var and local var)
    copyList = [] #Instructions to copy public var to local var
                  # to be used to call the original subroutine
    extra = [] #additional variables to read in module
    flush = [] #list of unit to flush
    docstringIN = ["Input arguments:"]
    docstringOUT = ["Output arguments:"]
    for N in scopeNode.findall('.//{*}dummy-arg-LT/{*}arg-N/{*}N'):
        var = scopeNode.varList.findVar(n2name(N), exactScope=True)
        if var is None:
            raise PybindingError(n2name(N) + ' seems to be undeclared')
        vartype = var['t'].replace(' ', '').upper()
        if vartype == 'TYPE(DIMPHYEX_T)':
            moduleList.append('USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t')
            declList.append('TYPE(DIMPHYEX_t) :: D')
            for v in ['NIT', 'NKT', 'NKL', 'NVEXT']:
                declList.append('INTEGER, INTENT(IN) :: ' + v)
                argList1.append((v, 'INTEGER', False, 'IN', None))
                docstringIN.append(f"    {v} (IN) to replace the FORTRAN D structure")
            copyList.append('D%NIT=NIT')
            copyList.append('D%NIB=1')
            copyList.append('D%NIE=NIT')
            copyList.append('D%NJT=1')
            copyList.append('D%NJB=1')
            copyList.append('D%NJE=1')
            copyList.append('D%NIJT=NIT')
            copyList.append('D%NIJB=1')
            copyList.append('D%NIJE=NIT')
            copyList.append('D%NKL=NKL')
            copyList.append('D%NKT=NKT')
            copyList.append('D%NKA=MERGE(1, NKT, NKL==1)')
            copyList.append('D%NKU=MERGE(NKT, 1, NKL==1)')
            copyList.append('D%NKB=D%NKA+NKL*NVEXT')
            copyList.append('D%NKE=D%NKU-NKL*NVEXT')
            copyList.append('D%NKTB=1+NVEXT')
            copyList.append('D%NKTE=NKT-NVEXT')
            copyList.append('D%NIBC=1')
            copyList.append('D%NJBC=1')
            copyList.append('D%NIEC=NIT')
            copyList.append('D%NJEC=1')
            argList2.append('D')
        elif vartype in ['TYPE(NSV_T)',]:
            moduleList.append('USE MODD_' + vartype[5:-3] + ', ONLY: T' + vartype[5:-3])
            declList.append('INTEGER, INTENT(IN) :: NSV')
            argList1.append(('NSV', 'INTEGER', False, 'IN', None))
            docstringIN.append(f"    NSV (IN) to replace the FORTRAN TNSV structure")
            copyList.append('TNSV%NSV=NSV')
            argList2.append('T'+vartype[5:-3])
        elif vartype in ['TYPE(CST_T)', 'TYPE(ELEC_PARAM_T)', 'TYPE(ELEC_DESCR_T)',
                         'TYPE(PARAM_LIMA_T)','TYPE(PARAM_LIMA_WARM_T)','TYPE(PARAM_LIMA_COLD_T)','TYPE(PARAM_LIMA_MIXED_T)',]:
            moduleList.append('USE MODD_' + vartype[5:-3] + ', ONLY: ' + vartype[5:-3])
            argList2.append(vartype[5:-3])
        elif vartype in ['TYPE(NEB_T)', 'TYPE(RAIN_ICE_PARAM_T)', 'TYPE(PARAM_ICE_T)',
                         'TYPE(RAIN_ICE_DESCR_T)', 'TYPE(TURB_T)', 'TYPE(PARAM_MFSHALL_T)']:
            moduleList.append('USE MODD_' + vartype[5:-3] + '_N, ONLY: ' + vartype[5:-3] + 'N')
            argList2.append(vartype[5:-3] + 'N')
        elif vartype == 'TYPE(CSTURB_T)':
            moduleList.append('USE MODD_CTURB, ONLY: CSTURB')
            argList2.append(vartype[5:-3])
        elif vartype == 'TYPE(TBUDGETCONF_T)':
            moduleList.append('USE MODD_BUDGET, ONLY: TBUCONF_ASSOCIATE, TBUCONF')
            argList2.append('TBUCONF')
            copyList.append('CALL TBUCONF_ASSOCIATE')
        elif vartype == 'TYPE(TBUDGETDATA_PTR)':
            moduleList.append('USE MODD_BUDGET, ONLY: NBUDGET_RH, TBUDGETDATA_PTR')
            moduleList.append('USE MODE_BUDGET_OFFLINE, ONLY: TBUDGETDATA_OFFLINE')
            declList.append('TYPE(TBUDGETDATA_PTR), DIMENSION(NBUDGET_RH) :: YLBUDGET')
            declList.append('TYPE(TBUDGETDATA_OFFLINE), DIMENSION(NBUDGET_RH), TARGET :: YLOFFBUD')
            declList.append('INTEGER :: JRR')
            copyList.append('DO JRR=1, NBUDGET_RH\n  YLBUDGET(JRR)%PTR => YLOFFBUD(JRR)\n  YLBUDGET(JRR)%PTR%NBUDGET=JRR\nENDDO')
            argList2.append('YLBUDGET')
        elif vartype == 'TYPE(TFILEDATA)':
            moduleList.append('USE MODD_IO, ONLY: TFILEDATA')
            declList.append('TYPE(TFILEDATA) :: TPFILE')
            declList.append('INTEGER, INTENT(IN) :: KUNIT')
            argList1.append(('KUNIT', 'INTEGER', False, 'IN', None))
            docstringIN.append("    KUNIT (IN) to replace the FORTRAN TFILEDATA structure")
            copyList.append('TPFILE%NLU=KUNIT')
            argList2.append('TPFILE')
            if var['n'] == 'TPFILE' and name == "INI_PHYEX":
                declList.append('CHARACTER(LEN=200) :: SLOCAL_NAMFILE')
                declList.append('INTEGER :: ILOOPVAR_NAMFILE')
                declList.append('CHARACTER, DIMENSION(200), INTENT(IN) :: NAMFILE')
                argList1.append(('NAMFILE', 'CHARACTER', False, 'IN', ((None, '200'),)))
                docstringIN.append("    NAMFILE (IN) to replace the FORTRAN TFILEDATA structure")
                copyList.append('DO ILOOPVAR_NAMFILE=1, 200\n' + \
                                '  SLOCAL_NAMFILE(ILOOPVAR_NAMFILE:ILOOPVAR_NAMFILE)=' + \
                                     'NAMFILE(ILOOPVAR_NAMFILE)\n' + \
                                'ENDDO')
                copyList.append('IF (KUNIT/=0) OPEN(KUNIT, FILE=SLOCAL_NAMFILE)')
        elif vartype == 'TYPE(TLES_T)':
            moduleList.append('USE MODD_LES, ONLY: TLES_t')
            declList.append('TYPE(TLES_t) :: TLES')
            copyList.append('TLES%LLES=.FALSE.')
            copyList.append('TLES%LLES_CALL=.FALSE.')
            argList2.append('TLES')
        elif vartype == 'TYPE(PHYEX_T)':
            pass
        elif vartype.startswith('TYPE('):
            raise NotImplementedError('Does not know how to deal with' +
                                      f'argument of type {vartype} ' +
                                      f'in {fortran_in}')
        elif var['n'] == 'KSPLITR':
            moduleList.append('USE MODD_CLOUDPAR_n, ONLY: KSPLITR => NSPLITR')
            argList2.append('KSPLITR')
        elif vartype.startswith('CHARACTER(LEN=') and vartype[14].split(')')[0] != '1':
            charlen = int(vartype[14:].split(')')[0])
            #We must use a one-character length arrays,
            #we add a dimension to replace the string length
            #This is the variable we will receive (given by ctypesforfortran)
            replvar = var.copy()
            replvar['as'] = [(None, str(charlen))] + var['as']
            replvar['t'] = 'CHARACTER(KIND=C_CHAR)'
            argList1.append((replvar['n'], replvar['t'],
                             replvar['opt'], replvar['i'].upper(), replvar['as']))
            docstr = ', OPTIONAL' if replvar['opt'] else ''
            docstr = f"    {var['n']} ({var['i']}{docstr})"
            if replvar['i'] in ('IN', 'INOUT'):
                docstringIN.append(docstr)
            if replvar['i'] in ('OUT', 'INOUT'):
                docstringOUT.append(docstr)
            declList.append(pftin.varSpec2stmt(replvar))
            #But the original routine wants the same characters but with a multi-character string
            #We create this local var, and pass it to the routine
            localvar = var.copy()
            localvar['n'] = 'SLOCAL_' + var['n']
            localvar['i'] = None
            localvar['arg'] = False
            declList.append(pftin.varSpec2stmt(localvar))
            argList2.append(localvar['n'])
            #And we need to copy the received characters into the local var
            if len(var['as']) == 0:
                declList.append(f"INTEGER :: ILOOPVAR_{var['n']}")
                copyList.append(('DO ILOOPVAR_{name}=1, {len}\n' + \
                                 '  SLOCAL_{name}(ILOOPVAR_{name}:ILOOPVAR_{name})=' + \
                                      '{name}(ILOOPVAR_{name})\n' + \
                                 'ENDDO').format(name=var['n'], len=str(charlen)))
            else:
                loopVars = [f"ILOOPVAR_{var['n']}_{i}" for i in range(len(var['as']) + 1)]
                declList.append('INTEGER :: ' + ', '.join(loopVars))
                copystring = ''
                for d in range(len(var['as'])):
                    copystring += ("DO ILOOPVAR_{varname}_{i}=LBOUND(SLOCAL_{varname}, {i}), " + \
                                                             "UBOUND(SLOCAL_{varname}, {i})\n"
                                  ).format(varname=var['n'], i=str(d + 1))
                copystring += ('DO ILOOPVAR_{name}_0=1, {len}\n' + \
                               '  SLOCAL_{name}({loopv})(ILOOPVAR_{name}_0:ILOOPVAR_{name}_0)=' + \
                                    '{name}(ILOOPVAR_{name}_0, {loopv})\n' + \
                               'ENDDO\n').format(name=var['n'], loopv=', '.join(loopVars[1:]),
                                                 len=str(charlen))
                copystring += 'ENDDO\n' * (len(loopVars) - 1)
                copyList.append(copystring)
        else:
            if var['n'] == 'KLUOUT' and name == "INI_PHYEX":
                flush.append('KLUOUT')
            replvar = var.copy()
            localvar = var.copy()
            if vartype.startswith('LOGICAL'):
                replvar['t'] = 'LOGICAL(KIND=C_BOOL)'
                localvar['n'] = 'SLOCAL_' + var['n']
                localvar['i'] = None
                localvar['arg'] = False
                localvar['opt'] = False
                # The folowing copy is done even if the argument is a missing optional one
                # It's not correct, but how can it be otherwise? Protecting the copy isn't
                # enough, we'd also have to call the routine without the argument
                copyList.append('SLOCAL_{varname}={varname}'.format(varname=var['n']))
            argList1.append((replvar['n'], replvar['t'].replace(' ', '').upper(),
                             replvar['opt'], replvar['i'].upper(), replvar['as']))
            docstr = ', OPTIONAL' if replvar['opt'] else ''
            docstr = f"    {var['n']} ({var['i']}{docstr})"
            if replvar['i'] in ('IN', 'INOUT'):
                docstringIN.append(docstr)
            if replvar['i'] in ('OUT', 'INOUT'):
                docstringOUT.append(docstr)
            argList2.append(localvar['n'])
            for d in replvar['as']: # The following code handles MERGE() in declaration
                if 'NIT' in d[1] or 'NIJT' in d[1]:
                    d[1] = 'NIT'
                elif 'NJT' in d[1]:
                    d[1] = '1'
                elif 'NKT' in d[1]:
                    d[1] = 'NKT'
                elif 'KSV' in d[1]:
                    d[1] = 'KSV'
                elif 'NSV' in d[1]:
                    d[1] = 'NSV'
            declList.append(pftin.varSpec2stmt(replvar))
            if vartype.startswith('LOGICAL'):
                declList.append(pftin.varSpec2stmt(localvar))
            if 'JPSVMAX' not in extra and \
               'JPSVMAX' in [elem for d in var['as'] for di in d for elem in findcomponents(di)]:
                extra.append('JPSVMAX')
                moduleList.append('USE MODD_PARAMETERS, ONLY: JPSVMAX')
                moduleListHelper.append('USE MODD_PARAMETERS, ONLY: JPSVMAX')
            #if 'NSV' not in extra and \
            #   'NSV' in [elem for d in var['as'] for di in d for elem in findcomponents(di)]:
            #    extra.append('NSV')
            #    moduleList.append('USE MODD_NSV, ONLY: NSV')
            #    moduleListHelper.append('USE MODD_NSV, ONLY: NSV')
            for vv in ('NSP', 'NCARB', 'NSOA'):
                if vv not in extra and \
                   vv in [elem for d in var['as'] for di in d for elem in findcomponents(di)]:
                    extra.append(vv)
                    moduleList.append(f'USE MODD_CH_AEROSOL, ONLY: {vv}')
                    moduleListHelper.append(f'USE MODD_CH_AEROSOL, ONLY: {vv}')

    pftin.close()

    #Write output FORTRAN file
    exists = os.path.exists(fortran_out)
    with open(fortran_out, 'a', encoding="utf-8") as f:
        if exists:
            f.write('\n')

        #Main subroutine
        argList1name = [v[0] for v in argList1]
        f.write(('{kind} PY{name}({argList}) ' + \
                 'BIND(C, name="PY{name}")\n').format(kind=kind,
                                                      argList=', '.join(argList1name),
                                                      name=name))
        if scope.split('/')[0].split(':')[0] == 'module':
            f.write(f"\nUSE {scope.split('/')[0].split(':')[1]}\n")
        else:
            f.write(f'\nUSE MODI_{name}\n')
        f.write('\n'.join(moduleList) + '\n')
        f.write('\nUSE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR, C_BOOL\n')
        f.write('\nIMPLICIT NONE\n')
        f.write('\n'.join(declList) + '\n')
        f.write('\n'.join(copyList) + '\n')
        f.write(f"CALL {name}({', '.join(argList2)})\n")
        for unit in flush:
            f.write(f"FLUSH({unit})\n")
        f.write(f'END {kind} PY{name}')

        #Helper subroutine
        f.write('\n')
        extraname = ['K_' + e for e in extra]
        f.write(('{kind} HPY{name}({argList}) ' + \
                 'BIND(C, name="HPY{name}")\n').format(kind=kind,
                                                       argList=', '.join(argList1name + extraname),
                                                       name=name))
        f.write('\n'.join(moduleListHelper) + '\n')
        f.write('\nUSE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR, C_BOOL\n')
        f.write('INTEGER(KIND=4), INTENT(OUT) :: ' + (', '.join(argList1name + extraname)) + '\n')
        for varname, vartype, _, _, _ in argList1:
            f.write(vartype + ' :: ' + 'K_' + varname + '\n')
        for varname, vartype, _, _, _ in argList1:
            if vartype.startswith('CHARACTER'):
                f.write(f'{varname}=0\n')
            else:
                f.write('{varname}=KIND(K_{varname})\n'.format(varname=varname))
        for varname in extra:
            f.write(f'K_{varname}={varname}\n')
        f.write(f'END {kind} HPY{name}')

    #Write output pyton file
    exists = os.path.exists(python_out)
    with open(python_out, 'a', encoding="utf-8") as f:
        if not exists:
            f.write('#!/usr/bin/env python3\n')
            f.write('import numpy\n')
            f.write('import ctypesForFortran\n')
            f.write('from ctypesForFortran import string2array, array2string\n')
            f.write('import sys\n')
            f.write('from functools import lru_cache\n')
            f.write('\n')
            f.write('IN = ctypesForFortran.IN\n')
            f.write('OUT = ctypesForFortran.OUT\n')
            f.write('INOUT = ctypesForFortran.INOUT\n')
            f.write('MISSING = ctypesForFortran.MISSING\n')
            f.write('MAO = ctypesForFortran.MANDATORY_AFTER_OPTIONAL\n')
            f.write('\n')
            libso = os.path.join(os.path.dirname(os.path.realpath(python_out)), libso)
            f.write(f'ctypesFF, handle = ctypesForFortran.ctypesForFortranFactory("{libso}")\n')
            f.write('\n')
            f.write('def close():\n')
            f.write('    ctypesForFortran.dlclose(handle)\n')
        f.write('\n@lru_cache\n')
        f.write("@ctypesFF(prefix='', suffix='')\n")
        f.write(f"def HPY{name}():\n")
        f.write(f"    return ([], [(numpy.int32, None, OUT)] * {len(argList1 + extra)}, None)\n")
        f.write('\n')
        f.write("@ctypesFF(prefix='', suffix='', castInput=True, " + \
                f"indexing='{'F' if   Findexing else 'C'}')\n")
        argList = []
        missing = False
        for varname, _, opt, intent, _ in argList1:
            if intent in ('IN', 'INOUT'):
                if opt is True:
                    argList.append(varname + '=MISSING')
                    missing = True
                elif missing:
                    argList.append(varname + '=MAO')
                else:
                    argList.append(varname)
        outvar = [varname for (varname, _, _, intent, _) in argList1 if intent in ('OUT', 'INOUT')]
        argList.append('missingOUT=None')
        f.write(f"def PY{name}({', '.join(argList)}):\n")
        f.write( '    """\n')
        f.write( "    Wrapping to provide a python binding using pybinding " + \
                 "from the PHYEX package.\n")
        f.write(f"    Original file is {fortran_in}\n\n")
        f.write( "    Special notes for OPTIONAL variables:\n")
        f.write( "     - To be absent on input, the variable must receive the MISSING\n")
        f.write( "       special value.\n")
        f.write( "     - For a variable with OUT intent, the variable name must be added in the\n")
        f.write( "       missingOUT list to not be returned by the python function.\n")
        f.write( "     - For a variable with INOUT intent, its presence or absence is\n")
        f.write( "       simultaneous at input and output\n\n")
        f.write( "    " + "\n    ".join(docstringIN) + "\n")
        f.write( "    " + "\n    ".join(docstringOUT) + "\n")
        f.write( '    """\n')
        f.write( "    missingOUT = [] if missingOUT is None else missingOUT\n")
        f.write(f"    kinds = HPY{name}()\n")
        if len(extra) > 0:
            f.write(f"    {', '.join(['_' + e for e in extra])}, = kinds[{len(argList1)}:]\n")
        inlist = []
        alllist = []
        for ivar, (varname, vartype, opt, intent, dim) in enumerate(argList1):
            if vartype.startswith('INTEGER'):
                npkind = {4: 'numpy.int32', 8: 'numpy.int64'}
            elif vartype.startswith('REAL'):
                npkind = {4: 'numpy.float32', 8: 'numpy.float64'}
            elif vartype.startswith('LOGICAL(KIND=C_BOOL)'):
                npkind = {1: 'bool'}
            elif vartype.startswith('CHARACTER'):
                npkind = {0: 'str'}
            else:
                raise NotImplementedError('Argument type not yet implemented: ' + vartype)
            if vartype.startswith('CHARACTER') and intent in ('OUT', 'INOUT'):
                # In this case, the array2string decorator must me added with
                # the indexes of output strings
                raise NotImplementedError('Character strings in OUT/INOUT')
            if intent in ('IN', 'INOUT'):
                if vartype.startswith('CHARACTER'):
                    inlist.append(f'string2array({varname}, {dim[0][1]}), #{varname}')
                else:
                    inlist.append(f'{varname}, #{varname}')
            strdim = []
            if not (dim is None or len(dim) == 0):
                for d in dim:
                    if d[0] is not None:
                        raise NotImplementedError('Array with non-zero indexing not implemented')
                    strdimpart = d[1]
                    for elem in findcomponents(d[1]):
                        if elem in ('NIT', 'NIJT'):
                            strdimpart = re.sub(rf'\b{elem}\b', 'NIT', strdimpart)
                        elif elem == 'NJT':
                            strdimpart = re.sub(rf'\b{elem}\b', '1', strdimpart)
                        elif elem == 'NKT':
                            pass
                        elif elem in argList1name:
                            pass  # This dimension is a dummy argument
                        elif elem in ('JPSVMAX', 'NSV', 'NSP', 'NCARB', 'NSOA'):
                            strdimpart = re.sub(rf'\b{elem}\b', f'_{elem}', strdimpart)
                        else:
                            try:
                                int(elem)
                            except ValueError as verr:
                                raise NotImplementedError(f'Dimension not implemeneted ({elem}) ' + \
                                                          f'for variable {varname}') from verr
                    strdim.append(strdimpart)
            if len(strdim) == 0:
                strdim = 'None'
            else:
                if Findexing:
                    strdim = '(' + ', '.join(strdim) + ',)'
                else:
                    if vartype.startswith('CHARACTER'):
                        strdim = '(' + ', '.join(['1'] + strdim[1:][::-1] + [strdim[0]]) + ',)'
                    else:
                        strdim = '(' + ', '.join(strdim[::-1]) + ',)'
            kind = repr(npkind).replace("'", "") + f'[kinds[{ivar}]]'
            if opt and intent == 'OUT':
                kind = f'MISSING if "{varname.upper()}" in missingOUT else {kind}'
            alllist.append(f'({kind}, {strdim}, {intent}), #{varname}')
        f.write(("    return ([{inlist}\n" + \
                 "            ],\n" + \
                 "            [{alllist}\n" + \
                 "            ],\n" + \
                 "            None)\n").format(inlist='\n             '.join(inlist),
                                               alllist='\n             '.join(alllist)))

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Python binding maker')
    parser.add_argument('INPUT', type=str,
                        help='FORTRAN input file')
    parser.add_argument('SCOPE', type=str,
                          help='scope to deal with')
    parser.add_argument('FORTRAN', type=str,
                        help='FORTRAN output file')
    parser.add_argument('PYTHON', type=str,
                        help='PYTHON output file')
    parser.add_argument('LIBSO', type=str,
                        help='Shared lib path')
    parser.add_argument('--Findexing', default=False, action='store_true',
                        help='If set, array indexes are in the same order as the FORTRAN routine')
    args = parser.parse_args()
    pybinding(args.INPUT, args.SCOPE, args.FORTRAN, args.PYTHON, args.LIBSO, args.Findexing)
