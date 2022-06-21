#!/usr/bin/env python3

import os
import glob
import logging

def verify_mnh_expand(path):
    """
    Verifies if source files are ready for expansion through mnh_expand
    :param path: directory to recursively check or file name

    Presently the folowing tests are performed:
      - starting and closing directives are conform
      - each instruction in the bloc is an effectation to an array with the right number of dimensions

    Limitation:
      - if the '=' sign is not on same line than the left hand side of the affectation instruction,
        the instruction will no be checked
      - one-line version of IF or WHERE statement are really on one line (no continuation line)
    """
    lhschar = b'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890_(:)%, \t'

    if os.path.isdir(path):
        logging.info("Enters directory: " + path)
        for filename in glob.glob(os.path.join(path, '*')):
            verify_mnh_expand(filename)
    else:
        logging.debug("Checks filename: " + path)
        with open(path, 'rb') as f: #read as byte because some files contain non UTF-8 characters
            lines = f.readlines()
        inside = False
        for iline, line in enumerate(lines):
            line = line.strip()
            if line[:13] == b'!$mnh_expand_':
                #New mnh_expand bloc
                logging.debug('Opening directive found. Line {line} of file {filename}'.format(line=iline + 1, filename=path))
                if inside:
                    logging.error('New mnh_expand bloc detected whereas we are already in a bloc. ' +
                                  'Line {line} of file {filename}'.format(line=iline + 1, filename=path))
                inside = True
                open_directive = line[13:].split(b'(')[0]
                open_args = line[13:].split(b'(')[1].split(b')')[0].replace(b' ', b'')
                dim = len(line.split(b'(')[1].split(b','))
                if line[-1:] != b')':
                    logging.error('Open directive must end with a closing bracket. ' +
                                  'Line {line} of file {filename}'.format(line=iline + 1, filename=path))
            elif line[:17] == b'!$mnh_end_expand_':
                #End of a mnh_expand bloc
                logging.debug('Closing directive found. Line {line} of file {filename}'.format(line=iline + 1, filename=path))
                if not inside:
                    logging.error('End of a mnh_expand bloc detected whereas we are not in a bloc. ' +
                                  'Line {line} of file {filename}'.format(line=iline + 1, filename=path))
                else:
                    inside = False
                    end_directive = line[17:].split(b'(')[0]
                    end_args = line[17:].split(b'(')[1].split(b')')[0].replace(b' ', b'')
                    if end_directive != open_directive:
                        logging.error('The end directive ({enddirect}) is not consistent with the opening directive ({opendirect}). '
                                      'Line {line} of file {filename}'.format(enddirect=end_directive.decode('UTF-8'),
                                                                              opendirect=open_directive.decode('UTF-8'),
                                                                              line=iline + 1, filename=path))
                    if end_args.upper() != open_args.upper():
                        logging.error('The end args ({endargs}) are not consistent with the opening args ({openargs}). '
                                      'Line {line} of file {filename}'.format(endargs=end_args.decode('UTF-8'),
                                                                              openargs=open_args.decode('UTF-8'),
                                                                              line=iline + 1, filename=path))
                    if line[-1:] != b')':
                        logging.error('Closing directive must end with a closing bracket. ' +
                                      'Line {line} of file {filename}'.format(line=iline + 1, filename=path))
            elif inside:
                #We do not want to implement a full fortran parser, we are only interested in the left hand side of
                #affectation instructions. If left hand side is correct (an array element) the right hand side
                #will be necessarily correct (otherwise a compilation error will be thrown).

                #Suppresion of condition in 'IF' and 'WHERE' one-line instructions
                if line[:3].upper() in(b'IF ', b'IF('):
                    line = line[line.index(b'(') + 1:]
                    nb = 1
                    while nb >= 1 and(len(line) > 0):
                        if line[:1] == b'(': nb += 1
                        elif line[:1] == b')': nb -= 1
                        line = line[1:].strip()
                    if line.upper()[:5] in (b'THEN ', b'THEN!'): line = line[5:]
                elif line[:6].upper() in(b'WHERE ', b'WHERE('):
                    if open_directive != b'where':
                        logging.error('There is a WHERE statement in a mnh_expand array bloc. '
                                      'Line {line} of file {filename}'.format(line=iline + 1, filename=path))
                    line = line[line.index(b'(') + 1:]
                    nb = 1
                    while nb >= 1 and(len(line)>0):
                        if line[:1] == b'(': nb += 1
                        elif line[:1] == b')': nb -= 1
                        line = line[1:].strip()

                #Check if it is the left hand side of an affectation
                if line[:3].upper() == b'DO ':
                    logging.warning('A DO loop is inside a mnh_expand bloc, is order correct?. ' +
                                    'Line {line} of file {filename}'.format(line=iline + 1, filename=path))
                elif line[:5].upper() == b'CALL ':
                    logging.warning('A CALL statement is inside a mnh_expand bloc, is it correct? ' +
                                    'Line {line} of file {filename}'.format(line=iline + 1, filename=path))
                elif b'=' in line and all([c in lhschar for c in line.split(b'=')[0]]):
                    lhs = line.split(b'=')[0]
                    if not b'(' in lhs:
                        logging.error('Array on the left hand side of an effectation instruction must be written ' +
                                      'with opening and closing brackets. ' + 
                                      'Line {line} of file {filename}'.format(line=iline + 1, filename=path))
                    if lhs.count(b':') != dim:
                        logging.error('Array on the left hand side of an effectation instruction must have the same ' +
                                      'number of :-dimensions as the number defined in the directive. ' +
                                      'Line {line} of file {filename}'.format(line=iline + 1, filename=path))
                    if line[line.index(b'(')-1] in b' \t':
                        logging.error('There must be no space wetween the array name and the opening bracket in ' +
                                      'the affectation instruction. ' +
                                      'Line {line} of file {filename}'.format(line=iline + 1, filename=path))

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='mnh_expand checker')
    parser.add_argument("-v", "--verbose", dest="verbose", action="count", default=0,
                        help="Show warning (-v), info (-v -v) or debug (-v -v -v) messages")
    parser.add_argument('PATH', help="directory to recursively check, or filename")
    args = parser.parse_args()
    level = {0:'ERROR',
             1:'WARNING',
             2:'INFO',
             3:'DEBUG'}[args.verbose]
    logging.basicConfig(level=getattr(logging, level, None))
    verify_mnh_expand(args.PATH)
