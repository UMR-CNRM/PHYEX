#!/usr/bin/env python3

import sys
import re
from pathlib import Path

reg_proc_name = re.compile(r"^\s*(SUBROUTINE|FUNCTION)\s+(\w+)", re.IGNORECASE)
reg_proc_args = re.compile(r"\(([\w\s,]+)\)", re.IGNORECASE)
reg_func_result = re.compile(r"RESULT\s*\(([\w\s]+)\)", re.IGNORECASE)
reg_use_module_only = re.compile(r"(?:^\s*use\s+\w+\s*,\s*only\s*:)(.*)", re.IGNORECASE)
reg_use_module_name = re.compile(r"^\s*use\s+(\w+)", re.IGNORECASE)
reg_variable_name_simple = re.compile(r"::\s*(\w+)", re.IGNORECASE)
reg_string = re.compile(r"'[^']*'|\"[^\"]*\"")
reg_component = re.compile(r"%\s*\w+")  # derived type components (D%NIJT)
reg_identifier = re.compile(r"[A-Z_][A-Z0-9_]*")

# Macro directives and Fortran include statements are case-sensitive so they must never be uppercased
reg_case_is_important = re.compile(r"\s*#|\s*INCLUDE", re.IGNORECASE)

reg_all_decl_variable = [
    r"^\s*INTEGER\s*[\(*,:\s]",
    r"^\s*REAL\s*[\(*,:\s]",
    r"^\s*DOUBLE PRECISION\s*[\(*,:\s]",
    r"^\s*CHARACTER\s*[\(*,:\s]",
    r"^\s*LOGICAL\s*[\(*,:\s]",
    r"^\s*TYPE\s*\(",
    r"^\s*CLASS\s*\(",
]
reg_all_decl = [
    r"^\s*$",
    r"^\s*!",
    r"^\s*SUBROUTINE ",
    r"^\s*FUNCTION ",
    r"^\s* ",
    r"^\s*USE[,\s]",
    r"^\s*#",
    r"^\s*IMPLICIT NONE",
    r"^\s*include ",
    *reg_all_decl_variable,
]

reg_is_still_decl = re.compile("|".join(x for x in reg_all_decl), re.IGNORECASE)
reg_is_variable_decl = re.compile(
    "|".join(x for x in reg_all_decl_variable), re.IGNORECASE
)


def is_modi_needed(lines):
    # If a module is declared in the file, then no interface is needed
    for line in lines:
        if line.startswith("MODULE "):
            return False
    return True


def simplify_code(src):
    # Remove comments, remove double spaces, join lines using continuation marker (&)
    lines = []
    buf = ""
    for line in src:
        line = line.strip()
        if line.startswith("!"):
            continue
        while "  " in line:
            line = line.replace("  ", " ")
        if not reg_case_is_important.match(line):
            line = line.upper()
        if "!" in line and not line.startswith("#"):
            # When '!' is in the middle of a line, most of the time what follow it is comment.
            # But sometimes, it is a negation in an ifdef, see: #if ! defined foo
            line = line.split("!")[0].strip()
        if line.endswith("&"):
            buf += line.replace("&", "")
        else:
            lines.append(buf + line.replace("&", "") + "\n")
            buf = ""
    return lines


def remove_contained_procedures(lines):
    to_rm = []
    depth = 0

    for idx, line in enumerate(lines):
        if line.startswith("SUBROUTINE ") or line.startswith("FUNCTION "):
            depth = depth + 1
        if depth > 1:
            to_rm.append(idx)
        if line.startswith("END SUBROUTINE") or line.startswith("END FUNCTION"):
            depth = depth - 1
    for idx in reversed(to_rm):
        lines.pop(idx)


def remove_dangling_ifdef(lines):
    ifdef = []
    for idx, line in enumerate(lines):
        if line.startswith("#if"):
            ifdef.append(idx)
        if line.startswith("#endif"):
            ifdef.pop()
    for idx in reversed(ifdef):
        lines.pop(idx)


def remove_unsed_modules(lines):
    use_only = {}  # line_idx -> "use only" values
    # list of lines with used module without an only clause,
    # or imported value that are not used in the procedure signature
    to_rm = []

    # Retrieve all explicitly imported values (eg. use foo, only: stuff)
    # If there is no 'only' clause, the line is mark for deletion
    for idx, line in enumerate(lines):
        if line.startswith("USE "):
            matches = reg_use_module_only.match(line)
            if not matches:
                to_rm.append(idx)
            else:
                use_only[idx] = [v.strip() for v in matches.group(1).split(",") if v.strip()]

    # Identifiers used in the variable declarations (whole words only, string
    # constants and derived type components removed: 'TM06' does not use TM06,
    # LRMC01 does not use RMC01)
    used = set()
    for line in lines:
        if reg_is_variable_decl.match(line):
            line = reg_component.sub("", reg_string.sub("", line))
            used.update(reg_identifier.findall(line))

    # Only the imported values used in a variable declaration are kept, a line
    # without any of them is removed
    for idx, values in use_only.items():
        kept = []
        for val in values:
            # renaming (local => used): the local name is the one used
            local = val.split("=>")[0].strip()
            # generic operators and assignments are always kept
            if local.startswith(("OPERATOR", "ASSIGNMENT")) or local in used:
                kept.append(val)
        if kept:
            module = reg_use_module_name.match(lines[idx]).group(1)
            lines[idx] = f"USE {module}, ONLY: {', '.join(kept)}\n"
        else:
            to_rm.append(idx)

    for idx in reversed(sorted(to_rm)):
        lines.pop(idx)


def get_proc_signature(line):
    # Retrieve type of procedure, name, arguments (with result if any)
    matches = reg_proc_name.match(line)
    proc_type = matches.group(1)
    proc_name = matches.group(2)

    proc_args = []
    matches = reg_proc_args.search(line)
    if matches:
        proc_args = [arg.strip() for arg in matches.group(1).split(",")]

    if proc_type == "FUNCTION":
        matches = reg_func_result.search(line)
        if matches:
            proc_args.append(matches.group(1).strip())

    return proc_type, proc_name, proc_args


def break_long_lines(lines):
    newlines = []
    for line in lines:
        while len(line) > 80:
            break_idx_space = line[0:80].rfind(" ")
            break_idx_comma = line[0:80].rfind(",")

            # Find which one let us have the longest line of less than 80 chars
            if break_idx_comma > break_idx_space:
                break_idx = break_idx_comma
            else:
                break_idx = break_idx_space

            # If no space or comma was found, force break the line at 80 chars
            if break_idx == -1:
                break_idx = 80

            newlines.append(line[:break_idx] + "&\n")
            line = line[break_idx:]
        newlines.append(line)
    return newlines


def gen_modi(f90):
    modi_path = Path(modi_dir / Path("modi_" + str(f90.name)))
    src = open(f90, "r")
    lines = simplify_code(src)
    if not is_modi_needed(lines):
        modi_path.touch()
        return

    remove_contained_procedures(lines)

    modi_lines = []
    proc_type = ""
    proc_name = ""
    proc_args = []
    for line in lines:
        if not proc_name:
            if line.startswith("SUBROUTINE ") or line.startswith("FUNCTION "):
                proc_type, proc_name, proc_args = get_proc_signature(line)
                modi_lines.append(f"MODULE MODI_{proc_name}\nINTERFACE\n")

        matches = reg_is_still_decl.match(line)
        if not matches and proc_name:
            modi_lines.append(f"END {proc_type} {proc_name}\n")
            modi_lines.append(f"END INTERFACE\n")
            modi_lines.append(f"END MODULE MODI_{proc_name}\n")
            proc_name = ""
            continue

        if proc_name and reg_is_variable_decl.match(line) and "INTENT" not in line:
            # if there are no intent, it might be a local variable
            # but sometimes the intent is missing, so we check also
            # if it's not in the list of arguments, it can also be the result of a function
            matches = reg_variable_name_simple.search(line)
            if matches:
                var_name = matches.group(1)
                if var_name not in proc_args and var_name != proc_name:
                    continue

        if proc_name:
            modi_lines.append(line)

    remove_dangling_ifdef(modi_lines)
    remove_unsed_modules(modi_lines)
    modi_lines = break_long_lines(modi_lines)
    modi = open(modi_path, "w")
    modi.write("".join(modi_lines))


modi_dir = Path("modi")
args = sys.argv
try:
    idx = args.index("-o")
except:
    idx = None
if idx:
    modi_dir = Path(args[idx + 1])
    args.pop(idx + 1)
    args.pop(idx)

for arg in sys.argv[1:]:
    file = Path(arg)
    gen_modi(file)
