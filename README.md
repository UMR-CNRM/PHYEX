# PHYEX\_tools

This package contains tools related to the PHYEX package (https://github.com/QuentinRodier/PHYEX).
Specifically, the prep\_code.sh scripts prepare the source code for inclusion in the compilation machinery
of the different models.

And, the check\_commit\_ial.sh script compiles, executes IAL test cases and compares the results againts a reference simulation.

Moreover, the check\_commit\_mesonh.sh script compiles, executes a test case and compares the results againts a reference simulation.

The check\_commit\_testprogs.sh does the same for the offline test programs.

verify\_mnh\_expand.py checks the conformace of the code regarding mnh\_expand directives.

## Installation

Instructions can be found in INSTALL file.

## Usage

Help on check\_commit\_ial.sh, check\_commit\_mesonh.sh, check\_commit\_testprogs.sh,
verify\_mnh\_expand.py and prep\_code.sh can be printed with the '-h' option.

For check\_commit\_mesonh.sh the following environment variables can be set:
* MNHPACK: directory in which MNH pack will be created (default is $HOME/MesoNH/PHYEX)
* REFDIR: directory in which reference pack can be found (default is the pack directory near the check\_commit\_mesonh.sh file)
* TARGZDIR: directory in which the tar.gz file can be found (default is the pack directory near the check\_commit\_mesonh.sh file)
