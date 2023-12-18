#!/bin/bash

function usage {
  echo "Usage: $0 [-h] output"
  echo "output         output file with .pdf or .html extension"
  echo "-h             to get this help"
  echo
  echo "This script builds a standalone documentation. The format of the documentation"
  echo "is determined by the output file extension."
  echo "The script has been tested with .pdf and .html extensions."
  echo
  echo "The script must be stored in the tools directory (to find the relevant md files."
  echo "All the md files must begins with a uniq title line of first level (eg: '# title')"
}

#Script must be called with one argument: the output file to produce
output=""
while [ -n "$1" ]; do
  case "$1" in
    '-h') usage;;
    *) if [ -z "${output-}" ]; then
         output="$1"
       else
         echo "Only one argument is allowed, type $0 -h for help"
         exit 97
       fi;;
  esac
  shift
done
if [ -z "${output-}" ]; then
  echo "Script must be called with the output file as first and only argument"
  exit 98
fi
output="$(cd "$(dirname "${output}")" && pwd)"/$(basename $output)
format=$(echo $output | rev | cut -s -d. -f1 | rev)

#Ordered list of md files
mdfiles="PHYEX.md Developer.md CodingNorms.md Integrator.md Offline.md Plugging.md Tools.md"

#Resources needed
resources="AROMEworkflow1.svg AROMEworkflow2.svg"

#Script is assumed to be in the tools directory of PHYEX
PHYEXTOOLSDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ $(basename $PHYEXTOOLSDIR) != 'tools' ]; then
  echo "The script must be put in the tools directory of the PHYEX package"
  exit 96
fi

#Create a temporary directory and set it to be deleted at the end
dir=$(mktemp -d ${TMP:-/tmp}/d.XXXXXX) || exit 99
trap 'rm -rf "$dir"' EXIT

#Copy files in this directory
for file in $mdfiles $resources; do
  cp $PHYEXTOOLSDIR/../docs/$file $dir/
done

#Move to temp dir
cd $dir

#Create helper files
cat > rsvg-convert <<EOF
#!/bin/bash
#Comments from https://github.com/miyako/console-rsvg-convert
#-d, --dpi-x=<float> pixels per inch [optional; defaults to 90dpi]
#-p, --dpi-y=<float> pixels per inch [optional; defaults to 90dpi]
#-x, --x-zoom=<float> x zoom factor [optional; defaults to 1.0]
#-y, --y-zoom=<float> y zoom factor [optional; defaults to 1.0]
#-z, --zoom=<float> zoom factor [optional; defaults to 1.0]
#-w, --width=<int> width [optional; defaults to the SVG's width]
#-h, --height=<int> height [optional; defaults to the SVG's height]
#-f, --format=[png, pdf, ps, svg] [optional; defaults to 'png']
#-o, --output=<path> output filename [optional; defaults to stdout]
#-b, --background-color=[black, white, #abccee, #aaa...] set the background color [optional; defaults to None]
#-u, --base-uri=<uri>
#-v, --version show version information
#
#-u, --unlimited
#-f, --format=[eps, xml, recording]
#-a, --keep-aspect-ratio whether to preserve the aspect ratio [optional; defaults to FALSE]
#--keep-image-data
#--no-keep-image-data

#We only deal with
# -f: ignored because convert detect output format with the extension
# -d, --dpi-x, -p, --dpi-y
# -a: ignore ratio is always kept
# -o

format=''
preserveratio=0
input=''
output=''
dpiX=''
dpiY=''
while [ -n "\$1" ]; do
  case "\$1" in
    '-f') format=\$2; shift;;
    '-d') dpiX=\$2; shift;;
    '--dpi-x') dpiX=\$2; shift;;
    '-p') dpiY=\$2; shift;;
    '--dpi-y') dpiY=\$2; shift;;
    '-a') preserveratio=1;;
    '-o') output="\$2"; shift;;
    *) input="\$1";;
  esac
  shift
done

dpi=""
if [ "\$dpiX" != "" -a "\$dpiY" != "" ]; then
  dpi="-density \$dpiXx\$dpiY"
fi

convert \$dpi \$input \$output
EOF
chmod +x rsvg-convert

if [ "$format" == 'pdf' ]; then
  cat > titlesec.tex <<EOF
\usepackage{sectsty} \sectionfont{\clearpage}

EOF
else
  #normaly useless but toc is not displayed in html
  #if option --include-in-header isn't set
  touch titlesec.tex
fi

cat > title.md <<EOF
---
title: PHYEX (PHYsique EXternalisée)
geometry: margin=2cm
...
EOF

#Links between files
#All the files must begin with a first level title, we get the correspondance
declare -A sections
for file in $mdfiles; do
  if [ "$(head -1 $file | cut -c 1)" != '#' -o  "$(head -1 $file | cut -c 1-2)" == '##' ]; then
    echo "All the md files must begin (first line) with a top level title (one and only one '#')"
    echo "Please check $file file."
    exit 95
  fi
  sections[$file]=$(echo $(head -1 $file | cut -c 2-) | sed -e 's/\(.*\)/\L\1/' | sed 's/ /-/g')
done
#We replace links to these files by links to anchors
for mdfile in $mdfiles; do
  for file in $mdfiles; do
    sed -i "s/](.\/$file)/](#${sections[$file]})/g" $mdfile
    sed -i "s/]($file)/](#${sections[$file]})/g" $mdfile
  done
done

#Generate output file
pandoc --toc --toc-depth=2 --number-sections \
       --include-in-header titlesec.tex \
       -o $output \
       --self-contained \
       title.md $mdfiles
