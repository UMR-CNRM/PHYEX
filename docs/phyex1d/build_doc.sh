#!/usr/bin/bash

set -e

CWD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd $CWD

# Remove old build
rm -rf html/
mkdir html

# Temporary virtualenv
TEMPDIR=$(mktemp -d)
trap "\rm -rf $TEMPDIR" EXIT
cd $TEMPDIR
python3 -m venv phyex1denv
. phyex1denv/bin/activate
cd $CWD/..
pip install .
cd $CWD

# Main page
cat ../README.md  > html/README.md
echo "" >> html/README.md
echo "## Command line arguments" >> html/README.md
echo "Help message obtained by \`\`\`phyex1d -h\`\`\`" >> html/README.md
echo "\`\`\`" >> html/README.md
phyex1d -h >> html/README.md
echo "\`\`\`" >> html/README.md
echo "" >> html/README.md
echo "## Description of the format for the .grid files" >> html/README.md
echo "Docstring of the Grid class:" >> html/README.md
echo "\`\`\`" >> html/README.md
python3 -c "import phyex1d.grid; print(phyex1d.grid.Grid.__doc__)" >> html/README.md
echo "\`\`\`" >> html/README.md
head -1 html/README.md > html/README.new
echo "" >> html/README.new
cat html/README.md | grep ^\#\# | while read line; do
  level=$(echo $line | cut -d\  -f1 | wc -m)
  title=$(echo $line | cut -d\  -f2-)
  printf ' %.0s' $(seq 1 $((($level-2)*2)))
  echo "- [$title](#$(echo ${title,,} | sed 's/ /-/g' | sed 's/\.//g'  | sed 's/,//g'))"
done >> html/README.new
echo "" >> html/README.new
tail -n +2 html/README.md >> html/README.new
mv html/README.new html/README.md

doxygen Doxyfile
