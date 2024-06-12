DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

export PATH=$DIR:$PATH
export PYTHONPATH=$PYTHONPATH:$DIR/../build/with_fcm/arch_gnu/build/bin

PYFT_DIR=$DIR/site/pyft
. $PYFT_DIR/bin/env.sh
