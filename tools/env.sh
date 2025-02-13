DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

export PATH=$DIR:$PATH
export PYTHONPATH=$PYTHONPATH:$DIR/../build/with_fcm/arch_gnu/build/bin

. $DIR/site/pyfortool/bin/env.sh
