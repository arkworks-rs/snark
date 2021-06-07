#!/bin/bash
# shellcheck disable=SC2086

set -xeo pipefail

retval=0
cd proof-systems
cargo $CARGOARGS check || retval="$?"
cargo $CARGOARGS build --all-features --tests --examples || retval="$?"
exit "$retval"
