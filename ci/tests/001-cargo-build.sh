#!/bin/bash
# shellcheck disable=SC2086

set -xeo pipefail

retval=0
cargo $CARGOARGS build || retval="$?"
cargo $CARGOARGS build --all-features --tests || retval="$?"
exit "$retval"
