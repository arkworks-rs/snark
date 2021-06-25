#!/bin/bash

set -xeo pipefail

if grep -q 'Cargo.lock' .gitignore &> /dev/null; then
  rm -f Cargo.lock
fi

# shellcheck disable=SC2086
cargo $CARGOARGS clean

