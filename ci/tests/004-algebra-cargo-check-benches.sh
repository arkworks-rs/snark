#!/bin/bash
# shellcheck disable=SC2086

set -xeo pipefail

cd algebra
cargo $CARGOARGS check --all-features --benches
