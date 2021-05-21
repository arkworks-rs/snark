#!/bin/bash
# shellcheck disable=SC2086

set -xeo pipefail

cargo $CARGOARGS check
cargo $CARGOARGS check --all-features --tests
