#!/bin/bash
# shellcheck disable=SC2086

set -xeo pipefail

RUSTFLAGS="-C target-feature=+bmi2,+adx --emit=asm" cargo $CARGOARGS test --all-features
