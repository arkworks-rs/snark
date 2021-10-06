#!/bin/bash

# start with a clean environment
[ -n "${RUSTFLAGS:-}" ] ||
[ -n "${CARGOARGS:-}" ] ||
[ -n "${RUST_CROSS_TARGETS:-}" ] ||
[ -n "${CARGO_AUDIT_EXIT_ON_ERROR:-}" ] ||
[ -n "${RUSTUP_TOOLCHAIN:-}" ] &&
exec -c "$0"

workdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../.." &> /dev/null && pwd )"
mapfile -t globals < <(docker run --rm -v "$workdir":/workdir mikefarah/yq e '.env.global.[]' .travis.yml)
mapfile -t runners < <(docker run --rm -v "$workdir":/workdir mikefarah/yq e '.jobs.include.[].env' .travis.yml)

for i in "${!runners[@]}"; do
  # shellcheck disable=SC1090
  ( source <(for var in "${globals[@]}"; do echo "export $var"; done; echo "export ${runners[$i]}"); \
  "$workdir"/ci/script.sh 2>&1 | tee -a "$workdir/ci/devtools/job-$i.log" )
done
