#!/bin/bash

set -xeuo pipefail

workdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../.." &> /dev/null && pwd )"

# travis-cli lint
cd "$workdir"/ci/devtools/travis-cli
docker build --pull -t travis-cli .
docker run --rm -t -v "$workdir:/mnt" travis-cli lint /mnt/.travis.yml

# run_shellcheck
docker pull koalaman/shellcheck
cd "$workdir"
find . -type f -name '*.sh' -print0 | xargs -0 docker run --rm -t -v "$workdir:/mnt" koalaman/shellcheck
