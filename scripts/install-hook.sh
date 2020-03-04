#!/bin/env bash
# This script will install the provided directory ../.hooks as the hook
# directory for the present repo. See there for hooks, including a pre-commit
# hook that runs rustfmt on files before a commit.

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
HOOKS_DIR="${DIR}/../.hooks"

git config core.hooksPath "$HOOKS_DIR"
