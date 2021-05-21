#!/bin/bash

set -euo pipefail

run_test () {
  local testscript="$1"
  local code=0
  local start=""
  local end=""
  start="$(date +%s%N)"
  eval "$testscript" && code=$? || code=$?
  end="$(date +%s%N)"
  seconds="$(printf "%010u" $((end-start)) | sed "s/.........$/.&/g")"
  return $code
}

# absolute path to project from relative location of this script
WORKDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." &> /dev/null && pwd )"
export WORKDIR
export LC_ALL=C.UTF-8
# TESTS="001,002,003" to tests=( 001 002 003 )
mapfile -t tests < <(tr "," "\n" <<< "${TESTS:-}")
scripts=()
returncodes=()
durations=()
exitcode=0

# determine what tests to run
if printf "%s\n" "${tests[@]}" | grep -q -P "^\*$"; then
  mapfile -t scripts < <(find "$WORKDIR/ci/tests" -name "*.sh" -type f | sort)
else
  for test in "${tests[@]}"; do
    mapfile -t -O "${#scripts[@]}" scripts < <(find "$WORKDIR/ci/tests" -name "${test}*.sh" -type f)
  done
fi

# run tests and store results and return codes
cd "$WORKDIR"
for i in "${!scripts[@]}"; do
  returncode=0
  seconds=0
  echo -e "\n\n$(date --utc +%FT%T.%3NZ)  Info: Running test ${scripts[$i]}\n\n"
  run_test "${scripts[$i]}" && returncode=$? || returncode=$?
  if [ $returncode -ne 0 ]; then
    exitcode=$returncode
    echo -e "\n\n$(date --utc +%FT%T.%3NZ) Error: Test ${scripts[$i]} failed with code $returncode"
  fi
  returncodes[$i]="$returncode"
  durations[$i]="$seconds"
done

# print short test summary
echo -e "\n\n############ TEST SUMMARY ############\n"
for i in "${!scripts[@]}"; do
  echo -e "Testscript: ${scripts[$i]} Status: $([ "${returncodes[$i]}" -eq 0 ] && echo -n passed || echo -n failed) Exitcode: ${returncodes[$i]} Duration: ${durations[$i]} seconds\n"
done

exit $exitcode
