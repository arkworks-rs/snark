#!/bin/bash
set -euo pipefail

# Add local zenbuilder user, either use LOCAL_USER_ID:LOCAL_GRP_ID
# if set via environment or fallback to 9001:9001
USER_ID="${LOCAL_USER_ID:-9001}"
GRP_ID="${LOCAL_GRP_ID:-9001}"
if [ "$USER_ID" != "0" ]; then
    export USERNAME=zenbuilder
    getent group "$GRP_ID" &> /dev/null || groupadd -g "$GRP_ID" "$USERNAME"
    id -u "$USERNAME" &> /dev/null || useradd --shell /bin/bash -u "$USER_ID" -g "$GRP_ID" -o -c "" -m "$USERNAME"
    CURRENT_UID="$(id -u "$USERNAME")"
    CURRENT_GID="$(id -g "$USERNAME")"
    export HOME=/home/"$USERNAME"
    if [ "$USER_ID" != "$CURRENT_UID" ] || [ "$GRP_ID" != "$CURRENT_GID" ]; then
        echo "WARNING: User with differing UID ${CURRENT_UID}/GID ${CURRENT_GID} already exists, most likely this container was started before with a different UID/GID. Re-create it to change UID/GID."
    fi
else
    export USERNAME=root
    export HOME=/root
    CURRENT_UID="$USER_ID"
    CURRENT_GID="$GRP_ID"
    echo "WARNING: Starting container processes as root. This has some security implications and goes against docker best practice."
fi

# Check if CARGOARGS or RUSTUP_TOOLCHAIN is set and is a toolchain of format <channel>[-<date>][-<host>]
# see https://rust-lang.github.io/rustup/concepts/toolchains.html
# and https://rust-lang.github.io/rustup/overrides.html for order.
# If we don't have the corresponding toolchain installed, get it with rustup and
# also install all targets from $RUST_CROSS_TARGETS for the newly added toolchain.
# $RUST_CROSS_TARGETS is set at container build time, but can be overwritten via Travis env.
export CARGO_HOME=/opt/rust/cargo
for toolchain in $(tr " " "\n" <<< "${CARGOARGS:-}" | sed 's/^+//g') "${RUSTUP_TOOLCHAIN:-}"; do
  if grep -qE '^(nightly|stable|beta|[0-9]+\.[0-9]+(\.[0-9]+)?)(-[0-9]{4}-[0-9]{2}-[0-9]{2})?(-[a-zA-Z0-9_-]+)?$' <<< "$toolchain"; then
    toolchain_have="$(rustup toolchain list | grep 'default' | cut -d ' ' -f 1)"
    if ! grep -q "$toolchain" <<< "$toolchain_have"; then
      rustup toolchain install "$toolchain"
      for target in $(tr "," "\n" <<< "${RUST_CROSS_TARGETS:-}"); do
        rustup target add --toolchain "$toolchain" "$target"
      done
    fi
    break
  fi
done
# this installs all targets from RUST_CROSS_TARGETS for the default toolchain
for target in $(tr "," "\n" <<< "${RUST_CROSS_TARGETS:-}"); do
  rustup target add "$target"
done

rustup show
cargo -vV
echo
lscpu
echo
free -h
echo
echo "Username: $USERNAME, HOME: $HOME, UID: $CURRENT_UID, GID: $CURRENT_GID"
echo "CARGOARGS: ${CARGOARGS:-unset}, RUSTFLAGS: ${RUSTFLAGS:-unset}, RUST_CROSS_TARGETS: ${RUST_CROSS_TARGETS:-unset}, RUSTUP_TOOLCHAIN: ${RUSTUP_TOOLCHAIN:-unset}"

# Fix ownership of everything in /build recursively
chown -fR "$CURRENT_UID":"$CURRENT_GID" /build

# We need to unfortunately do this as we can't guarantee that the unprivileged
# user doesn't need to install a different version of Rust.
# The needed Rust version could be overriden by the code in the git repo,
# see https://rust-lang.github.io/rustup/overrides.html#the-toolchain-file.
# Without doing this any attempt to install a new version of Rust would fail
# with permission errors.
# To guarantee proper regression tests with a static version number of Rust,
# the CARGOARGS env var can be set to the desired version.
chown -fR "$CURRENT_UID":"$CURRENT_GID" /opt/rust

if [ "$USERNAME" = "root" ]; then
  exec "$@"
else
  exec gosu "$USERNAME" "$@"
fi
