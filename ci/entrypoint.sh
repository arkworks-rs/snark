#!/bin/bash
set -euo pipefail

# check required vars are set
if [ -z "${CONTAINER_JAVA_VER+x}" ] && [ -z "${CONTAINER_RUST_VER+x}" ]; then
  echo "CONTAINER_JAVA_VER and CONTAINER_RUST_VER environment variables need to be set!"
  exit 1
fi

# Add local zenbuilder user
# Either use LOCAL_USER_ID:LOCAL_GRP_ID if set via environment
# or fallback to 9001:9001

USER_ID=${LOCAL_USER_ID:-2000}
GRP_ID=${LOCAL_GRP_ID:-2000}

getent group zenbuilder > /dev/null 2>&1 || groupadd -g $GRP_ID zenbuilder
id -u zenbuilder > /dev/null 2>&1 || useradd --shell /bin/bash -u $USER_ID -g $GRP_ID -o -c "" -m zenbuilder

LOCAL_UID=$(id -u zenbuilder)
LOCAL_GID=$(getent group zenbuilder | cut -d ":" -f 3)

if [ ! "$USER_ID" == "$LOCAL_UID" ] || [ ! "$GRP_ID" == "$LOCAL_GID" ]; then
    echo "Warning: User zenbuilder with differing UID $LOCAL_UID/GID $LOCAL_GID already exists, most likely this container was started before with a different UID/GID. Re-create it to change UID/GID."
fi

echo "Starting with UID/GID: $LOCAL_UID:$LOCAL_GID"

export HOME=/home/zenbuilder

# Get Java $CONTAINER_JAVA_VER
apt-get update
DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y "$CONTAINER_JAVA_VER" maven
apt-get -y clean
apt-get -y autoclean
rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*.deb

# Get Rust $CONTAINER_RUST_VER
curl https://sh.rustup.rs -sSf | gosu zenbuilder bash -s -- --default-toolchain none -y
gosu zenbuilder echo 'source $HOME/.cargo/env' >> $HOME/.bashrc
export PATH="/home/zenbuilder/.cargo/bin:${PATH}"
gosu zenbuilder rustup toolchain install "$CONTAINER_RUST_VER"
gosu zenbuilder rustup target add --toolchain "$CONTAINER_RUST_VER" x86_64-pc-windows-gnu
# fix "error: could not compile `api`." "/usr/bin/ld: unrecognized option '--nxcompat'"
# https://github.com/rust-lang/rust/issues/32859#issuecomment-284308455
# appears to be fixed in rust 1.42.0
gosu zenbuilder cat << EOF > $HOME/.cargo/config
[target.x86_64-pc-windows-gnu]
linker = "$(which x86_64-w64-mingw32-gcc)"
EOF

# Print version information
gosu zenbuilder java -version
gosu zenbuilder rustc --version

# Fix ownership recursively
chown -RH zenbuilder:zenbuilder /build

exec "$@"
