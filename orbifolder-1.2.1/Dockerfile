FROM ubuntu:22.04

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential \
    autoconf \
    automake \
    libtool \
    libgsl-dev \
    libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*

COPY . /orbifolder
WORKDIR /orbifolder

# Reconfigure and build
RUN autoreconf -fi || true
RUN ./configure || (cat config.log | tail -50 && false)
RUN make -j$(nproc)

# The binary should be in src/orbifolder/
CMD ["bash"]
