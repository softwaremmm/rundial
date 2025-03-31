# ---- Rust Build Stage ----
FROM rust:1.84 as chef

# RUN apt-get update && apt-get install -y clang llvm-dev
RUN cargo install cargo-chef
WORKDIR /app

FROM chef as planner
COPY Cargo.* .
COPY src ./src
RUN cargo chef prepare --recipe-path recipe.json

FROM chef as builder
COPY --from=planner /app/recipe.json recipe.json
RUN cargo chef cook --release --recipe-path recipe.json

# Copy the source code and build
COPY Cargo.* .
COPY src ./src
RUN cargo build --release

# ---- Conda Build Stage ----
FROM continuumio/miniconda3 as conda_builder
WORKDIR /app

COPY env.yml /app/env.yml
RUN conda env create --file env.yml


# Pack conda installation
RUN conda install -c conda-forge conda-pack
RUN conda-pack -n rundial -o /app/conda_env.tar.gz


# ---- Conda Build Stage ----
FROM debian:stable-slim as runtime
WORKDIR /app

RUN apt-get update && apt-get install -y libssl-dev procps && rm -rf /var/lib/apt/lists/*

COPY --from=conda_builder /app/conda_env.tar.gz /app/conda_env.tar.gz
RUN mkdir -p /opt/conda
RUN tar -xzf /app/conda_env.tar.gz -C /opt/conda && rm /app/conda_env.tar.gz

ENV PATH="/opt/conda/bin:$PATH"


# Install Python code if testing is true (default false)
COPY ./pyproject.toml /app/pyproject.toml

# Install pytest if TESTING is true (default false)
ARG TESTING=false
RUN if [ "$TESTING" = "true" ]; then \
        pip install .[dev]; \
    fi

# Copy the rundial rust build
COPY --from=builder /app/target/release/rundial /usr/local/bin/rundial
