# ---- Rust Build Stage ----
FROM rust:1.88-slim as chef

RUN apt-get update && apt-get upgrade -y && apt-get clean && rm -rf /var/lib/apt/lists/*
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
RUN conda env update -n base --file env.yml && conda clean -afy


# ---- Conda Build Stage ----
FROM conda_builder as runtime
WORKDIR /app

ENV PATH="/opt/conda/bin:$PATH"


# Install Python code
COPY ./pyproject.toml /app/pyproject.toml
COPY src/helper_scripts /app/src/helper_scripts

# Install pytest if TESTING is true (default false)
ARG TESTING=false
RUN if [ "$TESTING" = "true" ]; then \
        pip install .[dev]; \
    else \
        pip install .; \
    fi

# Copy the rundial rust build
COPY --from=builder /app/target/release/rundial /usr/local/bin/rundial
