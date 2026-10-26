# Reproducible build of the Clarabel WebAssembly module (no local Rust toolchain needed).
FROM rust:1.86-slim
RUN rustup target add wasm32-unknown-unknown \
 && cargo install wasm-bindgen-cli --version 0.2.100 --locked
WORKDIR /src
