.PHONY: build
build:
	cargo build --release

.PHONY: test
test:
	cargo test

.PHONY: benchmark
benchmark:
	snakemake -s Snakefile --cores 8

.PHONY: lint
lint:
	cargo clippy --all-targets --all-features -- -D warnings && cargo fmt --all -- --check