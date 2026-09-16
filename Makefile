clean:
	find . -type d -name .nextflow | xargs rm -rf
	find . -type d -name .nf-test | xargs sudo rm -rf
	find . -type d -name work | xargs sudo rm -rf
	find . -type f -regex '.*\.nextflow\.log.*' | xargs rm -f
	find . -type f -name trace*txt | xargs rm -f

test:
	cargo test
	nf-test test tests/nextflow/*.nf.test

test_local:
	cargo test
	docker build -t test_container_rundial .
	nf-test test tests/nextflow/*.nf.test --profile local_docker

run_example:
	nextflow run . --workflow clair3 \
		--input_dir test_data/assemblers \
		--publish_dir results \
		--ref-fasta data/h37rv_20231215.fa.gz \
		--clair3_models_dir data/clair3_models \
		--basecalling_model dna_r10.4.1_e8.2_400bps_sup@v4.3.0 \
		-profile local_docker -resume

container:
	docker build -t test_container_rundial .

pc:
	pre-commit run --all-files

clippy:
	cargo clippy --all --all-features --tests -- -D warnings

clippy_fix:
	cargo clippy --all --all-features --tests --fix --allow-dirty -- -D warnings

release:
	cargo build --release
