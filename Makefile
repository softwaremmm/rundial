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
