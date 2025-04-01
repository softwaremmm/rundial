"""Script to run bcftools mpileup in parallel on a sorted bam file."""

import argparse
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor


def main():
    """Main function to parse arguments and run bcftools mpileup in parallel."""
    parser = argparse.ArgumentParser(description="Run bcftools mpileup in parallel.")
    parser.add_argument(
        "-a",
        "--alignment",
        required=True,
        type=str,
        help="The sorted bam alignment file.",
    )
    parser.add_argument(
        "-r", "--reference", required=True, type=str, help="The reference fasta file."
    )
    parser.add_argument(
        "-o", "--output", required=True, type=str, help="The output vcf file."
    )
    parser.add_argument(
        "-t", "--num-threads", type=int, default=4, help="Number of threads to use."
    )
    parser.add_argument(
        "--settings", required=True, type=str, help="Settings for bcftools mpileup."
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=500000,
        help="Number of bases to process in each batch.",
    )

    args = parser.parse_args()
    alignment = args.alignment
    reference = args.reference
    output = args.output
    num_threads = int(args.num_threads)
    settings = args.settings
    batch_size = args.batch_size

    # first need to index the bam file and reference
    subprocess.run(["samtools", "index", alignment], check=True)
    subprocess.run(["samtools", "faidx", reference], check=True)

    # read the reference index as a tsv
    contigs = []
    with open(reference + ".fai", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            contigs.append((parts[0], int(parts[1])))

    print(contigs)
    jobs = []
    for contig, length in contigs:
        if length <= batch_size * 1.5:
            # Process contig in one go
            jobs.append(
                {
                    "contig": contig,
                    "start": 1,
                    "end": length,
                }
            )
            continue

        # Process contig in batches
        for start in range(1, length + 1, batch_size):
            end = min(start + batch_size - 1, length)
            jobs.append(
                {
                    "contig": contig,
                    "start": start,
                    "end": end,
                }
            )

    commands = [
        (job["contig"], job["start"], job["end"], alignment, reference, settings)
        for job in jobs
    ]

    with ProcessPoolExecutor(max_workers=num_threads) as executor:
        _results = list(executor.map(run_pileup, commands))

    command = f"ls pileups/*.bcf | sort -V | bcftools concat -Ou -f - -o {output}"
    print("Running command:", command)
    subprocess.run(command, shell=True, check=True)


def run_pileup(args: tuple[str, int, int, str, str, str]) -> str:
    """Runs mpileup on a single region."""
    contig, start, end, alignment, reference, settings = args
    os.makedirs("pileups", exist_ok=True)
    output = f"pileups/{contig}_{start}_{end}.bcf"

    # setting -d 10000 ensures no differences from chunking the reference
    command = f"bcftools mpileup -f {reference} -d 10000 -t {contig}:{start}-{end} {settings} -Ou -o {output} {alignment}"
    print("Running command:", command)
    subprocess.run(command, shell=True, check=True)

    return output


if __name__ == "__main__":
    main()
