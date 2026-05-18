# pylint: disable-msg=too-many-locals
# pylint: disable-msg=too-many-arguments
# pylint: disable-msg=too-many-positional-arguments
# pylint: disable-msg=no-else-return

"""Recreate genome creation report json but using mask to provide extra stats"""

import argparse
import json
from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass
class HetSettings:
    min_allele_dp: int
    min_allele_pc: float
    min_strand_pc: float


def read_fasta(path: str) -> dict[str, str]:
    """Read fasta file and return sequences as a dict of contig to sequence"""

    sequences = {}

    with open(path, "r", encoding="utf-8") as file:
        lines = file.readlines()
        current_contig = None
        current_sequence = []
        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                if current_contig is not None:
                    sequences[current_contig] = "".join(current_sequence)
                current_contig = line[1:].split()[0]
                current_sequence = []
            else:
                current_sequence.append(line)
        if current_contig is not None:
            sequences[current_contig] = "".join(current_sequence)

    return sequences


def read_vcf(vcf_path: str) -> pd.DataFrame:
    """Read VCF file into DataFrame, skipping header lines"""
    columns = [
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
        "SAMPLE",
    ]
    usecols = ["CHROM", "POS", "REF", "ALT", "FILTER", "INFO", "FORMAT", "SAMPLE"]
    df = pd.read_csv(vcf_path, sep="\t", comment="#", names=columns, usecols=usecols)
    return df


def read_mask(
    mask_file: str, convert_to_1_indexed: bool, fasta_file: str
) -> pd.DataFrame:
    """Read mask file and return dataframe with cols: contig, position"""

    # mask file may either have header with contig, position or just be a list of positions.
    with open(mask_file, "r", encoding="utf-8") as f:
        first_line = f.readline().strip()
        if "\t" in first_line:
            # assume tab-delimited with contig and position
            df = pd.read_csv(mask_file, sep="\t")
        else:
            # assume just a list of positions
            df = pd.read_csv(mask_file, header=None, names=["position"])
            # use first contig from fasta as default contig for all positions
            with open(fasta_file, "r", encoding="utf-8") as fasta_f:
                for line in fasta_f:
                    if line.startswith(">"):
                        default_contig = line[1:].split()[0]
                        break
            df["contig"] = default_contig

    if convert_to_1_indexed:
        df["position"] = df["position"] + 1

    return df[["contig", "position"]]


def make_format_lookup(format_str: str, field_names: list[str]) -> dict[str, int]:
    """Create a lookup from field name to index in the FORMAT column"""
    format_fields = format_str.split(":")
    return {
        field_name: format_fields.index(field_name)
        for field_name in field_names
        if field_name in format_fields
    }


def extract_info_fields(info_str: str, field_names: list[str]) -> dict[str, str]:
    """Convert INFO string to dict and extract specified fields"""
    info_fields = dict(
        field.split("=", 1) for field in info_str.split(";") if "=" in field
    )
    return {
        field_name: info_fields[field_name]
        for field_name in field_names
        if field_name in info_fields
    }


def to_int_list(s: str) -> list[int]:
    return list(map(int, s.split(",")))


def set_depth_columns(vcf: pd.DataFrame) -> pd.DataFrame:
    """Ideally want to set ADF, ADR, COV, and HQ_DP columns.
    This function is designed to work with a variety of VCFs"""

    vcf["n_alleles"] = (
        vcf["ALT"].fillna(".").apply(lambda s: 1 if s == "." else s.count(",") + 1)
    )

    format_lookup_cache: dict[str, dict[str, int]] = {}
    for format_str in vcf["FORMAT"].unique():
        format_lookup_cache[format_str] = make_format_lookup(
            format_str, ["ADF", "ADR", "COV", "DP4", "DP_ACGT", "AD"]
        )

    def get_depths(row) -> tuple[list[int] | None, list[int] | None, list[int], int]:
        if row["SAMPLE"].startswith("./."):
            return [0], [0], [0], 0

        format_lookup = format_lookup_cache[row["FORMAT"]]
        sample_values = row["SAMPLE"].split(":")
        format_fields = {
            field_name: sample_values[index]
            for field_name, index in format_lookup.items()
        }
        info_fields = extract_info_fields(
            row["INFO"], ["ADF", "ADR", "COV", "DP4", "DP_ACGT", "AD"]
        )
        fields = {**info_fields, **format_fields}

        if "ADF" in fields:
            adf = to_int_list(fields["ADF"])
            adr = to_int_list(fields["ADR"])
            assert len(adf) == len(adr)
            cov = [a + r for a, r in zip(adf, adr)]
            hq_df = sum(cov)
            return adf, adr, cov, hq_df
        elif row["n_alleles"] <= 2 and "DP4" in fields:
            dp4 = to_int_list(fields["DP4"])
            assert len(dp4) == 4
            if row["n_alleles"] == 1:
                adf = [dp4[0]]
                adr = [dp4[1]]
                cov = [dp4[0] + dp4[1]]
                return adf, adr, cov, sum(dp4)
            else:
                adf = [dp4[0], dp4[2]]
                adr = [dp4[1], dp4[3]]
                cov = [dp4[0] + dp4[1], dp4[2] + dp4[3]]
                return adf, adr, cov, sum(dp4)
        elif not row["indel_like"] and "DP_ACGT" in fields:
            dp_acgt = to_int_list(fields["DP_ACGT"])
            assert len(dp_acgt) == 8
            forward_dp = {
                "A": dp_acgt[0],
                "C": dp_acgt[2],
                "G": dp_acgt[4],
                "T": dp_acgt[6],
            }
            reverse_dp = {
                "A": dp_acgt[1],
                "C": dp_acgt[3],
                "G": dp_acgt[5],
                "T": dp_acgt[7],
            }
            alleles = (
                [row["REF"]] + row["ALT"].split(",")
                if row["ALT"] != "."
                else [row["REF"]]
            )
            adf = [forward_dp[allele] for allele in alleles]
            adr = [reverse_dp[allele] for allele in alleles]
            cov = [a + r for a, r in zip(adf, adr)]
            hq_df = sum(cov)
            return adf, adr, cov, hq_df
        # From here on we give up on strand-specific information
        elif "COV" in fields:
            cov = to_int_list(fields["COV"])
            return None, None, cov, sum(cov)
        elif "AD" in fields:
            cov = to_int_list(fields["AD"])
            return None, None, cov, sum(cov)

        raise ValueError(
            f"Could not find suitable depth fields for row with FORMAT {row['FORMAT']}.\nrow={row}"
        )

    if vcf.empty:
        vcf[["ADF", "ADR", "COV", "HQ_DP"]] = pd.DataFrame(
            [[None, None, [0], 0]], index=vcf.index
        )
    else:
        vcf[["ADF", "ADR", "COV", "HQ_DP"]] = vcf.apply(
            get_depths, axis=1, result_type="expand"
        )
    return vcf


def set_het_flag(vcf: pd.DataFrame, het_settings: HetSettings) -> pd.DataFrame:
    def is_het(row) -> bool:
        if len(row["COV"]) < 2:
            return False
        hq_dp = row["HQ_DP"]

        n_supported_alleles = 0
        if row["ADF"] is not None and row["ADR"] is not None:
            for adf, adr in zip(row["ADF"], row["ADR"]):
                total = adf + adr
                if (
                    total < het_settings.min_allele_dp
                    or total / hq_dp < het_settings.min_allele_pc
                ):
                    continue
                if min(adf, adr) / total < het_settings.min_strand_pc:
                    continue
                n_supported_alleles += 1
        else:
            n_supported_alleles = sum(
                1
                for dp in row["COV"]
                if dp >= het_settings.min_allele_dp
                and dp / hq_dp >= het_settings.min_allele_pc
            )

        if n_supported_alleles < 2:
            return False
        return True

    if vcf.empty:
        vcf["is_het"] = False
    else:
        vcf["is_het"] = vcf.apply(is_het, axis=1)
    return vcf


def set_is_filtered(vcf: pd.DataFrame) -> pd.DataFrame:
    allowed_filters = {
        "PASS",
        ".",
        "MIN_FRS",
        "MIN_GCP",
        "OVERLAP_BETTER_VARIANT",
        "SNP_IN_SUPPORT_VCF",
    }
    vcf["is_filtered"] = ~vcf["FILTER"].str.split(";").apply(
        lambda filters: all(f in allowed_filters for f in filters)
    )
    return vcf


def set_is_masked(vcf: pd.DataFrame, mask: pd.DataFrame) -> pd.DataFrame:
    def in_mask(contig, pos, mask_set, ref):
        for i in range(len(ref)):
            if (contig, pos + i) in mask_set:
                return True
        return False

    mask_set = set(zip(mask["contig"], mask["position"]))

    if mask_set and not vcf.empty:
        vcf["is_masked"] = vcf.apply(
            lambda row: in_mask(row["CHROM"], row["POS"], mask_set, row["REF"]), axis=1
        )
    else:
        vcf["is_masked"] = False
    return vcf


def report_het_stats(hets_df: pd.DataFrame, cluster_window: int) -> dict:
    """Report mixed calls, snps and snp clusters"""
    stats = {}
    stats["Mixed calls"] = len(hets_df)
    snps = hets_df[~hets_df["indel_like"]]
    stats["Mixed snps"] = len(snps)

    # each row is a new cluster if it's more than cluster_window away from the previous row
    snps = snps.sort_values("POS")
    stats["Mixed snp clusters"] = (
        snps["POS"].diff().fillna(cluster_window + 1) > cluster_window
    ).sum()
    stats["Isolated mixed snps"] = (
        (snps["POS"].diff().fillna(cluster_window + 1) > cluster_window)
        & (snps["POS"].diff(-1).fillna(cluster_window + 1).abs() > cluster_window)
    ).sum()

    return stats


def annotate_vcf(
    vcf: pd.DataFrame,
    mask: pd.DataFrame,
    het_settings: HetSettings,
) -> pd.DataFrame:
    vcf["indel_like"] = (vcf["REF"].str.len() != 1) | vcf["ALT"].fillna(".").apply(
        lambda s: any(len(a) != 1 for a in s.split(",")) if s != "." else False
    )
    vcf = set_depth_columns(vcf)
    vcf = set_is_filtered(vcf)
    vcf = set_het_flag(
        vcf,
        het_settings,
    )
    vcf = set_is_masked(vcf, mask)
    return vcf


def count_deletions(vcf: pd.DataFrame) -> tuple[int, int]:
    """Count the number of deleted bases as total and unmasked"""
    # Only count completely non-filtered indels
    indels = vcf[vcf["indel_like"] & vcf["FILTER"].isin([".", "PASS"])].copy()

    if indels.empty:
        return 0, 0

    def get_deletion_length(row):
        gt = row["SAMPLE"].split(":")[0]
        if "/" in gt:
            part1, part2 = gt.split("/")
            if part1 != part2 or part1 == "." or part1 == "0":
                return 0
            allele_index = int(part1) - 1
        else:
            if gt != "0" and gt != ".":
                allele_index = int(gt) - 1
            else:
                return 0

        alt = row["ALT"].split(",")[allele_index]
        if len(alt) >= len(row["REF"]):
            return 0
        return len(row["REF"]) - len(alt)

    indels["deletion_length"] = indels.apply(get_deletion_length, axis=1)
    total_deleted_bases = indels["deletion_length"].sum()
    unmasked_deleted_bases = indels[~indels["is_masked"]]["deletion_length"].sum()
    return total_deleted_bases, unmasked_deleted_bases


def count_nulls(fasta_file: str, mask: pd.DataFrame) -> tuple[int, int, int]:
    """Count the number of null bases as total and unmasked, and return genome length"""
    fasta_seqs = read_fasta(fasta_file)
    total_nulls = 0
    unmasked_nulls = 0
    genome_length = 0

    for contig, fasta_seq in fasta_seqs.items():
        contig_mask = set(mask[mask["contig"] == contig]["position"])
        null_positions = set(
            index + 1
            for index, base in enumerate(fasta_seq)
            if base.upper() == "N" or base == "-"
        )
        total_nulls += len(null_positions)
        unmasked_nulls += len(null_positions - contig_mask)
        genome_length += len(fasta_seq)

    return total_nulls, unmasked_nulls, genome_length


def vcf_file_to_stats(
    vcf_path: str,
    fasta_file: str,
    mask: pd.DataFrame,
    het_settings: HetSettings,
    cluster_window: int,
) -> dict:
    """From a vcf and fasta return df of het calls, set of indel positions and set of null positions"""
    vcf = read_vcf(vcf_path)

    total_nulls, unmasked_nulls, genome_length = count_nulls(fasta_file, mask)

    # Only actually focused on variants for the most part
    variants = vcf[vcf["ALT"] != "."].copy()
    variants = annotate_vcf(variants, mask, het_settings)

    total_dels, unmasked_dels = count_deletions(variants)

    hets = variants[
        variants["is_het"] & ~variants["is_filtered"] & ~variants["is_masked"]
    ].copy()

    het_stats = report_het_stats(hets, cluster_window)

    real_nulls = max(total_nulls - total_dels, 0)
    unmasked_real_nulls = max(unmasked_nulls - unmasked_dels, 0)
    fixed_coverage = 1 - (real_nulls / genome_length)
    unmasked_fixed_coverage = 1 - (unmasked_real_nulls / genome_length)

    all_stats = {
        "Reference genome length": genome_length,
        "Null bases": real_nulls,
        "Deleted bases": total_dels,
        "Fixed coverage": 100 * fixed_coverage,
        **het_stats,
    }
    if len(mask) > 0:
        all_stats["Unmasked null bases"] = unmasked_real_nulls
        all_stats["Unmasked deleted bases"] = unmasked_dels
        all_stats["Fixed coverage with mask"] = 100 * unmasked_fixed_coverage

    return all_stats


def process_sample(
    vcf_file: str,
    fasta_file: str,
    mask: pd.DataFrame,
    het_settings: HetSettings,
    cluster_window: int,
) -> dict:
    stats = vcf_file_to_stats(vcf_file, fasta_file, mask, het_settings, cluster_window)

    # need to convert any numpy types to native python types for json serialization
    stats = {
        k: (v.item() if isinstance(v, np.generic) else v) for k, v in stats.items()
    }

    return stats


def main():
    parser = argparse.ArgumentParser(description="Assess genome creation stats")
    parser.add_argument("--vcf", required=True, help="VCF file")
    parser.add_argument("--fasta", required=True, help="fasta file")
    parser.add_argument("--mask", help="file path to 0-indexed mask")
    parser.add_argument(
        "--min_allele_dp",
        type=int,
        default=3,
        help="minimum depth for an allele to be considered supported",
    )
    parser.add_argument(
        "--min_allele_pc",
        type=float,
        default=5,
        help="minimum fraction of reads for an allele to be considered supported",
    )
    parser.add_argument(
        "--min_strand_pc",
        type=float,
        default=0.0,
        help="minimum fraction of reads on each strand for an allele to be considered supported",
    )
    parser.add_argument(
        "--cluster_window",
        type=int,
        default=12,
        help="window size for clustering het snps (default: 12)",
    )
    parser.add_argument("-o", "--output", required=True, help="output json report")
    args = parser.parse_args()

    het_settings = HetSettings(
        min_allele_dp=args.min_allele_dp,
        min_allele_pc=args.min_allele_pc / 100
        if args.min_allele_pc >= 1
        else args.min_allele_pc,
        min_strand_pc=args.min_strand_pc,
    )

    # default contig is for TB, which only has one contig so historic mask is just a list of sites
    mask = read_mask(args.mask, convert_to_1_indexed=True, fasta_file=args.fasta)

    stats = process_sample(
        args.vcf, args.fasta, mask, het_settings, args.cluster_window
    )
    json_content = {
        "Sequencing Quality": stats,
    }
    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(json_content, f, indent=4)


if __name__ == "__main__":
    main()
