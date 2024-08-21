#![allow(clippy::needless_return)]

use std::path::PathBuf;

use clap::{Parser, Subcommand};

use rundial::consensus::make_consensus;
use rundial::filter_vcf::filter_vcf;

#[derive(Parser, Debug)]
#[command(name = "Rundial")]
#[command(version, about = "tool for processing vcfs", long_about = None)]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    #[command(
        name = "consensus",
        about = "Generates consensus sequence from VCF and reference"
    )]
    Consensus {
        #[arg(short = 'i', long, value_name = "FILE")]
        main_vcf: PathBuf,
        #[arg(short = 's', long, value_name = "FILE")]
        support_vcf: Option<PathBuf>,
        #[arg(short, long, value_name = "FILE")]
        ref_fasta: PathBuf,
        #[arg(short, long, value_name = "FILE")]
        output_root: PathBuf,
        #[arg(short, long, value_name = "FILE")]
        params: PathBuf,

        #[arg(short, long)]
        verbose: bool,
    },
    #[command(name = "filter", about = "Add FILTERs to VCF based on parameters")]
    Filter {
        #[arg(short, long, value_name = "FILE")]
        in_vcf: PathBuf,
        #[arg(short, long, value_name = "FILE")]
        out_vcf: PathBuf,
        #[arg(short, long, value_name = "FILE")]
        params: PathBuf,

        #[arg(short = 'x', long)]
        overwrite: bool,
        #[arg(short, long)]
        verbose: bool,
    },
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Consensus {
            main_vcf,
            support_vcf,
            ref_fasta,
            output_root,
            params,
            verbose,
        } => {
            let support_str = support_vcf
                .as_ref()
                .map(|p| p.to_str().expect("Invalid path"));
            make_consensus(
                main_vcf.to_str().expect("Invalid path"),
                support_str,
                ref_fasta.to_str().expect("Invalid path"),
                output_root.to_str().expect("Invalid path"),
                params.to_str().expect("Invalid path"),
                verbose,
            )?;
        }
        Commands::Filter {
            in_vcf,
            out_vcf,
            params,
            overwrite,
            verbose,
        } => {
            filter_vcf(in_vcf, out_vcf, params, overwrite, verbose)?;
        }
    }

    Ok(())
}
