#![allow(clippy::needless_return)]

use std::path::PathBuf;

use clap::Parser;

mod parameter_structs;
mod filter_vcf;
pub mod vcf;
use filter_vcf::filter_vcf;

use regex::Regex;

#[derive(Parser)]
#[command(name = "VCF-filtering")]
#[command(about = "Filters VCF based on yaml file parameters", long_about = None)]
#[command(version)]
struct Cli {
    #[arg(short, long, value_name = "FILE")]
    in_vcf: PathBuf,
    #[arg(short, long, value_name = "FILE")]
    out_vcf: PathBuf,
    #[arg(short, long, value_name = "FILE")]
    config: PathBuf,

    #[arg(short = 'x', long)]
    overwrite: bool,
    #[arg(short, long)]
    verbose: bool,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    filter_vcf(
        cli.in_vcf,
        cli.out_vcf,
        cli.config,
        cli.overwrite,
        cli.verbose,
    )?;

    Ok(())
}
