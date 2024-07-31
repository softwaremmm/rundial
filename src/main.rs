use std::io::{self, BufWriter};
use std::path::PathBuf;
use std::fs::File;
use std::io::Read;

use clap::Parser;

use noodles::vcf;
use noodles::vcf::header::record::value::{map::Filter, Map};

mod parameter_structs;
use parameter_structs::FilterParams;

mod filter_vcf;
use filter_vcf::filter_vcf;



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

    filter_vcf(cli.in_vcf, cli.out_vcf, cli.config, cli.overwrite, cli.verbose)?;

    Ok(())
}
