use std::io::Write;
use std::{collections::HashMap, io::BufReader};

use niffler;
use noodles::fasta::record::Sequence;

use noodles::fasta::{
    self as noodles_fasta,
    record::{Definition, Record},
};

use super::Result;
use super::{FILTERED, HET, MASKED, NULL};

pub fn read_fasta(fasta_file: &str) -> Result<HashMap<String, Vec<char>>> {
    let mut consensus: HashMap<String, Vec<char>> = HashMap::new();

    let (reader, _format) = niffler::from_path(fasta_file)
        .map_err(|e| format!("Failed to open fasta input file {fasta_file}. Error: {e}"))?;
    let buf_reader = BufReader::new(reader);
    let mut fasta_reader = noodles_fasta::Reader::new(buf_reader);
    for record in fasta_reader.records() {
        let record = record?;
        let chrom = String::from_utf8_lossy(record.definition().name()).to_string();
        let seq: Vec<char> = record
            .sequence()
            .as_ref()
            .iter()
            .map(|c| *c as char)
            .collect();
        consensus.insert(chrom, seq);
    }

    return Ok(consensus);
}

fn potentially_gzipped_writer(file: &str) -> Result<Box<dyn Write>> {
    let (nif_format, level) = if file.ends_with(".gz") {
        (
            niffler::compression::Format::Gzip,
            niffler::compression::Level::One,
        )
    } else {
        (
            niffler::compression::Format::No,
            niffler::compression::Level::Zero,
        )
    };

    let niffler_writer = niffler::to_path(file, nif_format, level)
        .map_err(|e| format!("Failed to open fasta output file {file}. Error: {e}"))?;

    return Ok(niffler_writer);
}

pub fn save_fasta(consensus: &HashMap<String, Vec<char>>, output_file: &str) -> Result<()> {
    let mut writer = noodles_fasta::Writer::new(potentially_gzipped_writer(output_file)?);

    for (chrom, seq) in consensus.iter() {
        let definition = Definition::new(chrom.clone(), None);
        let sequence = Sequence::from(seq.iter().collect::<String>().as_bytes().to_vec());
        let record = Record::new(definition, sequence);
        writer.write_record(&record)?;
    }
    return Ok(());
}

pub fn clean_fasta_characters(consensus: &mut HashMap<String, Vec<char>>) {
    for (_, seq) in consensus.iter_mut() {
        for base in seq.iter_mut() {
            if [FILTERED, HET, MASKED].contains(base) {
                *base = NULL;
            }
        }
    }
}
