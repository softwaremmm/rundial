//! Module for parsing vcf files with sensible types
//!
//! The [VCFReader] will take a file and provide access to the [VCFHeader] and
//! can be iterated to access [VariantRecord]s.
//!
//! The [VCFWriter] can then be used to write to files.

use std::{
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

pub mod variant_record;
pub use variant_record::VariantRecord;
pub mod vcf_header;
pub use vcf_header::VCFHeader;
pub mod vcf_error;
pub use vcf_error::VCFError;
pub mod record_value;
pub use record_value::RecordValue;
pub mod genotype;
pub use genotype::Genotype;

/// Reader for a vcf file. Can be iterated to access [VariantRecord]s
pub struct VCFReader<W>
where
    W: BufRead,
{
    reader: W,
    header: VCFHeader,
    curr_line: Option<String>,
}

impl<W> VCFReader<W>
where
    W: BufRead,
{
    /// Instantiate reader from buffered reader.
    ///
    /// Will return an error if the file does not fit the expected format
    pub fn new(mut reader: W) -> Result<Self> {
        let mut header_lines: Vec<String> = Vec::new();
        let mut curr_line = String::new();
        let mut final_line: Option<String> = None;
        loop {
            let size = reader.read_line(&mut curr_line)?;
            if size == 0 {
                // Reached end of line. Probably an empty vcf
                break;
            }
            if curr_line.starts_with('#') {
                header_lines.push(curr_line.trim().to_string());
                curr_line.clear();
            } else {
                final_line = Some(curr_line.trim().to_string());
                break;
            }
        }

        let header = VCFHeader::from_lines(header_lines);

        return Ok(VCFReader {
            reader,
            header,
            curr_line: final_line,
        });
    }

    /// Provides a reference to the header.
    pub fn header(&self) -> &VCFHeader {
        &self.header
    }
}

impl<W> Iterator for VCFReader<W>
where
    W: BufRead,
{
    type Item = Result<VariantRecord>;

    /// Provides the next [VariantRecord] in the vcf
    ///
    /// Will return None once all records have been read.
    /// Will return Some(Error) if any of the vcf lines cannot be parsed
    /// Otherwise returns Some([VariantRecord]).
    fn next(&mut self) -> Option<Self::Item> {
        match &mut self.curr_line {
            None => return None,
            Some(line) => {
                let next = match VariantRecord::from_string(&self.header, line.trim()) {
                    Err(err) => return Some(Err(err)),
                    Ok(record) => Some(Ok(record)),
                };

                // read next line
                line.clear();
                match self.reader.read_line(line) {
                    Err(err) => return Some(Err(Box::new(err))),
                    Ok(size) => {
                        if size == 0 {
                            self.curr_line = None;
                        }
                    }
                }

                return next;
            }
        }
    }
}

impl VCFReader<BufReader<Box<dyn std::io::Read>>> {
    pub fn from_path<P: AsRef<Path>>(file_path: P) -> Result<Self> {
        let (reader, _format) = niffler::from_path(file_path)
            .map_err(|e| format!("Failed to read input vcf file. Error: {e}"))?;
        let buf_reader = BufReader::new(reader);
        VCFReader::new(buf_reader)
    }
}

/// Struct for writing VCF file
///
/// Note: the file header cannot be changed once the writer is created.
pub struct VCFWriter<W>
where
    W: Write,
{
    pub writer: W,
    header: VCFHeader,
}

impl<W> VCFWriter<W>
where
    W: Write,
{
    /// Create new writer.
    pub fn new(mut writer: W, header: VCFHeader) -> Result<Self> {
        writer.write_all(header.to_string().as_bytes())?;

        Ok(VCFWriter { writer, header })
    }

    /// Write record to file.
    ///
    /// Note that the record is not validated against the [VCFWriter]'s header.
    pub fn write_record(&mut self, record: &VariantRecord) -> Result<()> {
        self.writer.write_all(record.to_string().as_bytes())?;
        self.writer.write_all("\n".as_bytes())?;
        Ok(())
    }

    /// Provides a reference to the header.
    pub fn header(&self) -> &VCFHeader {
        &self.header
    }
}

impl VCFWriter<BufWriter<Box<dyn Write>>> {
    pub fn to_path<P: AsRef<Path>>(file_path: P, header: VCFHeader) -> Result<Self> {
        let niffler_writer = potentially_gzipped_writer(file_path)?;
        let buf_writer = BufWriter::new(niffler_writer);
        VCFWriter::new(buf_writer, header)
    }
}

fn potentially_gzipped_writer<P: AsRef<Path>>(file: P) -> Result<Box<dyn Write>> {
    let (nif_format, level) = if file.as_ref().extension().unwrap_or_default() == "gz" {
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

    let file_str = file.as_ref().display().to_string();
    let niffler_writer = niffler::to_path(file, nif_format, level)
        .map_err(|e| format!("Failed to open fasta output file {file_str}. Error: {e}"))?;

    return Ok(niffler_writer);
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::path::Path;

    use super::*;
    use std::io;

    fn read_lines<P: AsRef<Path>>(path: P) -> io::Result<Vec<String>> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        reader.lines().collect()
    }

    #[test]
    fn test_read_and_write_file() -> Result<()> {
        let vcf_reader = VCFReader::new(BufReader::new(File::open("test_data/example.vcf")?))?;
        let header = vcf_reader.header();

        let new_header = header.clone();
        let output_file = "tests/test_outputs/test_read_and_write_file.vcf";
        let mut vcf_writer =
            VCFWriter::new(BufWriter::new(File::create(output_file)?), new_header)?;

        for record in vcf_reader {
            let record = record?;
            vcf_writer.write_record(&record)?;
        }

        vcf_writer.writer.flush()?;

        // check files are equal
        let initial_lines = read_lines("test_data/example.vcf")?;
        let new_lines = read_lines(output_file)?;

        assert_eq!(initial_lines, new_lines);
        Ok(())
    }

    #[test]
    fn test_read_and_write_gzipped() -> Result<()> {
        let vcf_reader = VCFReader::from_path("test_data/example.vcf.gz")?;
        let header = vcf_reader.header();

        let output_file = "tests/test_outputs/test_read_and_write_file.vcf.gz";
        let mut vcf_writer = VCFWriter::to_path(output_file, header.clone())?;
        for record in vcf_reader {
            let record = record?;
            vcf_writer.write_record(&record)?;
        }
        vcf_writer.writer.flush()?;
        drop(vcf_writer);
        print!("flushed");

        // check files are equal
        let initial_lines = read_lines("test_data/example.vcf")?;

        let (niffler_reader, format) = niffler::from_path(output_file)?;
        assert_eq!(format, niffler::compression::Format::Gzip);
        let buf_reader = BufReader::new(niffler_reader);
        let new_lines: Vec<String> = buf_reader.lines().map(|l| l.unwrap()).collect();

        assert_eq!(initial_lines, new_lines);
        Ok(())
    }
}
