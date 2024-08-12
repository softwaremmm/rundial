//! Module for parsing vcf files with sensible types
//!
//! The [VCFReader] will take a file and provide access to the [VCFHeader] and
//! can be iterated to access [VariantRecord]s.
//!
//! The [VCFWriter] can then be used to write to files.

use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

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
pub struct VCFReader {
    reader: BufReader<File>,
    header: VCFHeader,
    curr_line: Option<String>,
}

impl VCFReader {
    /// Instantiate reader from vcf file.
    ///
    /// Will return an error if the file does not fit the expected format
    pub fn new(mut reader: BufReader<File>) -> Result<Self, Box<dyn Error>> {
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

impl Iterator for VCFReader {
    type Item = Result<VariantRecord, Box<dyn Error>>;

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

/// Struct for writing VCF file
///
/// Note: the file header cannot be changed once the writer is created.
pub struct VCFWriter {
    pub writer: BufWriter<File>,
    header: VCFHeader,
}

impl VCFWriter {
    /// Create new writer.
    pub fn new(mut writer: BufWriter<File>, header: VCFHeader) -> Result<Self, Box<dyn Error>> {
        writer.write_all(header.to_string().as_bytes())?;

        Ok(VCFWriter { writer, header })
    }

    /// Write record to file.
    ///
    /// Note that the record is not validated against the [VCFWriter]'s header.
    pub fn write_record(&mut self, record: &VariantRecord) -> Result<(), Box<dyn Error>> {
        self.writer.write_all(record.to_string().as_bytes())?;
        self.writer.write_all("\n".as_bytes())?;
        Ok(())
    }

    /// Provides a reference to the header.
    pub fn header(&self) -> &VCFHeader {
        &self.header
    }
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use super::*;
    use std::io;
    use tempfile::NamedTempFile;

    fn read_lines<P: AsRef<Path>>(path: P) -> io::Result<Vec<String>> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        reader.lines().collect()
    }

    #[test]
    fn test_read_and_write_file() -> Result<(), Box<dyn Error>> {
        let vcf_reader = VCFReader::new(BufReader::new(File::open("test_data/example.vcf")?))?;
        let header = vcf_reader.header();

        let new_header = header.clone();
        let temp_file = NamedTempFile::new().expect("Failed to create temporary file");
        let mut vcf_writer =
            VCFWriter::new(BufWriter::new(File::create(temp_file.path())?), new_header)?;

        for record in vcf_reader {
            let record = record?;
            vcf_writer.write_record(&record)?;
        }

        vcf_writer.writer.flush()?;

        // check files are equal
        let initial_lines = read_lines("test_data/example.vcf")?;
        let new_lines = read_lines(temp_file.path())?;

        if initial_lines != new_lines {
            let mut file = File::create("tests/test_outputs/test_read_and_write_file.vcf")?;
            for l in new_lines.iter() {
                write!(file, "{}", l)?;
            }
            panic!("VCFWriter produced a different file to that read in.\nResult written to \"tests/test_outputs/test_read_and_write_file.vcf\"");
        }

        Ok(())
    }
}
