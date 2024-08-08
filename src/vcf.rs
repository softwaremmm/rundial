use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

pub mod variant_record;
pub use variant_record::{VariantRecord};
pub mod vcf_header;
pub use vcf_header::VCFHeader;
pub mod vcf_error;
pub use vcf_error::VCFError;
pub mod record_value;
pub use record_value::RecordValue;

pub struct VCFReader {
    pub reader: BufReader<File>,
    pub header: VCFHeader,
    pub curr_line: Option<String>,
}

impl VCFReader {
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

    pub fn header(&self) -> &VCFHeader {
        &self.header
    }
}

impl Iterator for VCFReader {
    type Item = Result<VariantRecord, Box<dyn Error>>;

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

pub struct VCFWriter {
    pub writer: BufWriter<File>,
    pub header: VCFHeader,
}

impl VCFWriter {
    pub fn new(mut writer: BufWriter<File>, header: &VCFHeader) -> Result<Self, Box<dyn Error>> {
        writer.write_all(header.to_string().as_bytes())?;

        Ok(VCFWriter {
            writer,
            header: header.clone(),
        })
    }

    pub fn write_record(&mut self, record: &VariantRecord) -> Result<(), Box<dyn Error>> {
        self.writer.write_all(record.to_string().as_bytes())?;
        self.writer.write_all("\n".as_bytes())?;
        Ok(())
    }
}