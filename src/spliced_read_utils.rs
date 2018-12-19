use rust_htslib::bam::record::{Cigar, CigarError, CigarStringView};

/// Special struct for dealing with spliced reads.
///
/// Modified to handle N CIGARs apppropriately for RNASeq data.
///
/// Finds the position within a spliced read.
pub trait SplicedReadCigarStringView {
    fn read_pos_spliced(
        &self,
        ref_pos: u32,
        include_softclips: bool,
        include_dels: bool,
        record_pos: u32,
    ) -> Result<Option<u32>, CigarError>;
}

impl SplicedReadCigarStringView for CigarStringView {
    /// Trait implementation for regular CigarStringView.
    fn read_pos_spliced(
        &self,
        ref_pos: u32,
        include_softclips: bool,
        include_dels: bool,
        record_pos: u32,
    ) -> Result<Option<u32>, CigarError> {
        let mut rpos = record_pos; // reference position
        let mut qpos = 0u32; // position within read
        let mut j = 0; // index into cigar operation vector
        for (i, c) in self.iter().enumerate() {
            match c {
                &Cigar::Match(_) | &Cigar::Diff(_) | &Cigar::Equal(_) | &Cigar::Ins(_) => {
                    j = i;
                    break;
                }
                &Cigar::SoftClip(l) => {
                    j = i;
                    if include_softclips {
                        rpos = rpos.saturating_sub(l);
                    }
                    break;
                }
                &Cigar::Del(_) => {
                    return Err(CigarError::UnexpectedOperation(
                        "'deletion' (D) found before any operation describing read sequence"
                            .to_owned(),
                    ));
                }
                &Cigar::RefSkip(_) => {
                    return Err(CigarError::UnexpectedOperation(
                        "'reference skip' (N) found before any operation describing read sequence"
                            .to_owned(),
                    ));
                }
                &Cigar::HardClip(_) if i > 0 && i < self.len() - 1 => {
                    return Err(CigarError::UnexpectedOperation(
                        "'hard clip' (H) found in between operations, contradicting SAMv1 spec that hard clips can only be at the ends of reads".to_owned()
                    ));
                }
                &Cigar::Pad(_) | &Cigar::HardClip(_) if i == self.len() - 1 => return Ok(None),
                &Cigar::Pad(_) | &Cigar::HardClip(_) => (),
            }
        }

        let contains_ref_pos = |cigar_op_start: u32, cigar_op_length: u32| {
            cigar_op_start <= ref_pos && cigar_op_start + cigar_op_length > ref_pos
        };

        while rpos <= ref_pos && j < self.len() {
            match &self[j] {
                // potential SNV evidence
                &Cigar::Match(l) | &Cigar::Diff(l) | &Cigar::Equal(l)
                    if contains_ref_pos(rpos, l) =>
                {
                    // difference between desired position and first position of current cigar
                    // operation
                    qpos += ref_pos - rpos;
                    return Ok(Some(qpos));
                }
                &Cigar::SoftClip(l) if include_softclips && contains_ref_pos(rpos, l) => {
                    qpos += ref_pos - rpos;
                    return Ok(Some(qpos));
                }
                &Cigar::Del(l) if include_dels && contains_ref_pos(rpos, l) => {
                    // qpos shall resemble the start of the deletion
                    return Ok(Some(qpos));
                }

                &Cigar::RefSkip(l) if contains_ref_pos(rpos, l) => {
                    // Added to handle introns. This translates as though we said:
                    // If the normal rust-htslib would have said that this position
                    // is contained by the skipped region (intron), then rather than
                    // letting the loop exit and return nothing,
                    // instead move through without incrementing anything
                    j += 1;
                }

                // for others, just increase pos and qpos as needed
                &Cigar::Match(l) | &Cigar::Diff(l) | &Cigar::Equal(l) => {
                    rpos += l;
                    qpos += l;
                    j += 1;
                }
                &Cigar::SoftClip(l) => {
                    qpos += l;
                    j += 1;
                    if include_softclips {
                        rpos += l;
                    }
                }
                &Cigar::Ins(l) => {
                    qpos += l;
                    j += 1;
                }
                &Cigar::RefSkip(l) | &Cigar::Del(l) => {
                    rpos += l;
                    j += 1;
                }
                &Cigar::Pad(_) => {
                    j += 1;
                }
                &Cigar::HardClip(_) if j < self.len() - 1 => {
                    return Err(CigarError::UnexpectedOperation(
                        "'hard clip' (H) found in between operations, contradicting SAMv1 spec that hard clips can only be at the ends of reads".to_owned()
                    ));
                }
                &Cigar::HardClip(_) => return Ok(None),
            }
        }

        Ok(None)
    }
}
