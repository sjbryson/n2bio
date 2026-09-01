//! n2bio/peat/src/pass.rs
//! 

use clap::Args;
use n2core::sam::{SamStr, SamFields, SamFlags, SamTags, AlignmentStats};
use crate::cli::ThresholdMetrics;

// ============================================================================
// SAM Filter Functions
// ============================================================================

pub(crate) fn threshold_args(args: &ThresholdMetrics) -> bool {
        args.align_prop.is_some() || args.align_ident.is_some() || args.align_score.is_some() || 
        args.align_length.is_some() || args.base_score.is_some() || args.mapq.is_some() 
    }

/// Low pass filter logic - must be unmapped or pass all defined thresholds
pub(crate) fn lowpass_filter(sam: &SamStr, args: &ThresholdMetrics, thresholds: bool) -> bool {
    // Keep all unmapped reads
    if !sam.is_mapped() {
        return true;
    }
    // If it is mapped and there are no defined thresholds - return false
    if !thresholds {
        return false;
    }
    // Otherwise evaluate all optional filters
    if args.align_prop.is_some_and(|max: f32| sam.calculate_alignment_proportion().ok().flatten().is_some_and(|val: f32| val > max)) {
        return false;
    }
    if args.align_ident.is_some_and(|max: f32| sam.calculate_alignment_identity().ok().flatten().is_some_and(|val: f32| val > max)) {
        return false;
    }
    if args.align_score.is_some_and(|max: i32| sam.get_int_tag("AS").is_some_and(|val: i32| val > max)) {
        return false;
    }
    if args.align_length.is_some_and(|max: u32| sam.calculate_alignment_length().ok().flatten().is_some_and(|val: u32| val > max)) {
        return false;
    }
    if args.base_score.is_some_and(|max: f32| sam.calculate_base_score().ok().flatten().is_some_and(|val: f32| val > max)) {
        return false;
    }
    if args.mapq.is_some_and(|max: u32| sam.mapq() > max) {
        return false;
    }
    // If it is mapped and not exceeding any of the specified max thresholds - return true
    true
}

/// Highpass filter logic - must be mapped or pass all defined thresholds
pub(crate) fn highpass_filter(sam: &SamStr, args: &ThresholdMetrics, thresholds: bool) -> bool {
    // Check if read is mapped
    if !sam.is_mapped() {
        return false;
    }
    // If it is mapped and there are no defined thresholds - return true
    if !thresholds {
        return true;
    }
    // Evaluate optional filters
    if args.align_prop.is_some_and(|min: f32| sam.calculate_alignment_proportion().ok().flatten().is_some_and(|val: f32| val < min)) {
        return false;
    }
    if args.align_ident.is_some_and(|min: f32| sam.calculate_alignment_identity().ok().flatten().is_some_and(|val: f32| val < min)) {
        return false;
    }
    if args.align_score.is_some_and(|min: i32| sam.get_int_tag("AS").is_some_and(|val: i32| val < min)) {
        return false;
    }
    if args.align_length.is_some_and(|min: u32| sam.calculate_alignment_length().ok().flatten().is_some_and(|val: u32| val < min)) {
        return false;
    }
    if args.base_score.is_some_and(|min: f32| sam.calculate_base_score().ok().flatten().is_some_and(|val: f32| val < min)) {
        return false;
    }
    if args.mapq.is_some_and(|min: u32| sam.mapq() < min) {
        return false;
    }
    // If it is mapped and passed all of the specified min thresholds - return true
    true
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // Mock SAM string helpers
    // Unmapped: FLAG = 4
    const UNMAPPED_SAM: &str = "read1\t4\t*\t0\t0\t*\t*\t0\t0\tATCG\tIIII";
    // Mapped: FLAG = 0, MAPQ = 30, CIGAR = 10M, AS:i:100
    const MAPPED_SAM: &str = "read2\t0\tchr1\t1\t30\t10M\t*\t0\t0\tATCGATCGAT\tIIIIIIIIII\tAS:i:100";

    fn empty_args() -> ThresholdMetrics {
        ThresholdMetrics {
            align_prop: None,
            align_ident: None,
            align_score: None,
            align_length: None,
            base_score: None,
            mapq: None,
        }
    }

    // --- Threshold Detection Tests ---

    #[test]
    fn test_threshold_args_detection() {
        let empty: ThresholdMetrics = empty_args();
        assert!(!threshold_args(&empty));

        let mut active: ThresholdMetrics = empty_args();
        active.mapq = Some(20);
        assert!(threshold_args(&active));
    }

    // --- Lowpass Filter Tests ---

    #[test]
    fn test_lowpass_unmapped_always_passes() {
        let sam: SamStr<'_> = SamStr::new(UNMAPPED_SAM);
        let args: ThresholdMetrics = empty_args();
        let thresholds: bool = threshold_args(&args);

        // Unmapped passes regardless of thresholds state
        assert!(lowpass_filter(&sam, &args, thresholds));
    }

    #[test]
    fn test_lowpass_mapped_no_thresholds_fails() {
        let sam: SamStr<'_> = SamStr::new(MAPPED_SAM);
        let args: ThresholdMetrics = empty_args();
        let thresholds: bool = threshold_args(&args);

        // Mapped read with no active thresholds fails lowpass
        assert!(!lowpass_filter(&sam, &args, thresholds));
    }

    #[test]
    fn test_lowpass_mapped_within_thresholds_passes() {
        let sam: SamStr<'_> = SamStr::new(MAPPED_SAM);
        let mut args: ThresholdMetrics = empty_args();
        args.mapq = Some(60);         // record is 30 <= 60
        args.align_score = Some(150); // record AS is 100 <= 150
        let thresholds: bool = threshold_args(&args);

        assert!(lowpass_filter(&sam, &args, thresholds));
    }

    #[test]
    fn test_lowpass_mapped_exceeding_threshold_fails() {
        let sam: SamStr<'_> = SamStr::new(MAPPED_SAM);
        let mut args: ThresholdMetrics = empty_args();
        args.mapq = Some(20);   // record is 30 > 20 (exceeds max)
        let thresholds: bool = threshold_args(&args);

        assert!(!lowpass_filter(&sam, &args, thresholds));
    }

    // --- Highpass Filter Tests ---

    #[test]
    fn test_highpass_unmapped_always_fails() {
        let sam: SamStr<'_> = SamStr::new(UNMAPPED_SAM);
        let mut args: ThresholdMetrics = empty_args();
        args.mapq = Some(10);
        let thresholds = threshold_args(&args);

        // Highpass rejects all unmapped reads upfront
        assert!(!highpass_filter(&sam, &args, thresholds));
    }

    #[test]
    fn test_highpass_mapped_no_thresholds_passes() {
        let sam: SamStr<'_> = SamStr::new(MAPPED_SAM);
        let args: ThresholdMetrics = empty_args();
        let thresholds: bool = threshold_args(&args);

        // Mapped read with no defined thresholds passes highpass
        assert!(highpass_filter(&sam, &args, thresholds));
    }

    #[test]
    fn test_highpass_mapped_meeting_thresholds_passes() {
        let sam: SamStr<'_> = SamStr::new(MAPPED_SAM);
        let mut args: ThresholdMetrics = empty_args();
        args.mapq = Some(20);        // record is 30 >= 20
        args.align_score = Some(50); // record AS is 100 >= 50
        let thresholds: bool = threshold_args(&args);

        assert!(highpass_filter(&sam, &args, thresholds));
    }

    #[test]
    fn test_highpass_mapped_failing_threshold_fails() {
        let sam: SamStr<'_> = SamStr::new(MAPPED_SAM);
        let mut args: ThresholdMetrics = empty_args();
        args.mapq = Some(40); // record is 30 < 40 (fails min requirement)
        let thresholds: bool = threshold_args(&args);

        assert!(!highpass_filter(&sam, &args, thresholds));
    }
}