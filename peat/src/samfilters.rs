//! n2bio/peat/src/pass.rs
//! 
#![allow(unused)]

use clap::Args;
use n2core::sam::{SamStr, SamFields, SamFlags, SamTags, AlignmentStats};
use crate::cli::ThresholdMetrics;


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
    if args.base_score.is_some_and(|max: f32| sam.calculate_as_al().ok().flatten().is_some_and(|val: f32| val > max)) {
        return false;
    }
    if args.mapq.is_some_and(|max: u32| sam.mapq() > max) {
        return false;
    }
    // If it is mapped and not exceeding any of the specified max thresholds
    true
}

/// Highpass filter logic - must be mapped or pass all defined thresholds
pub(crate) fn highpass_filter(sam: &SamStr, args: &ThresholdMetrics, thresholds: bool) -> bool {
    // Check if read is mapped
    if !sam.is_mapped() {
        return false;
    }
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
    if args.base_score.is_some_and(|min: f32| sam.calculate_as_al().ok().flatten().is_some_and(|val: f32| val < min)) {
        return false;
    }
    if args.mapq.is_some_and(|min: u32| sam.mapq() < min) {
        return false;
    }
    // If it is mapped and passed all of the specified min thresholds
    true
}