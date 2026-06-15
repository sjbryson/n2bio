//! n2bio/pfqsim/src/qualities.rs
//! 

use rand::Rng;
use rand_distr::{ Normal, Distribution };
use crate::simstats::{ QualityModel, NormalDistParams };

pub struct QualityScores {
    pub r1_qual: Vec<Normal<f64>>,
    pub r2_qual: Vec<Normal<f64>>,
}

impl QualityScores {
    pub fn new(model: &QualityModel) -> Result<Self, &'static str> {
        let to_normal = |params: &NormalDistParams| {
            let std_dev: f64 = if params.std_dev <= 0.0 { 0.1 } else { params.std_dev };
            Normal::new(params.mean, std_dev).unwrap()
        };

        let r1_qual: Vec<Normal<f64>> = model.r1_quals.iter().map(to_normal).collect();
        let r2_qual: Vec<Normal<f64>> = model.r2_quals.iter().map(to_normal).collect();

        Ok(Self { r1_qual, r2_qual })
    }

    /// Generates a quality string of ASCII characters for the specified read (1 or 2)
    pub fn generate<R: Rng>(&self, rng: &mut R, length: usize, read: u8) -> Vec<u8> {
        let distributions: &Vec<Normal<f64>> = match read {
            1 => &self.r1_qual,
            2 => &self.r2_qual,
            _ => panic!("Read argument must be 1 or 2"),
        };

        distributions.iter().take(length).map(|dist| {
            let q: f64 = dist.sample(rng).round();
            let clamped_q: u8 = q.clamp(0.0, 41.0) as u8;
            clamped_q + 33
        }).collect()
    }
}