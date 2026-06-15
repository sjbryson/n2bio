//! n2bio/pfqsim/src/inserts.rs
//! 

use rand::Rng;
use rand_distr::{Normal, Distribution};
use crate::simstats::InsertModel; 

pub struct InsertSize {
    dist: Normal<f64>,
}

impl InsertSize {
    pub fn new(params: &InsertModel) -> Result<Self, &'static str> {
        // Access the nested NormalDistParams fields within InsertModel
        let mean: f64 = params.insert_dist.mean;
        let std_dev: f64 = if params.insert_dist.std_dev <= 0.0 { 0.1 } else { params.insert_dist.std_dev };
        
        let dist: Normal<f64> = Normal::new(mean, std_dev)
            .map_err(|_| "Failed to create normal distribution for insert size")?;
        
        Ok(Self { dist })
    }

    /// Samples a random insert size, rounding to the nearest integer
    pub fn sample<R: Rng>(&self, rng: &mut R) -> usize {
        self.dist.sample(rng).round().max(0.0) as usize
    }
}