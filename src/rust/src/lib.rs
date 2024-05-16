use extendr_api::prelude::*;
use rand::{Rng, thread_rng};
use rand::distributions::Uniform;
//use itertools::Itertools;

/// Return a random sample from the d unit simplex
/// @export
#[extendr]
fn sample_unit_simplex(d: u32) -> Vec<f64> {

  let mut rng = thread_rng();
  let unif = Uniform::new(0.0, 1.0);
  let mut vals: Vec<f64> = (0..d-1).map(|_| rng.sample(&unif)).collect();
  vals.push(1.0);
  vals.push(0.0);
  vals.sort_by(|a, b| a.partial_cmp(b).unwrap());

  let mut valsret = vec![0.0; d as usize];
  for ii in 1..d+1 {
    let i = ii as usize;
    valsret[i-1] = vals[i] - vals[i - 1];
  }
  valsret

}


/// Return n random samples from the d unit simplex
/// @export
#[extendr]
fn sample_unit_simplexn(d: u32, n: u32) -> Vec<f64> {

  let mut valsret: Vec<f64> = Vec::with_capacity((d * n) as usize);

  for ii in 0..n {
    let _i = ii as usize;
    valsret.append(&mut sample_unit_simplex(d));
  }
  valsret

}

/// calculate multinomial probabilities
/// @export
#[extendr]
fn calc_probs_rust(sar: Vec<f64>, logt: Vec<f64>, logc: Vec<f64>, d: u32, n: u32, nt: u32) -> Vec<f64> {

  let du: usize = d as usize;
  let nu: usize = n as usize;
  let ntu: usize = nt as usize;

  let sar_ni: Vec<usize> = (0..nu).map(|i| du * i).collect();
  let the_ni: Vec<usize> = (0..ntu).map(|i| du * i).collect();

  let mut res: Vec<f64> = Vec::with_capacity(nt as usize);
  let mut restmp: f64 = 0.0;
  let mut i: usize = 0;

  for tj in the_ni.into_iter() {
    for vj in (*sar_ni).into_iter() {
    restmp = restmp + (
      (0..du).map(|j| sar[*vj + j] * logt[tj + j]).sum::<f64>() + logc[i]
      ).exp();
    i = i + 1;
  }
  i = 0;
  res.push(restmp);
  restmp = 0.0;
  }

  res

}




// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod xactonomial;
    fn sample_unit_simplex;
    fn sample_unit_simplexn;
    fn calc_probs_rust;
   // fn rust_sspace;
}
