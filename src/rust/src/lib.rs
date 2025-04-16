pub mod rbindings;
use extendr_api::prelude::*;
use crate::rbindings::*;

/// Enumerate the multinomial sample space
/// @param d The dimension
/// @param n The sample size
/// @returns A vector enumerating the sample space, to be converted to a matrix
/// with d columns and choose(n + d - 1, d - 1) rows
/// @export
/// @examples
/// matrix(sspace_multinom(3, 5), ncol = 3, byrow = TRUE)
#[extendr]
fn sspace_multinom(d: u32, n: u32) -> Vec<u32> {

  let di = d as usize;
  let mut bins = vec![0; di];
  bins[0] = n;
  let mut res = bins.clone();

  loop {
    if &bins[di - 1] == &n {
      break;
    }
    if bins[0] > 0 {
      bins[0] -= 1;
      bins[1] += 1;
    } else {
      let mut nz = 1usize;
      while bins[nz-1] == 0 {
        nz += 1;
      }
      bins[0] = &bins[nz-1] - 1;
      bins[nz] += 1;
      bins[nz-1] = 0;
      }
    res.extend(&bins);
  }
  res

}


/// Return a single random sample from the d unit simplex
/// @param d Dimension
/// @returns A vector of length d
/// @examples
/// sample_unit_simplex(3)
#[extendr(use_rng = true)]
fn sample_unit_simplex(d: u32) -> Vec<f64> {

  unsafe {
      GetRNGstate();

  let mut vals: Vec<f64> = (0..d-1).map(|_|
    Rf_runif(0.0, 1.0)
      ).collect();
  vals.push(1.0);
  vals.push(0.0);
  vals.sort_by(|a, b| a.partial_cmp(b).unwrap());

  PutRNGstate();
  let mut valsret = vec![0.0; d as usize];
  for ii in 1..d+1 {
    let i = ii as usize;
    valsret[i-1] = vals[i] - vals[i - 1];
  }
  valsret
  }
}

/// Sample n times from the unit simplex in d dimensions
/// @param d the dimension
/// @param n the number of samples to take uniformly in the d space
/// @returns The grid over Theta, the parameter space. To be converted to a matrix with d columns and nsamp rows
/// @export
/// @examples
/// matrix(sample_unit_simplexn(3, 10), ncol = 3, byrow = TRUE)
#[extendr]
fn sample_unit_simplexn(d: u32, n: u32) -> Vec<f64> {

  let mut valsret: Vec<f64> = Vec::with_capacity((d * n) as usize);

  for ii in 0..n {
    let _i = ii as usize;
    valsret.append(&mut sample_unit_simplex(d));
  }
  valsret

}

/// Calculate multinomial probabilities
/// @param sar The unrolled matrix containing the portion of the sample space to sum over
/// @param logt The vector of candidate theta values, as sampled from the null space
/// @param logc The vector of log multinomial coefficients see \link{log_multinom_coef}
/// @param d The total dimension, sum(d_j)
/// @param n The sample size
/// @param nt The number of candidate theta values
/// @returns A vector of probabilities
/// @export
/// @examples
/// sspace_3_5 <- sspace_multinom(3, 5)
/// calc_multinom_probs(sspace_3_5, sample_unit_simplexn(3, 10),
///   apply(matrix(sspace_3_5, ncol = 3, byrow = TRUE), 1, log_multinom_coef, sumx = 5), 3, 5, 10)
///
#[extendr]
fn calc_multinom_probs(sar: Vec<f64>, logt: Vec<f64>, logc: Vec<f64>, d: u32, n: u32, nt: u32) -> Vec<f64> {

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
    fn sample_unit_simplexn;
    fn calc_multinom_probs;
    fn sspace_multinom;
}
