* This is an update to fix the installation failures for 1.2.1
  - The `cargo test` step isn't part of rextendr's standard build recipe (and Rust-level tests belong in CI, not the CRAN install path), so I removed it (sorry about that)
  - `R CMD check --as-cran` passes locally with `checking Rust compilation ... OK` and no home-directory writes.
