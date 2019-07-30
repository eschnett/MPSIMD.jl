# MPSIMD

Arbitrary precision integer arithmetic with an efficient (?) SIMD
representation

[![Build Status](https://dev.azure.com/schnetter/MPSIMD.jl/_apis/build/status/eschnett.MPSIMD.jl?branchName=master)](https://dev.azure.com/schnetter/MPSIMD.jl/_build/latest?definitionId=7&branchName=master)
TODO [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3355421.svg)](https://doi.org/10.5281/zenodo.3355421)



After some experimentation I am giving up on this project for two reasons:

* Using SIMD code in Julia leads to crashes too often. In particular I
  am unable to run benchmarks.

* Intel CPUs have a scalar multiplication instruction that multiplies
  two 64-bit integers and produces a 128-bit integer. SIMD
  instructions can at most multiple 32-bit integers, producing 32-bit
  integers. Emulating the scalar multiplication instruction with SIMD
  code requires many SIMD operations, and likely won't lead to a
  benefit even with AVX512 instructions.
