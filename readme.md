

# fortran-gmm-EM

This repository contains an extremely efficient implementation of the [expectation-maximization (EM) algorithm](https://en.wikipedia.org/wiki/Expectation%E2%80%93maximization_algorithm) for Gaussian mixtures. 
BLAS is used wherever possible for maximum efficiency.
The code is also parallelized using OpenMP. For maximum efficiency a parallel version of BLAS should be used, such as for example [OpenBLAS](https://www.openblas.net/) or intel MKL.

## Compiling the code

To compile the example you need a fortran compiler as well as an installation of Blas/LAPACK. 
Using CMake the code can be compiled with the following commands.
Two CMake flags are provided, that allow to compile with or without OpenMP parallelization and with the intel or gnu compiler.

```bash
mkdir build
cd build
cmake -DOPENMP=ON -DINTEL=OFF .. # compile with OpenMP parallelization and gfortran
make
```


## Using the code

An example on how to use the code can be found in `src/emExample.f90`. 
With the derived type `gaussian` handling of Gaussian mixtures (GM) is made easy as a GM can be represented as an array of `type(gaussian)`. 
The `gaussian` type contains mean, covariance and also the respective weight of the Gaussian in the mixture. 
Further it also containes the inverse of the covariance matrix and also their Cholesky decomposition. These quantities are usually not required by the user directly but they are used to calculate the probability density or when the GM is sampled. 
It is therefore important to call `gaussianFromCovariance` after manually changing the covariance matrix, so that these values can be recomputed. 

A Gaussian mixture can be fit with the following command:
```fortran
call fit(dim, nSamples, samples, nGaussians, gmmFit, eps)
```
The arguments are the following:

| Argument      | Description        | dimension  |
| ------------- | :----------------- | :--------- |
| dim           | dimension of the data | scalar |
| nSamples      | number of training samples | scalar |
| samples       | training samples  | (dim, nSamples) | 
| nGaussians    | number of Gaussians in the Gaussian mixture | scalar |
| gmmFit        | array of `type(Gaussian)` containing initial guess on input and fitting result on output | (nGaussians) |
| eps           | change in log likelihood per sample until which fit is done | scalar |

The gaussian mixtures can be sampled using:
```fortran
call sampleGMM(dim, nGaussians, gmm, sample)
```

The log PDF of the Gaussian mixture at coordinates `x` can be computed using:
```fortran
logPDF = gmmLogPdf(dimension, nGaussians, gaussians, x) 
```

The code also contains subroutines to write and read Gaussian mixtures from disk:
```fortran
call writeGMM('example-gmm.txt', dimension, nGaussians, gaussians)
call readGMM('example-gmm.txt', dimension, nGaussians, gaussians)
```


Please contact me regarding any bugs you encounter using this code. 

