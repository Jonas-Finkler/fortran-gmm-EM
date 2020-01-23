! Copyright (C) 2020 Jonas A. Finkler
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

program emExample
    use precision
    use gmmEM
    use random

    implicit none

    integer, parameter :: dim = 10
    integer, parameter :: nSamples = 1000000
    integer, parameter :: nGaussians = 3

    real(dp), allocatable :: samples(:,:)
    type(gaussian) :: gmm(nGaussians), gmmFit(nGaussians)
    integer :: i, j

    allocate(samples(dim,nSamples))

    ! use this routine to initialize the gaussian and allocate all the variables in the gaussian type
    call initUnitGaussians(dim, nGaussians, gmm)

    ! now we creates some ranoom gaussians
    do i=1,nGaussians
        gmm(i)%weight = 1._dp / nGaussians
        do j=1,dim
            gmm(i)%mean(j) = random_normal() * 0.1_dp
            gmm(i)%covariance(j,j) = random_normal()**2
        end do
        ! gaussianFromCovariance should be called everytime the cov. matrix is modified. It recalculates the Cholesky decomposition and the precision matrix.
        call gaussianFromCovariance(dim, gmm(i))
    end do

    ! generate samples from these gaussians
    do i=1, nSamples
        call sampleGMM(dim, nGaussians, gmm, samples(:,i))
    end do

    call initRandomGaussians(dim, nGaussians, 1.0e-2_dp, gmmFit)
    ! fit gmm to the samples
    call fit(dim, nSamples, samples, nGaussians, gmmFit, 1.e-10_dp)

    ! comapre the result with the original parameters. The order of the Gaussians is random, we therefore compare all the pairs.
    write(*,*)
    do i=1,nGaussians
        do j=1,nGaussians
            write(*,'(A,I3.3,A,I3.3,A)') '    ---    ', i, '~', j, '    ---'
            write(*,*) 'weight     ', (gmm(i)%weight - gmmFit(j)%weight)**2
            write(*,*) 'mean       ', sum((gmm(i)%mean - gmmFit(j)%mean)**2)
            write(*,*) 'covariance ', sum((gmm(i)%covariance - gmmFit(j)%covariance)**2)
        end do
    end do

end program emExample
