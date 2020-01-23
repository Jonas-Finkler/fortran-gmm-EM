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


module gmmEM
    use precision
    use random

    implicit none

    real(dp) :: regCovar = 0._dp! 1.e-8_dp !todo: these thigs should be parameters.
    real(dp) :: initVar = 1.0e-2_dp
    real(dp), parameter :: PI = 4 * atan (1.0_dp)
    real(dp), parameter :: log2PI = log(2._dp * PI)

    type gaussian
        real(dp) :: weight
        real(dp), allocatable :: mean(:), covariance(:, :), precision(:, :), cholesky(:,:)
        real(dp) :: logDet
    end type gaussian

contains

    ! runs the EM algorithm until change in loglikelihood per sample is smaller than eps
    subroutine fit(dimension, nSamples, samples, nGaussians, gaussians, eps)
        integer, intent(in) :: dimension, nSamples, nGaussians
        real(dp), intent(in) :: samples(dimension, nSamples), eps
        type(gaussian), intent(inout) :: gaussians(nGaussians)
        real(dp) :: ll, lastll
        real(dp), allocatable :: sampleWeights(:, :), samplesT(:, :)
        integer :: i, j
        real(dp) :: st, et

        allocate(sampleWeights(nSamples, nGaussians))

        call E(dimension, nSamples, samples, nGaussians, gaussians, sampleWeights, ll)

        lastll = -huge(ll)

        do while((ll - lastll) / nSamples > eps)
    !        call cpu_time(st)
            call M(dimension, nSamples, samples, nGaussians, gaussians, sampleWeights)
    !        call cpu_time(et)
    !        write(*,*) "M", et-st
            lastll = ll
    !        call cpu_time(st)
            call E(dimension, nSamples, samples, nGaussians, gaussians, sampleWeights, ll)
    !        call cpu_time(et)
    !        write(*,*) "E", et-st
            write(*, *) ll, (ll - lastll) / nSamples
        end do
        !write(*, *) "weights", (gaussians(j)%weight, j = 1, nGaussians)

        deallocate(sampleWeights)
    end subroutine fit


    subroutine E(dimension, nSamples, samples, nGaussians, gaussians, sampleWeights, llh)
        integer, intent(in) :: dimension, nSamples, nGaussians
        real(dp), intent(in) :: samples(dimension, nSamples)
        type(gaussian), intent(inout) :: gaussians(nGaussians)
        real(dp), intent(out) :: sampleWeights(nSamples, nGaussians)
        real(dp), intent(out) :: llh
        real(dp) :: norm(nSamples)
        integer :: i, j
        real(dp), allocatable :: tmp(:,:)
        allocate(tmp(dimension, nSamples))
        do i = 1, nGaussians
            !$omp parallel do schedule(static, 1024)
            do j = 1, nSamples
                sampleWeights(j, i) = logPdf(dimension, gaussians(i), samples(:, j))
            end do
            !$omp end parallel do
        end do
        !$omp parallel do schedule(static, 1024)
        do i = 1, nSamples
            norm(i) = logSumExp(nGaussians, sampleWeights(i, :))
        end do
        !$omp end parallel do

        llh = sum(norm)

        do i = 1, nGaussians
            !$omp parallel do schedule(static, 1024)
            do j=1, nSamples
                sampleWeights(j, i) = exp(sampleWeights(j, i) - norm(j))
            end do
            !$omp end parallel do
        end do

    end subroutine

    subroutine M(dimension, nSamples, samples, nGaussians, gaussians, sampleWeights)
        integer, intent(in) :: dimension, nSamples, nGaussians
        real(dp), intent(in) :: samples(dimension, nSamples) ! using the transpose here improves vectorization and caching -> faster
        type(gaussian), intent(inout) :: gaussians(nGaussians)
        real(dp), intent(in) :: sampleWeights(nSamples, nGaussians)
        real(dp) :: norm(nSamples), v(dimension), A(dimension, dimension)
        integer :: i, j, c, k, l
        real(dp) :: tot
        real(dp), allocatable :: meanedSamples(:, :)
        real(dp), allocatable :: tmp(:)
        ! real(dp) :: tmp(nSamples)

        allocate(meanedSamples(dimension, nSamples))
        allocate(tmp(nSamples))

        do i = 1, nGaussians
            tot = sum(sampleWeights(:, i))
            gaussians(i)%weight = tot / nSamples

            call DGEMV('N', dimension, nSamples, 1.d0 / tot, samples, dimension, sampleWeights(:,i), &
                1, 0.d0, gaussians(i)%mean, 1)

            gaussians(i)%covariance(:,:) = 0._dp
            do j=1,dimension
                gaussians(i)%covariance(j,j) = regCovar
            end do
            tmp(:) = sqrt(sampleWeights(:,i) / tot)

            !omp parallel do schedule(static, 1024)
            do j=1,nSamples
                meanedSamples(:,j) = (samples(:,j) - gaussians(i)%mean(:)) * tmp(j)
            end do
            !omp end parallel do

            ! cov = cov + meanS @ meanS^T
            call DSYRK('L', 'N', dimension, nSamples, 1.d0, meanedSamples, dimension, 1.d0, gaussians(i)%covariance, dimension)

            ! todo: check if this is necessary
            do j=1,dimension
                do k=1,j-1
                    gaussians(i)%covariance(k,j) = gaussians(i)%covariance(j,k)
                end do
            end do

            call gaussianFromCovariance(dimension, gaussians(i))
        end do
        deallocate(meanedSamples)
    end subroutine

    subroutine logLikelihood(dimension, nSamples, samples, nGaussians, gaussians, llh)
        integer, intent(in) :: dimension, nSamples, nGaussians
        real(dp), intent(in) :: samples(dimension, nSamples)
        type(gaussian), intent(inout) :: gaussians(nGaussians)
        real(dp), intent(out) :: llh
        !real(dp) :: lPdf(nGaussians)
        integer :: i, j

        llh = 0._dp

        !$omp parallel do private(j), reduction(+:llh)
        do i = 1, nSamples
            llh = llh + gmmLogPdf(dimension, nGaussians, gaussians, samples(:,i))
        end do
        !$omp end parallel do

    end subroutine

    !
    subroutine initRandomGaussians(dimension, nGaussians, initVar, gaussians)
        integer, intent(in) :: dimension, nGaussians
        real(dp), intent(in) :: initVar
        type(gaussian), intent(out) :: gaussians(nGaussians)
        integer :: i, j

        do i = 1, nGaussians
            allocate(gaussians(i)%mean(dimension))
            allocate(gaussians(i)%covariance(dimension, dimension))
            allocate(gaussians(i)%precision(dimension, dimension))
            allocate(gaussians(i)%cholesky(dimension, dimension))
            gaussians(i)%weight = 1._dp / nGaussians
            gaussians(i)%covariance = 0._dp
            do j = 1, dimension
                gaussians(i)%mean(j) = random_normal() * initVar
                gaussians(i)%covariance(j, j) = 1._dp
            end do
            call gaussianFromCovariance(dimension, gaussians(i))
        end do
    end subroutine

    ! fills array with unit gaussians
    subroutine initUnitGaussians(dimension, nGaussians, gaussians)
        integer, intent(in) :: dimension, nGaussians
        type(gaussian), intent(out) :: gaussians(nGaussians)
        integer :: i, j

        do i = 1, nGaussians
            allocate(gaussians(i)%mean(dimension))
            allocate(gaussians(i)%covariance(dimension, dimension))
            allocate(gaussians(i)%precision(dimension, dimension))
            allocate(gaussians(i)%cholesky(dimension, dimension))
            gaussians(i)%weight = 1._dp / nGaussians
            gaussians(i)%covariance = 0._dp
            do j = 1, dimension
                gaussians(i)%covariance(j, j) = 1._dp
            end do
            call gaussianFromCovariance(dimension, gaussians(i))
        end do
    end subroutine

    subroutine gaussianFromCovariance(dimension, g)
        integer, intent(in) :: dimension
        type(gaussian), intent(inout) :: g
        integer :: info
        real(dp) :: invChol(dimension, dimension)
        integer :: i, j, k
        g%cholesky = g%covariance

        ! cholesky decomposition
        call dpotrf('L', dimension, g%cholesky, dimension, info)
        if(info /= 0) stop "cholesky factorization failed!"
        invChol = g%cholesky
        ! invert triangular matrix
        call dtrtri('L', 'N', dimension, invChol, dimension, info)
        if(info /= 0) stop "triangular inversion failed!"

        ! precision = invChol**T * invChol
        g%precision = 0._dp
        do i = 1, dimension
            do j = i, dimension
                g%precision(i, j) = sum(invChol(j:, i) * invChol(j:, j))
                g%precision(j, i) = g%precision(i, j) ! important, because full matrix is used in function logPdf
            end do
        end do

        g%logDet = 0._dp
        do i = 1, dimension
            g%logDet = g%logDet + log(g%cholesky(i, i))
        end do
        g%logDet = g%logDet * 2._dp
    end subroutine gaussianFromCovariance


    function gmmLogPdf(dimension, nGaussians, gaussians, x) result(p)
        integer, intent(in) :: dimension, nGaussians
        type(gaussian), intent(in) :: gaussians(nGaussians)
        real(dp), intent(in) :: x(dimension)
        real(dp) :: p
        real(dp) :: lPdf(nGaussians)
        integer :: i, n

        n = 0
        do i=1,nGaussians
            if (gaussians(i)%weight > 0._dp) then
                n = n+1
                lPdf(n) = logPdf(dimension,gaussians(i),x)
            end if
        end do
        p = logSumExp(n, lPdf(:n))

    end function gmmLogPdf


    function logPdf(dimension, g, x) result(p)
        integer, intent(in) :: dimension
        type(gaussian), intent(in) :: g
        real(dp), intent(in) :: x(dimension)
        real(dp) :: p
        real(dp) :: d(dimension), s
        integer :: i, j

        d = x - g%mean
        s = 0._dp
        do i = 1, dimension
            s = s + sum(d(:) * g%precision(:, i)) * d(i)
        end do

        p = log(g%weight) - 0.5_dp * (s + dimension * log2PI + g%logDet)

        if(p > huge(p)) then
            write(*, *) "info: ", s, g%logDet
            stop "inf at logPdf"
        end if
    end function logPdf

    function logSumExp(n, x) result(s)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(:) !implicit shape array to avoid creation of temporary array
        real(dp) :: s, m
        integer :: i

        m = maxval(x)
        s = 0._dp
        do i = 1, n
            s = s + exp(x(i) - m)
        end do
        s = m + log(s)
    end function logSumExp

    subroutine writeGMM(file, dimension, nGaussians, gaussians)
        character(len=*), intent(in) :: file
        integer, intent(in) :: nGaussians, dimension
        type(gaussian), intent(in) :: gaussians(nGaussians)
        integer :: i,j

        open(45, file = file, action = "write")
        write(45,*) nGaussians,dimension
        do i=1,nGaussians
            write(45,*) gaussians(i)%weight
            write(45,*) gaussians(i)%mean
            do j=1,dimension
                write(45,*) gaussians(i)%covariance(:,j)
            end do
        end do
        close(45)
    end subroutine writeGMM

    subroutine readGMM(file, dimension, nGaussians, gaussians)
        character(len=*), intent(in) :: file
        integer, intent(in) :: nGaussians, dimension
        type(gaussian), intent(inout) :: gaussians(nGaussians)
        integer :: i,j
        integer :: nG, dim, unit

        open(newUnit = unit, file = file, action = "read")
        read(unit,*) nG,dim

        if(nG /= nGaussians .or. dim /= dimension) stop "gmm does not fit required parameters"

        do i=1,nGaussians
            allocate(gaussians(i)%mean(dimension))
            allocate(gaussians(i)%covariance(dimension, dimension))
            allocate(gaussians(i)%precision(dimension, dimension))
            allocate(gaussians(i)%cholesky(dimension, dimension))
            read(unit,*) gaussians(i)%weight
            read(unit,*) gaussians(i)%mean
            do j=1,dimension
                read(unit,*) gaussians(i)%covariance(:,j)
            end do
            call gaussianFromCovariance(dimension, gaussians(i))
        end do

        close(unit)
    end subroutine readGMM

    subroutine sampleGMM(dimension, nGaussians, gaussians, sample)
        integer, intent(in) :: dimension, nGaussians
        type(gaussian), intent(in) :: gaussians(nGaussians)
        real(dp), intent(out) :: sample(dimension)
        integer :: i
        real(dp) :: rnd, c

        call random_number(rnd)
        c = 0._dp
        i = 0
        do while(c < rnd .and. i < nGaussians)
            i = i + 1
            c = c + gaussians(i)%weight
        end do
        call sampleGaussian(dimension, gaussians(i), sample)
    end subroutine sampleGMM

    subroutine sampleGaussian(dimension, g, sample)
        integer, intent(in) :: dimension
        type(gaussian), intent(in) :: g
        real(dp), intent(out) :: sample(dimension)
        integer :: i
        real(dp) :: rnd(dimension)

        do i=1,dimension
            rnd(i) = random_normal()
        end do

        do i=1,dimension
            sample(i) = g%mean(i) + sum(g%cholesky(i,:i) * rnd(:i))
        end do

    end subroutine sampleGaussian
end module gmmEM
