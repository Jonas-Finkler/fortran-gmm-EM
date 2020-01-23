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


module linalg
    use precision

    implicit none
    real(dp), parameter :: PI = 4 * atan (1.0_dp)

contains

    ! knuth shuffle from rosetta code
    subroutine shuffle(n, a)
        integer, intent(in) :: n
        integer, intent(inout) :: a(n)
        integer :: i, randpos, temp
        real :: r

        do i = n, 2, -1
            call random_number(r)
            randpos = int(r * i) + 1
            temp = a(randpos)
            a(randpos) = a(i)
            a(i) = temp
        end do

    end subroutine shuffle

    ! modified (numerically stable) Gram-Schmidt process
    subroutine gramSchmidt(d, n, nu, v, u)
        implicit none
        integer, intent(in) :: d ! dimensionality of the vectors
        integer, intent(in) :: n ! number of vectors that are created
        integer, intent(in) :: nu  ! number of vectors already orthogonal and in the first nu columns of u
        real(dp), intent(in) :: v(d, n) ! linearly independent vectors
        real(dp), intent(inout) :: u(d, n) ! the remaining n-nu cols will be filled with orthonormal vectors orthogonal to the first nu ones alredy present

        integer :: i, j

        do i=nu+1,n
          u(:,i) = v(:,i)
          do j=1,i-1
            u(:,i) = u(:,i) - u(:,j) * sum(u(:,i) * u(:,j)) / sum(u(:,j)**2)
          end do
          u(:,i) = u(:,i) / sqrt(sum(u(:,i)**2))
        end do

    end subroutine

    ! evecs(:,i) is eigenvector to eval(i), only needs lower half of matrix
    subroutine eigensystemDiag(n, m, evals, evecs)
        integer, intent(in) :: n
        real(dp), intent(in) :: m(n, n)
        real(dp), intent(out) :: evals(n), evecs(n, n)
        real(dp) :: A(n, n), W(n), LWORKQUERY(1)
        real(dp), allocatable :: WORK(:)
        integer :: LWORK, INFO, LDA

        A = m
        LDA = n

        call DSYEV('V', 'L', n, A, LDA, W, LWORKQUERY, -1, INFO)

        LWORK = int(LWORKQUERY(1))
        allocate(WORK(LWORK))
        call DSYEV('V', 'L', n, A, LDA, W, WORK, LWORK, INFO)
        deallocate(WORK)
        if(INFO/=0) stop "DSYEV failed"
        evecs = A
        evals = W
    end subroutine eigensystemDiag

    subroutine eigenvaluesDiag(n, m, evals)
        integer, intent(in) :: n
        real(dp), intent(in) :: m(n, n)
        real(dp), intent(out) :: evals(n)
        real(dp) :: A(n, n), W(n), LWORKQUERY(1)
        real(dp), allocatable :: WORK(:)
        integer :: LWORK, INFO, LDA

        A = m
        LDA = n

        call DSYEV('N', 'L', n, A, LDA, W, LWORKQUERY, -1, INFO)

        LWORK = int(LWORKQUERY(1))
        allocate(WORK(LWORK))
        call DSYEV('N', 'L', n, A, LDA, W, WORK, LWORK, INFO)
        deallocate(WORK)
        if(INFO/=0) stop "DSYEV failed"
        evals = W
    end subroutine eigenvaluesDiag

    function matMulVec3(mat, vec) result(res)
        implicit none
        real(dp), dimension(3, 3), intent(in) :: mat
        real(dp), dimension(3), intent(in) :: vec
        real(dp), dimension(3) :: res

        res(1) = sum(mat(1, :) * vec(:))
        res(2) = sum(mat(2, :) * vec(:))
        res(3) = sum(mat(3, :) * vec(:))
    end function

    function matMulVecNM(n, m, mat, vec) result(res)
        integer, intent(in) :: n, m
        real(dp), intent(in) :: mat(n,m), vec(m)
        real(dp) :: res(n)
        integer :: i

        do i=1,n
            res(i) = sum(mat(i,:) * vec(:))
        end do

    end function

    function cross(a, b) result(c)
        implicit none
        real(dp), dimension(3), intent(in) :: a, b
        real(dp), dimension(3) :: c

        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)
    end function

    subroutine inv3D(a1, a2, a3, b1, b2, b3)
        implicit none

        real(dp), intent(in), dimension(3) :: a1, a2, a3
        real(dp), intent(out), dimension(3) :: b1, b2, b3
        real(dp) :: det

        det = a1(1) * a2(2) * a3(3) + a1(2) * a2(3) * a3(1) + a1(3) * &
                a2(1) * a3(2) - a1(1) * a2(3) * a3(2) - a1(2) * a2(1) * a3(3) - a1(3) * a2(2) * a3(1)

        b1 = [a2(2) * a3(3) - a2(3) * a3(2), a1(3) * a3(2) - a1(2) * a3(3), a1(2) * a2(3) - a1(3) * a2(2)]
        b2 = [a2(3) * a3(1) - a2(1) * a3(3), a1(1) * a3(3) - a1(3) * a3(1), a1(3) * a2(1) - a1(1) * a2(3)]
        b3 = [a2(1) * a3(2) - a2(2) * a3(1), a1(2) * a3(1) - a1(1) * a3(2), a1(1) * a2(2) - a1(2) * a2(1)]

        b1 = b1 / det
        b2 = b2 / det
        b3 = b3 / det
    end subroutine

    subroutine inv3DM(A, B)
        implicit none
        real(dp), intent(in), dimension(3, 3) :: A
        real(dp), intent(out), dimension(3, 3) :: B
        real(dp), dimension(3) :: a1, a2, a3
        real(dp), dimension(3) :: b1, b2, b3
        real(dp) :: det

        a1 = A(1, :)
        a2 = A(2, :)
        a3 = A(3, :)

        det = a1(1) * a2(2) * a3(3) + a1(2) * a2(3) * a3(1) + a1(3) * &
                a2(1) * a3(2) - a1(1) * a2(3) * a3(2) - a1(2) * a2(1) * a3(3) - a1(3) * a2(2) * a3(1)

        b1 = [a2(2) * a3(3) - a2(3) * a3(2), a1(3) * a3(2) - a1(2) * a3(3), a1(2) * a2(3) - a1(3) * a2(2)]
        b2 = [a2(3) * a3(1) - a2(1) * a3(3), a1(1) * a3(3) - a1(3) * a3(1), a1(3) * a2(1) - a1(1) * a2(3)]
        b3 = [a2(1) * a3(2) - a2(2) * a3(1), a1(2) * a3(1) - a1(1) * a3(2), a1(1) * a2(2) - a1(2) * a2(1)]

        b1 = b1 / det
        b2 = b2 / det
        b3 = b3 / det

        B(1, :) = b1
        B(2, :) = b2
        B(3, :) = b3
    end subroutine


    function arrayCompare(a, b)
        logical :: arrayCompare
        integer, intent(in) :: a(:), b(:)
        integer :: i
        arrayCompare = .false.
        if(size(a) /= size(b)) return
        do i = 1, size(a)
            if(a(i)/=b(i)) return
        end do
        arrayCompare = .true.
    end function arrayCompare


end module linalg
