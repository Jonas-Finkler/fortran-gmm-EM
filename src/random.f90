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

module random
    use precision
    implicit none

contains

    ! uses the Box-Müller algorithm to generate normal distributed random numbers
    ! threadsafe -> can be used in omp parallelized sections
    function random_normal() result (rnor)
        real(dp), parameter :: PI = 4._dp * datan(1._dp)
        real(dp) :: rnor
        real(dp) :: u(2)
        real(dp), save :: cache
        logical, save :: has_cache = .false.
        ! each thread has its own cache
        !$omp threadprivate(cache, has_cache)

        if (has_cache) then
            has_cache = .false.
            rnor = cache
        else
            call random_number(u)
            rnor  = sqrt(-2.0 * log(u(1))) * sin(2 * PI * u(2))
            cache = sqrt(-2.0 * log(u(1))) * cos(2 * PI * u(2))
            has_cache = .true.
        end if
    end function random_normal

end module random
