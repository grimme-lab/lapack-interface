! This file is part of lapack-interface.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module lapack_trf
  use lapack_kinds, only : sp, dp, ik
  use lapack_utils, only : lapack_env
  implicit none
  private

  public :: lapack_sytrf, lapack_getrf, lapack_potrf


  interface lapack_sytrf
    pure subroutine ssytrf(uplo, n, a, lda, ipiv, work, lwork, info)
      import :: sp, ik
      integer(ik), intent(in) :: lda
      real(sp), intent(inout) :: a(lda, *)
      character(len=1), intent(in) :: uplo
      integer(ik), intent(out) :: ipiv(*)
      integer(ik), intent(out) :: info
      integer(ik), intent(in) :: n
      real(sp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
    end subroutine ssytrf
    pure subroutine dsytrf(uplo, n, a, lda, ipiv, work, lwork, info)
      import :: dp, ik
      integer(ik), intent(in) :: lda
      real(dp), intent(inout) :: a(lda, *)
      character(len=1), intent(in) :: uplo
      integer(ik), intent(out) :: ipiv(*)
      integer(ik), intent(out) :: info
      integer(ik), intent(in) :: n
      real(dp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
    end subroutine dsytrf

    module procedure :: lapack_ssytrf
    module procedure :: lapack_dsytrf
  end interface lapack_sytrf


  !> Computes the Cholesky factorization of a real symmetric
  !> positive definite matrix A.
  !>
  !> The factorization has the form
  !>    A = U**T * U,  if UPLO = 'U', or
  !>    A = L  * L**T,  if UPLO = 'L',
  !> where U is an upper triangular matrix and L is lower triangular.
  !>
  !> This is the block version of the algorithm, calling Level 3 BLAS.
  interface lapack_potrf
    pure subroutine spotrf(uplo, n, a, lda, info)
      import :: sp, ik
      real(sp), intent(inout) :: a(lda, *)
      character(len=1), intent(in) :: uplo
      integer(ik), intent(out) :: info
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
    end subroutine spotrf
    pure subroutine dpotrf(uplo, n, a, lda, info)
      import :: dp, ik
      real(dp), intent(inout) :: a(lda, *)
      character(len=1), intent(in) :: uplo
      integer(ik), intent(out) :: info
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
    end subroutine dpotrf

    module procedure :: lapack_spotrf
    module procedure :: lapack_dpotrf
  end interface lapack_potrf


  !> Computes an LU factorization of a general M-by-N matrix A
  !> using partial pivoting with row interchanges.
  !>
  !> The factorization has the form
  !>    A = P * L * U
  !> where P is a permutation matrix, L is lower triangular with unit
  !> diagonal elements (lower trapezoidal if m > n), and U is upper
  !> triangular (upper trapezoidal if m < n).
  interface lapack_getrf
    pure subroutine sgetrf(m, n, a, lda, ipiv, info)
      import :: sp, ik
      real(sp), intent(inout) :: a(lda, *)
      integer(ik), intent(out) :: ipiv(*)
      integer(ik), intent(out) :: info
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
    end subroutine sgetrf
    pure subroutine dgetrf(m, n, a, lda, ipiv, info)
      import :: dp, ik
      real(dp), intent(inout) :: a(lda, *)
      integer(ik), intent(out) :: ipiv(*)
      integer(ik), intent(out) :: info
      integer(ik), intent(in) :: m
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
    end subroutine dgetrf

    module procedure :: lapack_sgetrf
    module procedure :: lapack_dgetrf
  end interface lapack_getrf


contains


  subroutine lapack_ssytrf(amat, ipiv, uplo, info)
    real(sp), intent(inout) :: amat(:, :)
    integer(ik), intent(out) :: ipiv(:)
    character(len=1), intent(in), optional :: uplo
    integer(ik), intent(out), optional :: info
    character(len=1) :: ula
    integer(ik) :: stat, n, lda, lwork
    real(sp), allocatable :: work(:)

    ula = 'u'
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    n = size(amat, 2)

    lwork = n * lapack_env(1_ik, 'SSYTRF', ula, n, -1_ik, -1_ik, -1_ik)
    allocate(work(lwork), stat=stat)
    if (stat == 0) then
      call lapack_sytrf(ula, n, amat, lda, ipiv, work, lwork, stat)
    else
      stat = -1000_ik
    end if
    if (present(info)) then
      info = stat
    else
      if (stat /= 0) error stop "[lapack] sytrf failed"
    end if
  end subroutine lapack_ssytrf

  subroutine lapack_dsytrf(amat, ipiv, uplo, info)
    real(dp), intent(inout) :: amat(:, :)
    integer(ik), intent(out) :: ipiv(:)
    character(len=1), intent(in), optional :: uplo
    integer(ik), intent(out), optional :: info
    character(len=1) :: ula
    integer(ik) :: stat, n, lda, lwork
    real(dp), allocatable :: work(:)

    ula = 'u'
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    n = size(amat, 2)

    lwork = n * lapack_env(1_ik, 'DSYTRF', ula, n, -1_ik, -1_ik, -1_ik)
    allocate(work(lwork), stat=stat)
    if (stat == 0) then
      call lapack_sytrf(ula, n, amat, lda, ipiv, work, lwork, stat)
    else
      stat = -1000_ik
    end if
    if (present(info)) then
      info = stat
    else
      if (stat /= 0) error stop "[lapack] sytrf failed"
    end if
  end subroutine lapack_dsytrf


  subroutine lapack_spotrf(amat, uplo, info)
    real(sp), intent(inout) :: amat(:, :)
    integer(ik), intent(out), optional :: info
    character(len=1), intent(in), optional :: uplo
    integer(ik) :: m, n, lda, stat
    character(len=1) :: ula

    ula = 'u'
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call lapack_potrf(ula, n, amat, lda, stat)
    if (present(info)) then
      info = stat
    else
      if (stat /= 0) error stop "[lapack] potrf failed"
    end if
  end subroutine lapack_spotrf

  subroutine lapack_dpotrf(amat, uplo, info)
    real(dp), intent(inout) :: amat(:, :)
    integer(ik), intent(out), optional :: info
    character(len=1), intent(in), optional :: uplo
    integer(ik) :: m, n, lda, stat
    character(len=1) :: ula

    ula = 'u'
    if (present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    n = size(amat, 2)
    call lapack_potrf(ula, n, amat, lda, stat)
    if (present(info)) then
      info = stat
    else
      if (stat /= 0) error stop "[lapack] potrf failed"
    end if
  end subroutine lapack_dpotrf


  subroutine lapack_sgetrf(amat, ipiv, info)
    real(sp), intent(inout) :: amat(:, :)
    integer(ik), intent(out) :: ipiv(:)
    integer(ik), intent(out), optional :: info
    integer(ik) :: m, n, lda, stat

    lda = max(1, size(amat, 1))
    m = size(amat, 1)
    n = size(amat, 2)
    call lapack_getrf(m, n, amat, lda, ipiv, stat)
    if (present(info)) then
      info = stat
    else
      if (stat /= 0) error stop "[lapack] getrf failed"
    end if
  end subroutine lapack_sgetrf


  subroutine lapack_dgetrf(amat, ipiv, info)
    real(dp), intent(inout) :: amat(:, :)
    integer(ik), intent(out) :: ipiv(:)
    integer(ik), intent(out), optional :: info
    integer(ik) :: m, n, lda, stat

    lda = max(1, size(amat, 1))
    m = size(amat, 1)
    n = size(amat, 2)
    call lapack_getrf(m, n, amat, lda, ipiv, stat)
    if (present(info)) then
      info = stat
    else
      if (stat /= 0) error stop "[lapack] getrf failed"
    end if
  end subroutine lapack_dgetrf

end module lapack_trf
