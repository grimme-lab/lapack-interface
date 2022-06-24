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

!> Reduces a real symmetric-definite generalized eigenproblem to standard form.
module lapack_gst
  use lapack_kinds, only : sp, dp, ik
  implicit none
  private

  public :: lapack_sygst


  !> Reduces a real symmetric-definite generalized eigenproblem to standard form.
  !>
  !> If ITYPE = 1, the problem is A*x = lambda*B*x,
  !> and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
  !>
  !> If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
  !> B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.
  !>
  !> B must have been previously factorized as U**T*U or L*L**T by POTRF.
  interface lapack_sygst
    pure subroutine ssygst(itype, uplo, n, a, lda, b, ldb, info)
      import :: sp, ik
      real(sp), intent(inout) :: a(lda, *)
      real(sp), intent(in) :: b(ldb, *)
      integer(ik), intent(in) :: itype
      character(len=1), intent(in) :: uplo
      integer(ik), intent(out) :: info
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: ldb
    end subroutine ssygst
    pure subroutine dsygst(itype, uplo, n, a, lda, b, ldb, info)
      import :: dp, ik
      real(dp), intent(inout) :: a(lda, *)
      real(dp), intent(in) :: b(ldb, *)
      integer(ik), intent(in) :: itype
      character(len=1), intent(in) :: uplo
      integer(ik), intent(out) :: info
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: ldb
    end subroutine dsygst
    module procedure :: lapack_ssygst
    module procedure :: lapack_dsygst
  end interface lapack_sygst

contains

  pure subroutine lapack_ssygst(amat, bmat, info, itype, uplo)
    real(sp), intent(inout) :: amat(:, :)
    real(sp), intent(in) :: bmat(:, :)
    integer(ik), intent(in), optional :: itype
    character(len=1), intent(in), optional :: uplo
    integer(ik), intent(out) :: info
    character(len=1) :: ula
    integer(ik) :: ita, n, lda, ldb

    ita = 1
    if(present(itype)) ita = itype
    ula = 'u'
    if(present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    n = size(amat, 2)
    call lapack_sygst(ita, ula, n, amat, lda, bmat, ldb, info)
  end subroutine lapack_ssygst

  pure subroutine lapack_dsygst(amat, bmat, info, itype, uplo)
    real(dp), intent(inout) :: amat(:, :)
    real(dp), intent(in) :: bmat(:, :)
    integer(ik), intent(in), optional :: itype
    character(len=1), intent(in), optional :: uplo
    integer(ik), intent(out) :: info
    character(len=1) :: ula
    integer(ik) :: ita, n, lda, ldb

    ita = 1
    if(present(itype)) ita = itype
    ula = 'u'
    if(present(uplo)) ula = uplo
    lda = max(1, size(amat, 1))
    ldb = max(1, size(bmat, 1))
    n = size(amat, 2)
    call lapack_sygst(ita, ula, n, amat, lda, bmat, ldb, info)
  end subroutine lapack_dsygst

end module lapack_gst
