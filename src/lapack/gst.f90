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

  public :: lapack_sygst, lapack_spgst, lapack_hegst, lapack_hpgst


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

  !> reduces a complex Hermitian-definite generalized
  !> eigenproblem to standard form.
  !>
  !> If ITYPE = 1, the problem is A*x = lambda*B*x,
  !> and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H)
  !>
  !> If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
  !> B*A*x = lambda*x, and A is overwritten by U*A*U**H or L**H*A*L.
  !>
  !> B must have been previously factorized as U**H*U or L*L**H by CPOTRF.
  interface lapack_hegst
    pure subroutine chegst(itype, uplo, n, a, lda, b, ldb, info)
      import :: sp
      complex(sp), intent(inout) :: a(lda, *)
      complex(sp), intent(in) :: b(ldb, *)
      integer, intent(in) :: itype
      character(len=1), intent(in) :: uplo
      integer, intent(out) :: info
      integer, intent(in) :: n
      integer, intent(in) :: lda
      integer, intent(in) :: ldb
    end subroutine chegst
    pure subroutine zhegst(itype, uplo, n, a, lda, b, ldb, info)
      import :: dp
      complex(dp), intent(inout) :: a(lda, *)
      complex(dp), intent(in) :: b(ldb, *)
      integer, intent(in) :: itype
      character(len=1), intent(in) :: uplo
      integer, intent(out) :: info
      integer, intent(in) :: n
      integer, intent(in) :: lda
      integer, intent(in) :: ldb
    end subroutine zhegst
  end interface lapack_hegst

  !> reduces a real symmetric-definite generalized eigenproblem
  !> to standard form, using packed storage.
  !>
  !> If ITYPE = 1, the problem is A*x = lambda*B*x,
  !> and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
  !>
  !> If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
  !> B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.
  !>
  !> B must have been previously factorized as U**T*U or L*L**T by DPPTRF.
  interface lapack_spgst
    pure subroutine sspgst(itype, uplo, n, ap, bp, info)
      import :: sp
      real(sp), intent(inout) :: ap(*)
      real(sp), intent(in) :: bp(*)
      integer, intent(in) :: itype
      character(len=1), intent(in) :: uplo
      integer, intent(out) :: info
      integer, intent(in) :: n
    end subroutine sspgst
    pure subroutine dspgst(itype, uplo, n, ap, bp, info)
      import :: dp
      real(dp), intent(inout) :: ap(*)
      real(dp), intent(in) :: bp(*)
      integer, intent(in) :: itype
      character(len=1), intent(in) :: uplo
      integer, intent(out) :: info
      integer, intent(in) :: n
    end subroutine dspgst
  end interface lapack_spgst

  !> reduces a complex Hermitian-definite generalized
  !> eigenproblem to standard form, using packed storage.
  !>
  !> If ITYPE = 1, the problem is A*x = lambda*B*x,
  !> and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H)
  !>
  !> If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
  !> B*A*x = lambda*x, and A is overwritten by U*A*U**H or L**H*A*L.
  !>
  !> B must have been previously factorized as U**H*U or L*L**H by CPPTRF.
  interface lapack_hpgst
    pure subroutine chpgst(itype, uplo, n, ap, bp, info)
      import :: sp
      complex(sp), intent(inout) :: ap(*)
      complex(sp), intent(in) :: bp(*)
      integer, intent(in) :: itype
      character(len=1), intent(in) :: uplo
      integer, intent(out) :: info
      integer, intent(in) :: n
    end subroutine chpgst
    pure subroutine zhpgst(itype, uplo, n, ap, bp, info)
      import :: dp
      complex(dp), intent(inout) :: ap(*)
      complex(dp), intent(in) :: bp(*)
      integer, intent(in) :: itype
      character(len=1), intent(in) :: uplo
      integer, intent(out) :: info
      integer, intent(in) :: n
    end subroutine zhpgst
  end interface lapack_hpgst

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
