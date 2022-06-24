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

module lapack_gvd
  use lapack_kinds, only : sp, dp, ik
  implicit none
  private

  public :: lapack_sygvd, lapack_spgvd, lapack_hegvd, lapack_hpgvd

  !> Computes all the eigenvalues, and optionally, the eigenvectors
  !> of a real generalized symmetric-definite eigenproblem, of the form
  !> A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and
  !> B are assumed to be symmetric and B is also positive definite.
  !> If eigenvectors are desired, it uses a divide and conquer algorithm.
  interface lapack_sygvd
    pure subroutine ssygvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork,    &
        & iwork, liwork, info)
      import :: ik, sp
      integer, parameter :: wp = sp
      real(wp), intent(inout) :: a(lda, *)
      real(wp), intent(inout) :: b(ldb, *)
      real(wp), intent(out) :: w(*)
      integer(ik), intent(in) :: itype
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: uplo
      integer(ik), intent(out) :: info
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: ldb
      real(wp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
      integer(ik), intent(inout) :: iwork(*)
      integer(ik), intent(in) :: liwork
    end subroutine ssygvd
    pure subroutine dsygvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork,    &
        & iwork, liwork, info)
      import :: ik, dp
      integer, parameter :: wp = dp
      real(wp), intent(inout) :: a(lda, *)
      real(wp), intent(inout) :: b(ldb, *)
      real(wp), intent(out) :: w(*)
      integer(ik), intent(in) :: itype
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: uplo
      integer(ik), intent(out) :: info
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: ldb
      real(wp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
      integer(ik), intent(inout) :: iwork(*)
      integer(ik), intent(in) :: liwork
    end subroutine dsygvd
  end interface lapack_sygvd

  interface lapack_hegvd
    pure subroutine chegvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, &
        & rwork, lrwork, iwork, liwork, info)
      import :: sp
      complex(sp), intent(inout) :: a(lda, *)
      complex(sp), intent(inout) :: b(ldb, *)
      real(sp), intent(out) :: w(*)
      integer, intent(in) :: itype
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: uplo
      integer, intent(out) :: info
      integer, intent(in) :: n
      integer, intent(in) :: lda
      integer, intent(in) :: ldb
      complex(sp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      real(sp), intent(inout) :: rwork(*)
      integer, intent(in) :: lrwork
      integer, intent(inout) :: iwork(*)
      integer, intent(in) :: liwork
    end subroutine chegvd
    pure subroutine zhegvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, &
        & rwork, lrwork, iwork, liwork, info)
      import :: dp
      complex(dp), intent(inout) :: a(lda, *)
      complex(dp), intent(inout) :: b(ldb, *)
      real(dp), intent(out) :: w(*)
      integer, intent(in) :: itype
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: uplo
      integer, intent(out) :: info
      integer, intent(in) :: n
      integer, intent(in) :: lda
      integer, intent(in) :: ldb
      complex(dp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      real(dp), intent(inout) :: rwork(*)
      integer, intent(in) :: lrwork
      integer, intent(inout) :: iwork(*)
      integer, intent(in) :: liwork
    end subroutine zhegvd
  end interface lapack_hegvd

  interface lapack_spgvd
    pure subroutine sspgvd(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, lwork, &
        & iwork, liwork, info)
      import :: sp
      real(sp), intent(inout) :: ap(*)
      real(sp), intent(inout) :: bp(*)
      real(sp), intent(out) :: w(*)
      integer, intent(in) :: itype
      character(len=1), intent(in) :: uplo
      real(sp), intent(out) :: z(ldz, *)
      integer, intent(out) :: info
      character(len=1), intent(in) :: jobz
      integer, intent(in) :: n
      integer, intent(in) :: ldz
      real(sp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(inout) :: iwork(*)
      integer, intent(in) :: liwork
    end subroutine sspgvd
    pure subroutine dspgvd(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, lwork, &
        & iwork, liwork, info)
      import :: dp
      real(dp), intent(inout) :: ap(*)
      real(dp), intent(inout) :: bp(*)
      real(dp), intent(out) :: w(*)
      integer, intent(in) :: itype
      character(len=1), intent(in) :: uplo
      real(dp), intent(out) :: z(ldz, *)
      integer, intent(out) :: info
      character(len=1), intent(in) :: jobz
      integer, intent(in) :: n
      integer, intent(in) :: ldz
      real(dp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(inout) :: iwork(*)
      integer, intent(in) :: liwork
    end subroutine dspgvd
  end interface lapack_spgvd

  interface lapack_hpgvd
    pure subroutine chpgvd(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, lwork, &
        & rwork, lrwork, iwork, liwork, info)
      import :: sp
      complex(sp), intent(inout) :: ap(*)
      complex(sp), intent(inout) :: bp(*)
      real(sp), intent(out) :: w(*)
      integer, intent(in) :: itype
      character(len=1), intent(in) :: uplo
      complex(sp), intent(out) :: z(ldz, *)
      integer, intent(out) :: info
      character(len=1), intent(in) :: jobz
      integer, intent(in) :: n
      integer, intent(in) :: ldz
      complex(sp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      real(sp), intent(inout) :: rwork(*)
      integer, intent(in) :: lrwork
      integer, intent(inout) :: iwork(*)
      integer, intent(in) :: liwork
    end subroutine chpgvd
    pure subroutine zhpgvd(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, lwork, &
        & rwork, lrwork, iwork, liwork, info)
      import :: dp
      complex(dp), intent(inout) :: ap(*)
      complex(dp), intent(inout) :: bp(*)
      real(dp), intent(out) :: w(*)
      integer, intent(in) :: itype
      character(len=1), intent(in) :: uplo
      complex(dp), intent(out) :: z(ldz, *)
      integer, intent(out) :: info
      character(len=1), intent(in) :: jobz
      integer, intent(in) :: n
      integer, intent(in) :: ldz
      complex(dp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      real(dp), intent(inout) :: rwork(*)
      integer, intent(in) :: lrwork
      integer, intent(inout) :: iwork(*)
      integer, intent(in) :: liwork
    end subroutine zhpgvd
  end interface lapack_hpgvd

end module lapack_gvd
