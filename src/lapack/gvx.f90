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

module lapack_gvx
  use lapack_kinds, only : sp, dp, ik
  implicit none
  private

  public :: lapack_sygvx, lapack_spgvx, lapack_hegvx, lapack_hpgvx

  !> Computes selected eigenvalues, and optionally, eigenvectors
  !> of a real generalized symmetric-definite eigenproblem, of the form
  !> A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
  !> and B are assumed to be symmetric and B is also positive definite.
  !> Eigenvalues and eigenvectors can be selected by specifying either a
  !> range of values or a range of indices for the desired eigenvalues.
  interface lapack_sygvx
    pure subroutine ssygvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il,  &
        & iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info)
      import :: ik, sp
      integer, parameter :: wp = sp
      real(wp), intent(inout) :: a(lda, *)
      real(wp), intent(inout) :: b(ldb, *)
      real(wp), intent(out) :: w(*)
      integer(ik), intent(in) :: itype
      character(len=1), intent(in) :: uplo
      real(wp), intent(out) :: z(ldz, *)
      real(wp), intent(in) :: vl
      real(wp), intent(in) :: vu
      integer(ik), intent(in) :: il
      integer(ik), intent(in) :: iu
      integer(ik), intent(out) :: m
      integer(ik), intent(out) :: ifail(*)
      real(wp), intent(in) :: abstol
      integer(ik), intent(out) :: info
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: range
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: ldz
      real(wp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
      integer(ik), intent(in) :: iwork(*)
    end subroutine ssygvx
    pure subroutine dsygvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il,  &
        & iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info)
      import :: ik, dp
      integer, parameter :: wp = dp
      real(wp), intent(inout) :: a(lda, *)
      real(wp), intent(inout) :: b(ldb, *)
      real(wp), intent(out) :: w(*)
      integer(ik), intent(in) :: itype
      character(len=1), intent(in) :: uplo
      real(wp), intent(out) :: z(ldz, *)
      real(wp), intent(in) :: vl
      real(wp), intent(in) :: vu
      integer(ik), intent(in) :: il
      integer(ik), intent(in) :: iu
      integer(ik), intent(out) :: m
      integer(ik), intent(out) :: ifail(*)
      real(wp), intent(in) :: abstol
      integer(ik), intent(out) :: info
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: range
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: ldb
      integer(ik), intent(in) :: ldz
      real(wp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
      integer(ik), intent(in) :: iwork(*)
    end subroutine dsygvx
  end interface lapack_sygvx

  interface lapack_hegvx
    pure subroutine chegvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, &
        & il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info)
      import :: sp
      complex(sp), intent(inout) :: a(lda, *)
      complex(sp), intent(inout) :: b(ldb, *)
      real(sp), intent(out) :: w(*)
      integer, intent(in) :: itype
      character(len=1), intent(in) :: uplo
      complex(sp), intent(out) :: z(ldz, *)
      real(sp), intent(in) :: vl
      real(sp), intent(in) :: vu
      integer, intent(in) :: il
      integer, intent(in) :: iu
      integer, intent(out) :: m
      integer, intent(out) :: ifail(*)
      real(sp), intent(in) :: abstol
      integer, intent(out) :: info
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: range
      integer, intent(in) :: n
      integer, intent(in) :: lda
      integer, intent(in) :: ldb
      integer, intent(in) :: ldz
      complex(sp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      real(sp), intent(in) :: rwork(*)
      integer, intent(in) :: iwork(*)
    end subroutine chegvx
    pure subroutine zhegvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, &
        & il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info)
      import :: dp
      complex(dp), intent(inout) :: a(lda, *)
      complex(dp), intent(inout) :: b(ldb, *)
      real(dp), intent(out) :: w(*)
      integer, intent(in) :: itype
      character(len=1), intent(in) :: uplo
      complex(dp), intent(out) :: z(ldz, *)
      real(dp), intent(in) :: vl
      real(dp), intent(in) :: vu
      integer, intent(in) :: il
      integer, intent(in) :: iu
      integer, intent(out) :: m
      integer, intent(out) :: ifail(*)
      real(dp), intent(in) :: abstol
      integer, intent(out) :: info
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: range
      integer, intent(in) :: n
      integer, intent(in) :: lda
      integer, intent(in) :: ldb
      integer, intent(in) :: ldz
      complex(dp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      real(dp), intent(in) :: rwork(*)
      integer, intent(in) :: iwork(*)
    end subroutine zhegvx
  end interface lapack_hegvx

  interface lapack_spgvx
    pure subroutine sspgvx(itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, &
        & abstol, m, w, z, ldz, work, iwork, ifail, info)
      import :: sp
      real(sp), intent(inout) :: ap(*)
      real(sp), intent(inout) :: bp(*)
      real(sp), intent(out) :: w(*)
      integer, intent(in) :: itype
      character(len=1), intent(in) :: uplo
      real(sp), intent(out) :: z(ldz, *)
      real(sp), intent(in) :: vl
      real(sp), intent(in) :: vu
      integer, intent(in) :: il
      integer, intent(in) :: iu
      integer, intent(out) :: m
      integer, intent(out) :: ifail(*)
      real(sp), intent(in) :: abstol
      integer, intent(out) :: info
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: range
      integer, intent(in) :: n
      integer, intent(in) :: ldz
      real(sp), intent(in) :: work(*)
      integer, intent(in) :: iwork(*)
    end subroutine sspgvx
    pure subroutine dspgvx(itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, &
        & abstol, m, w, z, ldz, work, iwork, ifail, info)
      import :: dp
      real(dp), intent(inout) :: ap(*)
      real(dp), intent(inout) :: bp(*)
      real(dp), intent(out) :: w(*)
      integer, intent(in) :: itype
      character(len=1), intent(in) :: uplo
      real(dp), intent(out) :: z(ldz, *)
      real(dp), intent(in) :: vl
      real(dp), intent(in) :: vu
      integer, intent(in) :: il
      integer, intent(in) :: iu
      integer, intent(out) :: m
      integer, intent(out) :: ifail(*)
      real(dp), intent(in) :: abstol
      integer, intent(out) :: info
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: range
      integer, intent(in) :: n
      integer, intent(in) :: ldz
      real(dp), intent(in) :: work(*)
      integer, intent(in) :: iwork(*)
    end subroutine dspgvx
  end interface lapack_spgvx

  interface lapack_hpgvx
    pure subroutine chpgvx(itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, &
        & abstol, m, w, z, ldz, work, rwork, iwork, ifail, info)
      import :: sp
      complex(sp), intent(inout) :: ap(*)
      complex(sp), intent(inout) :: bp(*)
      real(sp), intent(out) :: w(*)
      integer, intent(in) :: itype
      character(len=1), intent(in) :: uplo
      complex(sp), intent(out) :: z(ldz, *)
      real(sp), intent(in) :: vl
      real(sp), intent(in) :: vu
      integer, intent(in) :: il
      integer, intent(in) :: iu
      integer, intent(out) :: m
      integer, intent(out) :: ifail(*)
      real(sp), intent(in) :: abstol
      integer, intent(out) :: info
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: range
      integer, intent(in) :: n
      integer, intent(in) :: ldz
      complex(sp), intent(in) :: work(*)
      real(sp), intent(in) :: rwork(*)
      integer, intent(in) :: iwork(*)
    end subroutine chpgvx
    pure subroutine zhpgvx(itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, &
        & abstol, m, w, z, ldz, work, rwork, iwork, ifail, info)
      import :: dp
      complex(dp), intent(inout) :: ap(*)
      complex(dp), intent(inout) :: bp(*)
      real(dp), intent(out) :: w(*)
      integer, intent(in) :: itype
      character(len=1), intent(in) :: uplo
      complex(dp), intent(out) :: z(ldz, *)
      real(dp), intent(in) :: vl
      real(dp), intent(in) :: vu
      integer, intent(in) :: il
      integer, intent(in) :: iu
      integer, intent(out) :: m
      integer, intent(out) :: ifail(*)
      real(dp), intent(in) :: abstol
      integer, intent(out) :: info
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: range
      integer, intent(in) :: n
      integer, intent(in) :: ldz
      complex(dp), intent(in) :: work(*)
      real(dp), intent(in) :: rwork(*)
      integer, intent(in) :: iwork(*)
    end subroutine zhpgvx
  end interface lapack_hpgvx

end module lapack_gvx
