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

module lapack_evx
  use lapack_kinds, only : sp, dp, ik
  implicit none
  private

  public :: lapack_syevx, lapack_spevx, lapack_heevx, lapack_hpevx


  !> Computes selected eigenvalues and, optionally, eigenvectors
  !> of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
  !> selected by specifying either a range of values or a range of indices
  !> for the desired eigenvalues.
  interface lapack_syevx
    pure subroutine ssyevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w,&
        & z, ldz, work, lwork, iwork, ifail, info)
      import :: sp, ik
      integer, parameter :: wp = sp
      real(wp), intent(inout) :: a(lda, *)
      real(wp), intent(out) :: w(*)
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
      integer(ik), intent(in) :: ldz
      real(wp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
      integer(ik), intent(in) :: iwork(*)
    end subroutine ssyevx
    pure subroutine dsyevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w,&
        & z, ldz, work, lwork, iwork, ifail, info)
      import :: dp, ik
      integer, parameter :: wp = dp
      real(wp), intent(inout) :: a(lda, *)
      real(wp), intent(out) :: w(*)
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
      integer(ik), intent(in) :: ldz
      real(wp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
      integer(ik), intent(in) :: iwork(*)
    end subroutine dsyevx
  end interface lapack_syevx

  interface lapack_spevx
    pure subroutine sspevx(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, &
        & m, w, z, ldz, work, iwork, ifail, info)
      import :: sp, ik
      real(sp), intent(inout) :: ap(*)
      real(sp), intent(out) :: w(*)
      character(len=1), intent(in) :: uplo
      real(sp), intent(out) :: z(ldz, *)
      real(sp), intent(in) :: vl
      real(sp), intent(in) :: vu
      integer(ik), intent(in) :: il
      integer(ik), intent(in) :: iu
      integer(ik), intent(out) :: m
      integer(ik), intent(out) :: ifail(*)
      real(sp), intent(in) :: abstol
      integer(ik), intent(out) :: info
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: range
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: ldz
      real(sp), intent(in) :: work(*)
      integer(ik), intent(in) :: iwork(*)
    end subroutine sspevx
    pure subroutine dspevx(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, &
        & m, w, z, ldz, work, iwork, ifail, info)
      import :: dp, ik
      real(dp), intent(inout) :: ap(*)
      real(dp), intent(out) :: w(*)
      character(len=1), intent(in) :: uplo
      real(dp), intent(out) :: z(ldz, *)
      real(dp), intent(in) :: vl
      real(dp), intent(in) :: vu
      integer(ik), intent(in) :: il
      integer(ik), intent(in) :: iu
      integer(ik), intent(out) :: m
      integer(ik), intent(out) :: ifail(*)
      real(dp), intent(in) :: abstol
      integer(ik), intent(out) :: info
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: range
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: ldz
      real(dp), intent(in) :: work(*)
      integer(ik), intent(in) :: iwork(*)
    end subroutine dspevx
  end interface lapack_spevx

  interface lapack_heevx
    pure subroutine cheevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, &
        & abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info)
      import :: sp, ik
      complex(sp), intent(inout) :: a(lda, *)
      real(sp), intent(out) :: w(*)
      character(len=1), intent(in) :: uplo
      complex(sp), intent(out) :: z(ldz, *)
      real(sp), intent(in) :: vl
      real(sp), intent(in) :: vu
      integer(ik), intent(in) :: il
      integer(ik), intent(in) :: iu
      integer(ik), intent(out) :: m
      integer(ik), intent(out) :: ifail(*)
      real(sp), intent(in) :: abstol
      integer(ik), intent(out) :: info
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: range
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: ldz
      complex(sp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
      real(sp), intent(in) :: rwork(*)
      integer(ik), intent(in) :: iwork(*)
    end subroutine cheevx
    pure subroutine zheevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, &
        & abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info)
      import :: dp, ik
      complex(dp), intent(inout) :: a(lda, *)
      real(dp), intent(out) :: w(*)
      character(len=1), intent(in) :: uplo
      complex(dp), intent(out) :: z(ldz, *)
      real(dp), intent(in) :: vl
      real(dp), intent(in) :: vu
      integer(ik), intent(in) :: il
      integer(ik), intent(in) :: iu
      integer(ik), intent(out) :: m
      integer(ik), intent(out) :: ifail(*)
      real(dp), intent(in) :: abstol
      integer(ik), intent(out) :: info
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: range
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: ldz
      complex(dp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
      real(dp), intent(in) :: rwork(*)
      integer(ik), intent(in) :: iwork(*)
    end subroutine zheevx
  end interface lapack_heevx


  interface lapack_hpevx
    pure subroutine chpevx(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, &
        & m, w, z, ldz, work, rwork, iwork, ifail, info)
      import :: sp, ik
      complex(sp), intent(inout) :: ap(*)
      real(sp), intent(out) :: w(*)
      character(len=1), intent(in) :: uplo
      complex(sp), intent(out) :: z(ldz, *)
      real(sp), intent(in) :: vl
      real(sp), intent(in) :: vu
      integer(ik), intent(in) :: il
      integer(ik), intent(in) :: iu
      integer(ik), intent(out) :: m
      integer(ik), intent(out) :: ifail(*)
      real(sp), intent(in) :: abstol
      integer(ik), intent(out) :: info
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: range
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: ldz
      complex(sp), intent(in) :: work(*)
      real(sp), intent(in) :: rwork(*)
      integer(ik), intent(in) :: iwork(*)
    end subroutine chpevx
    pure subroutine zhpevx(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, &
        & m, w, z, ldz, work, rwork, iwork, ifail, info)
      import :: dp, ik
      complex(dp), intent(inout) :: ap(*)
      real(dp), intent(out) :: w(*)
      character(len=1), intent(in) :: uplo
      complex(dp), intent(out) :: z(ldz, *)
      real(dp), intent(in) :: vl
      real(dp), intent(in) :: vu
      integer(ik), intent(in) :: il
      integer(ik), intent(in) :: iu
      integer(ik), intent(out) :: m
      integer(ik), intent(out) :: ifail(*)
      real(dp), intent(in) :: abstol
      integer(ik), intent(out) :: info
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: range
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: ldz
      complex(dp), intent(in) :: work(*)
      real(dp), intent(in) :: rwork(*)
      integer(ik), intent(in) :: iwork(*)
    end subroutine zhpevx
  end interface lapack_hpevx

end module lapack_evx
