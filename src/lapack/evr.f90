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

module lapack_evr
  use lapack_kinds, only : sp, dp, ik
  implicit none
  private

  public :: lapack_syevr, lapack_heevr


  interface lapack_syevr
    pure subroutine ssyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, &
        & abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info)
      import :: sp, ik
      real(sp), intent(inout) :: a(lda, *)
      real(sp), intent(out) :: w(*)
      character(len=1), intent(in) :: uplo
      real(sp), intent(out) :: z(ldz, *)
      real(sp), intent(in) :: vl
      real(sp), intent(in) :: vu
      integer(ik), intent(in) :: il
      integer(ik), intent(in) :: iu
      integer(ik), intent(out) :: m
      integer(ik), intent(out) :: isuppz(*)
      real(sp), intent(in) :: abstol
      integer(ik), intent(out) :: info
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: range
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: ldz
      real(sp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
      integer(ik), intent(inout) :: iwork(*)
      integer(ik), intent(in) :: liwork
    end subroutine ssyevr
    pure subroutine dsyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, &
        & abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info)
      import :: dp, ik
      real(dp), intent(inout) :: a(lda, *)
      real(dp), intent(out) :: w(*)
      character(len=1), intent(in) :: uplo
      real(dp), intent(out) :: z(ldz, *)
      real(dp), intent(in) :: vl
      real(dp), intent(in) :: vu
      integer(ik), intent(in) :: il
      integer(ik), intent(in) :: iu
      integer(ik), intent(out) :: m
      integer(ik), intent(out) :: isuppz(*)
      real(dp), intent(in) :: abstol
      integer(ik), intent(out) :: info
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: range
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: ldz
      real(dp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
      integer(ik), intent(inout) :: iwork(*)
      integer(ik), intent(in) :: liwork
    end subroutine dsyevr
  end interface lapack_syevr

  interface lapack_heevr
    pure subroutine cheevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, &
        & abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, &
        & liwork, info)
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
      integer(ik), intent(out) :: isuppz(*)
      real(sp), intent(in) :: abstol
      integer(ik), intent(out) :: info
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: range
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: ldz
      complex(sp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
      real(sp), intent(inout) :: rwork(*)
      integer(ik), intent(in) :: lrwork
      integer(ik), intent(inout) :: iwork(*)
      integer(ik), intent(in) :: liwork
    end subroutine cheevr
    pure subroutine zheevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, &
        & abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, &
        & liwork, info)
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
      integer(ik), intent(out) :: isuppz(*)
      real(dp), intent(in) :: abstol
      integer(ik), intent(out) :: info
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: range
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
      integer(ik), intent(in) :: ldz
      complex(dp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
      real(dp), intent(inout) :: rwork(*)
      integer(ik), intent(in) :: lrwork
      integer(ik), intent(inout) :: iwork(*)
      integer(ik), intent(in) :: liwork
    end subroutine zheevr
  end interface lapack_heevr

end module lapack_evr
