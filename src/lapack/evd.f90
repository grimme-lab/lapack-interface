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

module lapack_evd
  use lapack_kinds, only : sp, dp, ik
  implicit none
  private

  public :: lapack_syevd, lapack_spevd, lapack_heevd, lapack_hpevd


  interface lapack_syevd
    pure subroutine ssyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, &
        & liwork, info)
      import :: sp, ik
      real(sp), intent(inout) :: a(lda, *)
      real(sp), intent(out) :: w(*)
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: uplo
      integer(ik), intent(out) :: info
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
      real(sp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
      integer(ik), intent(inout) :: iwork(*)
      integer(ik), intent(in) :: liwork
    end subroutine ssyevd
    pure subroutine dsyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, &
        & liwork, info)
      import :: dp, ik
      real(dp), intent(inout) :: a(lda, *)
      real(dp), intent(out) :: w(*)
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: uplo
      integer(ik), intent(out) :: info
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
      real(dp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
      integer(ik), intent(inout) :: iwork(*)
      integer(ik), intent(in) :: liwork
    end subroutine dsyevd
  end interface lapack_syevd

  interface lapack_heevd
    pure subroutine cheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, &
        & lrwork, iwork, liwork, info)
      import :: sp, ik
      complex(sp), intent(inout) :: a(lda, *)
      real(sp), intent(out) :: w(*)
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: uplo
      integer(ik), intent(out) :: info
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
      complex(sp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
      real(sp), intent(inout) :: rwork(*)
      integer(ik), intent(in) :: lrwork
      integer(ik), intent(inout) :: iwork(*)
      integer(ik), intent(in) :: liwork
    end subroutine cheevd
    pure subroutine zheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, &
        & lrwork, iwork, liwork, info)
      import :: dp, ik
      complex(dp), intent(inout) :: a(lda, *)
      real(dp), intent(out) :: w(*)
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: uplo
      integer(ik), intent(out) :: info
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: lda
      complex(dp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
      real(dp), intent(inout) :: rwork(*)
      integer(ik), intent(in) :: lrwork
      integer(ik), intent(inout) :: iwork(*)
      integer(ik), intent(in) :: liwork
    end subroutine zheevd
  end interface lapack_heevd

  interface lapack_spevd
    pure subroutine sspevd(jobz, uplo, n, ap, w, z, ldz, work, lwork, iwork, &
        & liwork, info)
      import :: sp, ik
      real(sp), intent(inout) :: ap(*)
      real(sp), intent(out) :: w(*)
      character(len=1), intent(in) :: uplo
      real(sp), intent(out) :: z(ldz, *)
      integer(ik), intent(out) :: info
      character(len=1), intent(in) :: jobz
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: ldz
      real(sp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
      integer(ik), intent(inout) :: iwork(*)
      integer(ik), intent(in) :: liwork
    end subroutine sspevd
    pure subroutine dspevd(jobz, uplo, n, ap, w, z, ldz, work, lwork, iwork, &
        & liwork, info)
      import :: dp, ik
      real(dp), intent(inout) :: ap(*)
      real(dp), intent(out) :: w(*)
      character(len=1), intent(in) :: uplo
      real(dp), intent(out) :: z(ldz, *)
      integer(ik), intent(out) :: info
      character(len=1), intent(in) :: jobz
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: ldz
      real(dp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
      integer(ik), intent(inout) :: iwork(*)
      integer(ik), intent(in) :: liwork
    end subroutine dspevd
  end interface lapack_spevd

  interface lapack_hpevd
    pure subroutine chpevd(jobz, uplo, n, ap, w, z, ldz, work, lwork, rwork, &
        & lrwork, iwork, liwork, info)
      import :: sp, ik
      complex(sp), intent(inout) :: ap(*)
      real(sp), intent(out) :: w(*)
      character(len=1), intent(in) :: uplo
      complex(sp), intent(out) :: z(ldz, *)
      integer(ik), intent(out) :: info
      character(len=1), intent(in) :: jobz
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: ldz
      complex(sp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
      real(sp), intent(inout) :: rwork(*)
      integer(ik), intent(in) :: lrwork
      integer(ik), intent(inout) :: iwork(*)
      integer(ik), intent(in) :: liwork
    end subroutine chpevd
    pure subroutine zhpevd(jobz, uplo, n, ap, w, z, ldz, work, lwork, rwork, &
        & lrwork, iwork, liwork, info)
      import :: dp, ik
      complex(dp), intent(inout) :: ap(*)
      real(dp), intent(out) :: w(*)
      character(len=1), intent(in) :: uplo
      complex(dp), intent(out) :: z(ldz, *)
      integer(ik), intent(out) :: info
      character(len=1), intent(in) :: jobz
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: ldz
      complex(dp), intent(inout) :: work(*)
      integer(ik), intent(in) :: lwork
      real(dp), intent(inout) :: rwork(*)
      integer(ik), intent(in) :: lrwork
      integer(ik), intent(inout) :: iwork(*)
      integer(ik), intent(in) :: liwork
    end subroutine zhpevd
  end interface lapack_hpevd

end module lapack_evd
