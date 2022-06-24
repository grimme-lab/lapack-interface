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

module lapack_ev
  use lapack_kinds, only : sp, dp, ik
  implicit none
  private

  public :: lapack_syev, lapack_spev, lapack_heev, lapack_hpev


  interface lapack_syev
    pure subroutine ssyev(jobz, uplo, n, a, lda, w, work, lwork, info)
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
    end subroutine ssyev
    pure subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
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
    end subroutine dsyev
  end interface lapack_syev

  interface lapack_heev
    pure subroutine cheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
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
      real(sp), intent(in) :: rwork(*)
    end subroutine cheev
    pure subroutine zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
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
      real(dp), intent(in) :: rwork(*)
    end subroutine zheev
  end interface lapack_heev

  interface lapack_spev
    pure subroutine sspev(jobz, uplo, n, ap, w, z, ldz, work, info)
      import :: sp, ik
      real(sp), intent(inout) :: ap(*)
      real(sp), intent(out) :: w(*)
      character(len=1), intent(in) :: uplo
      real(sp), intent(out) :: z(ldz, *)
      integer(ik), intent(out) :: info
      character(len=1), intent(in) :: jobz
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: ldz
      real(sp), intent(in) :: work(*)
    end subroutine sspev
    pure subroutine dspev(jobz, uplo, n, ap, w, z, ldz, work, info)
      import :: dp, ik
      real(dp), intent(inout) :: ap(*)
      real(dp), intent(out) :: w(*)
      character(len=1), intent(in) :: uplo
      real(dp), intent(out) :: z(ldz, *)
      integer(ik), intent(out) :: info
      character(len=1), intent(in) :: jobz
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: ldz
      real(dp), intent(in) :: work(*)
    end subroutine dspev
  end interface lapack_spev

  interface lapack_hpev
    pure subroutine chpev(jobz, uplo, n, ap, w, z, ldz, work, rwork, info)
      import :: sp, ik
      complex(sp), intent(inout) :: ap(*)
      real(sp), intent(out) :: w(*)
      character(len=1), intent(in) :: uplo
      complex(sp), intent(out) :: z(ldz, *)
      integer(ik), intent(out) :: info
      character(len=1), intent(in) :: jobz
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: ldz
      complex(sp), intent(in) :: work(*)
      real(sp), intent(in) :: rwork(*)
    end subroutine chpev
    pure subroutine zhpev(jobz, uplo, n, ap, w, z, ldz, work, rwork, info)
      import :: dp, ik
      complex(dp), intent(inout) :: ap(*)
      real(dp), intent(out) :: w(*)
      character(len=1), intent(in) :: uplo
      complex(dp), intent(out) :: z(ldz, *)
      integer(ik), intent(out) :: info
      character(len=1), intent(in) :: jobz
      integer(ik), intent(in) :: n
      integer(ik), intent(in) :: ldz
      complex(dp), intent(in) :: work(*)
      real(dp), intent(in) :: rwork(*)
    end subroutine zhpev
  end interface lapack_hpev

end module lapack_ev
