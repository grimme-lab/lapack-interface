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

module lapack_gv
   use lapack_kinds, only : sp, dp, ik
   implicit none
   private

   public :: lapack_sygv, lapack_spgv, lapack_hegv, lapack_hpgv


   interface lapack_sygv
     pure subroutine ssygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, &
         & info)
       import :: sp
       real(sp), intent(inout) :: a(lda, *)
       real(sp), intent(inout) :: b(ldb, *)
       real(sp), intent(out) :: w(*)
       integer, intent(in) :: itype
       character(len=1), intent(in) :: jobz
       character(len=1), intent(in) :: uplo
       integer, intent(out) :: info
       integer, intent(in) :: n
       integer, intent(in) :: lda
       integer, intent(in) :: ldb
       real(sp), intent(inout) :: work(*)
       integer, intent(in) :: lwork
     end subroutine ssygv
     pure subroutine dsygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, &
         & info)
       import :: dp
       real(dp), intent(inout) :: a(lda, *)
       real(dp), intent(inout) :: b(ldb, *)
       real(dp), intent(out) :: w(*)
       integer, intent(in) :: itype
       character(len=1), intent(in) :: jobz
       character(len=1), intent(in) :: uplo
       integer, intent(out) :: info
       integer, intent(in) :: n
       integer, intent(in) :: lda
       integer, intent(in) :: ldb
       real(dp), intent(inout) :: work(*)
       integer, intent(in) :: lwork
     end subroutine dsygv
   end interface lapack_sygv

   interface lapack_hegv
     pure subroutine chegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, &
         & rwork, info)
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
       real(sp), intent(in) :: rwork(*)
     end subroutine chegv
     pure subroutine zhegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, &
         & rwork, info)
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
       real(dp), intent(in) :: rwork(*)
     end subroutine zhegv
   end interface lapack_hegv

   interface lapack_spgv
     pure subroutine sspgv(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, info)
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
       real(sp), intent(in) :: work(*)
     end subroutine sspgv
     pure subroutine dspgv(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, info)
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
       real(dp), intent(in) :: work(*)
     end subroutine dspgv
   end interface lapack_spgv

   interface lapack_hpgv
     pure subroutine chpgv(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, rwork, &
         & info)
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
       complex(sp), intent(in) :: work(*)
       real(sp), intent(in) :: rwork(*)
     end subroutine chpgv
     pure subroutine zhpgv(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, rwork, &
         & info)
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
       complex(dp), intent(in) :: work(*)
       real(dp), intent(in) :: rwork(*)
     end subroutine zhpgv
   end interface lapack_hpgv

 end module lapack_gv
