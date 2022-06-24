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

module lapack_trs
   use lapack_kinds, only : sp, dp, ik
   implicit none
   private

   public :: lapack_sytrs, lapack_getrs


   interface lapack_sytrs
      pure subroutine ssytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: sp, ik
         integer(ik), intent(in) :: lda
         integer(ik), intent(in) :: ldb
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(inout) :: b(ldb, *)
         integer(ik), intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
         integer(ik), intent(in) :: nrhs
      end subroutine ssytrs
      pure subroutine dsytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp, ik
         integer(ik), intent(in) :: lda
         integer(ik), intent(in) :: ldb
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(inout) :: b(ldb, *)
         integer(ik), intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
         integer(ik), intent(in) :: nrhs
      end subroutine dsytrs

      module procedure :: lapack_ssytrs
      module procedure :: lapack_dsytrs
   end interface lapack_sytrs


   !> Solves a system of linear equations
   !>    A * X = B  or  A**T * X = B
   !> with a general N-by-N matrix A using the LU factorization computed
   !> by ?GETRF.
   interface lapack_getrs
      pure subroutine sgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: sp, ik
         real(sp), intent(in) :: a(lda, *)
         integer(ik), intent(in) :: ipiv(*)
         real(sp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: trans
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
         integer(ik), intent(in) :: nrhs
         integer(ik), intent(in) :: lda
         integer(ik), intent(in) :: ldb
      end subroutine sgetrs
      pure subroutine dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp, ik
         real(dp), intent(in) :: a(lda, *)
         integer(ik), intent(in) :: ipiv(*)
         real(dp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: trans
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
         integer(ik), intent(in) :: nrhs
         integer(ik), intent(in) :: lda
         integer(ik), intent(in) :: ldb
      end subroutine dgetrs

      module procedure :: lapack_sgetrs
      module procedure :: lapack_dgetrs
   end interface lapack_getrs


contains


subroutine lapack_ssytrs(amat, bmat, ipiv, uplo, info)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout) :: bmat(:, :)
   integer(ik), intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer(ik), intent(out), optional :: info
   character(len=1) :: ula
   integer(ik) :: stat, n, nrhs, lda, ldb

   ula = 'u'
   if (present(uplo)) ula = uplo
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2)
   nrhs = size(bmat, 2)
   call lapack_sytrs(ula, n, nrhs, amat, lda, ipiv, bmat, ldb, stat)
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[lapack] sytrs failed"
   end if
end subroutine lapack_ssytrs


subroutine lapack_dsytrs(amat, bmat, ipiv, uplo, info)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout) :: bmat(:, :)
   integer(ik), intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer(ik), intent(out), optional :: info
   character(len=1) :: ula
   integer(ik) :: stat, n, nrhs, lda, ldb

   ula = 'u'
   if (present(uplo)) ula = uplo
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2)
   nrhs = size(bmat, 2)
   call lapack_sytrs(ula, n, nrhs, amat, lda, ipiv, bmat, ldb, stat)
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[lapack] sytrs failed"
   end if
end subroutine lapack_dsytrs


subroutine lapack_sgetrs(amat, bmat, ipiv, trans, info)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout) :: bmat(:, :)
   integer(ik), intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: trans
   integer(ik), intent(out), optional :: info
   character(len=1) :: tra
   integer(ik) :: n, nrhs, lda, ldb, stat

   tra = 'n'
   if (present(trans)) tra = trans
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2)
   nrhs = size(bmat, 2)
   call lapack_getrs(tra, n, nrhs, amat, lda, ipiv, bmat, ldb, stat)
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[lapack] getrs failed"
   end if
end subroutine lapack_sgetrs


subroutine lapack_dgetrs(amat, bmat, ipiv, trans, info)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout) :: bmat(:, :)
   integer(ik), intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: trans
   integer(ik), intent(out), optional :: info
   character(len=1) :: tra
   integer(ik) :: n, nrhs, lda, ldb, stat

   tra = 'n'
   if (present(trans)) tra = trans
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2)
   nrhs = size(bmat, 2)
   call lapack_getrs(tra, n, nrhs, amat, lda, ipiv, bmat, ldb, stat)
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[lapack] getrs failed"
   end if
end subroutine lapack_dgetrs

end module lapack_trs
