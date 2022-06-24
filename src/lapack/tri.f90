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

module lapack_tri
   use lapack_kinds, only : sp, dp, ik
   use lapack_utils, only : lapack_env
   implicit none
   private

   public :: lapack_sytri, lapack_getri


   interface lapack_sytri
      pure subroutine ssytri(uplo, n, a, lda, ipiv, work, info)
         import :: sp, ik
         integer(ik), intent(in) :: lda
         real(sp), intent(inout) :: a(lda, *)
         integer(ik), intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
         real(sp), intent(in) :: work(*)
      end subroutine ssytri
      pure subroutine dsytri(uplo, n, a, lda, ipiv, work, info)
         import :: dp, ik
         integer(ik), intent(in) :: lda
         real(dp), intent(inout) :: a(lda, *)
         integer(ik), intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
         real(dp), intent(in) :: work(*)
      end subroutine dsytri

      module procedure :: lapack_ssytri
      module procedure :: lapack_dsytri
   end interface lapack_sytri


   !> Computes the inverse of a matrix using the LU factorization
   !> computed by ?GETRF.
   !>
   !> This method inverts U and then computes inv(A) by solving the system
   !> inv(A)*L = inv(U) for inv(A).
   interface lapack_getri
      pure subroutine sgetri(n, a, lda, ipiv, work, lwork, info)
         import :: sp, ik
         real(sp), intent(inout) :: a(lda, *)
         integer(ik), intent(in) :: ipiv(*)
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
         integer(ik), intent(in) :: lda
         real(sp), intent(inout) :: work(*)
         integer(ik), intent(in) :: lwork
      end subroutine sgetri
      pure subroutine dgetri(n, a, lda, ipiv, work, lwork, info)
         import :: dp, ik
         real(dp), intent(inout) :: a(lda, *)
         integer(ik), intent(in) :: ipiv(*)
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
         integer(ik), intent(in) :: lda
         real(dp), intent(inout) :: work(*)
         integer(ik), intent(in) :: lwork
      end subroutine dgetri

      module procedure :: lapack_sgetri
      module procedure :: lapack_dgetri
   end interface lapack_getri


contains


subroutine lapack_ssytri(amat, ipiv, uplo, info)
   real(sp), intent(inout) :: amat(:, :)
   integer(ik), intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer(ik), intent(out), optional :: info
   character(len=1) :: ula
   integer(ik) :: stat, n, lda
   real(sp), allocatable :: work(:)

   ula = 'u'
   if (present(uplo)) ula = uplo
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   allocate(work(n), stat=stat)
   if (stat==0) then
      call lapack_sytri(ula, n, amat, lda, ipiv, work, stat)
   else
      stat = -1000_ik
   end if
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[lapack] sytri failed"
   end if
end subroutine lapack_ssytri


subroutine lapack_dsytri(amat, ipiv, uplo, info)
   real(dp), intent(inout) :: amat(:, :)
   integer(ik), intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer(ik), intent(out), optional :: info
   character(len=1) :: ula
   integer(ik) :: stat, n, lda
   real(dp), allocatable :: work(:)

   ula = 'u'
   if (present(uplo)) ula = uplo
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   allocate(work(n), stat=stat)
   if (stat==0) then
      call lapack_sytri(ula, n, amat, lda, ipiv, work, stat)
   else
      stat = -1000_ik
   end if
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[lapack] sytri failed"
   end if
end subroutine lapack_dsytri


subroutine lapack_sgetri(amat, ipiv, info)
   real(sp), intent(inout) :: amat(:, :)
   integer(ik), intent(in) :: ipiv(:)
   integer(ik), intent(out), optional :: info
   integer(ik) :: n, lda, lwork, stat
   real(sp), allocatable :: work(:)

   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   lwork = n * lapack_env(1_ik, 'SGETRI', ' ', n, -1_ik, -1_ik, -1_ik)
   allocate(work(lwork), stat=stat)
   if (stat == 0) then
      call lapack_getri(n, amat, lda, ipiv, work, lwork, stat)
   end if
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[lapack] getri failed"
   end if
end subroutine lapack_sgetri


subroutine lapack_dgetri(amat, ipiv, info)
   real(dp), intent(inout) :: amat(:, :)
   integer(ik), intent(in) :: ipiv(:)
   integer(ik), intent(out), optional :: info
   integer(ik) :: n, lda, lwork, stat
   real(dp), allocatable :: work(:)

   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   lwork = n * lapack_env(1_ik, 'DGETRI', ' ', n, -1_ik, -1_ik, -1_ik)
   allocate(work(lwork), stat=stat)
   if (stat == 0) then
      call lapack_getri(n, amat, lda, ipiv, work, lwork, stat)
   end if
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[lapack] getri failed"
   end if
end subroutine lapack_dgetri

end module lapack_tri
