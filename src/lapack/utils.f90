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

module lapack_utils
   implicit none
   private

   public :: lapack_env

   interface lapack_env
      pure integer function ilaenv(ispec, name, opts, n1, n2, n3, n4)
         integer, intent(in) :: ispec
         character(len=1), intent(in) :: name
         character(len=1), intent(in) :: opts
         integer, intent(in) :: n1
         integer, intent(in) :: n2
         integer, intent(in) :: n3
         integer, intent(in) :: n4
      end function ilaenv
   end interface lapack_env

end module lapack_utils
