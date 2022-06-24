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

module lapack
  use lapack_ev, only : lapack_syev, lapack_spev, lapack_heev, lapack_hpev
  use lapack_evd, only : lapack_syevd, lapack_spevd, lapack_heevd, lapack_hpevd
  use lapack_evr, only : lapack_syevr, lapack_heevr
  use lapack_evx, only : lapack_syevx, lapack_spevx, lapack_heevx, lapack_hpevx
  use lapack_gst, only : lapack_sygst
  use lapack_gv, only : lapack_sygv, lapack_spgv, lapack_hegv, lapack_hpgv
  use lapack_gvd, only : lapack_sygvd, lapack_spgvd, lapack_hegvd, lapack_hpgvd
  use lapack_gvx, only : lapack_sygvx, lapack_spgvx, lapack_hegvx, lapack_hpgvx
  use lapack_trf, only : lapack_sytrf, lapack_potrf, lapack_getrf
  use lapack_tri, only : lapack_sytri, lapack_getri
  use lapack_trs, only : lapack_sytrs, lapack_getrs
  implicit none
  public
end module lapack

