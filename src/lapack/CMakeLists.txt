# This file is part of lapack-interface.
# SPDX-Identifier: Apache-2.0
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
  APPEND srcs
  "${dir}/ev.f90"
  "${dir}/evd.f90"
  "${dir}/evr.f90"
  "${dir}/evx.f90"
  "${dir}/gst.f90"
  "${dir}/gv.f90"
  "${dir}/gvd.f90"
  "${dir}/gvx.f90"
  "${dir}/kinds.F90"
  "${dir}/trf.f90"
  "${dir}/tri.f90"
  "${dir}/trs.f90"
  "${dir}/utils.f90"
)

set(srcs "${srcs}" PARENT_SCOPE)
