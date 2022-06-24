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

project(
  'lapack-interface',
  'fortran',
  version: '1.0.0',
  license: 'Apache-2.0',
  meson_version: '>=0.55',
  default_options: [
    'buildtype=debugoptimized',
    'default_library=both',
  ],
)
install = not (meson.is_subproject() and get_option('default_library') == 'static')

# General configuration information
lib_deps = []
subdir('config')

# Collect source of the project
srcs = []
subdir('src')

# LAPACK interface library target
lapack_interface_lib = library(
  meson.project_name(),
  sources: srcs,
  version: meson.project_version(),
  dependencies: lib_deps,
  install: install,
)

# Export dependency for other projects and test suite
lapack_interface_inc = lapack_interface_lib.private_dir_include()
lapack_interface_dep = declare_dependency(
  link_with: lapack_interface_lib,
  include_directories: lapack_interface_inc,
  dependencies: lib_deps,
)

# Package the license files
lapack_interface_lic = files(
  'LICENSE',
)

if install
  # Distribute the license files in share/licenses/<name>
  install_data(
    lapack_interface_lic,
    install_dir: get_option('datadir')/'licenses'/meson.project_name()
  )

  module_id = meson.project_name() / fc_id + '-' + fc.version()
  meson.add_install_script(
    find_program(files('config'/'install-mod.py')),
    get_option('includedir') / module_id,
  )

  pkg = import('pkgconfig')
  pkg.generate(
    lapack_interface_lib,
    description: 'Interfaces and wrappers to LAPACK',
    subdirs: ['', module_id],
  )
endif

