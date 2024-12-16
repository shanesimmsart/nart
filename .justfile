# list available recipes
@default:
  just --list --unsorted

# setup build config, defaults to Ninja
@setup generator='Ninja':
  cmake --preset unix -G '{{generator}}'

# set build option to value
@set option value:
  cmake --preset unix -D {{option}}={{value}}

# configure all build options via TUI
@config:
  ccmake --preset unix

# build target, defaults to all tools
@build target='all':
  cmake --build --preset unix --target {{target}}

# remove build and reset
@clean:
  rm -rf build
