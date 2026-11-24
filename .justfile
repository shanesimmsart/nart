# list available recipes
@default:
  just --list --unsorted

# setup build config, defaults to Ninja
@setup mode='release' generator='Ninja':
  cmake --preset {{mode}} -G '{{generator}}'

# set build option to value
@set option value mode='release':
  cmake --preset {{mode}} -D {{option}}={{value}}

# configure all build options via TUI
@config mode='release':
  ccmake --preset {{mode}}

# build target, defaults to all tools
@build target='all' mode='release':
  cmake --preset {{mode}}
  cmake --build --preset {{mode}} --target {{target}} --log-level=ERROR

# convenience aliases
@debug:
  just build --mode debug

@release:
  just build --mode release

# remove build and reset
@clean:
  rm -rf build
