Requirements
------------
- [Qt 6.2 or later](https://www.qt.io/)
- [CMake 3.16+](https://cmake.org/)
- [Python 3](https://www.python.org/)
- [Boost.Python](https://www.boost.org/) (linked against Python 3)
- [`libpng`](http://www.libpng.org/pub/png/libpng.html)
- [Lemon](https://www.hwaci.com/sw/lemon/)
- [Flex](https://github.com/westes/flex)
- [ninja](https://ninja-build.org/) (recommended)

--------------------------------------------------------------------------------

Arch Linux
----------
```
sudo pacman -S --needed git base-devel cmake ninja qt6-base python \
                        boost libpng lemon flex

git clone https://github.com/nsfm/stibium
cd stibium
mkdir build
cd build
cmake -GNinja ..
ninja
./app/antimony
```

--------------------------------------------------------------------------------

Debian / Ubuntu
---------------
```
sudo apt install git build-essential cmake ninja-build \
                 qt6-base-dev libgl1-mesa-dev \
                 python3-dev libboost-python-dev \
                 libpng-dev lemon flex

git clone https://github.com/nsfm/stibium
cd stibium
mkdir build
cd build
cmake -GNinja ..
ninja
./app/antimony
```

--------------------------------------------------------------------------------

macOS
-----
With [homebrew](https://brew.sh/):
```
brew install cmake ninja qt python boost-python3 libpng lemon flex

git clone https://github.com/nsfm/stibium
cd stibium
mkdir build
cd build
cmake -DCMAKE_PREFIX_PATH=$(brew --prefix qt) -GNinja ..
ninja

open app/Antimony.app
```
(`brew install qt` is Qt 6; do not use the legacy `qt@5` keg.)

--------------------------------------------------------------------------------

Tests
-----
The graph and fab test suites are wired into CTest:
```
cd build
ctest
```

Install
-------
`sudo ninja install` copies the executable to `/usr/local/bin` and the
Python node/shape libraries to `/usr/local/share/antimony`.

Packaging
---------
The legacy packaging under `deploy/` (macOS bundle script, Debian
source packaging) predates Qt6 and the Stibium rename; it is retained
for reference only and will be rebuilt (CI + AppImage/Flatpak + macOS
bundle) at the first-release milestone — see TODO.md.

Windows
-------
Upstream Antimony had an experimental MSYS2/Qt5 recipe; it has not been
attempted against Qt6/Stibium. The Qt6 port removes some old obstacles,
so a Windows build is plausible — contributions welcome.
