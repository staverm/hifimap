# hifimap

Hifimap is a mapping tool optimized for HiFi reads. Please see the original paper `/hifimap_paper.pdf` for more details.

# Usage

To build hifimap run the following commands:

```bash
git clone https://github.com/staverm/hifimap && cd hifimap && mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release .. && make
```

which will create hifimap library, executable and unit tests. 

Running `make install` will install the executable. In order to install the library, both biosoup and thread_pool (see Dependencies) need to be installed beforehand, and option `hifimap_install` used while configuring the build. Once the library is installed, a package will be copied to your system that can be searched and linked with:

```cmake
find_package(hifimap)
target_link_libraries(<target> hifimap::hifimap)
```

#### Build options

- `hifimap_install`: generate library install target
- `hifimap_build_exe`: build executable
- `hifimap_build_tests`: build unit tests

#### Dependencies

- gcc 4.8+ | clang 3.5+
- cmake 3.11+
- pthread
- (hifimap_exe)(hifimap_test) zlib 1.2.8+

###### Hidden

- rvaser/biosoup 0.10.0
- rvaser/thread_pool 3.0.3
- (hifimap_exe)(hifimap_test) rvaser/bioparser 3.0.13
- (hifimap_test) google/googletest 1.10.0
