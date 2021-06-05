cmake_minimum_required(VERSION 3.11)

project(hifimap VERSION 1.0.0
            LANGUAGES CXX
            DESCRIPTION "Hifimap is a mapping tool optimized for HiFi reads.")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

include(FetchContent)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(hifimap_main_project ON)
endif ()
option(hifimap_install "Generate library install target" OFF)
option(hifimap_build_exe "Build executable" ${hifimap_main_project})
option(hifimap_build_tests "Build unit tests" ${hifimap_main_project})

find_package(biosoup 0.10.0 QUIET)
find_package(thread_pool 3.0.3 QUIET)

if (NOT biosoup_FOUND)
  if (hifimap_install)
    message(FATAL_ERROR "Missing package biosoup!")
  endif ()

  FetchContent_Declare(
    biosoup
    GIT_REPOSITORY https://github.com/rvaser/biosoup
    GIT_TAG 0.10.0)

  FetchContent_GetProperties(biosoup)
  if (NOT biosoup_POPULATED)
    FetchContent_Populate(biosoup)
    add_subdirectory(
      ${biosoup_SOURCE_DIR}
      ${biosoup_BINARY_DIR}
      EXCLUDE_FROM_ALL)
  endif ()
endif ()

if (NOT thread_pool_FOUND)
  if (hifimap_install)
    message(FATAL_ERROR "Missing package thread_pool!")
  endif ()

  FetchContent_Declare(
    thread_pool
    GIT_REPOSITORY https://github.com/rvaser/thread_pool
    GIT_TAG 3.0.3)

  FetchContent_GetProperties(thread_pool)
  if (NOT thread_pool_POPULATED)
    FetchContent_Populate(thread_pool)
    add_subdirectory(
      ${thread_pool_SOURCE_DIR}
      ${thread_pool_BINARY_DIR}
      EXCLUDE_FROM_ALL)
  endif ()
endif ()

if (hifimap_build_exe OR hifimap_build_tests)
  find_package(bioparser 3.0.13 QUIET)
  if (NOT bioparser_FOUND)
    FetchContent_Declare(
      bioparser
      GIT_REPOSITORY https://github.com/rvaser/bioparser
      GIT_TAG 3.0.13)

    FetchContent_GetProperties(bioparser)
    if (NOT bioparser_POPULATED)
      FetchContent_Populate(bioparser)
      add_subdirectory(
        ${bioparser_SOURCE_DIR}
        ${bioparser_BINARY_DIR}
        EXCLUDE_FROM_ALL)
    endif ()
  endif ()
endif ()

if (hifimap_build_tests)
  find_package(GTest 1.10.0 QUIET)
  if (NOT GTest_FOUND)
    FetchContent_Declare(
      googletest
      GIT_REPOSITORY https://github.com/google/googletest
      GIT_TAG release-1.10.0)

    FetchContent_GetProperties(googletest)
    if (NOT googletest_POPULATED)
      FetchContent_Populate(googletest)
      add_subdirectory(
        ${googletest_SOURCE_DIR}
        ${googletest_BINARY_DIR}
        EXCLUDE_FROM_ALL)
      add_library(GTest::Main ALIAS gtest_main)
    endif ()
  endif ()
endif ()

add_library(hifimap
  src/minimizer_engine.cpp)
add_library(hifimap::hifimap ALIAS hifimap)

target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)

target_link_libraries(hifimap
  biosoup::biosoup
  thread_pool::thread_pool)

if (hifimap_install)
  configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    COMPATIBILITY SameMajorVersion)

  install(
    TARGETS hifimap
    EXPORT ${PROJECT_NAME}Targets
    DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(
    DIRECTORY include/hifimap
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  install(
    EXPORT ${PROJECT_NAME}Targets
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
  install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
          ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
endif ()

if (hifimap_build_exe)
  add_executable(hifimap_exe
    src/main.cpp)

  target_link_libraries(hifimap_exe
    hifimap
    bioparser::bioparser)

  target_compile_definitions(hifimap_exe PRIVATE VERSION="${PROJECT_VERSION}")
  set_property(TARGET hifimap_exe PROPERTY OUTPUT_NAME hifimap)

  install(TARGETS hifimap_exe DESTINATION ${CMAKE_INSTALL_BINDIR})
endif ()

if (hifimap_build_tests)
  add_executable(hifimap_test
    test/minimizer_engine_test.cpp)

  target_link_libraries(hifimap_test
    hifimap
    bioparser::bioparser
    GTest::Main)

  target_compile_definitions(hifimap_test
    PRIVATE TEST_DATA="${PROJECT_SOURCE_DIR}/test/data/sample.fasta.gz")
endif ()
