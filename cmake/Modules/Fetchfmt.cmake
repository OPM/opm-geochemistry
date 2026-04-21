include(FetchContent)

if(NOT TARGET fmt::fmt AND NOT TARGET fmt::fmt-header-only)
  FetchContent_Declare(
      fmt
      DOWNLOAD_EXTRACT_TIMESTAMP ON
      URL https://github.com/fmtlib/fmt/archive/refs/tags/10.1.1.tar.gz
      URL_HASH SHA512=288c349baac5f96f527d5b1bed0fa5f031aa509b4526560c684281388e91909a280c3262a2474d963b5d1bf7064b1c9930c6677fe54a0d8f86982d063296a54c)
  FetchContent_MakeAvailable(fmt)
endif()

set(fmt_FOUND TRUE)
