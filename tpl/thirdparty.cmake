# This maintains the configuration for third-party packages that can be installed 
# with the Coek Project.  Simply update this file to change the revision.
# One can use different revision on different platforms.
# e.g.
# if (UNIX)
#   ..
# else (APPLE)
#   ..
# endif()

add_revision(catch2
  SRC Catch2-2.13.6
  URL "https://github.com/catchorg/Catch2/archive/refs/tags/v2.13.6.tar.gz"
  URL_MD5 c7c7ef181b9e08418fd9f2ef8159d03f
  )

#add_revision(cppad
#  SRC CppAD-20210000.6
#  URL "https://github.com/coin-or/CppAD/archive/refs/tags/20210000.6.tar.gz"
#  URL_MD5 d63b03bce0417c420e610cb1cfb64d33
#  )
add_revision(cppad
  SRC CppAD-20220000.4
  URL "https://github.com/coin-or/CppAD/archive/refs/tags/20220000.4.tar.gz"
  URL_MD5 ddf5459514b3435ff6c4d537d39608dc
  )

add_revision(rapidjson
  SRC rapidjson-1.1.0
  URL "https://github.com/Tencent/rapidjson/archive/refs/tags/v1.1.0.tar.gz"
  URL_MD5 badd12c511e081fec6c89c43a7027bce
  )

add_revision(pybind11
  SRC "pybind11-2.6.2"
  URL "https://github.com/pybind/pybind11/archive/refs/tags/v2.6.2.tar.gz"
  URL_MD5 c5ea9c4c57082e05efe14e4b34323bfd
  )

add_revision(fmtlib
  SRC "fmt-8.0.0"
  URL "https://github.com/fmtlib/fmt/archive/refs/tags/8.0.0.tar.gz"
  URL_MD5 001d59c967115682aed0e836ba5753a8
  )

