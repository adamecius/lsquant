#!/bin/bash
# compiler launcher for the legacy-final provenance build on Aron: the only Eigen >= 3.4 there (lsq-env, 3.5.0) needs C++14
a=(); for x in "$@"; do [ "$x" = -std=c++11 ] && x=-std=c++14; [ "$x" = -std=gnu++11 ] && x=-std=gnu++14; a+=("$x"); done; exec "${a[@]}"
