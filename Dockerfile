# Dockerfile for building htslib dependency + 
# all three SIMD versions of the vargas binary

FROM rikorose/gcc-cmake:gcc-6
ADD . /
RUN chmod +x /docker_compile.sh
RUN ./docker_compile.sh
