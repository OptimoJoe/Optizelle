# Combine all the license information into a single file
project(licenses)

# Concatenates two files
function(cat in out)
    file(READ ${in} contents)
    file(APPEND ${out} "${contents}")
endfunction()

# Includes a license into the overall license
set(out "${CMAKE_CURRENT_BINARY_DIR}/LICENSE.txt")
file(REMOVE ${out})
function(add_license package interface text)
    file(APPEND ${out} "-~-~${package} (${interface})-~-~\n")
    cat(${text} ${out})
    file(APPEND ${out} "\n")
endfunction()

# Combine all the licenses
add_license("Optizelle" "all interfaces" optizelle.txt)
add_license("JsonCpp" "all interfaces" jsoncpp.txt)
add_license("BLAS/LAPACK" "all interfaces when building default BLAS/LAPACK"
    blas_lapack.txt)
add_license("CMake" "all interfaces" cmake.txt)
add_license("WiX" "all interfaces when building Windows installer" wix.txt)
add_license("GCC" "All interfaces when using GCC.  Note, we do *not* modify the compiler.  See the GCC Runtime Library Exception." gcc.txt)
add_license("GCC" "all interfaces when using GCC" gcc_runtime.txt)
add_license("TeX Live" "documentation when using TeX Live for pdflatex"
    texlive.txt)
add_license("Python" "Python interface" python.txt)
add_license("NumPy" "Python interface" numpy.txt)
add_license("MATLAB" "MATLAB interface" matlab.txt)
add_license("JSONlab" "MATLAB and/or Octave interface" jsonlab.txt)
add_license("Octave" "Octave interface" octave.txt)

# Install the license
install(FILES ${out} DESTINATION share/optizelle)
