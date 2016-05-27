project(preprocess)

# Preprocessing sources 
set(preprocess_srcs "exception.cpp" "stream.cpp")
add_library(preprocess SHARED ${preprocess_srcs})
