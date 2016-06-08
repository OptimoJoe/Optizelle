project(exceptions)

# Set what interfaces we're using
#set(interfaces "cpp" "python" "matlab" "octave")
set(interfaces "cpp" "python")

# Copy in files for our unit tests 
is_unit("${interfaces}" enable_unit)

# Compile and add our unit tests 
compile_add_unit(objective "${interfaces}")
