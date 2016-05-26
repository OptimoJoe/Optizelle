# Tests our build setup flags.  Since we include these flags in our
# documentation, this verifies that our flags mostly work correctly.
project(config)

# Set the directory name for the configuration
set(config_dir "config_unit")

# Read in the configuration flags
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/flags.txt" flags)

# Create a place for the build 
add_test(
    NAME prepare_configuration_check
    COMMAND ${CMAKE_COMMAND}
        "-E"
        "make_directory"
        "${config_dir}"
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")

# Configure a build with all the flags that we read in above 
add_test(
    NAME configuration_check 
    COMMAND ${CMAKE_COMMAND}
        ${flags}
        "${CMAKE_SOURCE_DIR}"
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/${config_dir}")

# Remove the result since it could interfere with tests later
add_test(
    NAME cleanup_configuration_check
    COMMAND ${CMAKE_COMMAND}
        "-E"
        "remove_directory"
        "${config_dir}"
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
