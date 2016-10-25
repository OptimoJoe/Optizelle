# See if we can find loadjson.m.  In truth, I'm not really sure where to
# look for this, so it'll likely fail, which is ok.
find_path(JSONLAB_DIR loadjson.m)
if(JSONLAB_DIR STREQUAL JSONLAB_DIR-NOTFOUND) 
    message(FATAL_ERROR
        "The base directory for jsonlab must be defined in the variable JSONLAB_DIR.")
endif()
mark_as_advanced(JSONLAB_DIR FORCE)
