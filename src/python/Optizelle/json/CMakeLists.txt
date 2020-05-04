project(optizelle_python_json)

install(
    FILES
        __init__.py
        Unconstrained.py
        EqualityConstrained.py
        InequalityConstrained.py
        Constrained.py
        Serialization.py
    DESTINATION ${PYTHON_SITE_PACKAGES}/Optizelle/json)
