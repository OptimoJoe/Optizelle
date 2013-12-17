project(optizelle_python_json)

install(
    FILES
        __init__.py
        Unconstrained.py
        EqualityConstrained.py
        InequalityConstrained.py
        Constrained.py
    DESTINATION share/optizelle/python/Optizelle/json)
