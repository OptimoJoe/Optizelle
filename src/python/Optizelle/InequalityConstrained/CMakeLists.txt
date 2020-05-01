project(optizelle_python_inequality_constrained)

install(
    FILES __init__.py State.py Functions.py Algorithms.py Restart.py
    DESTINATION ${PYTHON_SITE_PACKAGES}/Optizelle/InequalityConstrained)
