project(optizelle_python_unconstrained)

install(
    FILES __init__.py State.py Functions.py Algorithms.py Restart.py
    DESTINATION ${PYTHON_SITE_PACKAGES}/Optizelle/Unconstrained)
