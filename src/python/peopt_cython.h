#ifndef PEOPT_CYTHON_H
#define PEOPT_CYTHON_H
PyObject* pypeopt(
    PyObject* opt_type_,
    PyObject* vs_,
    PyObject* fns_,
    PyObject* pts_,
    PyObject* fname_
);
#endif
