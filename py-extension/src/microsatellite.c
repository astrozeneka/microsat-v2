#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <string.h>

static PyObject *microsatellite_searchssr(PyObject *self, PyObject *args){
    char *seq;
    int mono = 0;
    int di = 0;
    int tri = 0;
    int tetra = 0;
    int penta = 0;
    int hexa = 0;
    int rep[6];

    size_t len;
    int start;
    int length;
    int repeat;
    int i;
    int j;
    char motif[7] = "\0";

    PyObject *result = PyList_New(0);
    PyObject *tmp;

    return result;
}

static PyMethodDef MicrosatelliteMethods[] = {
        {"searchssr", microsatellite_searchssr, METH_VARARGS, "Search for microsatellites in a sequence"},
        {NULL, NULL, 0, NULL}
};
static struct PyModuleDef microsatellite_definition = {
        PyModuleDef_HEAD_INIT,
        "microsatellite",
        "Search for microsatellites in a genome",
        -1,
        MicrosatelliteMethods
};
PyMODINIT_FUNC PyInit_microsatellite(void){
    Py_Initialize();
    return PyModule_Create(&microsatellite_definition);
}