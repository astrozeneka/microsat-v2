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

    if(!PyArg_ParseTuple(args, "s(iiiiii)", &seq, &mono, &di, &tri, &tetra, &penta, &hexa))
        return NULL;
    rep[0] = mono;
    rep[1] = di;
    rep[2] = tri;
    rep[3] = tetra;
    rep[4] = penta;
    rep[5] = hexa;

    len = strlen(seq);
    for (i=0; i<len; i++) {
        if (seq[i] == 78)
            continue;

        for (j=1; j<=6; j++) {
            start = i;
            length = j;
            while(start+length<len && seq[i]==seq[i+j] && seq[i]!=78){
                i++;
                length++;
            }
            repeat = length/j;
            if(repeat>=rep[j-1]) {
                strncpy(motif, seq+start, j);
                motif[j] = '\0';
                length = repeat*j;
                tmp = Py_BuildValue("(siiiii)", motif, j, repeat, start+1, start+length, length);
                PyList_Append(result, tmp);
                Py_DECREF(tmp);
                i = start + length;
                j = 0;
            } else {
                i = start;
            }
        }
    }
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