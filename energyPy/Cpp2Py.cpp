#define MPICH_SKIP_MPICXX 1
#define OMPI_SKIP_MPICXX  1
#include <Python.h>
#include <mpi4py/mpi4py.h>
#include <numpy/arrayobject.h>
#include "energy.h"

#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <unordered_map>
/* -------------------------------------------------------------------------- */





/*
This function reads the arguments that is called in python
and stored initially as PyObject. The PyObjects are then used
to derive a corresponding C-Object from it. In the case of the
Communicator, we use the PyMPIComm_Get interface function that
generates C-datatype MPI_Comm from a PyObject.
*/

static PyObject *
Cpp2Py_CalculateEnergy(PyObject *self, PyObject *args)
{
  int argc = PyTuple_GET_SIZE(args);
  
  PyObject *py_comm = NULL;
  PyArrayObject *B;
  PyArrayObject *I;
  int nB, nI;
  double En = 0.0;
  // Reads the python MPIComm, numpy array, integer and string.
  // The O refers to a PyObject or PyArrayObject as an argument.
  // i refers to an integer.
  // s refers to a string.
  
  if (!PyArg_ParseTuple(args, "iOiO:Cpp2Py", &nB, &B, &nI, &I))
    return NULL;

  // Grab the data from the PyArrayObject.
  
  double *B_c                  = (double *)    B->data;
  double *I_c                  = (double *)    I->data;
  int i;
  
  En = ComputeEnergy(nB, B_c, nI, I_c);
  
  return PyFloat_FromDouble(En);
}








static struct PyMethodDef methods[] = {
  {"CalculateEnergy", (PyCFunction)Cpp2Py_CalculateEnergy, METH_VARARGS, NULL},


    
  {NULL, NULL, 0, NULL} /* sentinel */
};


/* --- Python 3 --- */
// This is the default structure of a python module.
// For now we only add methods which is the PartitionData function.
static struct PyModuleDef module = {
  PyModuleDef_HEAD_INIT,
  "Cpp2Py",       // m_name
  NULL,         // m_doc
  -1,           // m_size
  methods,      // m_methods
  NULL,         // m_reload
  NULL,         // m_traverse
  NULL,         // m_clear
  NULL          // m_free
};

// This constructs the PyModuleDef that is defined above.
PyMODINIT_FUNC
PyInit_Cpp2Py(void)
{
  PyObject *m = NULL;

  /* Initialize mpi4py's C-API */
  if (import_mpi4py() < 0) goto bad;

  /* Module initialization  */
  m = PyModule_Create(&module);
  if (m == NULL) goto bad;

  return m;

 bad:
  return NULL;
}
