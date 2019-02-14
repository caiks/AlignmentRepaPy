#include "AlignmentForeign.c"
#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdio.h>

void printArrayInt(PyObject* arr)
{
    npy_intp i;

    printf("NDIM(arr)=%d\n", PyArray_NDIM(arr));
    printf("SIZE(arr)=%d\n", PyArray_SIZE(arr));
    printf("ITEMSIZE(arr)=%d\n", PyArray_ITEMSIZE(arr));
    printf("STRIDES[0](arr)=%d\n", PyArray_STRIDES(arr)[0]);
    for (i = 0; i<PyArray_SIZE(arr); i++)
        printf("arr[%d]=%d\n", i, *(int*)((char*)(PyArray_DATA(arr)) + i * sizeof(int)));
}

void printArrayDouble(PyObject* arr)
{
    npy_intp i;

    printf("NDIM(arr)=%d\n", PyArray_NDIM(arr));
    printf("SIZE(arr)=%d\n", PyArray_SIZE(arr));
    printf("ITEMSIZE(arr)=%d\n", PyArray_ITEMSIZE(arr));
    printf("STRIDES[0](arr)=%d\n", PyArray_STRIDES(arr)[0]);
    for (i = 0; i<PyArray_SIZE(arr); i++)
        printf("arr[%d]=%.2f\n", i, *(double*)((char*)(PyArray_DATA(arr)) + i * sizeof(double)));
}

static PyObject *
ctest(PyObject *self, PyObject *args)
{
    printf("test ok\n");
    Py_RETURN_NONE;
}

static PyObject *
nptest(PyObject *dummy, PyObject *args)
{
    PyObject *arg1 = NULL, *arg2 = NULL, *out = NULL;
    PyObject *arr1 = NULL, *arr2 = NULL, *oarr = NULL;

    npy_intp i;

    if (!PyArg_ParseTuple(args, "OOO!", &arg1, &arg2,
        &PyArray_Type, &out)) return NULL;

    arr1 = PyArray_FROM_OTF(arg1, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
//    arr1 = PyArray_FROM_OTF(arg1, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (arr1 == NULL) return NULL;
    arr2 = PyArray_FROM_OTF(arg2, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (arr2 == NULL) goto fail;
#if NPY_API_VERSION >= 0x0000000c
    oarr = PyArray_FROM_OTF(out, NPY_DOUBLE, NPY_ARRAY_INOUT_ARRAY2);
#else
    oarr = PyArray_FROM_OTF(out, NPY_DOUBLE, NPY_ARRAY_INOUT_ARRAY);
#endif
    if (oarr == NULL) goto fail;

    printf("NDIM(arr1)=%d\n", PyArray_NDIM(arr1));
    printf("SIZE(arr1)=%d\n", PyArray_SIZE(arr1));
    printf("ITEMSIZE(arr1)=%d\n", PyArray_ITEMSIZE(arr1));
    printf("STRIDES[0](arr1)=%d\n", PyArray_STRIDES(arr1)[0]);
    printf("sizeof(double)=%d\n", sizeof(double));
    for (i=0;i<PyArray_SIZE(arr1);i++)
    {
        printf("arr1[i]=%.2f\n", *(double*)((char*)(PyArray_DATA(arr1)) + i * sizeof(double)));
        printf("arr1[i]=%.2f\n", *(double*)(PyArray_GETPTR1(arr1,i)));
    }

    printf("NDIM(arr1)=%d\n", PyArray_NDIM(oarr));
    printf("SIZE(arr1)=%d\n", PyArray_SIZE(oarr));
    printf("ITEMSIZE(arr1)=%d\n", PyArray_ITEMSIZE(oarr));
    printf("STRIDES[0](arr1)=%d\n", PyArray_STRIDES(oarr)[0]);
    printf("sizeof(double)=%d\n", sizeof(double));
    *(double*)((char*)(PyArray_DATA(oarr)) + 2 * sizeof(double)) = 4.567;
    *(double*)(PyArray_GETPTR1(oarr, 3)) = 1.23;
    for (i = 0; i<PyArray_SIZE(oarr); i++)
    {
        printf("oarr[i]=%.2f\n", *(double*)((char*)(PyArray_DATA(oarr)) + i * sizeof(double)));
        printf("oarr[i]=%.2f\n", *(double*)(PyArray_GETPTR1(oarr, i)));
    }


    /* code that makes use of arguments */
    /* You will probably need at least
    nd = PyArray_NDIM(<..>)    -- number of dimensions
    dims = PyArray_DIMS(<..>)  -- npy_intp array of length nd
    showing length in each dim.
    dptr = (double *)PyArray_DATA(<..>) -- pointer to data.

    If an error occurs goto fail.
    */

    Py_DECREF(arr1);
    Py_DECREF(arr2);
#if NPY_API_VERSION >= 0x0000000c
    PyArray_ResolveWritebackIfCopy(oarr);
#endif
    Py_DECREF(oarr);
    Py_INCREF(Py_None);
    return Py_None;

fail:
    Py_XDECREF(arr1);
    Py_XDECREF(arr2);
#if NPY_API_VERSION >= 0x0000000c
    PyArray_DiscardWritebackIfCopy(oarr);
#endif
    Py_XDECREF(oarr);
    return NULL;
}

/* 
    void listVarsArrayHistoriesReduce_u(double f, int n, int* ppkk, int* pskk, int z, int* prr, double* pmv)
*/
static PyObject *
listVarsArrayHistoriesReduce_u_py(PyObject *dummy, PyObject *args)
{
    double f;
    int n;
    int z;
    PyObject *arg_ppkk = NULL, *arg_pskk = NULL, *arg_prr = NULL, *arg_pmv = NULL;
    PyObject *ppkk = NULL, *pskk = NULL, *prr = NULL, *pmv = NULL;

    if (!PyArg_ParseTuple(args, "diO!O!iO!O!", 
            &f, 
            &n, 
            &PyArray_Type, &arg_ppkk, 
            &PyArray_Type, &arg_pskk, 
            &z, 
            &PyArray_Type, &arg_prr, 
            &PyArray_Type, &arg_pmv)) 
        return NULL;

    ppkk = PyArray_FROM_OTF(arg_ppkk, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (ppkk == NULL) return NULL;
    pskk = PyArray_FROM_OTF(arg_pskk, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (pskk == NULL) return NULL;
    prr = PyArray_FROM_OTF(arg_prr, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (prr == NULL) return NULL;
#if NPY_API_VERSION >= 0x0000000c
    pmv = PyArray_FROM_OTF(arg_pmv, NPY_FLOAT64, NPY_ARRAY_INOUT_ARRAY2);
#else
    pmv = PyArray_FROM_OTF(arg_pmv, NPY_FLOAT64, NPY_ARRAY_INOUT_ARRAY);
#endif
    if (pmv == NULL) goto fail;
/*
    printf("f=%g\n", f);
    printf("n=%d\n", n);
    printf("z=%d\n", z);
    printArrayInt(ppkk);
    printArrayInt(pskk);
    printArrayInt(prr);
*/
    listVarsArrayHistoriesReduce_u(
        f, 
        n, 
        (int*)(PyArray_DATA(ppkk)), 
        (int*)(PyArray_DATA(pskk)), 
        z, 
        (int*)(PyArray_DATA(prr)), 
        (double*)(PyArray_DATA(pmv)));
/*
    printArrayDouble(pmv);
*/
    Py_DECREF(ppkk);
    Py_DECREF(pskk);
    Py_DECREF(prr);
#if NPY_API_VERSION >= 0x0000000c
    PyArray_ResolveWritebackIfCopy(pmv);
#endif
    Py_DECREF(pmv);
    Py_INCREF(Py_None);
    return Py_None;

fail:
    Py_XDECREF(ppkk);
    Py_XDECREF(pskk);
    Py_XDECREF(prr);
#if NPY_API_VERSION >= 0x0000000c
    PyArray_DiscardWritebackIfCopy(pmv);
#endif
    Py_XDECREF(pmv);
    return NULL;
}

/*
int listVarsArrayHistoriesAlignedTop_u(
    int xmax, int omax, int n, int* svv, int m, int z1, int z2,
    int* ppww, int* phh1, double* pxx1, int* phh2, double* pxx2,
    int* tww1, int* tww2, double* ts1, double* ts2, int* ts3, int* s)
*/
static PyObject *
listVarsArrayHistoriesAlignedTop_u_py(PyObject *dummy, PyObject *args)
{
    int xmax;
    int omax;
    int n;
    PyObject *arg_svv = NULL, *svv = NULL;
    int m;
    int z1;
    int z2;
    PyObject *arg_ppww = NULL, *ppww = NULL;
    PyObject *arg_phh1 = NULL, *phh1 = NULL;
    PyObject *arg_pxx1 = NULL, *pxx1 = NULL;
    PyObject *arg_phh2 = NULL, *phh2 = NULL;
    PyObject *arg_pxx2 = NULL, *pxx2 = NULL;
    PyObject *arg_tww1 = NULL, *tww1 = NULL;
    PyObject *arg_tww2 = NULL, *tww2 = NULL;
    PyObject *arg_ts1 = NULL, *ts1 = NULL;
    PyObject *arg_ts2 = NULL, *ts2 = NULL;
    PyObject *arg_ts3 = NULL, *ts3 = NULL;
    PyObject *arg_s = NULL, *s = NULL;
    PyObject *arg_t = NULL, *t = NULL;
    int t1 = 0;

    if (!PyArg_ParseTuple(args, "iiiO!iiiO!O!O!O!O!O!O!O!O!O!O!O!",
        &xmax,
        &omax,
        &n,
        &PyArray_Type, &arg_svv,
        &m,
        &z1,
        &z2,
        &PyArray_Type, &arg_ppww,
        &PyArray_Type, &arg_phh1,
        &PyArray_Type, &arg_pxx1,
        &PyArray_Type, &arg_phh2,
        &PyArray_Type, &arg_pxx2,
        &PyArray_Type, &arg_tww1,
        &PyArray_Type, &arg_tww2,
        &PyArray_Type, &arg_ts1,
        &PyArray_Type, &arg_ts2,
        &PyArray_Type, &arg_ts3,
        &PyArray_Type, &arg_s,
        &PyArray_Type, &arg_t))
        return NULL;

    svv = PyArray_FROM_OTF(arg_svv, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (svv == NULL) return NULL;
    ppww = PyArray_FROM_OTF(arg_ppww, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (ppww == NULL) return NULL;
    phh1 = PyArray_FROM_OTF(arg_phh1, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (phh1 == NULL) return NULL;
    pxx1 = PyArray_FROM_OTF(arg_pxx1, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    if (pxx1 == NULL) return NULL;
    phh2 = PyArray_FROM_OTF(arg_phh2, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (phh2 == NULL) return NULL;
    pxx2 = PyArray_FROM_OTF(arg_pxx2, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    if (pxx2 == NULL) return NULL;
#if NPY_API_VERSION >= 0x0000000c
    tww1 = PyArray_FROM_OTF(arg_tww1, NPY_INT32, NPY_ARRAY_INOUT_ARRAY2);
    tww2 = PyArray_FROM_OTF(arg_tww2, NPY_INT32, NPY_ARRAY_INOUT_ARRAY2);
    ts1 = PyArray_FROM_OTF(arg_ts1, NPY_FLOAT64, NPY_ARRAY_INOUT_ARRAY2);
    ts2 = PyArray_FROM_OTF(arg_ts2, NPY_FLOAT64, NPY_ARRAY_INOUT_ARRAY2);
    ts3 = PyArray_FROM_OTF(arg_ts3, NPY_INT32, NPY_ARRAY_INOUT_ARRAY2);
    s = PyArray_FROM_OTF(arg_s, NPY_INT32, NPY_ARRAY_INOUT_ARRAY2);
    t = PyArray_FROM_OTF(arg_t, NPY_INT32, NPY_ARRAY_INOUT_ARRAY2);
#else
    tww1 = PyArray_FROM_OTF(arg_tww1, NPY_INT32, NPY_ARRAY_INOUT_ARRAY);
    tww2 = PyArray_FROM_OTF(arg_tww2, NPY_INT32, NPY_ARRAY_INOUT_ARRAY);
    ts1 = PyArray_FROM_OTF(arg_ts1, NPY_FLOAT64, NPY_ARRAY_INOUT_ARRAY);
    ts2 = PyArray_FROM_OTF(arg_ts2, NPY_FLOAT64, NPY_ARRAY_INOUT_ARRAY);
    ts3 = PyArray_FROM_OTF(arg_ts3, NPY_INT32, NPY_ARRAY_INOUT_ARRAY);
    s = PyArray_FROM_OTF(arg_s, NPY_INT32, NPY_ARRAY_INOUT_ARRAY);
    t = PyArray_FROM_OTF(arg_t, NPY_INT32, NPY_ARRAY_INOUT_ARRAY);
#endif
    if (tww1 == NULL) goto fail;
    if (tww2 == NULL) goto fail;
    if (ts1 == NULL) goto fail;
    if (ts2 == NULL) goto fail;
    if (ts3 == NULL) goto fail;
    if (s == NULL) goto fail;
    if (t == NULL) goto fail;
/*
    printf("xmax=%d\n", xmax);
    printf("omax=%d\n", omax);
    printf("n=%d\n", n);
    printArrayInt(svv);
    printf("m=%d\n", m);
    printf("z1=%d\n", z1);
    printf("z2=%d\n", z2);
    printArrayInt(ppww);
    printArrayInt(phh1);
    printArrayDouble(pxx1);
    printArrayInt(phh2);
    printArrayDouble(pxx2);
*/
    t1 = listVarsArrayHistoriesAlignedTop_u(
        xmax,
        omax,
        n,
        (int*)(PyArray_DATA(svv)),
        m,
        z1,
        z2,
        (int*)(PyArray_DATA(ppww)),
        (int*)(PyArray_DATA(phh1)),
        (double*)(PyArray_DATA(pxx1)),
        (int*)(PyArray_DATA(phh2)),
        (double*)(PyArray_DATA(pxx2)),
        (int*)(PyArray_DATA(tww1)),
        (int*)(PyArray_DATA(tww2)),
        (double*)(PyArray_DATA(ts1)),
        (double*)(PyArray_DATA(ts2)),
        (int*)(PyArray_DATA(ts3)),
        (int*)(PyArray_DATA(s)));

        *(int*)(PyArray_GETPTR1(t,0)) = t1;
/*
    printArrayInt(tww1);
    printArrayInt(tww2);
    printArrayDouble(ts1);
    printArrayDouble(ts2);
    printArrayInt(ts3);
    printArrayInt(s);
    printArrayInt(t);
*/
    Py_DECREF(svv);
    Py_DECREF(ppww);
    Py_DECREF(phh1);
    Py_DECREF(pxx1);
    Py_DECREF(phh2);
    Py_DECREF(pxx2);
#if NPY_API_VERSION >= 0x0000000c
    PyArray_ResolveWritebackIfCopy(tww1);
    PyArray_ResolveWritebackIfCopy(tww2);
    PyArray_ResolveWritebackIfCopy(ts1);
    PyArray_ResolveWritebackIfCopy(ts2);
    PyArray_ResolveWritebackIfCopy(ts3);
    PyArray_ResolveWritebackIfCopy(s);
    PyArray_ResolveWritebackIfCopy(t);
#endif
    Py_DECREF(tww1);
    Py_DECREF(tww2);
    Py_DECREF(ts1);
    Py_DECREF(ts2);
    Py_DECREF(ts3);
    Py_DECREF(s);
    Py_DECREF(t);
    Py_INCREF(Py_None);
    return Py_None;

fail:
    Py_XDECREF(svv);
    Py_XDECREF(ppww);
    Py_XDECREF(phh1);
    Py_XDECREF(pxx1);
    Py_XDECREF(phh2);
    Py_XDECREF(pxx2);
#if NPY_API_VERSION >= 0x0000000c
    PyArray_DiscardWritebackIfCopy(tww1);
    PyArray_DiscardWritebackIfCopy(tww2);
    PyArray_DiscardWritebackIfCopy(ts1);
    PyArray_DiscardWritebackIfCopy(ts2);
    PyArray_DiscardWritebackIfCopy(ts3);
    PyArray_DiscardWritebackIfCopy(s);
    PyArray_DiscardWritebackIfCopy(t);
#endif
    Py_XDECREF(tww1);
    Py_XDECREF(tww2);
    Py_XDECREF(ts1);
    Py_XDECREF(ts2);
    Py_XDECREF(ts3);
    Py_XDECREF(s);
    Py_XDECREF(t);
    return NULL;
}

/*
int listVarsListTuplesArrayHistoriesAlignedTop_u(
    int dense,
    int xmax, int omax, int n, int* svv, int m, int d, int e,
    int z1, int z2,
    int* ppww, int* ppdd,
    int* phh1, double* pxx1, int* phh2, double* pxx2,
    int* tww1, int* tww2, double* ts1, double* ts2, int* ts3, int* s)
*/
static PyObject *
listVarsListTuplesArrayHistoriesAlignedTop_u_py(PyObject *dummy, PyObject *args)
{
    int dense;
    int xmax;
    int omax;
    int n;
    PyObject *arg_svv = NULL, *svv = NULL;
    int m;
    int d;
    int e;
    int z1;
    int z2;
    PyObject *arg_ppww = NULL, *ppww = NULL;
    PyObject *arg_ppdd = NULL, *ppdd = NULL;
    PyObject *arg_phh1 = NULL, *phh1 = NULL;
    PyObject *arg_pxx1 = NULL, *pxx1 = NULL;
    PyObject *arg_phh2 = NULL, *phh2 = NULL;
    PyObject *arg_pxx2 = NULL, *pxx2 = NULL;
    PyObject *arg_tww1 = NULL, *tww1 = NULL;
    PyObject *arg_tww2 = NULL, *tww2 = NULL;
    PyObject *arg_ts1 = NULL, *ts1 = NULL;
    PyObject *arg_ts2 = NULL, *ts2 = NULL;
    PyObject *arg_ts3 = NULL, *ts3 = NULL;
    PyObject *arg_s = NULL, *s = NULL;
    PyObject *arg_t = NULL, *t = NULL;
    int t1 = 0;

    if (!PyArg_ParseTuple(args, "iiiiO!iiiiiO!O!O!O!O!O!O!O!O!O!O!O!O!",
        &dense,
        &xmax,
        &omax,
        &n,
        &PyArray_Type, &arg_svv,
        &m,
        &d,
        &e,
        &z1,
        &z2,
        &PyArray_Type, &arg_ppww,
        &PyArray_Type, &arg_ppdd,
        &PyArray_Type, &arg_phh1,
        &PyArray_Type, &arg_pxx1,
        &PyArray_Type, &arg_phh2,
        &PyArray_Type, &arg_pxx2,
        &PyArray_Type, &arg_tww1,
        &PyArray_Type, &arg_tww2,
        &PyArray_Type, &arg_ts1,
        &PyArray_Type, &arg_ts2,
        &PyArray_Type, &arg_ts3,
        &PyArray_Type, &arg_s,
        &PyArray_Type, &arg_t))
        return NULL;

    svv = PyArray_FROM_OTF(arg_svv, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (svv == NULL) return NULL;
    ppww = PyArray_FROM_OTF(arg_ppww, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (ppww == NULL) return NULL;
    ppdd = PyArray_FROM_OTF(arg_ppdd, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (ppdd == NULL) return NULL;
    phh1 = PyArray_FROM_OTF(arg_phh1, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (phh1 == NULL) return NULL;
    pxx1 = PyArray_FROM_OTF(arg_pxx1, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    if (pxx1 == NULL) return NULL;
    phh2 = PyArray_FROM_OTF(arg_phh2, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (phh2 == NULL) return NULL;
    pxx2 = PyArray_FROM_OTF(arg_pxx2, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    if (pxx2 == NULL) return NULL;
#if NPY_API_VERSION >= 0x0000000c
    tww1 = PyArray_FROM_OTF(arg_tww1, NPY_INT32, NPY_ARRAY_INOUT_ARRAY2);
    tww2 = PyArray_FROM_OTF(arg_tww2, NPY_INT32, NPY_ARRAY_INOUT_ARRAY2);
    ts1 = PyArray_FROM_OTF(arg_ts1, NPY_FLOAT64, NPY_ARRAY_INOUT_ARRAY2);
    ts2 = PyArray_FROM_OTF(arg_ts2, NPY_FLOAT64, NPY_ARRAY_INOUT_ARRAY2);
    ts3 = PyArray_FROM_OTF(arg_ts3, NPY_INT32, NPY_ARRAY_INOUT_ARRAY2);
    s = PyArray_FROM_OTF(arg_s, NPY_INT32, NPY_ARRAY_INOUT_ARRAY2);
    t = PyArray_FROM_OTF(arg_t, NPY_INT32, NPY_ARRAY_INOUT_ARRAY2);
#else
    tww1 = PyArray_FROM_OTF(arg_tww1, NPY_INT32, NPY_ARRAY_INOUT_ARRAY);
    tww2 = PyArray_FROM_OTF(arg_tww2, NPY_INT32, NPY_ARRAY_INOUT_ARRAY);
    ts1 = PyArray_FROM_OTF(arg_ts1, NPY_FLOAT64, NPY_ARRAY_INOUT_ARRAY);
    ts2 = PyArray_FROM_OTF(arg_ts2, NPY_FLOAT64, NPY_ARRAY_INOUT_ARRAY);
    ts3 = PyArray_FROM_OTF(arg_ts3, NPY_INT32, NPY_ARRAY_INOUT_ARRAY);
    s = PyArray_FROM_OTF(arg_s, NPY_INT32, NPY_ARRAY_INOUT_ARRAY);
    t = PyArray_FROM_OTF(arg_t, NPY_INT32, NPY_ARRAY_INOUT_ARRAY);
#endif
    if (tww1 == NULL) goto fail;
    if (tww2 == NULL) goto fail;
    if (ts1 == NULL) goto fail;
    if (ts2 == NULL) goto fail;
    if (ts3 == NULL) goto fail;
    if (s == NULL) goto fail;
    if (t == NULL) goto fail;
    /*
    printf("dense=%d\n", dense);
    printf("xmax=%d\n", xmax);
    printf("omax=%d\n", omax);
    printf("n=%d\n", n);
    printArrayInt(svv);
    printf("m=%d\n", m);
    printf("d=%d\n", d);
    printf("e=%d\n", e);
    printf("z1=%d\n", z1);
    printf("z2=%d\n", z2);
    printArrayInt(ppww);
    printArrayInt(ppdd);
    printArrayInt(phh1);
    printArrayDouble(pxx1);
    printArrayInt(phh2);
    printArrayDouble(pxx2);
    */
    t1 = listVarsListTuplesArrayHistoriesAlignedTop_u(
        dense,
        xmax,
        omax,
        n,
        (int*)(PyArray_DATA(svv)),
        m,
        d,
        e,
        z1,
        z2,
        (int*)(PyArray_DATA(ppww)),
        (int*)(PyArray_DATA(ppdd)),
        (int*)(PyArray_DATA(phh1)),
        (double*)(PyArray_DATA(pxx1)),
        (int*)(PyArray_DATA(phh2)),
        (double*)(PyArray_DATA(pxx2)),
        (int*)(PyArray_DATA(tww1)),
        (int*)(PyArray_DATA(tww2)),
        (double*)(PyArray_DATA(ts1)),
        (double*)(PyArray_DATA(ts2)),
        (int*)(PyArray_DATA(ts3)),
        (int*)(PyArray_DATA(s)));

    *(int*)(PyArray_GETPTR1(t, 0)) = t1;
    /*
    printArrayInt(tww1);
    printArrayInt(tww2);
    printArrayDouble(ts1);
    printArrayDouble(ts2);
    printArrayInt(ts3);
    printArrayInt(s);
    printArrayInt(t);
    */
    Py_DECREF(svv);
    Py_DECREF(ppww);
    Py_DECREF(ppdd);
    Py_DECREF(phh1);
    Py_DECREF(pxx1);
    Py_DECREF(phh2);
    Py_DECREF(pxx2);
#if NPY_API_VERSION >= 0x0000000c
    PyArray_ResolveWritebackIfCopy(tww1);
    PyArray_ResolveWritebackIfCopy(tww2);
    PyArray_ResolveWritebackIfCopy(ts1);
    PyArray_ResolveWritebackIfCopy(ts2);
    PyArray_ResolveWritebackIfCopy(ts3);
    PyArray_ResolveWritebackIfCopy(s);
    PyArray_ResolveWritebackIfCopy(t);
#endif
    Py_DECREF(tww1);
    Py_DECREF(tww2);
    Py_DECREF(ts1);
    Py_DECREF(ts2);
    Py_DECREF(ts3);
    Py_DECREF(s);
    Py_DECREF(t);
    Py_INCREF(Py_None);
    return Py_None;

fail:
    Py_XDECREF(svv);
    Py_XDECREF(ppww);
    Py_XDECREF(ppdd);
    Py_XDECREF(phh1);
    Py_XDECREF(pxx1);
    Py_XDECREF(phh2);
    Py_XDECREF(pxx2);
#if NPY_API_VERSION >= 0x0000000c
    PyArray_DiscardWritebackIfCopy(tww1);
    PyArray_DiscardWritebackIfCopy(tww2);
    PyArray_DiscardWritebackIfCopy(ts1);
    PyArray_DiscardWritebackIfCopy(ts2);
    PyArray_DiscardWritebackIfCopy(ts3);
    PyArray_DiscardWritebackIfCopy(s);
    PyArray_DiscardWritebackIfCopy(t);
#endif
    Py_XDECREF(tww1);
    Py_XDECREF(tww2);
    Py_XDECREF(ts1);
    Py_XDECREF(ts2);
    Py_XDECREF(ts3);
    Py_XDECREF(s);
    Py_XDECREF(t);
    return NULL;
}


/*
int listVarsListTuplesArrayHistoriesAlignedExcludeHiddenTop_u(
    int dense,
    int xmax, int omax, int n, int* svv, int m, int d, int e,
    int z1, int z2,
    int ccl, int* ppccd, int* ppccu,
    int* ppww, int* ppdd,
    int* phh1, double* pxx1, int* phh2, double* pxx2,
    int* tww1, int* tww2, double* ts1, double* ts2, int* ts3, int* s)
*/
static PyObject *
listVarsListTuplesArrayHistoriesAlignedExcludeHiddenTop_u_py(PyObject *dummy, PyObject *args)
{
    int dense;
    int xmax;
    int omax;
    int n;
    PyObject *arg_svv = NULL, *svv = NULL;
    int m;
    int d;
    int e;
    int z1;
    int z2;
    int ccl;
    PyObject *arg_ppccd = NULL, *ppccd = NULL;
    PyObject *arg_ppccu = NULL, *ppccu = NULL;
    PyObject *arg_ppww = NULL, *ppww = NULL;
    PyObject *arg_ppdd = NULL, *ppdd = NULL;
    PyObject *arg_phh1 = NULL, *phh1 = NULL;
    PyObject *arg_pxx1 = NULL, *pxx1 = NULL;
    PyObject *arg_phh2 = NULL, *phh2 = NULL;
    PyObject *arg_pxx2 = NULL, *pxx2 = NULL;
    PyObject *arg_tww1 = NULL, *tww1 = NULL;
    PyObject *arg_tww2 = NULL, *tww2 = NULL;
    PyObject *arg_ts1 = NULL, *ts1 = NULL;
    PyObject *arg_ts2 = NULL, *ts2 = NULL;
    PyObject *arg_ts3 = NULL, *ts3 = NULL;
    PyObject *arg_s = NULL, *s = NULL;
    PyObject *arg_t = NULL, *t = NULL;
    int t1 = 0;

    if (!PyArg_ParseTuple(args, "iiiiO!iiiiiiO!O!O!O!O!O!O!O!O!O!O!O!O!O!O!",
        &dense,
        &xmax,
        &omax,
        &n,
        &PyArray_Type, &arg_svv,
        &m,
        &d,
        &e,
        &z1,
        &z2,
        &ccl,
        &PyArray_Type, &arg_ppccd,
        &PyArray_Type, &arg_ppccu,
        &PyArray_Type, &arg_ppww,
        &PyArray_Type, &arg_ppdd,
        &PyArray_Type, &arg_phh1,
        &PyArray_Type, &arg_pxx1,
        &PyArray_Type, &arg_phh2,
        &PyArray_Type, &arg_pxx2,
        &PyArray_Type, &arg_tww1,
        &PyArray_Type, &arg_tww2,
        &PyArray_Type, &arg_ts1,
        &PyArray_Type, &arg_ts2,
        &PyArray_Type, &arg_ts3,
        &PyArray_Type, &arg_s,
        &PyArray_Type, &arg_t))
        return NULL;

    svv = PyArray_FROM_OTF(arg_svv, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (svv == NULL) return NULL;
    ppccd = PyArray_FROM_OTF(arg_ppccd, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (ppccd == NULL) return NULL;
    ppccu = PyArray_FROM_OTF(arg_ppccu, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (ppccu == NULL) return NULL;
    ppww = PyArray_FROM_OTF(arg_ppww, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (ppww == NULL) return NULL;
    ppdd = PyArray_FROM_OTF(arg_ppdd, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (ppdd == NULL) return NULL;
    phh1 = PyArray_FROM_OTF(arg_phh1, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (phh1 == NULL) return NULL;
    pxx1 = PyArray_FROM_OTF(arg_pxx1, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    if (pxx1 == NULL) return NULL;
    phh2 = PyArray_FROM_OTF(arg_phh2, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (phh2 == NULL) return NULL;
    pxx2 = PyArray_FROM_OTF(arg_pxx2, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    if (pxx2 == NULL) return NULL;
#if NPY_API_VERSION >= 0x0000000c
    tww1 = PyArray_FROM_OTF(arg_tww1, NPY_INT32, NPY_ARRAY_INOUT_ARRAY2);
    tww2 = PyArray_FROM_OTF(arg_tww2, NPY_INT32, NPY_ARRAY_INOUT_ARRAY2);
    ts1 = PyArray_FROM_OTF(arg_ts1, NPY_FLOAT64, NPY_ARRAY_INOUT_ARRAY2);
    ts2 = PyArray_FROM_OTF(arg_ts2, NPY_FLOAT64, NPY_ARRAY_INOUT_ARRAY2);
    ts3 = PyArray_FROM_OTF(arg_ts3, NPY_INT32, NPY_ARRAY_INOUT_ARRAY2);
    s = PyArray_FROM_OTF(arg_s, NPY_INT32, NPY_ARRAY_INOUT_ARRAY2);
    t = PyArray_FROM_OTF(arg_t, NPY_INT32, NPY_ARRAY_INOUT_ARRAY2);
#else
    tww1 = PyArray_FROM_OTF(arg_tww1, NPY_INT32, NPY_ARRAY_INOUT_ARRAY);
    tww2 = PyArray_FROM_OTF(arg_tww2, NPY_INT32, NPY_ARRAY_INOUT_ARRAY);
    ts1 = PyArray_FROM_OTF(arg_ts1, NPY_FLOAT64, NPY_ARRAY_INOUT_ARRAY);
    ts2 = PyArray_FROM_OTF(arg_ts2, NPY_FLOAT64, NPY_ARRAY_INOUT_ARRAY);
    ts3 = PyArray_FROM_OTF(arg_ts3, NPY_INT32, NPY_ARRAY_INOUT_ARRAY);
    s = PyArray_FROM_OTF(arg_s, NPY_INT32, NPY_ARRAY_INOUT_ARRAY);
    t = PyArray_FROM_OTF(arg_t, NPY_INT32, NPY_ARRAY_INOUT_ARRAY);
#endif
    if (tww1 == NULL) goto fail;
    if (tww2 == NULL) goto fail;
    if (ts1 == NULL) goto fail;
    if (ts2 == NULL) goto fail;
    if (ts3 == NULL) goto fail;
    if (s == NULL) goto fail;
    if (t == NULL) goto fail;
    /*
    printf("dense=%d\n", dense);
    printf("xmax=%d\n", xmax);
    printf("omax=%d\n", omax);
    printf("n=%d\n", n);
    printArrayInt(svv);
    printf("m=%d\n", m);
    printf("d=%d\n", d);
    printf("e=%d\n", e);
    printf("z1=%d\n", z1);
    printf("z2=%d\n", z2);
    printf("ccl=%d\n", ccl);
    printArrayInt(ppccd);
    printArrayInt(ppccu);
    printArrayInt(ppww);
    printArrayInt(ppdd);
    printArrayInt(phh1);
    printArrayDouble(pxx1);
    printArrayInt(phh2);
    printArrayDouble(pxx2);
    */
    t1 = listVarsListTuplesArrayHistoriesAlignedExcludeHiddenTop_u(
        dense,
        xmax,
        omax,
        n,
        (int*)(PyArray_DATA(svv)),
        m,
        d,
        e,
        z1,
        z2,
        ccl,
        (int*)(PyArray_DATA(ppccd)),
        (int*)(PyArray_DATA(ppccu)),
        (int*)(PyArray_DATA(ppww)),
        (int*)(PyArray_DATA(ppdd)),
        (int*)(PyArray_DATA(phh1)),
        (double*)(PyArray_DATA(pxx1)),
        (int*)(PyArray_DATA(phh2)),
        (double*)(PyArray_DATA(pxx2)),
        (int*)(PyArray_DATA(tww1)),
        (int*)(PyArray_DATA(tww2)),
        (double*)(PyArray_DATA(ts1)),
        (double*)(PyArray_DATA(ts2)),
        (int*)(PyArray_DATA(ts3)),
        (int*)(PyArray_DATA(s)));

    *(int*)(PyArray_GETPTR1(t, 0)) = t1;
    /*
    printArrayInt(tww1);
    printArrayInt(tww2);
    printArrayDouble(ts1);
    printArrayDouble(ts2);
    printArrayInt(ts3);
    printArrayInt(s);
    printArrayInt(t);
    */
    Py_DECREF(svv);
    Py_DECREF(ppccd);
    Py_DECREF(ppccu);
    Py_DECREF(ppww);
    Py_DECREF(ppdd);
    Py_DECREF(phh1);
    Py_DECREF(pxx1);
    Py_DECREF(phh2);
    Py_DECREF(pxx2);
#if NPY_API_VERSION >= 0x0000000c
    PyArray_ResolveWritebackIfCopy(tww1);
    PyArray_ResolveWritebackIfCopy(tww2);
    PyArray_ResolveWritebackIfCopy(ts1);
    PyArray_ResolveWritebackIfCopy(ts2);
    PyArray_ResolveWritebackIfCopy(ts3);
    PyArray_ResolveWritebackIfCopy(s);
    PyArray_ResolveWritebackIfCopy(t);
#endif
    Py_DECREF(tww1);
    Py_DECREF(tww2);
    Py_DECREF(ts1);
    Py_DECREF(ts2);
    Py_DECREF(ts3);
    Py_DECREF(s);
    Py_DECREF(t);
    Py_INCREF(Py_None);
    return Py_None;

fail:
    Py_XDECREF(svv);
    Py_XDECREF(ppccd);
    Py_XDECREF(ppccu);
    Py_XDECREF(ppww);
    Py_XDECREF(ppdd);
    Py_XDECREF(phh1);
    Py_XDECREF(pxx1);
    Py_XDECREF(phh2);
    Py_XDECREF(pxx2);
#if NPY_API_VERSION >= 0x0000000c
    PyArray_DiscardWritebackIfCopy(tww1);
    PyArray_DiscardWritebackIfCopy(tww2);
    PyArray_DiscardWritebackIfCopy(ts1);
    PyArray_DiscardWritebackIfCopy(ts2);
    PyArray_DiscardWritebackIfCopy(ts3);
    PyArray_DiscardWritebackIfCopy(s);
    PyArray_DiscardWritebackIfCopy(t);
#endif
    Py_XDECREF(tww1);
    Py_XDECREF(tww2);
    Py_XDECREF(ts1);
    Py_XDECREF(ts2);
    Py_XDECREF(ts3);
    Py_XDECREF(s);
    Py_XDECREF(t);
    return NULL;
}


/*
void listListVarsArrayHistoryPairsPartitionIndependent_u(
    double z, int v, int n, int* svv, int m, int r,
    int* lyy, int* syy, int* pppp, double* aa1, double* aa2,
    double* bb1, double* bb2)
*/
static PyObject *
listListVarsArrayHistoryPairsPartitionIndependent_u_py(PyObject *dummy, PyObject *args)
{
    double z;
    int v;
    int n;
    PyObject *arg_svv = NULL, *svv = NULL;
    int m;
    int r;
    PyObject *arg_lyy = NULL, *lyy = NULL;
    PyObject *arg_syy = NULL, *syy = NULL;
    PyObject *arg_pppp = NULL, *pppp = NULL;
    PyObject *arg_aa1 = NULL, *aa1 = NULL;
    PyObject *arg_aa2 = NULL, *aa2 = NULL;
    PyObject *arg_bb1 = NULL, *bb1 = NULL;
    PyObject *arg_bb2 = NULL, *bb2 = NULL;

    if (!PyArg_ParseTuple(args, "diiO!iiO!O!O!O!O!O!O!",
        &z,
        &v,
        &n,
        &PyArray_Type, &arg_svv,
        &m,
        &r,
        &PyArray_Type, &arg_lyy,
        &PyArray_Type, &arg_syy,
        &PyArray_Type, &arg_pppp,
        &PyArray_Type, &arg_aa1,
        &PyArray_Type, &arg_aa2,
        &PyArray_Type, &arg_bb1,
        &PyArray_Type, &arg_bb2))
        return NULL;

    svv = PyArray_FROM_OTF(arg_svv, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (svv == NULL) return NULL;
    lyy = PyArray_FROM_OTF(arg_lyy, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (lyy == NULL) return NULL;
    syy = PyArray_FROM_OTF(arg_syy, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (syy == NULL) return NULL;
    pppp = PyArray_FROM_OTF(arg_pppp, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (pppp == NULL) return NULL;
    aa1 = PyArray_FROM_OTF(arg_aa1, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    if (aa1 == NULL) return NULL;
    aa2 = PyArray_FROM_OTF(arg_aa2, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    if (aa2 == NULL) return NULL;
#if NPY_API_VERSION >= 0x0000000c
    bb1 = PyArray_FROM_OTF(arg_bb1, NPY_FLOAT64, NPY_ARRAY_INOUT_ARRAY2);
    bb2 = PyArray_FROM_OTF(arg_bb2, NPY_FLOAT64, NPY_ARRAY_INOUT_ARRAY2);
#else
    bb1 = PyArray_FROM_OTF(arg_bb1, NPY_FLOAT64, NPY_ARRAY_INOUT_ARRAY);
    bb2 = PyArray_FROM_OTF(arg_bb2, NPY_FLOAT64, NPY_ARRAY_INOUT_ARRAY);
#endif
    if (bb1 == NULL) goto fail;
    if (bb2 == NULL) goto fail;
 /* 
    printf("z=%g\n", z);
    printf("v=%d\n", v);
    printf("n=%d\n", n);
    printArrayInt(svv);
    printf("m=%d\n", m);
    printf("r=%d\n", r);
    printArrayInt(lyy);
    printArrayInt(syy);
    printArrayInt(pppp);
    printArrayDouble(aa1);
    printArrayDouble(aa2);
*/   
    listListVarsArrayHistoryPairsPartitionIndependent_u(
        z,
        v,
        n,
        (int*)(PyArray_DATA(svv)),
        m,
        r,
        (int*)(PyArray_DATA(lyy)),
        (int*)(PyArray_DATA(syy)),
        (int*)(PyArray_DATA(pppp)),
        (double*)(PyArray_DATA(aa1)),
        (double*)(PyArray_DATA(aa2)),
        (double*)(PyArray_DATA(bb1)),
        (double*)(PyArray_DATA(bb2)));

 /* 
    printArrayDouble(bb1);
    printArrayDouble(bb2);
 */     
    Py_DECREF(svv);
    Py_DECREF(lyy);
    Py_DECREF(syy);
    Py_DECREF(pppp);
    Py_DECREF(aa1);
    Py_DECREF(aa2);
#if NPY_API_VERSION >= 0x0000000c
    PyArray_ResolveWritebackIfCopy(bb1);
    PyArray_ResolveWritebackIfCopy(bb2);
#endif
    Py_DECREF(bb1);
    Py_DECREF(bb2);
    Py_INCREF(Py_None);
    return Py_None;

fail:
    Py_XDECREF(svv);
    Py_XDECREF(lyy);
    Py_XDECREF(syy);
    Py_XDECREF(pppp);
    Py_XDECREF(aa1);
    Py_XDECREF(aa2);
#if NPY_API_VERSION >= 0x0000000c
    PyArray_DiscardWritebackIfCopy(bb1);
    PyArray_DiscardWritebackIfCopy(bb2);
#endif
    Py_XDECREF(bb1);
    Py_XDECREF(bb2);
    return NULL;
}


/*
int listListVarsArrayHistoryPairsSetTuplePartitionTop_u(
    int pmax, double z, int v, int n, int* svv, int q, double y1,
    int* qm, int* ql, int* qs, int* qp, double* aa1, double* aa2,
    int* tt)
*/
static PyObject *
listListVarsArrayHistoryPairsSetTuplePartitionTop_u_py(PyObject *dummy, PyObject *args)
{
    int pmax;
    double z;
    int v;
    int n;
    PyObject *arg_svv = NULL, *svv = NULL;
    int q;
    double y1;
    PyObject *arg_qm = NULL, *qm = NULL;
    PyObject *arg_ql = NULL, *ql = NULL;
    PyObject *arg_qs = NULL, *qs = NULL;
    PyObject *arg_qp = NULL, *qp = NULL;
    PyObject *arg_aa1 = NULL, *aa1 = NULL;
    PyObject *arg_aa2 = NULL, *aa2 = NULL;
    PyObject *arg_tt = NULL, *tt = NULL;
    PyObject *arg_t = NULL, *t = NULL;
    int t1 = 0;

    if (!PyArg_ParseTuple(args, "idiiO!idO!O!O!O!O!O!O!O!",
        &pmax,
        &z,
        &v,
        &n,
        &PyArray_Type, &arg_svv,
        &q,
        &y1,
        &PyArray_Type, &arg_qm,
        &PyArray_Type, &arg_ql,
        &PyArray_Type, &arg_qs,
        &PyArray_Type, &arg_qp,
        &PyArray_Type, &arg_aa1,
        &PyArray_Type, &arg_aa2,
        &PyArray_Type, &arg_tt,
        &PyArray_Type, &arg_t))
        return NULL;

    svv = PyArray_FROM_OTF(arg_svv, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (svv == NULL) return NULL;
    qm = PyArray_FROM_OTF(arg_qm, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (qm == NULL) return NULL;
    ql = PyArray_FROM_OTF(arg_ql, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (ql == NULL) return NULL;
    qs = PyArray_FROM_OTF(arg_qs, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (qs == NULL) return NULL;
    qp = PyArray_FROM_OTF(arg_qp, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (qp == NULL) return NULL;
    aa1 = PyArray_FROM_OTF(arg_aa1, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    if (aa1 == NULL) return NULL;
    aa2 = PyArray_FROM_OTF(arg_aa2, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    if (aa2 == NULL) return NULL;
#if NPY_API_VERSION >= 0x0000000c
    tt = PyArray_FROM_OTF(arg_tt, NPY_INT32, NPY_ARRAY_INOUT_ARRAY2);
    t = PyArray_FROM_OTF(arg_t, NPY_INT32, NPY_ARRAY_INOUT_ARRAY2);
#else
    tt = PyArray_FROM_OTF(arg_tt, NPY_INT32, NPY_ARRAY_INOUT_ARRAY);
    t = PyArray_FROM_OTF(arg_t, NPY_INT32, NPY_ARRAY_INOUT_ARRAY);
#endif
    if (tt == NULL) goto fail;
    if (t == NULL) goto fail;
    /*
    printf("pmax=%d\n", dense);
    printf("z=%g\n", xmax);
    printf("v=%d\n", v);
    printf("n=%d\n", n);
    printArrayInt(svv);
    printf("y1=%g\n", y1);
    printArrayInt(qm);
    printArrayInt(ql);
    printArrayInt(qs);
    printArrayInt(qp);
    printArrayDouble(aa1);
    printArrayDouble(aa2);
    */
    t1 = listListVarsArrayHistoryPairsSetTuplePartitionTop_u(
        pmax,
        z,
        v,
        n,
        (int*)(PyArray_DATA(svv)),
        q,
        y1,
        (int*)(PyArray_DATA(qm)),
        (int*)(PyArray_DATA(ql)),
        (int*)(PyArray_DATA(qs)),
        (int*)(PyArray_DATA(qp)),
        (double*)(PyArray_DATA(aa1)),
        (double*)(PyArray_DATA(aa2)),
        (int*)(PyArray_DATA(tt)));

    *(int*)(PyArray_GETPTR1(t, 0)) = t1;
    /*
    printArrayInt(tt);
    printArrayInt(t);
    */
    Py_DECREF(svv);
    Py_DECREF(qm);
    Py_DECREF(ql);
    Py_DECREF(qs);
    Py_DECREF(qp);
    Py_DECREF(aa1);
    Py_DECREF(aa2);
#if NPY_API_VERSION >= 0x0000000c
    PyArray_ResolveWritebackIfCopy(tt);
    PyArray_ResolveWritebackIfCopy(t);
#endif
    Py_DECREF(tt);
    Py_DECREF(t);
    Py_INCREF(Py_None);
    return Py_None;

fail:
    Py_XDECREF(svv);
    Py_XDECREF(qm);
    Py_XDECREF(ql);
    Py_XDECREF(qs);
    Py_XDECREF(qp);
    Py_XDECREF(aa1);
    Py_XDECREF(aa2);
#if NPY_API_VERSION >= 0x0000000c
    PyArray_DiscardWritebackIfCopy(tt);
    PyArray_DiscardWritebackIfCopy(t);
#endif
    Py_XDECREF(tt);
    Py_XDECREF(t);
    return NULL;
}


/*
int arrayHistoryPairsRollMax_u(
    int v, int n, int* svvy, int d, int nd,
    double* aay, double* aaxy, double* bby, double* bbxy,
    int* ppm)
*/
static PyObject *
arrayHistoryPairsRollMax_u_py(PyObject *dummy, PyObject *args)
{
    int v;
    int n;
    PyObject *arg_svvy = NULL, *svvy = NULL;
    int d;
    int nd;
    PyObject *arg_aay = NULL, *aay = NULL;
    PyObject *arg_aaxy = NULL, *aaxy = NULL;
    PyObject *arg_bby = NULL, *bby = NULL;
    PyObject *arg_bbxy = NULL, *bbxy = NULL;
    PyObject *arg_ppm = NULL, *ppm = NULL;
    PyObject *arg_t = NULL, *t = NULL;
    int t1 = 0;

    if (!PyArg_ParseTuple(args, "iiO!iiO!O!O!O!O!O!",
        &v,
        &n,
        &PyArray_Type, &arg_svvy,
        &d,
        &nd,
        &PyArray_Type, &arg_aay,
        &PyArray_Type, &arg_aaxy,
        &PyArray_Type, &arg_bby,
        &PyArray_Type, &arg_bbxy,
        &PyArray_Type, &arg_ppm,
        &PyArray_Type, &arg_t))
        return NULL;

    svvy = PyArray_FROM_OTF(arg_svvy, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (svvy == NULL) return NULL;
    aay = PyArray_FROM_OTF(arg_aay, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    if (aay == NULL) return NULL;
    aaxy = PyArray_FROM_OTF(arg_aaxy, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    if (aaxy == NULL) return NULL;
    bby = PyArray_FROM_OTF(arg_bby, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    if (bby == NULL) return NULL;
    bbxy = PyArray_FROM_OTF(arg_bbxy, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    if (bbxy == NULL) return NULL;
#if NPY_API_VERSION >= 0x0000000c
    ppm = PyArray_FROM_OTF(arg_ppm, NPY_INT32, NPY_ARRAY_INOUT_ARRAY2);
    t = PyArray_FROM_OTF(arg_t, NPY_INT32, NPY_ARRAY_INOUT_ARRAY2);
#else
    ppm = PyArray_FROM_OTF(arg_ppm, NPY_INT32, NPY_ARRAY_INOUT_ARRAY);
    t = PyArray_FROM_OTF(arg_t, NPY_INT32, NPY_ARRAY_INOUT_ARRAY);
#endif
    if (ppm == NULL) goto fail;
    if (t == NULL) goto fail;
    /*
    printf("v=%d\n", v);
    printf("n=%d\n", n);
    printArrayInt(svvy);
    printf("d=%d\n", d);
    printf("nd=%d\n", nd);
    printArrayDouble(aay);
    printArrayDouble(aaxy);
    printArrayDouble(bby);
    printArrayDouble(bbxy);
    */
    t1 = arrayHistoryPairsRollMax_u(
        v,
        n,
        (int*)(PyArray_DATA(svvy)),
        d,
        nd,
        (double*)(PyArray_DATA(aay)),
        (double*)(PyArray_DATA(aaxy)),
        (double*)(PyArray_DATA(bby)),
        (double*)(PyArray_DATA(bbxy)),
        (int*)(PyArray_DATA(ppm)));

    *(int*)(PyArray_GETPTR1(t, 0)) = t1;
    /*
    printArrayInt(ppm);
    printArrayInt(t);
    */
    Py_DECREF(svvy);
    Py_DECREF(aay);
    Py_DECREF(aaxy);
    Py_DECREF(bby);
    Py_DECREF(bbxy);
#if NPY_API_VERSION >= 0x0000000c
    PyArray_ResolveWritebackIfCopy(ppm);
    PyArray_ResolveWritebackIfCopy(t);
#endif
    Py_DECREF(ppm);
    Py_DECREF(t);
    Py_INCREF(Py_None);
    return Py_None;

fail:
    Py_XDECREF(svvy);
    Py_XDECREF(aay);
    Py_XDECREF(aaxy);
    Py_XDECREF(bby);
    Py_XDECREF(bbxy);
#if NPY_API_VERSION >= 0x0000000c
    PyArray_DiscardWritebackIfCopy(ppm);
    PyArray_DiscardWritebackIfCopy(t);
#endif
    Py_XDECREF(ppm);
    Py_XDECREF(t);
    return NULL;
}

static PyMethodDef AlignmentForeignPyMethods[] = {
    { "ctest",  ctest, METH_VARARGS, "c test" },
    { "nptest",  nptest, METH_VARARGS, "numpy test" },
    { "listVarsArrayHistoriesReduce_u_py",  listVarsArrayHistoriesReduce_u_py, METH_VARARGS, "listVarsArrayHistoriesReduce_u_py" },
    { "listVarsArrayHistoriesAlignedTop_u_py",  listVarsArrayHistoriesAlignedTop_u_py, METH_VARARGS, "listVarsArrayHistoriesAlignedTop_u_py" },
    { "listVarsListTuplesArrayHistoriesAlignedTop_u_py",  listVarsListTuplesArrayHistoriesAlignedTop_u_py, METH_VARARGS, "listVarsListTuplesArrayHistoriesAlignedTop_u_py" },
    { "listVarsListTuplesArrayHistoriesAlignedExcludeHiddenTop_u_py",  listVarsListTuplesArrayHistoriesAlignedExcludeHiddenTop_u_py, METH_VARARGS, "listVarsListTuplesArrayHistoriesAlignedExcludeHiddenTop_u_py" },
    { "listListVarsArrayHistoryPairsPartitionIndependent_u_py",  listListVarsArrayHistoryPairsPartitionIndependent_u_py, METH_VARARGS, "listListVarsArrayHistoryPairsPartitionIndependent_u_py" },
    { "listListVarsArrayHistoryPairsSetTuplePartitionTop_u_py",  listListVarsArrayHistoryPairsSetTuplePartitionTop_u_py, METH_VARARGS, "listListVarsArrayHistoryPairsSetTuplePartitionTop_u_py" },
    { "arrayHistoryPairsRollMax_u_py",  arrayHistoryPairsRollMax_u_py, METH_VARARGS, "arrayHistoryPairsRollMax_u_py" },
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef AlignmentForeignPy = {
    PyModuleDef_HEAD_INIT,
    "AlignmentForeignPy",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    AlignmentForeignPyMethods
};

PyMODINIT_FUNC PyInit_AlignmentForeignPy(void) {
    import_array();
    return PyModule_Create(&AlignmentForeignPy);
}

