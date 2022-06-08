#include <Python.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#define NPY_NO_DEPRECATED_API NPY_1_9_API_VERSION
#include <numpy/arrayobject.h>   /* NumPy  as seen from C */
//#if defined(__MACH__)
//    #include <stdlib.h>
//#else
//    #include <malloc.h>
//#endif
#define T uint64_t
#define BYTE uint8_t // unsigned char

static PyObject* bitSet(PyObject* self, PyObject *args) {
   T v;
   unsigned int b;
   if (!PyArg_ParseTuple(args, "Ki", &v, &b)) {
      return NULL;
   }
   bool c = ((v >> b) & 1) == 1;
   PyObject *return_val = Py_BuildValue("b", c);
   return return_val;
}

static inline int bitsSet_internal(T v) {
   v = v - ((v >> 1) & (T)~(T)0/3);                           // temp
   v = (v & (T)~(T)0/15*3) + ((v >> 2) & (T)~(T)0/15*3);      // temp
   v = (v + (v >> 4)) & (T)~(T)0/255*15;                      // temp
   int c = (T)(v * ((T)~(T)0/255)) >> (sizeof(T) - 1) * CHAR_BIT; // count
   return c;
}

static PyObject* bitsSet(PyObject* self, PyObject *args) {
   T v;
   if (!PyArg_ParseTuple(args, "K", &v)) {
      return NULL;
   }
   int c = bitsSet_internal(v);
   PyObject *return_val = Py_BuildValue("i", c);
   return return_val;
}

static inline T grayCode_internal(T x) {
   T c = (x >> 1) ^ x;
   return c;
}

static PyObject* grayCode(PyObject* self, PyObject *args) {
   uint64_t x;
   if (!PyArg_ParseTuple(args, "K", &x)) {
      return NULL;
   }
   T c = grayCode_internal(x);
   PyObject *return_val = Py_BuildValue("l", c);
   return return_val;
}

static PyObject* buildGraySequence(PyObject* self, PyObject *args) {
    int i = 0;
    T x = 0;

    int length, b;
    if (!PyArg_ParseTuple(args, "ii", &length, &b)) {
       return NULL;
    }
    npy_intp dims = length;
    //dims[0] = length;
    PyObject *result = PyArray_SimpleNew(1, &dims, NPY_ULONGLONG); // PyArray_LONG);
    T* resultDataPtr = (T*)(PyArray_DATA((PyArrayObject*)result));
    while (1==1) {
        T g = grayCode_internal(x);
        if (bitsSet_internal(g) == b) {
            resultDataPtr[i++] = g;
            if (i >= length) {
                break;
            }
        }
        x += 1;
    }
    return result;
}

static void do_xor_bool(bool* arr_a, bool* arr_b, int length, bool* outArr) {
    for (int i = 0; i < length; i++){ // we might have to change to non-continuous assumption...
        outArr[i] = arr_a[i] ^ arr_b[i];
    }
}


static void printBin(BYTE in) {
    for (int i = 0; i < 8 ; i++) {
        PySys_WriteStdout("%d ", (bool)(in & (1 << i)));
    }
}

static void do_xor_byte(BYTE* arr_a, BYTE* arr_b, npy_intp length, BYTE* outArr) {
    for (unsigned int i = 0; i < length; i++){
        outArr[i] = arr_a[i] ^ arr_b[i];
    }
}

static void byte_printmatrix(PyArrayObject* A){
   npy_intp dims_a_0 = PyArray_DIM(A,0);
   npy_intp dims_a_1 = PyArray_DIM(A,1);
   for (int i=0;i<dims_a_0;i++) {
      PySys_WriteStdout("| ");
      for (int j=0;j<dims_a_1;j++) {
         PySys_WriteStdout("%i ", ((BYTE*)PyArray_GETPTR2(A,i,j))[0]);
      }
      PySys_WriteStdout(" |\n");
   }
}

static void printmatrix(PyArrayObject* A){
   npy_intp dims_a_0 = PyArray_DIM(A,0);
   npy_intp dims_a_1 = PyArray_DIM(A,1);
   for (int i=0;i<dims_a_0;i++) {
      PySys_WriteStdout("| ");
      for (int j=0;j<dims_a_1;j++) {
         PySys_WriteStdout("%d ", ((bool*)PyArray_GETPTR2(A,i,j))[0]);
      }
      PySys_WriteStdout(" |\n");
   }
}

static bool isSolvable(PyArrayObject* A) {
   npy_intp dims_a_0 = PyArray_DIM(A,0);
   npy_intp dims_a_1 = PyArray_DIM(A,1);
   return dims_a_0 >= dims_a_1;
}

static PyObject* elimination(PyObject *self, PyObject *args)
{
   PyArrayObject *A, *b, *packet_mapping;
   if (!PyArg_ParseTuple(args, "O!O!O!", &PyArray_Type, &A, &PyArray_Type, &b, &PyArray_Type, &packet_mapping)) {
      PyErr_BadArgument();
      return NULL;
   }
   if (!isSolvable(A)) {
      Py_BuildValue("O", Py_False);
   }
   npy_intp dims_a_0 = PyArray_DIM(A,0); // rows
   npy_intp dims_a_1 = PyArray_DIM(A,1); // columns
   npy_intp dims_b_1 = PyArray_DIM(b,1);
   bool *dirty_rows = PyMem_RawMalloc(((unsigned int)dims_a_0) * sizeof(char));
   if (dirty_rows == NULL)
       return PyErr_NoMemory();
   for (int i = 0; i < dims_a_0; i++) {
       dirty_rows[i] = false;
   }
   bool dirty = false;
   uint8_t num_dirty_rows = 0;
    for (uint32_t i = 0; i < dims_a_1; i++) {
        for (uint32_t j = i - num_dirty_rows; j < dims_a_0; j++) {
            if (*((bool*)PyArray_GETPTR2(A, j, i)) ) {
                if (i == j)
                    //A[i,i] is true, no need to swap rows
                    break;
                // we found the first row j with A[j,i] = True. char_xor it with row i: A[i] = A[i] ^ A[j]
                // we go for a xor-swap
                do_xor_bool((bool*)PyArray_GETPTR2(A, i,0), (bool*)PyArray_GETPTR2(A, j,0),
                         dims_a_1, (bool*)PyArray_GETPTR2(A, i,0));
                do_xor_bool((bool*)PyArray_GETPTR2(A,i,0),  (bool*)PyArray_GETPTR2(A, j,0),
                         dims_a_1, (bool*)PyArray_GETPTR2(A,j,0));
                do_xor_bool((bool*)PyArray_GETPTR2(A, i,0), (bool*)PyArray_GETPTR2(A, j,0),
                         dims_a_1, (bool*)PyArray_GETPTR2(A, i,0));
                do_xor_byte((BYTE*)PyArray_GETPTR2(b,i,0), (BYTE*)PyArray_GETPTR2(b,j,0),
                         dims_b_1, (BYTE*)PyArray_GETPTR2(b,i,0));
                do_xor_byte((BYTE*)PyArray_GETPTR2(b,j,0), (BYTE*)PyArray_GETPTR2(b,i,0),
                         dims_b_1, (BYTE*)PyArray_GETPTR2(b,j,0));
                do_xor_byte((BYTE*)PyArray_GETPTR2(b,i,0), (BYTE*)PyArray_GETPTR2(b,j,0),
                         dims_b_1, (BYTE*)PyArray_GETPTR2(b,i,0));
                // swap packet_mapping...
                unsigned long tmp = (*(unsigned long*)PyArray_GETPTR1(packet_mapping, i));
                //if (i >= dims_mapping || j >= dims_mapping)
                //    PySys_WriteStdout("Packet Mappings ( dim= %lu ): i = %u, j = %u\n", dims_mapping, i, j);

                //npy_intp stride_i = PyArray_STRIDE(packet_mapping, i);
                //npy_intp stride_j = PyArray_STRIDE(packet_mapping, j);
                //const char * dataptr = PyArray_BYTES(packet_mapping);
                //PyObject * p1 = PyArray_GETITEM(packet_mapping, dataptr);
                PyObject * old_i = PyArray_GETITEM(packet_mapping, PyArray_GETPTR1(packet_mapping,i));
                PyObject * old_j = PyArray_GETITEM(packet_mapping, PyArray_GETPTR1(packet_mapping,j));

                //PySys_WriteStdout("mapping.dtype: %i", PyArray_TYPE(packet_mapping));
                //PySys_WriteStdout("mapping[%u] = %lu", i, old_i);
                //PySys_WriteStdout("mapping[%u] = %lu", j, old_j);

                //(*(unsigned long*)PyArray_GETPTR1(packet_mapping, i)) = (*(unsigned long*)PyArray_GETITEM(packet_mapping, j));
                //(*(unsigned long*)PyArray_GETPTR1(packet_mapping, j)) = tmp;
                PyArray_SETITEM(packet_mapping,PyArray_GETPTR1(packet_mapping, j), old_i);
                PyArray_SETITEM(packet_mapping,PyArray_GETPTR1(packet_mapping, i), old_j);
                //PySys_WriteStdout("after swap: mapping[%u] = %lu", i, (*(unsigned long*)PyArray_GETPTR1(packet_mapping, i)));
                //PySys_WriteStdout("after swap: mapping[%u] = %lun\n", j, (*(unsigned long*)PyArray_GETPTR1(packet_mapping, j)));
                break;
            }

        }
        // after this step we should have A[i][i] == true
        // ( IF the Matrix is singular we might have no "true" in column i. )
        // but we might be able to retrieve as many blocks as possible
        if (!*((bool*)PyArray_GETPTR2(A,i,i))) {
            PySys_WriteStdout("Could not decode Chunk %u\n", i);
            dirty_rows[i] = true;
            dirty = true;
            num_dirty_rows++;
            continue;
        }
        // eliminate from top to bottom
        for (uint32_t j = i + 1; j < dims_a_0; j++) {
            if (dirty_rows[j]) {
                continue;
            }
            if (*((bool*)PyArray_GETPTR2(A,j,i))) {
                do_xor_bool((bool*)PyArray_GETPTR2(A,j,0),
                            (bool*)PyArray_GETPTR2(A,i,0), dims_a_1, (bool*)PyArray_GETPTR2(A,j,0));
                do_xor_byte((BYTE*)PyArray_GETPTR2(b,j,0), (BYTE*)PyArray_GETPTR2(b,i,0),
                          dims_b_1, (BYTE*)PyArray_GETPTR2(b,j,0));
            }
        }

    }
    // eliminate from bottom to top
    for (int32_t col = dims_a_1 - 1; col >= 0; col--) {
        if (dirty_rows[col]) {
            continue; //skip this column if it was marked as dirty previously
        }
        for (int32_t row = dims_a_0 - 1; row >= 0; row--) {
            if (dirty_rows[row]) {
                continue; //skip this column if it was marked as dirty previously
            }
            /*int32_t lim = col;
            if (dirty) {
                lim = 0;
            }*/
            if (*((bool*)PyArray_GETPTR2(A,row,col)) && col != row) {
                do_xor_bool((bool*)PyArray_GETPTR2(A,row,0),
                             (bool*)PyArray_GETPTR2(A,col,0),
                             dims_a_1, (bool*)PyArray_GETPTR2(A,row,0));
                // char_xor the content
                do_xor_byte((BYTE*)PyArray_GETPTR2(b,row,0), (BYTE*)PyArray_GETPTR2(b,col,0),
                            dims_b_1, (BYTE*)PyArray_GETPTR2(b,row,0));
            }
        }
   }
   PyMem_RawFree(dirty_rows);
   if (dirty) {
       return Py_BuildValue("O", Py_False);
   } else {
       return Py_BuildValue("O", Py_True);
   }
}

static PyObject* xor_array(PyObject *self, PyObject *args)
{
   PyArrayObject *X, *Y;
   PyObject *out;
   int i;

   if (!PyArg_ParseTuple(args, "O!O!", &PyArray_Type, &X, &PyArray_Type, &Y)) {
      PyErr_BadArgument();
      return NULL;
   }
   npy_intp dims_x = PyArray_DIM(X,0);
   if (dims_x != PyArray_DIM(Y,0)) { // arrays have to be same size
      PyErr_SetString(PyExc_TypeError, "input dimensions differ");
      return NULL;
   }
   BYTE* x_DataPtr = (BYTE*)(PyArray_DATA((PyArrayObject*)X));
   BYTE* y_DataPtr = (BYTE*)(PyArray_DATA((PyArrayObject*)Y));


   npy_intp dims[1];
   dims[0] = dims_x;
   out = PyArray_SimpleNew(1, dims, NPY_BYTE); // we can treat any input as byte for bitwise xor...
    BYTE* out_DataPtr = (BYTE*)(PyArray_DATA((PyArrayObject*)out));
   for (i=0; i<dims_x;i++) {
      out_DataPtr[i] = x_DataPtr[i] ^ y_DataPtr[i];
   }
   return out;
}

static PyObject* microsatellite(PyObject* self,  PyObject *args)
{
   char *text;
   int lengthToLookFor;
   if (!PyArg_ParseTuple(args, "si", &text, &lengthToLookFor)) {
      return NULL;
   }
   int i = 0;
   int n = strlen(text);
   int res = 1;
   char *resChars = PyMem_RawMalloc((lengthToLookFor + 1) * sizeof(char));
   if (resChars == NULL)
       return PyErr_NoMemory();
   strncpy( resChars, text, lengthToLookFor );
   resChars[lengthToLookFor] = '\0';
   //text[:lengthToLookFor];
   int maxLength = 0;
   int ret;
   while (i <= n - 2 * lengthToLookFor) {
       ret = memcmp(&text[i],&text[i+lengthToLookFor], lengthToLookFor);
       if (ret == 0) {
           // we found one
            res += 1;
        }else{
            if (maxLength < res) {
                maxLength = res;
                //resChars = text[i: i + lengthToLookFor];
                strncpy( resChars, &text[i], lengthToLookFor );
            }
            res = 1;
        }
        i += lengthToLookFor;
   }
   if (maxLength < res) {
        maxLength = res;
        //resChars = text[i: i + lengthToLookFor];
        strncpy( resChars, &text[i], lengthToLookFor );
    }
   PyObject *return_val = Py_BuildValue("(is)", maxLength, resChars);
   PyMem_RawFree(resChars);
   return return_val;
}

static PyObject* longestSequenceOfChar(PyObject* self,  PyObject *args)
{
   char *text; //dna_seq
   char *char_x; // dna base to check (* for all)
   if (!PyArg_ParseTuple(args, "ss", &text, &char_x)) {
      return NULL;
   }
   if (strlen(char_x) != 1) {
      return NULL;
   }
   int c = 0; // max homopolymer length found so far
   char res = char_x[0];  // text[0]
   int n = strlen(text);
   int curr = 1;
   int i = 0;
   while (i < n-1) {
        if (i < n - 1 && text[i] == text[i+1]) {
           curr += 1;
        } else {
            if (curr > c) {
                if (char_x[0] == '*' || text[i] == char_x[0]) {
                    c = curr;
                    res = text[i];
                }
            }
            curr = 1;
        }
        i += 1;
   }
   if (curr > c && (char_x[0] == '*' || text[i] == char_x[0])) {
       c = curr;
       res = text[i];
   }
   PyObject *return_val = Py_BuildValue("(ci)", res, c);
   return return_val;
}

static PyObject* repeatRegion(PyObject* self,  PyObject *args)
{
   char *text;
   int res = 0;
   int lengthToLookFor;
   if (!PyArg_ParseTuple(args, "si", &text, &lengthToLookFor)) {
      return NULL;
   }
   int len = strlen(text);
   char *subseq = PyMem_RawMalloc((lengthToLookFor+1) * sizeof(char));
   if (subseq == NULL)
       return PyErr_NoMemory();
   strncpy( subseq, text, lengthToLookFor );
   subseq[lengthToLookFor] = '\0';
   for (int i = 0; i < len-lengthToLookFor;i++) {
        strncpy( subseq, &text[i], lengthToLookFor);
        if (strstr(&text[i+1], subseq) != NULL) {
            res = 1;
            break;
        }
   }
   PyObject *return_val = Py_BuildValue("i", res);
   PyMem_RawFree(subseq);
   return return_val;
}

static PyObject* smallRepeatRegion(PyObject* self,  PyObject *args)
{
   char *text;
   float res = 1.0;
   int lengthToLookFor;
   if (!PyArg_ParseTuple(args, "si", &text, &lengthToLookFor)) {
      return NULL;
   }
   int len = strlen(text);
   char *subseq = PyMem_RawMalloc((lengthToLookFor+1) * sizeof(char));
   if (subseq == NULL)
       return PyErr_NoMemory();
   strncpy(subseq, text, lengthToLookFor );
   subseq[lengthToLookFor] = '\0';
   for (int i = 0; i <= len-lengthToLookFor;i++) {
        strncpy(subseq, &text[i], lengthToLookFor);
        if (strstr(&text[i+1], subseq) != NULL) {
            res += 1.0;
        }
   }
   if (res * lengthToLookFor / strlen(text) > 0.44) {
       res = 1.0;
   } else {
       res = res * lengthToLookFor / strlen(text) * 0.5;
   }
   PyObject *return_val = Py_BuildValue("d", res);
   PyMem_RawFree(subseq);
   return return_val;
}

static PyObject* getQUAT(PyObject* self,  PyObject *args)
{
   int bit1, bit2;
   char *res = PyMem_RawMalloc(2 * sizeof(char));
   if (res == NULL)
       return PyErr_NoMemory();
   res[0] = 'A';
   res[1] = '\0';
   if (!PyArg_ParseTuple(args, "pp", &bit1, &bit2)) {
      return NULL;
   }
   if (bit1 && bit2 ) {
       res[0] = 'T';
   } else if (!bit1 && bit2) {
       res[0] = 'C';
   } else if (bit1 && !bit2) {
       res[0] = 'G';
   }
   PyObject *return_val = Py_BuildValue("s", res);
   PyMem_RawFree(res);
   return return_val;
}

static PyObject* byte2QUATS(PyObject* self,  PyObject *args)
{
   int byte;
   int bit1, bit2;
   char *res = PyMem_RawMalloc(5 * sizeof(char));
   if (res == NULL)
       return PyErr_NoMemory();
   res[0] = 'A';
   res[1] = 'A';
   res[2] = 'A';
   res[3] = 'A';
   res[4] = '\0';
   if (!PyArg_ParseTuple(args, "i", &byte)) {
      return NULL;
   }
   int pos = 6;
   for (int i=0; i < 4; i++) {
        bit1 = ((byte >> (pos+1))  & 0x01);
        bit2 = ((byte >> pos)  & 0x01);
        if (bit1 && bit2 ) {
           res[i] = 'T';
        } else if (!bit1 && bit2) {
           res[i] = 'C';
        } else if (bit1 && !bit2) {
           res[i] = 'G';
        }
        pos -= 2;
   }
   PyObject *return_val = Py_BuildValue("s", res);
   PyMem_RawFree(res);
   return return_val;
}

static PyObject* strContainsSub(PyObject* self,  PyObject *args)
{
   char *text;
   char *substr;
   int res = 0;
   if (!PyArg_ParseTuple(args, "ss", &text, &substr)) {
      return NULL;
   }
   if (strstr(text, substr) != NULL) {
      res = 1;
   }
   PyObject *return_val = Py_BuildValue("O", res ? Py_True: Py_False);
   return return_val;
}

static char cdnarules_sat_docs[] =
   "microsatellite(text, lengthToLookFor): Finds the maximum microsatellite with the given lenght!\n";
static char cdnarules_lseq_docs[] =
   "longestSequenceOfChar(text, character_to_look_for): Finds the maximum sequence of a char (or any char if given '*')!\n";
static char cdnarules_rreg_docs[] =
   "repeatRegion(text, lengthToLookFor): returns true if a region of 'lengthToLookFor' is repeated within text (includes overlapping texts)";
static char cdnarules_srreg_docs[] =
   "smallRepeatRegion(text, lengthToLookFor): returns a error-value based on the number repeats the text contains";
static char cdnarules_getq_docs[] =
   "getQUAT(bit1,bit2): returns the DNA base for the given bits";
static char cdnarules_byte2quats_docs[] =
   "byte2QUATS(byte): Converts a given byte to DNA representation";
static char cdnarules_strcsubstr[] =
   "strContainsSub(text,substr): Returns true if the substr is present in text";
static char bitsSet_docs[] =
    "bitsSet(integer): returns the number of bits set int given integer";
static char grayCode_docs[] =
    "grayCode(integer): create a new int for graycode construction";
static char buildGraySequence_docs[] =
"buildGraySequence(length, b): build up a grey sequence of length <length> with bit b set";
static char bitSet_docs[] =
"bitSet(X,b): returns if bit b is set in X";
static char xorarray_docs[] =
"xor_array(X,Y): returns the xor of the two input arrays";
static char elimination_docs[] =
"elimination(A,b,packet_mapping): performs gaussian elimination on A and b. returns true (for now)";
static PyMethodDef cdnarules_funcs[] = {
   //{"microsatellite", (PyCFunction)microsatellite, METH_NOARGS, cdnarules_docs},
   {"bitsSet", bitsSet, METH_VARARGS, bitsSet_docs},
   {"microsatellite", microsatellite, METH_VARARGS, cdnarules_sat_docs},
   {"longestSequenceOfChar", longestSequenceOfChar, METH_VARARGS, cdnarules_lseq_docs},
   {"repeatRegion", repeatRegion, METH_VARARGS, cdnarules_rreg_docs},
   {"smallRepeatRegion", smallRepeatRegion, METH_VARARGS, cdnarules_srreg_docs},
   {"getQUAT", getQUAT, METH_VARARGS, cdnarules_getq_docs},
   {"byte2QUATS", byte2QUATS, METH_VARARGS, cdnarules_byte2quats_docs},
   {"strContainsSub", strContainsSub, METH_VARARGS, cdnarules_strcsubstr},
   {"bitsSet", bitsSet, METH_VARARGS, bitsSet_docs},
   {"grayCode", grayCode, METH_VARARGS, grayCode_docs},
   {"buildGraySequence", buildGraySequence, METH_VARARGS, buildGraySequence_docs},
   {"bitSet", bitSet, METH_VARARGS, bitSet_docs},
   {"xorArray", xor_array, METH_VARARGS, xorarray_docs},
   {"elimination", elimination, METH_VARARGS, elimination_docs},
   {NULL, NULL, 0, NULL}
};

static struct PyModuleDef cdnarules =
{
    PyModuleDef_HEAD_INIT,
    "cdnarules", /* name of module */
    "Extension module for fast DNARules processing!",
    //"usage: Combinations.uniqueCombinations(lstSortableItems, comboSize)\n", /* module documentation, may be NULL */
    -1,   /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    cdnarules_funcs
};

PyMODINIT_FUNC PyInit_cdnarules(void)
{
    //PyObject* module =  PyModule_Create(&cdnarules);
    import_array();
    //return module;
    return PyModule_Create(&cdnarules);
}

/*void initcdnarules(void)
{
   Py_Initialize3("cdnarules", cdnarules_funcs, "Extension module for fast DNARules processing!");

}*/