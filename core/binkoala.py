# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 07:42:46 2014

@author: naspert
"""

import numpy
import struct
import sys
"""
Bin files I/O
=============

Koala uses two different file formats having the same '.bin' extension. One is
for storing (real) floating-point data (e.g. phase) and the other is for storing
complex-valued floating point data (e.g. reference hologram). The routines in this 
module allow to read and write these formats.

"""
def read_mat_bin(fname):
    """
    Reads a .bin file containing floating-point values (real) saved by Koala

    Parameters
    ----------
    fname : string
        Path to the file
    
    Returns
    -------
    buffer : ndarray
        An array containing the floating-point values read from the file
    
    header : struct
        A structure containing the fields read from the .bin header

    See Also
    --------
    write_mat_bin



    Example
    --------
    >>> (buf, hdr) = read_mat_bin('test/file.bin')
    >>> hdr
    array([(1, 0, 23, 650, 650, 3.5365880535209726e-07, 3.312875662686565e-07, 1)], 
      dtype=[('version', 'u1'), ('endian', 'u1'), ('head_size', '<i4'), ('width', '<i4'), 
      ('height', '<i4'), ('px_size', '<f4'), ('hconv', '<f4'), ('unit_code', 'u1')])
    >>> buf
    array([[-3.129637  , -2.98734117, -3.01134109, ...,  0.60778421,
         1.0016439 ,  1.54053307],
       [-3.02219605, -2.72388673, -3.08404136, ...,  1.25938392,
         1.62311745,  1.75243735],
       [ 0.38411593, -0.56527954, -0.0919028 , ...,  2.616081  ,
         2.04891896,  1.47056615],
       ..., 
       [ 0.91004229,  0.41167739,  0.41661718, ..., -0.24841429,
         0.43015432,  1.45658135],
       [ 1.08976161,  0.99631542,  1.09343278, ..., -0.17306508,
         0.90412992,  1.96907508],
       [ 1.50583076,  1.60309029,  1.74559569, ..., -1.14475477,
        -2.74203372, -2.70048261]], dtype=float32)


    """
    kbin_header_dtype = numpy.dtype([
     ("version", "u1"), 
     ("endian", "u1"),
     ("head_size", "i4"),
     ("width", "i4"),
     ("height", "i4"),
     ("px_size", "f4"), # pixel size in [m]    
     ("hconv", "f4"),  # height conversion factor (-> m)
     ("unit_code", "u1")   # 1=rad 2=m
    ])
    f = open(fname, 'rb')
    kbin_header = numpy.fromfile(f, dtype=kbin_header_dtype, count=1)
    #print kbin_header
    tmp = numpy.fromfile(f, dtype='float32')
    #print 'array = {}'.format(len(tmp))
    f.close
    shape = (kbin_header['height'], kbin_header['width'])
    return (tmp.reshape(shape), kbin_header)

def write_mat_bin(fname, buf, width, height, px_size=1., hconv=1., unit_code=0):
    """
    Writes a .bin file containing floating-point values (real) using the same
    format as Koala. It consists of header, and the floating-point values of
    the input array, written as their little-endian 32 bits binary representation.

    Header format
    -------------
    The header consists of the following fields:
        + Version of the file format (1 byte, current is 1)
        + Endianness of the file (1 byte, 0 for little-endian)
        + Header size (4 bytes integer) in number of bytes (current = 23)
        + Width of the array (4 bytes integer)
        + Height of the array (4 bytes integer)
        + Pixel size (4 bytes float)
        + Height conversion factore (4 bytes float)
        + Unit code (1 byte)
    
    Parameters
    ----------
    fname : string
        Path to the file to be written
    buf : ndarray
        Buffer to be written
    width : int
        width of the input buffer    
    height : int
        height of the input buffer
    px_size: float, optional
        physical size of a pixel (in meters)
    hconv: float, optional
        height conversion factor (multiplying the values contained in the buffer
        by this value should result in having a buffer containing values in meters)
    unit_code: byte, optional
        Optional indication of the physical data stored in the buffer
        (0 -> no unit, 1 -> radians, 2 -> meters)

    See Also
    --------
    read_mat_bin



    Example
    --------
    >>> a = numpy.eye(5)
    >>> write_mat_bin('test/out.bin', a, 5, 5)


    """
    fd = open(fname, "wb")
    hdr = numpy.zeros((1, 23), dtype='B')
    if sys.byteorder == 'little':
        endian = 0
    else:
        endian = 1
    struct.pack_into('=BBiiiffB', hdr, 0, 1, endian, 23, width, height, px_size, hconv, unit_code)
    fd.write(hdr)
    fd.write(buf.astype(dtype=numpy.float32, order='C'))
    fd.close()

def read_mat_cplx_bin(fname):
    """
    Reads a .bin file containing floating-point values (complex) saved by Koala

    Parameters
    ----------
    fname : string
        Path to the file
    
    Returns
    -------
    buffer : ndarray
        An array containing the complex floating-point values read from the file
    
    
    See Also
    --------
    write_mat_cplx_bin



    Example
    -------
    >>> buf = read_mat_cplx_bin('test/file_cplx.bin')
    >>> buf
    array([[  0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j, ...,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j],
       [  0.00000000e+00 +0.00000000e+00j,
          4.97599517e-09 +9.14632536e-10j,
          5.99623329e-09 -1.52811275e-08j, ...,
          1.17636354e-07 -1.01500063e-07j,
          6.33714581e-10 +5.61812996e-09j,
          0.00000000e+00 +0.00000000e+00j],
       ..., 
       [  0.00000000e+00 +0.00000000e+00j,
         -1.26479121e-08 -2.92324431e-09j,
         -4.59448168e-09 +9.28236474e-08j, ...,
         -4.15031316e-08 +1.48466597e-07j,
          4.41099779e-08 -1.27046489e-08j,
          0.00000000e+00 +0.00000000e+00j],
       [ -0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j, ...,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j,
          0.00000000e+00 +0.00000000e+00j]], dtype=complex64)



    """
    kcplx_header_dtype = numpy.dtype([
    ("width", "i4"),
    ("height", "i4")
    ])
    f = open(fname, 'rb')
    kcplx_header = numpy.fromfile(f, dtype=kcplx_header_dtype, count=1)
    shape = (kcplx_header['height'], kcplx_header['width']) 
    #print kcplx_header    
    tmp = numpy.fromfile(f, dtype='float32')
    f.close
    real_tmp = (tmp[0:kcplx_header['height']*kcplx_header['width']]).reshape(shape)
    imag_tmp = (tmp[kcplx_header['height']*kcplx_header['width']:]).reshape(shape)
    #print tmp
    #print 'array = {}'.format(len(tmp)) 
    return real_tmp + 1j*imag_tmp
    
def write_mat_cplx_bin(fname, buf, width, height):
    """
    Writes a .bin file containing floating-point values (complex) using the same
    format as Koala. It consists of a basic header (just the width and height of
    the array), followed by the real part of the buffer, and finally the imaginary
    part of the input array, written as their little-endian 32 bits binary representation.

    Header format
    -------------
    The header consists of the following fields:
        + Width of the array (4 bytes integer)
        + Height of the array (4 bytes integer)
    
    Parameters
    ----------
    fname : string
        Path to the file to be written
    buf : ndarray
        Buffer to be written
    width : int
        width of the input buffer    
    height : int
        height of the input buffer

    See Also
    --------
    read_mat_cplx_bin



    Example
    --------
    >>> a = numpy.eye(5) + 1j*numpy.ones((5,5))
    >>> write_mat_bin('test/out.bin', a, 5, 5)


    """
    fd = open(fname, "wb")
    hdr = numpy.zeros((1, 8), dtype='B')
    
    struct.pack_into('=ii', hdr, 0, width, height)
    fd.write(hdr)
    r = numpy.real(buf).astype(dtype=numpy.float32, order='C')
    fd.write(r)
    r = numpy.imag(buf).astype(dtype=numpy.float32, order='C')
    fd.write(r)
    fd.close()
