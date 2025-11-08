# Python Code to read images
"""
Author: Bhargav K
Compress images using SVD, QR ITERATION, EIGENVALUE DECOMPOSITION.
"""


import matplotlib.pyplot as plt
import numpy as np
import ctypes
import time

# froberr = []
# kvals = []

lib = ctypes.CDLL("../c_backend/libsvd.so")

lib.truncatedsvd.argtypes = [
    ctypes.POINTER(ctypes.POINTER(ctypes.c_double)),  
    ctypes.c_int,  
    ctypes.c_int,  
    ctypes.c_int, 
    ctypes.POINTER(ctypes.POINTER(ctypes.c_double))   
]
lib.truncatedsvd.restype = None

def grayscale(image, k):

    A = np.array(image, dtype = np.float64)
    m, n = image.shape
    A_out = np.zeros_like(A)

    # Pass A to convert to double **
    A_ptr = numpy_to_ptr2d(A)
    Aout_ptr = numpy_to_ptr2d(A_out)
    
    # measure time
    start = time.perf_counter() #starting time
    lib.truncatedsvd(A_ptr, m, n, k, Aout_ptr) # Call C function
    end = time.perf_counter() #time end
    print(f"k = {k}, Time elapsed = {end-start: .5f} seconds")
    finalimage = A_out


    error = np.linalg.norm(A - finalimage, ord='fro')
    normA = np.linalg.norm(A, ord='fro')
    print("Frobenius error:", error)
    print("Normalized Frobenius error:", error / normA)  
    # froberr.append(error / normA)
    return finalimage  

def rgb(image, k):

    A = np.array(image, dtype = np.float64)
    m = A.shape[0]
    n = A.shape[1]
    c = A.shape[2]
    finalimage = np.zeros_like(A, dtype=np.float64)

    start = 0
    start = time.perf_counter() #start time
    for x in range(c):
        A_ch = np.ascontiguousarray(A[:, :, x], dtype=np.float64)
        A_out = np.zeros_like(A_ch)

        A_ptr = numpy_to_ptr2d(A_ch)
        Aout_ptr = numpy_to_ptr2d(A_out)

        
        lib.truncatedsvd(A_ptr, m, n, k, Aout_ptr)
        finalimage[:, :, x] = A_out
    end = time.perf_counter() #edning time

    print(f"For k = {k}, Time elapsed = {end-start: .5f} seconds")
    #frob error
    error = np.linalg.norm((A - finalimage).ravel(), 2)
    normA = np.linalg.norm(A.ravel(), 2)

    print("Frobenius error:", error)
    print("Normalized Frobenius error:", error /   normA)
    # froberr.append(error/ normA)

    return finalimage


def numpy_to_ptr2d(arr):
    m = arr.shape[0]  
    arr_ptr = (ctypes.POINTER(ctypes.c_double) * m)()


    for i in range(m):
        arr_ptr[i] = arr[i].ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    return arr_ptr

num = int(input("Enter the number of images: "))



for i in range(num):
    s = input("Enter image name(image.png/image.jpg/image.jpeg) without spaces: ")
    a = "../../../figs/" + s

    # err handling
    try:
        image = plt.imread(a)

    except FileNotFoundError:
        print("File not found!!")

        continue


    k = int(input("Enter number of singular values to keep: "))
    #kvals.append(k)
    if image.max() <= 1.0:
        image = (image * 255).astype(np.float64)


    if image.ndim == 2:  # Grayscale
        finalimage = grayscale(image, k)

    elif image.ndim == 3:  # RGB
        finalimage = rgb(image, k)

    else:
        print("Invalid input!!")
        continue

    finalimage = np.clip(finalimage, 0, 255).astype(np.uint8)
    
    output = f"../../../figs/k{k}_final_{s}"

    plt.imsave(output, finalimage)

    if finalimage.ndim == 2:
        plt.imshow(finalimage, cmap = "gray")
    elif finalimage.ndim == 3:
        plt.imshow(finalimage)

    plt.axis("off")
    plt.title(f"finalimage_k_{k}_{s}")
    plt.show()

# plt.plot(kvals, froberr, marker='o')
# plt.xlabel("Number of singular values (k)")
# plt.ylabel("Normalized Frobenius Error")

# plt.title("Image Quality vs Number of Singular Values")
# plt.grid(True)
# plt.savefig("../figs/image_quality_graph.png")

