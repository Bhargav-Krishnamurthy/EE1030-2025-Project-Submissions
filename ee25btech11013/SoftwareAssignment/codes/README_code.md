\section*{README}

\subsection*{Image Compression using Truncated SVD \\ QR-Based Eigenvalue Decomposition (Hybrid C + Python)}

This project implements image compression using Truncated Singular Value Decomposition (SVD).
The core numerical engine is written in C, performing eigenvalue decomposition using QR iteration
on \(A^{\top}A\). The C backend is loaded into Python using \texttt{ctypes}, while Python handles image
I/O and visualization.

\subsection*{Project Structure}

\begin{verbatim}
project_root/
│
├── report.tex / report.pdf         # LaTeX project report
├── tables/                         # Tables used in report
├── figs/                           # Plots / images
│
├── codes/
│   ├── c_libs/                     # Additional C utility sources (optional)
│   ├── c_main/                     # Standalone C test programs
│   ├── python_driver/              # Optional python drivers
│   │
│   └── hybrid_c_python/
│       ├── c_backend/              # C implementation + shared library
│       │   ├── main.c
│       │   ├── matrix.c
│       │   ├── matrix.h
│       │   ├── matrix.o
│       │   ├── svd.c
│       │   ├── svd.h
│       │   ├── svd.o
│       │   └── libsvd.so           # Compiled shared object
│       │
│       └── python_frontend/
│           └── main.py             # Python driver
│
├── gvv.sty
└── gvv-book.sty
\end{verbatim}

\subsection*{Building the Shared C Library (\texttt{libsvd.so})}

Move into the C backend directory:
\begin{verbatim}
cd codes/hybrid_c_python/c_backend/
\end{verbatim}

Compile object files:
\begin{verbatim}
gcc -O2 -fPIC -c main.c matrix.c svd.c
\end{verbatim}

Link into shared library:
\begin{verbatim}
gcc -shared -o libsvd.so main.o matrix.o svd.o -lm
\end{verbatim}

\noindent
Output:
\begin{verbatim}
libsvd.so
\end{verbatim}

\noindent
Flags:
\begin{itemize}
\item \verb|-O2| -- optimization
\item \verb|-fPIC| -- position independent code
\item \verb|-lm| -- link math library
\end{itemize}

\subsection*{Running the Python Frontend}

Navigate to:
\begin{verbatim}
cd codes/hybrid_c_python/python_frontend/
\end{verbatim}

Run:
\begin{verbatim}
python3 main.py
\end{verbatim}



\subsection*{End-to-End Workflow}

\begin{enumerate}
\item Build shared object:
\begin{verbatim}
cd codes/hybrid_c_python/c_backend/
gcc -O2 -fPIC -c main.c matrix.c svd.c
gcc -shared -o libsvd.so main.o matrix.o svd.o -lm
\end{verbatim}

\item Run Python:
\begin{verbatim}
cd ../python_frontend/
python3 main.py
\end{verbatim}
\end{enumerate}

The script:
\begin{itemize}
\item Loads the image into Python
\item Sends matrix to C backend
\item Performs truncated SVD via QR iteration
\item Reconstructs compressed image
\item Saves/displays result
\end{itemize}


\subsection*{Dependencies}

\textbf{C:}
\begin{itemize}
\item gcc
\item math.h
\end{itemize}

\textbf{Python:}
\begin{itemize}
\item numpy
\item matplotlib
\item ctypes (built-in)
\end{itemize}

Install via:
\begin{verbatim}
pip install numpy matplotlib
\end{verbatim}








