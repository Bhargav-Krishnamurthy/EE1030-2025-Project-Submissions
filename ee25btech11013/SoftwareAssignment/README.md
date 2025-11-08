\section*{Project Overview}

This project implements image compression using Truncated Singular Value Decomposition (SVD).
The core numerical engine is written in~C, performing eigenvalue decomposition using QR iteration
on the symmetric matrix \( \vec{A}^{\top}\vec{A} \). A Python front end (via \texttt{ctypes}) handles
image loading, communication with the C backend, and visualization.

\subsection*{Key Features}
\begin{itemize}
    \item Hybrid C + Python design
    \item SVD computed using QR iteration on \( \vec{A}^{\top}\vec{A} \)
    \item Works for grayscale and RGB images
    \item Adjustable rank \(k\)
    \item Quantitative error evaluation using Frobenius and Normalised Frobenius norms
\end{itemize}

\subsection*{Method Summary}

\begin{align}
    &\text{1) Convert input image to matrix form } \vec{A}.\\[6pt]
    &\text{2) Form the symmetric matrix }
        \vec{M} = \vec{A}^{\top}\vec{A}. \\[6pt]
    &\text{3) Apply QR iteration to obtain }
        \vec{M} \approx \vec{V}\Lambda\vec{V}^{\top}. \\[6pt]
    &\text{4) Singular values: }
        \sigma_i = \sqrt{\lambda_i}. \\[6pt]
    &\text{5) Left singular vectors: }
        \vec{u}_i = \frac{1}{\sigma_i}\,\vec{A}\vec{v}_i. \\[6pt]
    &\text{6) Truncated reconstruction: }
        \vec{A}_k = \vec{U}_k\,\Sigma_k\,\vec{V}_k^{\top}.
\end{align}

An example of matrix notation used in this work is:
\begin{align}
    \vec{A} = \myvec{
        1 & 0 \\
        0 & 1
    }.
\end{align}

\subsection*{Architecture Overview}

\textbf{C Backend}
\begin{itemize}
    \item Performs QR-based eigenvalue decomposition of \( \vec{A}^{\top}\vec{A} \)
    \item Produces eigenvalues (singular values squared) and right singular vectors
    \item Compiled into a shared library \texttt{libsvd.so}
\end{itemize}

\textbf{Python Frontend}
\begin{itemize}
    \item Loads the image, converts it to \(\vec{A}\)
    \item Calls the C backend using \texttt{ctypes}
    \item Computes \( \vec{U} \) and reconstructs \( \vec{A}_k \)
    \item Displays and saves results
\end{itemize}

\subsection*{Performance Notes}
Increasing \(k\) reduces compression but improves quality.
Error metrics such as Frobenius error and Normalised Frobenius error are reported.

\subsection*{Goals}
\begin{itemize}
    \item Implement SVD from first principles
    \item Study numerical behaviour of QR iteration in truncated SVD
    \item Demonstrate C--Python interoperability
    \item Analyse compression results quantitatively
\end{itemize}

\subsection*{Motivation}
Standard SVD routines in libraries (NumPy/LAPACK) act as black boxes.
This implementation allows direct study of the underlying numerical method.

\subsection*{Deliverables}
\begin{itemize}
    \item C backend \texttt{libsvd.so} implementing truncated SVD
    \item Python driver for image compression
    \item Results, plots, and analysis
    \item Full technical report
\end{itemize}

\subsection*{Future Work}
\begin{itemize}
    \item Bidiagonalisation (Golub--Kahan)
    \item Lanczos iteration for faster truncated SVD
    \item GUI for interactive image compression
\end{itemize}

