% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes_triangular_mat.R
\name{upper_triangular_cor_mat}
\alias{upper_triangular_cor_mat}
\title{Class for symmetric correlation matrices}
\description{
The class allows to store the upper triangular matrix of a symmetric
correlation matrix in an memory-efficient form and return a data.table or
dense R matrix if need be.
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-upper_triangular_cor_mat-new}{\code{upper_triangular_cor_mat$new()}}
\item \href{#method-upper_triangular_cor_mat-print}{\code{upper_triangular_cor_mat$print()}}
\item \href{#method-upper_triangular_cor_mat-get_data_table}{\code{upper_triangular_cor_mat$get_data_table()}}
\item \href{#method-upper_triangular_cor_mat-get_cor_matrix}{\code{upper_triangular_cor_mat$get_cor_matrix()}}
\item \href{#method-upper_triangular_cor_mat-get_cor_vector}{\code{upper_triangular_cor_mat$get_cor_vector()}}
\item \href{#method-upper_triangular_cor_mat-clone}{\code{upper_triangular_cor_mat$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-upper_triangular_cor_mat-new"></a>}}
\if{latex}{\out{\hypertarget{method-upper_triangular_cor_mat-new}{}}}
\subsection{Method \code{new()}}{
Initialises the R6 class.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{upper_triangular_cor_mat$new(cor_coef, features, shift)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{cor_coef}}{Numerical vector. The correlation coefficients of the
upper triangular correlation matrix stored as a row-major vector}

\item{\code{features}}{String vector. The features of the correlation matrix.}

\item{\code{shift}}{Integer. Was a shift applied during the calculation of the
upper triangular matrix. Typically 0 (diagonal included) or 1 (diagonal
not included)}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Returns the initialised class.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-upper_triangular_cor_mat-print"></a>}}
\if{latex}{\out{\hypertarget{method-upper_triangular_cor_mat-print}{}}}
\subsection{Method \code{print()}}{
Print the class
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{upper_triangular_cor_mat$print()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
Returns the initialised class
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-upper_triangular_cor_mat-get_data_table"></a>}}
\if{latex}{\out{\hypertarget{method-upper_triangular_cor_mat-get_data_table}{}}}
\subsection{Method \code{get_data_table()}}{
Returns the data in form of a data.table.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{upper_triangular_cor_mat$get_data_table(factor = FALSE, .verbose = TRUE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{factor}}{Boolean. Shall the string columns be transformed into
factors. Reduces size of the object; however, takes longer to generate.}

\item{\code{.verbose}}{Boolean. Controls verbosity.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A data.table with three columns:
\itemize{
\item feature_a: The name of the first feature in the correlation matrix.
\item feature_b: The name of the second feature in the correlation
matrix.
\item cor: The correlation coefficients between these two features.
}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-upper_triangular_cor_mat-get_cor_matrix"></a>}}
\if{latex}{\out{\hypertarget{method-upper_triangular_cor_mat-get_cor_matrix}{}}}
\subsection{Method \code{get_cor_matrix()}}{
Return the full correlation matrix.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{upper_triangular_cor_mat$get_cor_matrix(.verbose = TRUE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{.verbose}}{Boolean. Controls verbosity.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Returns the correlation matrix as a dense R matrix.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-upper_triangular_cor_mat-get_cor_vector"></a>}}
\if{latex}{\out{\hypertarget{method-upper_triangular_cor_mat-get_cor_vector}{}}}
\subsection{Method \code{get_cor_vector()}}{
Return the correlation data and shift
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{upper_triangular_cor_mat$get_cor_vector()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A list with
\itemize{
\item cor_data - Numeric vector. The correlations.
\item n_features - Integer. Number of initial features.
\item shift - Integer. The applied shift.
}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-upper_triangular_cor_mat-clone"></a>}}
\if{latex}{\out{\hypertarget{method-upper_triangular_cor_mat-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{upper_triangular_cor_mat$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
