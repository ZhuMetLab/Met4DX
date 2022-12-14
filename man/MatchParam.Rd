% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions-Params.R
\name{MatchParam}
\alias{MatchParam}
\title{Matching parameter setup}
\usage{
MatchParam(
  ppm = 35,
  cutoff = 0.8,
  methodMatch = c("direct", "bootstrap"),
  methodScore = c("dp", "msdial", "nist", "bootstrap", "combined", "ratio", "gnps",
    "bonanza", "hybrid"),
  weightMZ = 0,
  weightIntensity = 1,
  includePrecursor = TRUE,
  ppmPrecursorFilter = 20,
  ms2range = NULL,
  thrIntensityAbs = 0,
  thrIntensityRel = 0.01,
  intensityExpNormed = TRUE,
  intensityLibNormed = TRUE,
  tuneLibSpectra = FALSE,
  useMS1ResDefine = TRUE,
  resDefineAt = 200,
  normIntensity = TRUE,
  intensityNormedMethod = c("maximum", "bonanza", "gnps", "hybrid")
)
}
\arguments{
\item{ppm}{\code{numeric} ppm tolerance for MS2 m/z match}

\item{cutoff}{\code{numeric} threshold for a avaliable match}

\item{methodMatch}{\code{character} method for matching with librarial spectra
(eithor "direct" or "bootstrapping")}

\item{methodScore}{\code{character} method for scoring the MSMS match ("dp",
"msdial", "nist", "bs", "combined" are avaliable)}

\item{weightMZ}{\code{numeric} weight of m/z when scoring}

\item{weightIntensity}{\code{numeric} weight of intensity when scoring}

\item{includePrecursor}{\code{logical} if consider the precursor fragments when
matching}

\item{ms2range}{\code{numeric} mass range setup when acquiring MSMS data}

\item{thrIntensityAbs}{\code{numeric} absolute intensity threshold to be removed}

\item{thrIntensityRel}{\code{numeric} relative intensity threshold to be removed}

\item{intensityExpNormed}{\code{logical} if the spectral intensity is normalized
in experiment spectra}

\item{intensityLibNormed}{\code{logical} if the spectral intensity is normalized
in library spectra}

\item{tuneLibSpectra}{\code{logical} if apply thrIntensityAbs or thrIntensityRel
to reference}

\item{useMS1ResDefine}{\code{logical} if use resDefineAt when calculating precursor m/z tolerance}

\item{resDefineAt}{\code{numeric} m/z threshold for using ppm tolerance for MS1 or
MS2 match, for smaller m/z values, using the tolerance of ppm.ms1 * res.defineat
to mathing experimental and librarial fragments.}
}
\value{
an \code{ParamMatch} object
}
\description{
Matching parameter setup
}
