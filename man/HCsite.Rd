% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HCsite.R
\name{HCsite}
\alias{HCsite}
\title{Extract high confident transcript start or end site}
\usage{
HCsite(
  site_stats,
  cut_off = 0.05,
  cut_at = c("fdr", "pval", "prob"),
  site_region = 50,
  min_read = 2,
  type = c("TSS", "TES")
)
}
\arguments{
\item{site_stats}{A DataFrame of TSS or TES statistics generated by \code{SummaryTranSite}.}

\item{cut_off}{A numeric value of BH adjsuted p-value (fdr), p-value or probablity cut-off.}

\item{cut_at}{The \code{cut_off} is applied to "fdr" (adjustd p-value), "pval" (p-value) or "prob" (probability).}

\item{site_region}{An integer of enriched site region. TSSs/TESs locate within enriched
site+/-site_region are all treated as high confident sites.}

\item{min_read}{The minimum reads in a bin to support high confident sites.}

\item{type}{"TSS" for transcript start site and "TES" for transcript end site.}
}
\value{
A DataFrame of enriched TSS or TES sites.
}
\description{
Extract high confident transcript start or end site
}
