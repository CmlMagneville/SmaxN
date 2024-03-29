% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_Pseudo_SmaxN_and_checks.R
\name{next.possible}
\alias{next.possible}
\title{Tell if the next value is possible or not given time between cameras}
\usage{
next.possible(max_bloc, next_possible, time_df)
}
\arguments{
\item{max_bloc}{a dataframe containing in rows the number already chosen on
the path and their coordinates ie. on which column (camera name) and row 
(time) they are. Values are ranged fom the biggest (first rows) to the 
lowest (last rows). Column names: value, cam_nm, timestep}

\item{next_possible}{a dataframe with column names: value, cam_nm, timestep 
with a row for the information about the next possible value. If several cells 
have the same value (ex-aequo), the dataframe has one row for each cell. This
dataframe is obtained through the \code{search.next.value} function.}

\item{time_df}{a numerical dataframe containing the minimal time 
needed for an individual of the studied species to go from a camera to 
another camera.There are as many rows as there are cameras and there are
as many columns as there are cameras, thus the dataframe is symmetrical
and the diagonal is filled with 0. This dataframe is the output of the 
\code{compute.cam.time} function.}
}
\value{
the next value to be chosen that is to say the biggest value of the
bloc after the last value in max_bloc which is possible according to the 
distance between cameras. Note: For one possible value, check if possible
with already on the path values and stop when one is not possible. This
function tests for all possible values (all ex-aequo) so we have a df 
showing TRUE or FALSE for all ex-aequo given the values already on the path.
}
\description{
This function a TRUE/FALSE value if the next value is possible. If there
are ex-aequo values, begin with the first one in the \code{next_possible} 
dataframe, use the other until it works.
}
