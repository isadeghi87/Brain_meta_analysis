################################################################
#
# This files defines the default configuration for the course
# Statistics Applied to Bioinformatics It should be loaded before
# starting the practicals. All the rest will be done via internet.  #
#
# If you install the course on your computer, you only need to change
# the variable dir.main

#### Load files from the web server
## dir.course <- 'http://jacques.van-helden.perso.luminy.univmed.fr/statistics_bioinformatics'
#dir.course <- 'http://pedagogix-tagc.univ-mrs.fr/courses/statistics_bioinformatics'
dir.course <- 'http://pedagogix-tagc.univ-mrs.fr/courses/statistics_bioinformatics'

## This is the local configuration on my computer
## You can uncomment and adapt this row if you work with a local copy of the scripts
## dir.course <- 'http://localhost/courses/statistics_bioinformatics'
#dir.course <- '~/statistics_bioinformatics'

## Directory from which the demo data sets will be downloaded
dir.data <- file.path(dir.course, 'data')
print(paste("Data repository", dir.data))

## Directory with the R scripts for the course
dir.R.files <- file.path(dir.course, "R-files")
print(paste("R scripts source", dir.R.files))
dir.util <- file.path(dir.R.files, 'util')

################################################################
# Load utilities
source(file.path(dir.util, 'util.R'))

################################################################
# Output directories
dir.home = Sys.getenv('HOME');
if (dir.home=="") {
  stop("The environment variable $HOME is not defined in your operating system. You need to define it for a proper configuration of the scripts.")
}

dir.main <- file.path(dir.home, 'course_stats_bioinfo')

## Define directory for storing your results
dir.results <- file.path(dir.main, "results")
print(paste("Results will be saved to", dir.results))

## Define directory for saving your figures
dir.figures <- file.path(dir.main, "figures")
print(paste("Figures will be saved to", dir.figures))

for (dir in c(dir.main, dir.figures, dir.results)) {
  if (!file.exists(dir)) {
    dir.create(dir, recursive=T,showWarnings=T)
  }
}

################################################################
# Global parameters
verbosity <- 1
export.formats.plots <- c("eps","pdf", "png")
export.formats.obj <- c("table")

################################################################
## Graphic device Special fix for a bug in Mac OSX 10.5: X11() has
## different fonts from the export devices (pdf, png, eps), which
## creates problem with the legends (legends are cut bu the box
## limit). The problem is fixed by opening a Mac-specific graphical
## window using quartz() instead of X11().
sysinfo <- Sys.info()
if (sysinfo["sysname"] == "Darwin") {
  X11 <- function (...) {
    quartz(...)
  }
}

################################################################
## Specify if drawings should be done in colors or not.  For the book
## I cannot use colors, for the presentations, I prefer to use them.
in.colors <- F

