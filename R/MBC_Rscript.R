# MBC_Rscript.R
# ------------------------------------------------------------------------------
# This script is called by a python and does bias correction of gcm climate 
# output using the MBC package. Specifically, it reads in numpy arrays that  
# represent climate data for a single grid cell in a gcm, does the bias 
# correction, and then exports it in the same numpy format.
# 
# Run on macOS monterey 12.4
# R-version: 4.2.0
# ------------------------------------------------------------------------------

# Import required libraries
require(MBC,quietly = TRUE) # Multivariate Bias Correction of Climate Model
                            # Outputs (v0.10-5)
require(reticulate,quietly = TRUE) # Interface to 'Python' (v1.25)
require(abind,quietly = TRUE) # Combine Multidimensional Arrays (v1.4-5)

# Initialize workspace
rm(list=ls())
# gc()
defaultW <- getOption("warn")
options(warn = -1) # Turn warnings off to minimize output to screen

# Load conda evnironemnt and call to numpy for load/export
reticulate::use_condaenv("/Users/adam/miniconda3/envs/pygeo38")
numpy <- reticulate::import("numpy")

# Import commands from Python subprocess call. Arrays imported have numpy .npy
# extension.
args <- commandArgs(trailingOnly = TRUE)

ref_arr <- numpy$load(args[1]) # Reference array (ERA5 regridded)
hst_arr <- numpy$load(args[2]) # Historical "raw" GCM output overlapping with 
                               # ERA5
sim_arr <- numpy$load(args[3]) # Future/projected "raw" GCM output to be 
                               # bias corrected
bc_filename <- args[4] # File name to use for exported bias corrected datasets
time <- numpy$load(args[5]) # Time id by which to group and bias correct 
                            # (generally months).

# Default bias correction method is "QDM", otherwise "MBCn" can be chosen
if (length(args) == 5){
  
  method <- "QDM"
  
} else {
 
  method <- args[6] 
  
}

# Parameters for QDM, see ?'MBC-package' for details
ratio_seq <- c(FALSE,TRUE,TRUE,TRUE)
trace <- c(Inf,0.0,0.0,0.05)

unique_time <- unique(time) # how many unique time units/designations are there?

# Empty matrices to fill with bias-corrected output.
hst_bc <- matrix(NA,nrow = nrow(ref_arr),ncol = ncol(ref_arr))
sim_bc <- matrix(NA,nrow = nrow(ref_arr),ncol = ncol(ref_arr))

for (t in unique_time){ # For each unique time component ...
  
  id <- time == t # Find rows for the given time unit (e.g. month 1 "January")

  # Subset imported numpy arrays for just these rows
  o.c <- ref_arr[id,] # Reference ERA5 data
  m.c <- hst_arr[id,] # GCM historical data that overlaps with ERA5
  m.p <- sim_arr[id,] # Projected data to bias correct
  
  # If any of these three variables has all observations missing then skip.
  if (any((sum(is.nan(o.c)) == length(o.c)) | 
            (sum(is.nan(m.c)) == length(m.c)) | 
            (sum(is.nan(m.p)) == length(m.p)))){
    
    next
    
  }
  
  if (method == "QDM"){ # Quantile Delta Mapping
  
    for (k in 1:ncol(o.c)){ # For each variable (i.e. columns) ...
      
      # ... do the bias correction
      bc = MBC::QDM(o.c = o.c[,k], 
                    m.c = m.c[,k],
                    m.p = m.p[,k],
                    ratio = ratio_seq[k],
                    trace = trace[k])
      
      # Record bias-corrected data in empty arrays
      hst_bc[id,k] <- bc$mhat.c
      sim_bc[id,k] <- bc$mhat.p
    
    }
    
  } else if (method == "MBCn"){ # Same as above but using MBCn
    
    bc <- MBC::MBCn(o.c = o.c,
                    m.c = m.c,
                    m.p = m.p,
                    ratio.seq = ratio_seq,
                    trace = trace,
                    iter = 30,
                    n.escore = round(nrow(o.c)/2),
                    silent = TRUE,
                    return.all = TRUE)
    
    # Select ouput for iteration where escore was minimized
    min_escore_id <- which(bc$escore.iter==min(bc$escore.iter))[[1]] - 2
    
    # Record output
    hst_bc[id,] <- bc$m.iter[[min_escore_id]]$m.c
    sim_bc[id,] <- bc$m.iter[[min_escore_id]]$m.p
    
  }


}

# Set up array for export
export_arr <- abind::abind(hst_bc,
                           sim_bc,
                           along = 0)

# Save bias-corrected ouput as numpy array
numpy$save(bc_filename,export_arr)

# Turn warnings back on
options(warn = defaultW)

# End of script ----------------------------------------------------------------