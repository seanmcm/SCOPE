# Libraries
library(matlabr)
library(readr)
library(stringr)
library(TMB) # needed for fishmethods
library(fishmethods)
library(lubridate)
library(dplyr)
library(ggplot2)
setwd("/Users/Sean/Dropbox/GIT_REPOSITORIES/bcipanama/")
run_matlab_script("../SCOPE_test/SCOPE.m", # folder location of SCOPE
                  verbose = TRUE, 
                  desktop = FALSE,
                  splash = FALSE, 
                  display = FALSE, 
                  wait = TRUE,
                  single_thread = FALSE)

# Step 1: Load all functions, libraries and run verification of SCOPE
source(file = "7_scope/scope_functions.R")

# Step 2: Turn on/off verification mode and choose simulation protocol
# verify = 0:on; 1:off
#simulate = 0:individual runs; 1:time series; 2:look-up table
set.options(verify = 0, simulate = 1)

# Step 3: Change input data inside set_parameter_filenames.csv
set.param_file("input_data_default.csv")

# Step 4: (Only for Time Series) Create input directory
dir.name <- c("test_ts_plot")
input.filename <- "test_ts_plot_input"
dir.create(sprintf("../SCOPE/input/dataset %s", dir.name))

# Step 5: Modify filenames.csv
# Note: Here you can specify the folder name of the output, specify the directory of the input for time series, specify the filename of input hypercube (here as input.filename), and remove the parameters (e.g. tts) that are already specified as column in the input.filename

set.filenames(simulation_name = dir.name, 
          dataset_dir = dir.name,
          meteo_ec_csv = input.filename,
          tts = "" # when empty, this automatically computes solar zenith angle (tts) from the tz, lat and long,
          )

# Step 6: Create input file for time series

    # Step 6.1: Create duration of time series and resolution
    start.time <- as.character("2012-07-03 00:00")
    end.time <- as.character("2012-07-4 00:00")
    span <- "30 min"
    
    # Step 6.2: Create constant parameters (this can be linked to some ground truth data)
    par.constant <- c(80, 20, 0, 0.012, 0.01, 0, 1.4, 3, -0.35, -0.15, 70, 12, 70, 0, 90)
    names(par.constant) <- c("Cab", "Cca", "Cant", "Cdm", "Cw", "Cs", "N", "LAI", "LIDFa", "LIDFb", "Vcmax25", "BallBerrySlope", "tts", "tto", "psi")
    
    # Step 6.3: Create variable parameters (still need to work on the bci meteo data)
    load(file="../bcipanama/3_output_clean_data/bci_ec_flux.Rda")
    par.variable.names <- c("Rs", "tair", "RH")
    par.variable <-bci_ec_flux[(bci_ec_flux$date >= start.time & bci_ec_flux$date <= end.time & !is.na(bci_ec_flux$date)), par.variable.names]
    
    colnames(par.variable) <- c("Rin", "Ta", "RH")
    
    # Step 6.4: Create a input variables for time series (Note: solar zenith and azimuth angle automatically computes based on time, lat and long)
    make.ts.input(start.time, end.time, span = span, 
              par.constant, par.variable,
              dir.name = paste0("../SCOPE/input/dataset ", dir.name), 
              file.name = input.filename)
    
# Step 7: Modify input values (link the leaf param from PCA result)

input_data_default(Cdm = NULL, 
                   Cab = NULL, 
                   Cca = NULL, 
                   Vcmax25 = NULL, 
                   Cw = NULL, 
                   Cant = NULL, 
                   LAT = 9.1, 
                   LON = -79.8, 
                   timezn = -5)    

# Run SCOPE.m
run_matlab_script("../SCOPE/SCOPE.m", verbose = TRUE, desktop = FALSE,
                  splash = FALSE, display = FALSE, wait = TRUE,
                  single_thread = FALSE)

# Call the last result
run <- 4 # number of runs in the same filename
type <- "fluxes"

out <- call.output(dir.name, run, type)

output <- read.csv(out)
output <- output[-1,]
output <- output %>% mutate_if(is.character, as.numeric)

# Plot
plot(output$DoY, output$Actot, type="l")

bci_ec_flux[(bci_ec_flux$date >= start.time & bci_ec_flux$date <= end.time & !is.na(bci_ec_flux$date)),] %>%
  #filter(FLAG == 1) %>%
  #filter(ustar >= 0.2) %>%
  #filter(between(date, as.POSIXct(start.time), as.POSIXct(end.time))) %>%
  ggplot(aes(x=doy.frac, y=gpp)) + 
  geom_line(aes(color="EC")) + 
  geom_line(data=output, aes(x=DoY, y=Actot, color="SCOPE")) +
  theme_classic() + 
  ylab(bquote('GPP,' ~mu~'mol' ~m^-2 ~s-1)) +
  xlab("DOY, fraction") 
