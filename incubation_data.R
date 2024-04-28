
#load packages
library(PoolDilutionR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)

theme_set(theme_minimal() + theme(text = element_text(size = 11)))

setwd("C:/Users/morr497/OneDrive - PNNL/Documents/methane/Summer '23/")
raw_data <- read.csv("PNNL-KM-CH4-d13C_ALL-DATA.csv")
#timestamps <- read.csv("picarro_times.csv")

# 'Pretty n' function to round a numeric value and print that # of digits
pn <- function(x, n) {
  formatC(round(unlist(x), n),
          digits = n, format = "f")
}

# 'Clean p value' function to pretty-print p value(s), specifically
pclean <- function(x, digits = 3, printP = TRUE) {
  x <- as.vector(x)
  ltstring <- paste0("< 0.", paste(rep("0", digits - 1), 
                                   collapse = ""), "1")
  valstring <- ifelse(x < 10 ^ -digits, 
                      ltstring, pn(x, digits))
  if(printP) {
    paste("P", ifelse(x < 10 ^ -digits, 
                      valstring, paste("=", valstring)))
  } else {
    valstring
  }
}

# # Get names of data files
# files <- list.files("picarro/", pattern = "*.csv", full.names = TRUE)
# # Helper function
# read_file <- function(f) {
#   message("Reading ", f)
#   read_csv(f, col_types = "ccdcddcddddddddddddddddddddddccccc") %>%
#     mutate(File = basename(f))
# }
# Read in and pre-process data

raw_data %>%
  mutate(Depth = recode(Depth, "45585" = "10-20")) %>%
  # fix obnoxious excel error for 10 to 20 cm depth
  filter(Sample != "bkgd") %>%
  # select only pool dilution data
  select(Sample, Depth, Core,
         X12CH4_ppm, X13CH4_ppm) %>%
  # select only the data needed
  mutate(id = factor(paste(Core, Depth, sep = "_"))) -> incdat_raw

summary(incdat_raw)

VOL_ML <- 130  # total volume of jar

incdat_raw %>%
  # Volume of jar = 130 ml or 0.130 L, 1 ppm = 0.001 ml/L
  # 10 mL from jar injected into 12 mL evacuated exetainer
  # need to check standards to see if this dilution has an effect
  # assuming for now that it doesn't and using calculations from 2022 work
  # therefore ppm to ml = ppm * 0.001 * VOL_ML/1000 * 2
  mutate(cal12CH4ml = X12CH4_ppm * 0.001 * VOL_ML/1000 * 2 * 1000,
         cal13CH4ml = X13CH4_ppm * 0.001 * VOL_ML/1000 * 2 * 1000,
         # for each 10 ml sample from 130 ml jar,
         # 10 ml of zero air injected
         # 12 parts sample to 1 part dilutant (13 parts total)
         # remaining gas in jar is diluted 12:1
         cal12CH4ml = if_else(Sample != "T0", cal12CH4ml * 1.083, cal12CH4ml),
         cal13CH4ml = if_else(Sample != "T0", cal13CH4ml * 1.083, cal13CH4ml),
         # calculate atom percent (AP) of 13C methane in sample over time
         AP_obs = cal13CH4ml / (cal12CH4ml + cal13CH4ml) * 100) -> incdat

###
# NOTE ALL T0 and T4 data are in separate files!!!
# need to bring those in!!!
# for now pretending that T1 = T0
###

#add duration between samples
incdat$time_days <- NA
incdat[incdat$Sample == "T1",]$time_days <- 0
incdat[incdat$Sample == "T2",]$time_days <- 35
incdat[incdat$Sample == "T3",]$time_days <- 70
incdat$time_days <- as.numeric(incdat$time_days)

incdat %>%
  ungroup() %>%
  group_by(id) %>%
  arrange(time_days) -> incdat

ggplot(data = incdat, aes(time_days, cal12CH4ml + cal13CH4ml)) +
  geom_point(aes(group = id, color = as.factor(Core))) +
  geom_line(aes(group= id)) +
  facet_wrap(~Core)


pk_results <- list()
incdat_out <- list()
all_predictions <- list()

for(i in unique(incdat$id)) {
  message("------------------- ", i)
  # Isolate this sample's data
  incdat %>%
    filter(id == i) %>%
    select(id, Sample, time_days,
           cal12CH4ml, cal13CH4ml,
           AP_obs) ->
    dat
  
  # Let optim() try different values for P and k until it finds best fit to data
  # Here we provide starting values for P (0.1), k (k0), and P_frac and k_frac are set to the methane default
  result <- pdr_optimize(time = dat$time_days,
                         m = dat$cal12CH4ml + dat$cal13CH4ml,
                         n = dat$cal13CH4ml,
                         P = 0.1,
                         pool = "CH4",
                         m_prec = 1,
                         ap_prec = 1,
                         include_progress = TRUE)
  
  # Save progress details separately so they don't print below
  progress_detail <- result$progress
  result$progress <- NULL
  
  #message("Optimizer solution:")
  #print(result)
  P <- result$par["P"]
  id <- dat$id[1]
  pk_results[[i]] <- tibble(id = id,
                            P = P,
                            k = result$par["k"],
                            k0 = result$initial_par["k"],
                            convergence = result$convergence,
                            message = result$message)
  
  # Predict based on the optimized parameters
  pred <- pdr_predict(time = dat$time_days,
                      m0 = dat$cal12CH4ml[1] + dat$cal13CH4ml[1],
                      n0 = dat$cal13CH4ml[1],
                      P = P,
                      k = result$par["k"],
                      pool = "CH4")
  dat <- bind_cols(dat, pred)
  
  # Predict based on ALL the models that were tried
  
  x <- split(progress_detail, seq_len(nrow(progress_detail)))
  all_preds <- lapply(x, FUN = function(x) {
    y1<- data.frame(P = x$P[1],
                    k = x$k[1],
                    time = seq(min(dat$time_days), max(dat$time_days), length.out = 5))
    y2 <- pdr_predict(time = y1$time,
                      m0 = dat$cal12CH4ml[1] + dat$cal13CH4ml[1],
                      n0 = dat$cal13CH4ml[1],
                      P = x$P[1],
                      k = x$k[1],
                      pool = "CH4")
    cbind(y1, y2)
  })
  all_predictions[[i]] <- bind_rows(all_preds)
  
  # Calculate implied consumption (ml) based on predictions
  # Equation 4: dm/dt = P - C, so C = P - dm/dt
  total_methane <- dat$cal12CH4ml + dat$cal13CH4ml
  change_methane <- c(0, diff(total_methane))
  change_time <- c(0, diff(dat$time_days))
  dat$Pt <- P * change_time #P is ml/day
  #amount of methane produced at time (t) of this incubation, a volume in mL
  dat$Ct <- dat$Pt - change_methane
  #amount of methane consumed at time (t) of this incubation, a volume in mL
  
  incdat_out[[i]] <- dat
}

pk_results <- do.call(rbind, pk_results)
incdat_out <- do.call(rbind, incdat_out)
all_predictions <- do.call(rbind, all_predictions)

ggplot(incdat_out, aes(time_days, AP_obs, color = id)) +
  geom_point(aes(shape = ""), size = 4) +
  geom_line(
    data = all_predictions,
    aes(time, AP_pred, group = paste(id, P, k)), color = "grey", linetype = 2
  ) +
  geom_line(aes(y = AP_pred, group = id, linetype = ""),
            size = 1.5
  ) +
  scale_linetype_manual(
    name = "Prediction",
    values = "dotted"
  ) +
  scale_shape_manual(
    name = "Observations",
    values = 20
  ) +
  scale_color_discrete(guide = "none") +
  facet_wrap(~id, scales="free_y") +
  xlab("\n Timestep \n") +
  ylab("\n (13C-CH4/Total CH4) x 100 \n") +
  ggtitle("\n Atom% 13C \n") +
  theme(legend.position = "bottom")

## ----Multisample Fit Total Pool, echo = FALSE, fig.height = 5.5, fig.width = 8.25----
ggplot(incdat_out, aes(time_days, cal12CH4ml + cal13CH4ml, color = id)) +
  geom_point(aes(shape = ""), size = 4) +
  geom_line(
    data = all_predictions,
    aes(time, mt, group = paste(id, P, k)), color = "grey", linetype = 2
  ) +
  geom_line(aes(y = mt, group = id, linetype = ""),
            size = 1.5
  ) +
  scale_linetype_manual(
    name = "Prediction",
    values = "dotted"
  ) +
  scale_shape_manual(
    name = "Observations",
    values = 20
  ) +
  scale_color_discrete(guide = "none") +
  facet_wrap(~id, scales="free_y") +
  xlab("\n Timestep \n") +
  ylab("\n Volume (mL) \n") +
  ggtitle("\n Total Methane \n") +
  theme(legend.position = "bottom")


core_ids <- read.csv("treatments.csv")

incdat_out %>%
  ungroup() %>%
  filter(Sample == "T2") %>%
  mutate(Core = as.integer(substr(id, 1, 1)),
         Depth = substr(id, 3, 7)) %>%
  left_join(core_ids) -> core_data

core_data %>%
  pivot_longer(cols = c("Pt", "Ct")) -> long_core_data

ggplot(core_data, aes(Depth, Pt, group = Depth)) +
  geom_boxplot(aes(fill = Depth)) + facet_grid(Atmosphere ~ Water)

ggplot(core_data, aes(Depth, Ct, group = Depth)) +
  geom_boxplot(aes(fill = Depth)) + facet_grid(Atmosphere ~ Water)
