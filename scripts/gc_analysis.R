# Script to ananlyse growth curves, producing growth rate, lag time, AUC etc
# https://rfortherestofus.com/2021/02/how-to-use-git-github-with-r great guide
# for using github and R!

# Dependencies -----------------------------------------------------------------
library(tidyverse)
library(drc)              # DRC curve fit functions

# Functions --------------------------------------------------------------------
## Imlpement the growthcurver fit in drc, from following 
## ~ https://www.r-bloggers.com/2020/02/self-starting-routines-for-nonlinear-regression-models/
growth_curver_model <- function(k, n0, r, t){
  
  # k = max growth
  # N0 = min growth
  # r = growth rate
  # t = time (predictor?)
  
  k / (1 + ( (k - n0) / n0) * exp(-r * t))
  
}

growth_curver_fit <- function(){
  
  fct <- function(x, parm) {
    growth_curver_model(k  = parm[, 1], # max growth
                        n0 = parm[, 2], # min growth
                        r  = parm[, 3], # growth rate
                        t  = x)
  }
  
  ssfct <- function(data){
    x   <- data[, 1]  #time
    y   <- data[, 2]  #OD600
    
    k  <- max(y)
    n0 <- min(y[y > 0])
    
    # make an initial estimate for r
    glm_mod <- stats::glm(y / k ~ x,
                          family = stats::quasibinomial("logit"))
    
    r <- stats::coef(glm_mod)[[2]]   # slope
    if (r <= 0) {
      # the slope should only be positive for a growing culture, so default
      # to something small
      r <- 0.001
    }
    
    return(c(k, n0, r))
    
  }
  names <- c("max_growth", "min_growth", "rate")
  text <- "Growth Curver Equation"
  
  ## Returning the function with self starter and names
  returnList <- list(fct = fct, ssfct = ssfct, names = names, text = text)
  class(returnList) <- "drcMean"
  invisible(returnList)
  
}

# Options ----------------------------------------------------------------------
r2_OD600_file <- "./raw_data/20250228_strains_GC_r2.csv"

r2_plate_plan_file <- "./raw_data/strains_GC_r2_well_key.csv"

# Load Data --------------------------------------------------------------------
r2_plate_plan <- read.csv(r2_plate_plan_file) %>%
  mutate(well = paste0(row, sprintf("%02d", column)))

r2_OD600_raw <- read.csv(r2_OD600_file,
                         skip = 6) %>%
  filter(row_number() != 1) %>%
  rename(time_h = X) %>%
  dplyr::select(-Well) %>%
  pivot_longer(cols = -time_h,
               names_to = "well",
               values_to = "OD600") %>%
  left_join(r2_plate_plan,
            by = c("well" = "well")) %>%
  mutate(time_h = as.numeric(time_h),
         OD600  = as.numeric(OD600),
         biological_rep = 1) %>%
  filter(time_h <= 24)

## add in other repeats here
all_raw_data <- r2_OD600_raw

# Process data -----------------------------------------------------------------
blanks <- all_raw_data %>%
  filter(time_h == 0) %>%
  transmute(well = well,
            blank_OD600 = OD600)

all_OD600_data <- all_raw_data %>%
  left_join(blanks,
            by = c("well" = "well")) %>%
  mutate(corrected_OD600 = case_when((OD600 - blank_OD600) <  0 ~ 0,
                                     (OD600 - blank_OD600) >= 0 ~ (OD600 - blank_OD600)),
         Strain = factor(x      = Strain, 
                         levels = c("PAO1",
                                    "PA14-N",
                                    "PA14-N  ∆sid",
                                    "PA14-O",
                                    "PA14-O ∆rhlA",
                                    "PA14-O ∆rhlB",
                                    "PA14-O ∆rhlC",
                                    "59.37",
                                    "59.38",
                                    "59.40",
                                    "60.92",
                                    "Blank")))

# Examine techincal repeats ----------------------------------------------------
techincal_rep_plot <- ggplot(data    = all_OD600_data,
                             mapping = aes(x = time_h,
                                           y = corrected_OD600)) +
  stat_summary(geom        = "errorbar",
               fun.data    = "mean_cl_boot",
               fun.args    = list(conf.int = 0.95),
               show.legend = FALSE,
               alpha       = 0.5) +
  stat_summary(geom = "point",
               fun  = mean,
               color = "black") +
  facet_grid(rows = vars(biological_rep),
             cols = vars(Strain)) +
  scale_x_continuous(breaks = seq(0, 24, 12),
                     name = "Time (h)") +
  scale_y_continuous(breaks = seq(0, 1.5, 0.5),
                     limits = c(0, 1.6),
                     name   = "OD600") +
  theme_bw()
  
techincal_rep_plot

# Examine biological repeats ---------------------------------------------------
bio_rep_data <- all_OD600_data %>%
  group_by(biological_rep,
           Strain,
           time_h) %>%
  summarise(mean_OD600 = mean(corrected_OD600)) 

## Plot curves 
biological_rep_plot <- ggplot(data =  bio_rep_data,
                              mapping = aes(x = time_h,
                                            y = mean_OD600)) +
  stat_summary(geom        = "errorbar",
               fun.data    = "mean_cl_boot",
               fun.args    = list(conf.int = 0.95),
               show.legend = FALSE,
               alpha       = 0.5) +
  stat_summary(geom = "point",
               fun  = mean,
               color = "black") +
  stat_smooth(method      = drm,
              method.args = list(fct = growth_curver_fit()),
              formula     = y ~ x,
              se          = FALSE,
              fullrange   = TRUE) +
  facet_grid(rows = vars(biological_rep),
             cols = vars(Strain)) +
  scale_x_continuous(breaks = seq(0, 24, 12),
                     name = "Time (h)") +
  scale_y_continuous(breaks = seq(0, 1.5, 0.5),
                     limits = c(0, 1.6),
                     name   = "OD600") +
  theme_bw()

biological_rep_plot

# Plot summary data ------------------------------------------------------------
## fit models to get parameters 
curve_models <- drm(formula = mean_OD600 ~ time_h,
                    curveid = Strain,
                    data    = bio_rep_data,
                    fct     = growth_curver_fit())

curve_params <- as.data.frame(summary(curve_models)[["coefficients"]]) %>%
  rownames_to_column(var = "rowname") %>%
  separate(col  = rowname,
           into = c("parameter", "strain"),
           sep  = ":",
           extra = "merge")

params_plot <- ggplot(data = curve_params %>%
                        filter(parameter != "min_growth") %>%
                        filter(strain != "Blank"),
                      mapping = aes(x = strain,
                                    y = Estimate,
                                    ymin = Estimate - `Std. Error`,
                                    ymax = Estimate + `Std. Error`)) +
  geom_pointrange() +
  facet_grid(rows = vars(parameter),
             scale = "free_y")








