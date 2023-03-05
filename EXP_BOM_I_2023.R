# //////////////////////////////////////////////////////////////////////////////////
# INSTITUTO TECNOLOGICO DE COSTA RICA
# Escuela de Ingenieria en Construccion
# https://www.tec.ac.cr
# Session: BOMBAS CENTRIFUGAS

# M.Sc. Eng. Maikel Mendez M
# Water Resources + GIS + DataScience
# Instituto Tecnologico de Costa Rica
# https://www.tec.ac.cr
# https://orcid.org/0000-0003-1919-141X
# https://www.scopus.com/authid/detail.uri?authorId=51665581300
# https://scholar.google.com/citations?user=JnmSVFYAAAAJ&hl=en
# https://www.youtube.com/c/maikelmendez
# https://github.com/maikelonu
# //////////////////////////////////////////////////////////////////////////////////

# INFO:
# Modelos_polinomiales no-lineales 
# Diagnostico_numerico_modelos
# Prediccion_(CI + PI)
# Analisis_grafico_ggplot2
# Shapiro_test
# //////////////////////////////////////////////////////////////////////////////////

# Workspace is cleared
rm(list = ls())

# Working directory is selected
# setwd("/media/maikel/Trabajo/R_ITC/R_LABHYD/EXP_BOM")
setwd("C:/DATOS/R_ITC/R_LABHYD/EXP_BOM")

# CRAN libraries are loaded
# require(Agreement)
require(DescTools)
require(effects)
require(ggplot2)
require(MASS)
require(nls2)
require(nlstools)
require(pastecs)
require(reshape)
require(visreg)

# /////////////////////////////////////////////////////////////
# BLOCK: Custom function, round data.frame to specif digits
# /////////////////////////////////////////////////////////////
round_df <- function(df, digits) {
  options(scipen = 0)
  options(scipen = -3)
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] <- round(df[,nums], digits = digits)
  return(df)
}

# ////////////////////////////////////////////////////////
# BLOCK: Data input
# ////////////////////////////////////////////////////////

# Input data is loaded and a data.frame is created
df.base <- read.table("bombas3.txt", header = TRUE)

# Desc {DescTools} function is requested
Desc(df.base, plotit = TRUE)

# names {base} function is requested
names(df.base)

# theoretical pump-data full capacity
q_lps_teo <- c(0.00,0.95,1.89,2.52,3.15,4.42,5.68,6.31,6.62) # LPS
 hm_m_teo <- c(38.87,39.63,38.11,37.35,35.82,31.25,24.39,19.05,16.01) # m

# a data.frame is created
df.pump <- data.frame(q_lps_teo,hm_m_teo)

# ////////////////////////////////////////////////////////
# BLOCK: Regression Models
# ////////////////////////////////////////////////////////

# ------------------------------------------------
# lm 2nd degree polynomial
# ------------------------------------------------

# A lm {stats}  Least Squares model is fitted using a 2nd degree polynomial
mod.lm2 <- lm(hm_m_teo ~ poly(q_lps_teo, 2, raw=TRUE), data = df.pump)

# A model summary is requested
summary(mod.lm2)

#The equation of polynomial of degree 2 of the model is:
# f(x)=  38.68824 + 1.35979(x) - 0.70814(x^2)

# /////////////////////////////////////////////////
# lm SECTIONS:
# /////////////////////////////////////////////////
# 
# Residuals section:
# it provides a quick summary (min, 1Q, median, 3Q, max) of the distribution.
# 
# Coefficients section: each coefficient is a Gaussian random variable
# Estimate represents the mean distribution of the variable
# Std.Error displays the standard error of the variable
# the t value is Estimate divided by Std.Error
# the p value indicates the probability of getting a value larger than the t value
#
# Residual standard error outputs the standard deviation of residuals
#
# The degree of freedom indicates the differences between the observation 
# in training samples and the number used in the model
#
# Multiple R-squared is obtained by dividing the sum of squares.
#
# Adjusted R-squared uses an unbiased estimate, and will 
# be slightly less than multiple R-squared

# /////////////////////////////////////////////////

# shapiro.test {stats} Normality Test is applied
# if p-value > 0.05 then normality stands true, meaning that
# the variable is parametric
shapiro.test(mod.lm2$residuals)

# A full lm-model diagnostic plot is requested
oldpar.lm2 <- par(oma=c(0,0,3,0), mfrow=c(2,2))
plot(mod.lm2)
par(oldpar.lm2)

# flowrate observations are isolated to be used in lm models
# in order to harmonize observed and theoretical values
q_lps_teo <- df.base$q_lps

# An observations data.frame is created
q.observed <- data.frame(q_lps_teo)

# Observations data.frame is used to predict theoretical values 
df.predCI <- as.data.frame(predict(mod.lm2, q.observed, interval = "confidence"))

# df.predCI data.frame is rounded to 3 decimals
df.predCI <- round(df.predCI,3)

# cbind {base} function is used to join data.frames
df.output <- cbind(df.base,df.predCI)

# A simple deviation (SD) column is created
df.output$SD <- (df.output$fit - df.output$hm_m)

# shapiro.test {stats} Normality Test is applied to df.base$SD
# if p-value > 0.05 then normality stands true, meaning that
# the variable is parametric
shapiro.test(df.output$SD)

# Descriptive statistics are requested and rounded to 5 decimals
df.output.desc <- as.data.frame(round(stat.desc(df.output$SD), 3))

# data.frame variables are renamed
names(df.output.desc) <- c("SD")

# ////////////////////////////////////////////////////////
# BLOCK: Graphical Analysis
# ////////////////////////////////////////////////////////

# theretical pump-data NPSH-bias
q_lps_npsh <- c(3.15,4.57,5.20,5.84,6.94) # LPS
 hm_m_npsh <- c(7.62,6.10,4.57,3.05,1.52) # m

# a data.frame is created
df.npsh <- data.frame(q_lps_npsh,hm_m_npsh)

# A ggplot object is created
fg01 <- ggplot() +
  geom_smooth(aes(x = q_lps_teo,y = hm_m_teo),data=df.pump, size = 1.5) +
  geom_point(aes(x = q_lps_teo,y = hm_m_teo),data=df.pump,size = 2.5, shape = 15) +
  geom_vline(aes(xintercept = q_lps_npsh,colour = hm_m_npsh),data=df.npsh,size = 0.55,linetype = 2) +
  geom_text(aes(x = q_lps_npsh,y = hm_m_npsh,label = hm_m_npsh,colour = hm_m_npsh),data=df.npsh,size = 7.5,angle = 41.940000000000005,parse = FALSE) +
  scale_colour_gradient(guide = guide_legend(),low = '#0000ff',high = '#ff0000') +
  geom_point(aes(x = q_lps,y = hm_m),data=df.base,shape = 13,colour = '#cc0000',size = 3.5) +
  geom_line(aes(x = q_lps,y = hm_m),data=df.base,colour = '#cc0033',size = 1.25,linetype = 2) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  ggtitle("Desempeno de bombas centrifugas experimental vs teorico (marcado por NPSH)") +
  xlab("Caudal (LPS)") +
  ylab("Hm (m)") +
  theme_bw(base_size = 22.0)

# A ggplot object is requested
fg01

# round_df function is applied to relevant data.frames
df.output <- round_df(df=df.output, digits=3)
df.output.desc <- round_df(df=df.output.desc, digits=3)

# Objects to export
# df.output, df.output.desc, summary(mod.lm2)
# Shapiro.test
# oldpar.lm2, fg01
write.csv(df.output, file = "df.output.csv")
write.csv(df.output.desc, file = "df.output.desc.csv")

# /////////////////////////////////////////////////////////////
# END OF SCRIPT
# /////////////////////////////////////////////////////////////
