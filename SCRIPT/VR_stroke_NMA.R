#################
# VR Stroke NMA #
#################

setwd("/Users/zhengtinghe/Library/CloudStorage/OneDrive-JohnsHopkins/Research/rehabilitation/VR stroke NMA/VR_stroke_NMA")

require(tidyverse)
require(metafor)
require(netmeta)
require(dmetar)


#################
# Prepare input #
#################

data <- read_csv("./INPUT/VR_stroke_NMA_data.csv")

data <- escalc(measure = "MD", vtype =  "HO",
       m1i = mean1, sd1i = sd1, n1i = n1,
       m2i = mean2, sd2i = sd2, n2i = n2, data = data)   # assumes homoscedasticity
data <- data %>% mutate(TE = yi, seTE = sqrt(vi))

as.matrix(table(data$source))   # no multi-arm study


#################
# Fit NMA model #
#################

netmeta.fixed.fit <- netmeta(TE = TE,
                             seTE = seTE,
                             treat1 = treat1,
                             treat2 = treat2,
                             studlab = source,
                             data = data,
                             sm = "MD",
                             fixed = TRUE,
                             random = FALSE,
                             reference.group = "control",
                             details.chkmultiarm = TRUE,
                             sep.trts = " vs ")
summary(netmeta.fixed.fit)

decomp.design(netmeta.fixed.fit)


#################
# Visualization #
#################

# The Network Graph
netgraph(netmeta.fixed.fit)

# Direct and Indirect Evidence
d.evidence <- direct.evidence.plot(netmeta.fixed.fit)
plot(d.evidence)

# Effect estimate table
netleague(netmeta.fixed.fit, 
          bracket = "(", # use round brackets
          digits = 2)      # round to two digits

# Treatment ranking
netrank(netmeta.fixed.fit, small.values = "good")

# Forest plot
forest(netmeta.fixed.fit, 
       reference.group = "control",
       sortvar = TE,
       xlim = c(-2.5, 15),
       smlab = "Comparison: VR vs control",
       drop.reference.group = FALSE#,
       #label.left = "Favors Intervention",
       #label.right = "Favors Care As Usual"
       )

       