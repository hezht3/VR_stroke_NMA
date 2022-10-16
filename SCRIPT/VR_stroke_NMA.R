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

data <- escalc(measure = "SMD", vtype =  "UB",
       m1i = mean1, sd1i = sd1, n1i = n1,
       m2i = mean2, sd2i = sd2, n2i = n2, data = data)
data <- data %>% mutate(TE = yi, seTE = sqrt(vi))

as.matrix(table(data$source))   # no multi-arm study


##############################
# Fit random effect MA model #
##############################

rma(data = data %>% filter(treat1 == "HMD"),
    yi = yi,
    vi = vi,
    method = "SJ",
    slab = source,
    digits = 2,
    level = 95)

rma(data = data %>% filter(treat1 == "Kinect"),
    yi = yi,
    vi = vi,
    method = "SJ",
    slab = source,
    digits = 2,
    level = 95)

rma(data = data %>% filter(treat1 == "Wii"),
    yi = yi,
    vi = vi,
    method = "SJ",
    slab = source,
    digits = 2,
    level = 95)

#################
# Fit NMA model #
#################

netmeta.fixed.fit <- netmeta(TE = TE,
                             seTE = seTE,
                             treat1 = treat1,
                             treat2 = treat2,
                             studlab = source,
                             data = data,
                             sm = "SMD",
                             fixed = TRUE,
                             random = FALSE,
                             reference.group = "control",
                             details.chkmultiarm = TRUE,
                             sep.trts = " vs ")
summary(netmeta.fixed.fit)

decomp.design(netmeta.fixed.fit)

netmeta.random.fit <- netmeta(TE = TE,
                              seTE = seTE,
                              treat1 = treat1,
                              treat2 = treat2,
                              studlab = source,
                              data = data,
                              sm = "SMD",
                              fixed = FALSE,
                              random = TRUE,
                              reference.group = "control",
                              details.chkmultiarm = TRUE,
                              sep.trts = " vs ")
summary(netmeta.random.fit)


#################
# Visualization #
#################

# The Network Graph
netgraph(netmeta.random.fit)

# Direct and Indirect Evidence
d.evidence <- direct.evidence.plot(netmeta.random.fit)
plot(d.evidence)

# Effect estimate table
netleague(netmeta.random.fit, 
          bracket = "(", # use round brackets
          digits = 2)      # round to two digits

# Treatment ranking
netrank(netmeta.fixed.fit, small.values = "good")

# Forest plot
forest(netmeta.random.fit, 
       reference.group = "control",
       sortvar = TE,
       xlim = c(-2.5, 3),
       smlab = "Comparison: VR vs control",
       drop.reference.group = FALSE#,
       #label.left = "Favors Intervention",
       #label.right = "Favors Care As Usual"
       )

       