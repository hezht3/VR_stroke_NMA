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

data <- data %>% 
    mutate(chronic = ifelse(source == "Mekbib, 2021" | source == "Kim, 2018" |
                                source == "Choi, 2014" | source == "Kong, 2016",
                            "subacute", "chronic"))

data <- escalc(measure = "SMD", vtype =  "UB",
       m1i = mean1, sd1i = sd1, n1i = n1,
       m2i = mean2, sd2i = sd2, n2i = n2, data = data)
data <- data %>% mutate(TE = yi, seTE = sqrt(vi))

as.matrix(table(data$source))   # no multi-arm study


##############################
# Fit random effect MA model #
##############################

data$order <- seq(nrow(data), 1, -1)

# Overall
vr_vs_control <- rma(data = data,
                     yi = yi,
                     vi = vi,
                     method = "DL",
                     slab = source,
                     digits = 2,
                     level = 95)

# Subgroup
hmd_vs_control <- rma(data = data,
                      yi = yi,
                      vi = vi,
                      subset = (treat1 == "HMD"),
                      method = "DL",
                      slab = source,
                      digits = 2,
                      level = 95)

kinect_vs_control <- rma(data = data,
                         yi = yi,
                         vi = vi,
                         subset = (treat1 == "Kinect"),
                         method = "DL",
                         slab = source,
                         digits = 2,
                         level = 95)

wii_vs_control <- rma(data = data,
                      yi = yi,
                      vi = vi,
                      subset = (treat1 == "Wii"),
                      method = "DL",
                      slab = source,
                      digits = 2,
                      level = 95)

tiff("./OUTPUT/Figure 1. Random effect meta analysis forest plot.tiff",
     width = 3000, height = 2300, pointsize = 10, res = 300)
forest(vr_vs_control, xlim = c(-8, 10), at = seq(-2, 8, 2), 
       ilab = cbind(data$n1, data$n2),
       ilab.xpos = c(-4,-2), cex = 0.75, ylim = c(-1, 28), 
       order = order(data$order), rows = c(3:7, 12:15, 20:24),
       xlab = "Standard Mean Difference", mlab = "RE Model for All VR", psize = 1)
text(-8, 27, "Author(s) and Year", pos = 4, cex = 0.8, font = 2)
text(-4, 27, "VR", cex = 0.8, font = 2)
text(-2, 27, "Control", cex = 0.8, font = 2)
text(10, 27, "Standard mean difference [95% CI]", pos = 2, cex = 0.8, font = 2)
text(-2.8, 28, "No. of participants", cex = 0.8, font = 2)
text(-8, c(25,16,8), pos = 4, c("Immersive virtual reality",
                                "Microsoft Kinect",
                                "Nintendo Wii"), cex = 0.8, font = 4)
addpoly(hmd_vs_control, row = 18.5, cex = 0.75, mlab = "RE Model for Subgroup")
addpoly(kinect_vs_control, row = 10.5, cex =.75, mlab = "RE Model for Subgroup")
addpoly(wii_vs_control, row = 1.5, cex =.75, mlab = "RE Model for Subgroup")
dev.off()

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

forest(netmeta.random.fit,
       ref = "control",
       pooled = "random",
       digits = 2,
       smlab = "Random effects model",
       xlab = "HbA1c difference",
       leftlabs = "Contrast to placebo")


#################
# Visualization #
#################

# The Network Graph
tiff("./OUTPUT/Figure 3. Network graph.tiff",
     width = 3000, height = 2000, pointsize = 10, res = 300)
netgraph(netmeta.random.fit, start = "random", iterate = TRUE, col.multiarm = "purple",
         points = TRUE, cex.points = 1.5, cex = 1, 
         labels = c("Conventional rehabilitation", "Head mounted devices (immersive)",
                    "Microsoft Kinect (non-immersive)", "Nintendo Wii (non-immersive)"))
dev.off()

# Direct and Indirect Evidence
d.evidence <- direct.evidence.plot(netmeta.random.fit)
plot(d.evidence)

# Effect estimate table
netleague(netmeta.random.fit, 
          random = netmeta.random.fit$random,
          seq = c("HMD", "Kinect", "Wii", "control"),
          bracket = "(", # use round brackets
          digits = 2)      # round to two digits

# Treatment ranking
netrank(netmeta.random.fit, small.values = "good")

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

comparison <- netsplit(netmeta.random.fit)
comparison$comparison <- c("Immersive virtual reality vs Conventional rehabilitation",
                           "Microsoft Kinect vs Conventional rehabilitation",
                           "Nintendo Wii vs Conventional rehabilitation",
                           "Immersive virtual reality vs Microsoft Kinect",
                           "Immersive virtual reality vs Nintendo Wii",
                           "Microsoft Kinect vs Nintendo Wii")
tiff("./OUTPUT/Figure 4. Network meta analysis (spllit).tiff",
     width = 2500, height = 1800, pointsize = 10, res = 300)
comparison %>% forest(show = "all")
dev.off()

# Pairwise forest plot
contrasts <- data.frame(control = rep(NA, 4), HMD = rep(NA, 4), Kinect = rep(NA, 4), Wii = rep(NA, 4))
rownames(contrasts) <- c("control", "HMD", "Kinect", "Wii")
for(i in  1:nrow(contrasts)) {
    for(j in 1:ncol(contrasts)) {
        contrasts[i, j] <- paste(rownames(contrasts)[i], colnames(contrasts)[j], sep = " vs ")
    }
}

pairwise_es <- list(contrasts,
                    as.data.frame(netmeta.random.fit$TE.random),
                    as.data.frame(netmeta.random.fit$lower.random),
                    as.data.frame(netmeta.random.fit$upper.random)) %>% 
    map(~ pivot_longer(.x, everything(), names_to = "names", values_to = "values"))

pairwise_es <- tibble(contrasts = pairwise_es[[1]]$values, SMD = pairwise_es[[2]]$values,
                      LCI = pairwise_es[[3]]$values, UCI = pairwise_es[[4]]$values) %>% 
    filter(SMD != 0) %>% 
    mutate(labels = paste(round(SMD, 2), " (", format(round(LCI, 2), digits = 2), ", ", format(round(UCI, 2), digits = 2), ")", sep = "")) %>% 
    arrange(SMD)

pairwise_es$contrasts <- factor(pairwise_es$contrasts, levels = pairwise_es$contrasts)
tiff("./OUTPUT/Figure 2. Network meta analysis forest plot.tiff",
     width = 3000, height = 1500, pointsize = 10, res = 300)
pairwise_es %>% 
    ggplot(aes(x = SMD, y = contrasts)) +
    geom_point() +
    geom_linerange(aes(xmin = LCI, xmax = UCI)) +
    geom_vline(xintercept = 0) +
    xlab("Standard Mean Difference (95% Confidence Interval))") +
    ylab("") +
    theme(axis.ticks = element_blank(), legend.background = element_blank(), 
          legend.key = element_blank(), panel.background = element_blank(), 
          panel.border = element_rect(fill = NA, 
                                      colour = "grey20"), strip.background = element_blank(), 
          plot.background = element_blank(), complete = TRUE) +
    guides(y.sec = ggh4x::guide_axis_manual(labels = c(pairwise_es$labels)))
dev.off()

# Funnel plot
tiff("./OUTPUT/Figure 5. Network meta analysis funnel plot.tiff",
     width = 2000, height = 1200, pointsize = 10, res = 300)
funnel(netmeta.random.fit,
       order = c("HMD", "Kinect", "Wii", "control"),
       col = c("blue", "red", "purple"),
       linreg = TRUE)
dev.off()
