
library(readr)
library(tidyr)

# Read data ====

EOG_PO60 <- read_csv("../../R data/EOG_PO60", col_names = FALSE)
EOG_PO70 <- read_csv("../../R data/EOG_PO70", col_names = FALSE)

EOG_GP60 <- read_csv("../../R data/EOG_GP60", col_names = FALSE)
EOG_GP70 <- read_csv("../../R data/EOG_GP70", col_names = FALSE)

# Rename cols
colnames(EOG_PO60) <- c("70", "75", "80", "85", "90", "95")
colnames(EOG_PO70) <- c("75", "80", "85", "90", "95")

colnames(EOG_GP60) <- c("ISI0", "ISI60", "ISI120", "ISI240")
colnames(EOG_GP70) <- c("ISI0", "ISI60", "ISI120", "ISI240")

# Pivot to long format ====
EOG_PO60_long <- pivot_longer(EOG_PO60, cols = c("70", "75", "80", "85", "90", "95"))
EOG_PO70_long <- pivot_longer(EOG_PO70, cols = c("75", "80", "85", "90", "95"))

EOG_GP60_long <- pivot_longer(EOG_GP60, cols = c("ISI0", "ISI60", "ISI120", "ISI240"))
EOG_GP70_long <- pivot_longer(EOG_GP70, cols = c("ISI0", "ISI60", "ISI120", "ISI240"))

# Run ANOVAs ====
PO60_aov <- aov(data = EOG_PO60_long, value ~ name)
PO70_aov <- aov(data = EOG_PO70_long, value ~ name)

GP60_aov <- aov(data = EOG_GP60_long, value ~ name)
GP70_aov <- aov(data = EOG_GP70_long, value ~ name)

# ANOVA summaries ====

summary(PO60_aov)
summary(PO70_aov)
summary(GP60_aov)
summary(GP70_aov)

# Inhibition ====

GP60_ISI0_inhib <- 1 - (mean(EOG_GP60$ISI0)/(mean(EOG_PO60$`90`)))
GP60_ISI60_inhib <- 1 - (mean(EOG_GP60$ISI60)/(mean(EOG_PO60$`90`)))
GP60_ISI120_inhib <- 1 - (mean(EOG_GP60$ISI120)/(mean(EOG_PO60$`90`)))
GP60_ISI240_inhib <- 1 - (mean(EOG_GP60$ISI240)/(mean(EOG_PO60$`90`)))

GP70_ISI0_inhib <- 1 - (mean(EOG_GP70$ISI0)/(mean(EOG_PO70$`90`)))
GP70_ISI60_inhib <- 1 - (mean(EOG_GP70$ISI60)/(mean(EOG_PO70$`90`)))
GP70_ISI120_inhib <- 1 - (mean(EOG_GP70$ISI120)/(mean(EOG_PO70$`90`)))
GP70_ISI240_inhib <- 1 - (mean(EOG_GP70$ISI240)/(mean(EOG_PO70$`90`)))

# t-tests

t_GP60_ISI0 <- t.test(EOG_GP60$ISI0, EOG_PO60$`90`, paired = T)
t_GP60_ISI60 <- t.test(EOG_GP60$ISI60, EOG_PO60$`90`, paired = T)
t_GP60_ISI120 <- t.test(EOG_GP60$ISI120, EOG_PO60$`90`, paired = T)
t_GP60_ISI240 <- t.test(EOG_GP60$ISI240, EOG_PO60$`90`, paired = T)

t_GP70_ISI0 <- t.test(EOG_GP70$ISI0, EOG_PO70$`90`, paired = T)
t_GP70_ISI60 <- t.test(EOG_GP70$ISI60, EOG_PO70$`90`, paired = T)
t_GP70_ISI120 <- t.test(EOG_GP70$ISI120, EOG_PO70$`90`, paired = T)
t_GP70_ISI240 <- t.test(EOG_GP70$ISI240, EOG_PO70$`90`, paired = T)

t_GP60_ISI0
t_GP60_ISI60
t_GP60_ISI120
t_GP60_ISI240

t_GP70_ISI0
t_GP70_ISI60
t_GP70_ISI120
t_GP70_ISI240

# list and correct p-values

p_values = c(t_GP60_ISI0$p.value, t_GP60_ISI60$p.value, t_GP60_ISI120$p.value, t_GP60_ISI240$p.value,
             t_GP70_ISI0$p.value, t_GP70_ISI60$p.value, t_GP70_ISI120$p.value, t_GP70_ISI240$p.value)

p_names = c("60_0", "60_60", "60_120", "60_240",
            "70_0", "70_60", "70_120", "70_240")

tab <- as.table(rbind(p.adjust(p_values, method = "BH")))
colnames(tab) <- p_names
