
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
