
# Libraries, load and arrange data ====
library(readr)
library(tidyr)

#Read data from CSV
amp_PO60_L <- read_csv("../../R data/amp_PO60_L", col_names = FALSE)
amp_PO60_R <- read_csv("../../R data/amp_PO60_R", col_names = FALSE)

amp_PO70_L <- read_csv("../../R data/amp_PO70_L", col_names = FALSE)
amp_PO70_R <- read_csv("../../R data/amp_PO70_R", col_names = FALSE)

amp_GP60_L <- read_csv("../../R data/amp_GP60_L", col_names = FALSE)
amp_GP60_R <- read_csv("../../R data/amp_GP60_R", col_names = FALSE)

amp_GP70_L <- read_csv("../../R data/amp_GP70_L", col_names = FALSE)
amp_GP70_R <- read_csv("../../R data/amp_GP70_R", col_names = FALSE)

#Rename cols
colnames(amp_PO60_L) <- c("70", "75", "80", "85", "90", "95")
colnames(amp_PO60_R) <- c("70", "75", "80", "85", "90", "95")

colnames(amp_PO70_L) <- c("75", "80", "85", "90", "95")
colnames(amp_PO70_R) <- c("75", "80", "85", "90", "95")

colnames(amp_GP60_L) <- c("ISI0", "ISI60", "ISI120", "ISI240")
colnames(amp_GP60_R) <- c("ISI0", "ISI60", "ISI120", "ISI240")

colnames(amp_GP70_L) <- c("ISI0", "ISI60", "ISI120", "ISI240")
colnames(amp_GP70_R) <- c("ISI0", "ISI60", "ISI120", "ISI240")

#Make to long format for aov
PO60_L_long <- pivot_longer(amp_PO60_L, cols = c("70", "75", "80", "85", "90", "95"))
PO60_R_long <- pivot_longer(amp_PO60_R, cols = c("70", "75", "80", "85", "90", "95"))

PO70_L_long <- pivot_longer(amp_PO60_L, cols = c("75", "80", "85", "90", "95"))
PO70_R_long <- pivot_longer(amp_PO60_R, cols = c("75", "80", "85", "90", "95"))

GP60_L_long <- pivot_longer(amp_PO60_L, cols = c("ISI0", "ISI60", "ISI120", "ISI240"))
GP60_R_long <- pivot_longer(amp_PO60_R, cols = c("ISI0", "ISI60", "ISI120", "ISI240"))

GP70_L_long <- pivot_longer(amp_PO60_L, cols = c("ISI0", "ISI60", "ISI120", "ISI240"))
GP70_R_long <- pivot_longer(amp_PO60_R, cols = c("ISI0", "ISI60", "ISI120", "ISI240"))

#Run anovas ====
PO60_L_aov <- aov(data = PO60_L_long, value ~ name)
PO60_R_aov <- aov(data = PO60_R_long, value ~ name)

PO70_L_aov <- aov(data = PO70_L_long, value ~ name)
PO70_R_aov <- aov(data = PO70_R_long, value ~ name)

GP60_L_aov <- aov(data = GP60_L_long, value ~ name)
GP60_R_aov <- aov(data = GP60_R_long, value ~ name)

GP70_L_aov <- aov(data = GP70_L_long, value ~ name)
GP70_R_aov <- aov(data = GP70_R_long, value ~ name)

#Print summaries of aovs ====

summary(PO60_L_aov)
summary(PO60_R_aov)

summary(PO70_L_aov)
summary(PO70_R_aov)

summary(GP60_L_aov)
summary(GP60_R_aov)

summary(GP70_L_aov)
summary(GP70_R_aov)

#t-test PO_90 vs ISI240 ====

t_60L_240 <- t.test(amp_PO60_L$`90`, amp_GP60_L$ISI240)
t_60L_120 <- t.test(amp_PO60_L$`90`, amp_GP60_L$ISI120)
t_60L_60 <- t.test(amp_PO60_L$`90`, amp_GP60_L$ISI60)
t_60L_0 <- t.test(amp_PO60_L$`90`, amp_GP60_L$ISI0)

t_60R_240 <- t.test(amp_PO60_R$`90`, amp_GP60_R$ISI240)
t_60R_120 <- t.test(amp_PO60_R$`90`, amp_GP60_R$ISI120)
t_60R_60 <- t.test(amp_PO60_R$`90`, amp_GP60_R$ISI60)
t_60R_0 <- t.test(amp_PO60_R$`90`, amp_GP60_R$ISI0)

t_70L_240 <- t.test(amp_PO70_L$`90`, amp_GP70_L$ISI240)
t_70L_120 <- t.test(amp_PO70_L$`90`, amp_GP70_L$ISI120)
t_70L_60 <- t.test(amp_PO70_L$`90`, amp_GP70_L$ISI60)
t_70L_0 <- t.test(amp_PO70_L$`90`, amp_GP70_L$ISI0)

t_70R_240 <- t.test(amp_PO70_R$`90`, amp_GP70_R$ISI240)
t_70R_120 <- t.test(amp_PO70_R$`90`, amp_GP70_R$ISI120)
t_70R_60 <- t.test(amp_PO70_R$`90`, amp_GP70_R$ISI60)
t_70R_0 <- t.test(amp_PO70_R$`90`, amp_GP70_R$ISI0)

#P.adjust and print table ====

p_values = c(t_60L_240$p.value, t_60L_120$p.value, t_60L_60$p.value, t_60L_0$p.value,
             t_60R_240$p.value, t_60R_120$p.value, t_60R_60$p.value, t_60R_0$p.value,
             t_70L_240$p.value, t_70L_120$p.value, t_70L_60$p.value, t_70L_0$p.value,
             t_70R_240$p.value, t_70R_120$p.value, t_70R_60$p.value, t_70R_0$p.value)

p_names = c("60L240", "60L120", "60L60", "60L0",
          "60R240", "60R120", "60R60", "60R0",
          "70L240", "70L120", "70L60", "70R0",
          "70R240", "70R120", "70R60", "70RL0")

tab <- as.table(rbind(p.adjust(p_values, method = "BH")))
colnames(tab) <- p_names








