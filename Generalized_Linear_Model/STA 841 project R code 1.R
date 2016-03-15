library(lme4) 

setwd("D:/Course material/Statistics and population genetics/STA 841/project/proposal")
mouse.data <- read.table('tumor data.txt', header = T)
mouse.data$ECM1xDay <- mouse.data$ECM1_level * mouse.data$Day
mouse.data$Weight_0xDay <- mouse.data$Weight_0 * mouse.data$Day
mouse.data$ECM1xWeight_0xDay <- mouse.data$ECM1_level * mouse.data$Weight_0 * mouse.data$Day

gauss.m1 <- glm(formula = Tumor_volume ~ Day + ECM1xDay + Weight_0xDay, family = gaussian(link = log), data = mouse.data)
summary(gauss.m1)
jpeg("gauss_m1_diagnostic_plot.jpg")
par(mfrow=c(2,2))
plot(gauss.m1)
dev.off()
gamma.m3 <- glm(formula = Tumor_volume ~ Day + ECM1xDay + Weight_0xDay, family = Gamma(link = 'log'), data = mouse.data)
summary(gamma.m3)
jpeg("gamma_m3_diagnostic_plot.jpg")
par(mfrow=c(2,2))
plot(gamma.m3)
dev.off()

gamma.m1 <- glm(formula = Tumor_volume ~ Day, family = Gamma(link = 'log'), data = mouse.data)
gamma.m2 <- glm(formula = Tumor_volume ~ Day + ECM1xDay, family = Gamma(link = 'log'), data = mouse.data)
gamma.m4 <- glm(formula = Tumor_volume ~ Day + ECM1xDay + Weight_0xDay + ECM1xWeight_0xDay, family = Gamma(link = 'log'), data = mouse.data)

anova(gamma.m1, gamma.m2, gamma.m3, gamma.m4, test = 'Chi')

gamma.m5 <- glm(formula = Tumor_volume ~ Day + ECM1xDay + ECM1_squaredxDay, family = Gamma(link = 'log'), data = mouse.data)
gamma.m6 <- glm(formula = Tumor_volume ~ Day + ECM1xDay + ECM1_squaredxDay + Weight_0xDay, family = Gamma(link = 'log'), data = mouse.data)
anova(gamma.m1, gamma.m2, gamma.m5, gamma.m6, test = 'Chi')
plot(gamma.m5)
mouse.data$ECM1_squaredxDay <- mouse.data$ECM1_level^2 * mouse.data$Day
mouse.data$ECM1_cubedxDay <- mouse.data$ECM1_level^3 * mouse.data$Day
mouse.data$logECM1xDay <- log(mouse.data$ECM1_level) * mouse.data$Day
mouse.data$logECM1xWeight_0xDay <- log(mouse.data$ECM1_level) * mouse.data$Weight_0 * mouse.data$Day

gamma.m7 <- glm(formula = Tumor_volume ~ Day + ECM1xDay + ECM1_squaredxDay + ECM1_cubedxDay, family = Gamma(link = 'log'), data = mouse.data)
gamma.m8 <- glm(formula = Tumor_volume ~ Day + ECM1xDay + ECM1_squaredxDay + ECM1_cubedxDay + Weight_0xDay, family = Gamma(link = 'log'), data = mouse.data)
anova(gamma.m1, gamma.m2, gamma.m5, gamma.m7, gamma.m8, test = 'Chi')

gamma.m9 <- glm(formula = Tumor_volume ~ Day + logECM1xDay, family = Gamma(link = 'log'), data = mouse.data)
gamma.m10 <- glm(formula = Tumor_volume ~ Day + logECM1xDay + Weight_0xDay, family = Gamma(link = 'log'), data = mouse.data)
gamma.m11 <- glm(formula = Tumor_volume ~ Day + logECM1xDay + Weight_0xDay + logECM1xWeight_0xDay, family = Gamma(link = 'log'), data = mouse.data)
anova(gamma.m1, gamma.m9, gamma.m10, gamma.m11, test = 'Chi')

gamma.m12 <- glm(formula = Tumor_volume ~ poly(Day,2) + ECM1xDay, family = Gamma(link = 'log'), data = mouse.data)
gamma.m13 <- glm(formula = Tumor_volume ~ poly(Day,2) + ECM1xDay + Weight_0xDay, family = Gamma(link = 'log'), data = mouse.data)
anova(gamma.m1, gamma.m2, gamma.m12, gamma.m13, test = 'Chi')

gamma.m14 <- glm(formula = Tumor_volume ~ Day + ECM1xDay + Weight_0xDay + ECM1_squaredxDay, family = Gamma(link = 'log'), data = mouse.data)
gamma.m15 <- glm(formula = Tumor_volume ~ Day + ECM1xDay + Weight_0xDay + ECM1_squaredxDay + ECM1_cubedxDay, family = Gamma(link = 'log'), data = mouse.data)
anova(gamma.m1, gamma.m2, gamma.m3, gamma.m14, gamma.m15, test = 'Chi')

gamma.m4 <- glm(formula = Tumor_volume ~ Day + ECM1xDay + Weight_0xDay + ECM1_squaredxDay, family = Gamma(link = 'log'), data = mouse.data)
gamma.m5 <- glm(formula = Tumor_volume ~ Day + ECM1xDay + Weight_0xDay + ECM1_squaredxDay + ECM1_cubedxDay, family = Gamma(link = 'log'), data = mouse.data)
anova(gamma.m1, gamma.m2, gamma.m3, gamma.m4, gamma.m5, test = 'Chi')
summary(gamma.m5)
gauss.mm1 <- glmer(formula = Tumor_volume ~ Day + (1|Cage), family = gaussian(link = log), data = mouse.data)

gamma.m6 <- glmer(formula = Tumor_volume ~ Day + ECM1xDay + Weight_0xDay + ECM1_squaredxDay + ECM1_cubedxDay + (Day|Cage), family = Gamma(link = 'log'), data = mouse.data)
gamma.m6 <- glmer(formula = Tumor_volume ~ Day + ECM1xDay + Weight_0xDay + (1|Cage), family = Gamma(link = 'log'), data = mouse.data)
gamma.m6 <- glmer(formula = Tumor_volume ~ Day + ECM1xDay + Weight_0xDay + ECM1_squaredxDay + (1|Cage), family = Gamma(link = 'log'), data = mouse.data)