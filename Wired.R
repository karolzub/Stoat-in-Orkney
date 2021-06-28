library(lmerTest)
library(popbio)
library(logistf)
library(spdep)
library(multcomp)
library(blme)
library(adehabitatHS)

trap<-read.table("C:/Orkney/Model/Wired.txt", header=TRUE, na.strings = "NA")
head(trap)

trap<-read.table("C:/Orkney/Model/Wired_SR.txt", header=TRUE, na.strings = "NA")
head(trap)

lm1 <- lmer(BM ~ Sex + Age + Wired + Month + (1|Part), data = trap)
summary(lm1)

sub <- subset (sub2, Wired2 == 0, data = sub2)
sub

sub1 <- subset (sub2, Age2 == 1, data = trap)
sub1

sub2 <- subset (trap, Sex == "M", data = trap)
sub2

sub3 <- subset (trap, Year == 1, data = trap)
sub3

sub3 <- subset (trap, Year == 2, data = trap)
sub3

males21 <- subset (sub3, Sex == "F", data = sub3)
males21

sub4 <- subset (trap, Location == "SR", data = trap)
sub4

sub5 <- subset (sub4, Sex == "M", data = sub4)
sub5
Age3 <- as.factor(trap$Age3)
lm1 <- glmer(Age3 ~ Sex + Wired + (1|Month), family = binomial, data = trap)
summary(lm1)

lm1 <- lmer(Length ~ Age + Wired + Month + (1|Part), data = sub2)
summary(lm1)

lm1 <- lmer(Zyg ~ Age + Wired + Month + (1|Location), data = sub2)
summary(lm1)

lm1 <- glm(BM ~ Age + Wired + Month + Part, data = sub2)
summary(lm1)

lm1 <- glm(BM ~ Age + Wired + Month, data = sub2)
summary(lm1)

lm1 <- glm(Zyg ~ Sex + Age + Wired + Month + Location, data = trap)
summary(lm1)

lm1 <- glm(Zyg ~ Age + Wired + Month + Location, data = sub)
summary(lm1)

lm1 <- glm(Length ~ Sex + Age + Wired + Month + Location, data = trap)
summary(lm1)

lm1 <- glm(Length ~ Age + Wired + Month + Location, data = sub)
summary(lm1)

plot(trap$Month, trap$BM)
boxplot(BM ~ Location, sub)
boxplot(BM ~ Month, sub)
boxplot(BM ~ Month, trap)

boxplot(Zyg ~ Part, sub2)
boxplot(Zyg ~ Wired, trap)
boxplot(BM ~ Year, sub5)

boxplot(BM ~ Sex, trap)

lm1 <- glm(BM ~ Age + Wired + Month, data = sub4)
summary(lm1)

dat1 <- xtabs(formula = Wired2 ~ Sex2, data=sub2)
dat1

write.table(dat1, file="C:/Orkney/Model/Habitat selection/dat1.txt", append = FALSE, sep = ",", dec = ".", row.names = TRUE, col.names = TRUE)

lm1 <- glm(Sex2 ~ Age + Wired + Month, data = trap)
summary(lm1)

plot(trap$BM, trap$Length)
hist(sub$BM)

