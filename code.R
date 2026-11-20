library(ThresholdROC)
library(GsymPoint)
library(pROC)
library(ThresholdROCsurvival)
library(CatPredi)
library(survival)



# 2-state example: AD
data(AD)

# descriptives
pdf("descriptive_AD.pdf", width = 5, height = 5.5)
boxplot(Tau~Status, AD, col = c("green", "red"), names = c("Control", "AD"), ylim = c(0, 2000),
        ylab = "Tau protein (pg/mL)")
dev.off()

# gsym.point (GsymPoint)
res <- gsym.point("GPQ", AD, "Tau", "Status", tag.healthy = 0)
summary(res)

# coords (pROC)
roc.obj <- roc(Status~Tau, AD)
coords(roc.obj, "best", ret = c("threshold", "sens", "spec"), best.method = "youden")
coords(roc.obj, "best", ret = c("threshold", "sens", "spec"), best.method = "closest.topleft")

# thres2 (ThresholdROC)
k1 <- subset(AD, Status == 0)$Tau
k2 <- subset(AD, Status == 1)$Tau
thres <- thres2(k1, k2, 0.2, matrix(c(0, 0, 1, 1), byr = TRUE, nc = 2, nr = 2),
                method = "empirical", ci.method = "boot")
thres
pdf("plot_thres2.pdf", width = 8, height = 6)
plot(thres, legend = FALSE, col = c(3, 2, 4), lwd = 3, ylim = c(0, 0.006),
     xlab = "Tau protein (pg/mL)")
legend("topright", c("Control", "AD", "Threshold", "CI"), col = c(3, 2, 4, 4), lwd = 3,
       lty = c(1, 1, 1, 2))
dev.off()




# 3-state example: chemo
data(chemo)
chemo$resp <- factor(chemo$resp, c("Stable-disease", "Partial-response", "Complete-response"))

# descriptives
pdf("descriptive_chemo.pdf", width = 5, height = 6)
boxplot(SUV~resp, chemo, col = c("red", "green", "blue"), ylim = c(0, 100), xlab = "",
        ylab = "SUV (%)", names = c("SD", "PR", "CR"))
dev.off()

# thres3 (ThresholdROC)
SD <- subset(chemo, resp == "Stable-disease")$SUV
PR <- subset(chemo, resp == "Partial-response")$SUV
CR <- subset(chemo, resp == "Complete-response")$SUV

thres3 <- thres3(SD, PR, CR, rho = prop.table(table(chemo$resp)), start = c(mean(SD), mean(CR)), ci.method = "param")
thres3

pdf("plot_thres3.pdf", height = 7, width = 9)
plot(thres3, legend = F, bw = rep(10, 3), col = 1:4, main = "", xlab = "SUV change (%)",
     xlim = c(-10, 130), ylim = c(0, 0.03), lwd = 3)
legend("topleft", c("Stable d.", "Partial r.", "Complete r.", "Thresholds", "CIs"),
       col = c(1:4, 4), lty = c(1, 1, 1, 1, 2), lwd = 3)
dev.off()




# right-censored outcome example: NSCLC
data(NSCLC)

# catpredi.survival (CatPredi)
res <- catpredi.survival(Surv(OS, ST)~1, cat.var = "COL", cat.points = 1,
                         data = NSCLC)
res
summary(res)
plot(res)

# diagnostic_assessment_continuous (ThresholdROCsurvival)
with(NSCLC, diagnostic_assessment_continuous(log(COL), OS,
                                             ST, 1095, method = "ICT", method.thres = "normal",
                                             var.equal = FALSE, m = 50, ci = TRUE))
