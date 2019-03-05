library(survival)
library(survminer)
kmsurvo<-Surv(combined$Rectime,combined$Censor.Status)
sfit <- survfit(Surv(combined$Rectime,combined$Censor.Status)~., data=combined)
ggsurvplot(cox_fit, pval=TRUE,
           xlab = "Time in days",
           legend.title = " ",
           legend.labs=c("Normal", "Altered"), 
           font.main = c(16, "plain", "darkblue"),
           font.tickslab = c(14, "plain", "darkgreen"),
           font.x = c(16, "plain", "red"),
           font.y = c(16, "plain", "darkred"),
           legend=c(.85,.85),
           title="WFDC2")
ggsurvplot(sfit, legend = c(0.2, 0.2))
plot(sfit)
cox1<-aareg(Surv(combined$Rectime,combined$Censor.Status)~.,data=combined[,11:ncol(combined)]) 
cox<-coxph(kmsurvo~.,data=combined[,11:ncol(combined)])
cox_fit <- survfit(cox1)
autoplot(cox1)

autoplot(cox_fit)
ggsurv<-ggsurvplot(sfit, data=combined,  pval=TRUE,
           
           #xlim = c(0,2000), # present narrower X axis, but not affect
           # survival estimates.
           xlab = "Time in days",   # customize X axis label.
           title="BRCA1",
           title=c(1,2000),
           main = "Survival curve",
           font.main = c(16, "plain", "darkblue"),
           font.tickslab = c(14, "plain", "darkgreen"),
           legend.title = " ",
           legend.labs=c("Normal", "Altered"), 
           font.x = c(16, "plain", "red"),
           font.y = c(16, "plain", "darkred"),
           legend=c(.85,.85))

ggsurv$plot <- ggsurv$plot + 
  theme(legend.text = element_text(size = 14, color = "black", face = "plain"))
ggsurv
