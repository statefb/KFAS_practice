library(KFAS) # R �R���\�[���ւ̃p�b�P�[�W�̓ǂݍ���(R �̋N�����ƂɕK�v)

### ��3.1 ###

Weight <- ts(scan("Weight.dat")) # �̏d�f�[�^�̓ǂݍ���

# ���[�J�����x�����f��
modLocallevel <- SSModel(Weight ~ SSMtrend(1, Q = NA), H = NA)
fitLocallevel <- fitSSM(modLocallevel, numeric(2), method = "BFGS")
kfsLocallevel <- KFS(fitLocallevel$model)

# 2 ���̃g�����h���f��
modTrend <- SSModel(Weight ~ SSMtrend(2, Q = c(list(0),list(NA))), H = NA)
fitTrend <- fitSSM(modTrend, numeric(2), method = "BFGS")
kfsTrend <- KFS(fitTrend$model)

# ���[�J�����`�g�����h���f��
modLocaltrend <- SSModel(Weight ~ SSMtrend(2, Q = c(list(NA),list(NA))), H = NA)
fitLocaltrend <- fitSSM(modLocaltrend, numeric(3), method = "BFGS")
kfsLocaltrend <- KFS(fitLocaltrend$model)

# ���������ƌX�������̕�������Ԃ̕`��i�}3.1�j
plot(Weight, lty = 3, type = "o", ylab = "��������")
lines(kfsLocallevel$alphahat[,"level"], lwd = 2, col = 8)
lines(kfsTrend$alphahat[,"level"], lwd = 2)
lines(kfsLocaltrend$alphahat[,"level"], lwd = 2, lty = 2)
plot(kfsTrend$alphahat[,"slope"], lwd = 2, ylab="�X������")
lines(kfsLocaltrend$alphahat[,"slope"], lwd = 2, lty = 2)

# �ő�ΐ��ޓx
likLocallevel <- kfsLocallevel$logLik - sum(kfsLocallevel$Finf>0) * log(2*pi)/2
likTrend <- kfsTrend$logLik - sum(kfsTrend$Finf>0) * log(2*pi)/2
likLocaltrend <- kfsLocaltrend$logLik - sum(kfsLocaltrend$Finf>0) * log(2*pi)/2

# AIC (�Ԓr���ʋK��)
aicLocallevel <- -2*likLocallevel + 2*(2+1)
aicTrend <- -2*likTrend + 2*(2+2)
aicLocaltrend <- -2*likLocaltrend + 2*(3+2)

# 1 ����\���̕��ϓ��덷
mseLocallevel <- sum(kfsLocallevel$v[3:60]^2) / 58
mseTrend <- sum(kfsTrend$v[3:60]^2) / 58
mseLocaltrend <- sum(kfsLocaltrend$v[3:60]^2) / 58

# �}3.1�̕`��
par(mfrow=c(2,1))   # �`��̈���Q����
par(mar=c(2,4,1,1)) # �`��̈�̗]���ݒ�
plot(Weight,type="o",lty=3,xlab="",ylab="��������")
lines(kfsLocallevel$alphahat[,"level"],lwd=2,col=8)
lines(kfsTrend$alphahat[,"level"],lwd=2)
lines(kfsLocaltrend$alphahat[,"level"],lwd=2,lty=5)

plot(kfsTrend$alphahat[,"slope"],type="l",lwd=2,xlab="",ylab="�X������")
lines(kfsLocaltrend$alphahat[,"slope"],lwd=2,lty=5)

### ��3.2 ###

sales <- read.csv("sales.csv")

#�_�~�[�ϐ��^(�Œ�) �̋G�ߐ������f��(�G�ߕϓ����Œ�̏ꍇ)
modSeasDummy0 <- SSModel(sales$Fabric ~ SSMtrend(2, Q = c(list(0), list(NA)))
  + SSMseasonal(12, sea.type="dummy"), H = NA)
fitSeasDummy0 <- fitSSM(modSeasDummy0, numeric(2), method = "BFGS")
kfsSeasDummy0 <- KFS(fitSeasDummy0$model)

#�_�~�[�ϐ��^(�ω�) �̋G�ߐ������f��(�G�ߕϓ����ω�����ꍇ)
modSeasDummy <- SSModel(sales$Fabric ~ SSMtrend(2, Q = c(list(0), list(NA)))
  + SSMseasonal(12, sea.type="dummy", Q = NA), H = NA)
fitSeasDummy <- fitSSM(modSeasDummy, numeric(3), method = "BFGS")
kfsSeasDummy <- KFS(fitSeasDummy$model)

#�O�p�֐��^(�Œ�) �̋G�ߐ������f��(�G�ߕϓ����Œ�̏ꍇ)
modSeasTri0 <- SSModel(sales$Fabric ~ SSMtrend(2, Q = c(list(0), list(NA)))
  + SSMseasonal(12, sea.type="trigonometric"), H = NA)
fitSeasTri0 <- fitSSM(modSeasTri0, numeric(2), method = "BFGS")
kfsSeasTri0 <- KFS(fitSeasTri0$model)

#�O�p�֐��^(�ω�) �̋G�ߐ������f��(�G�ߕϓ����ω�����ꍇ)
modSeasTri <- SSModel(sales$Fabric ~ SSMtrend(2, Q = c(list(0), list(NA)))
  + SSMseasonal(12, sea.type="trigonometric", Q = NA), H = NA)
updatefn <- function(pars, model){
model$H[] <- exp(pars[1])
diag(model$Q[,,1]) <- c(0,exp(pars[2]),rep(exp(pars[3:8]),c(rep(2,5),1)))
return(model)
}
fitSeasTri <- fitSSM(modSeasTri, c(6,0,1,2,0,0,0,0), updatefn, method="BFGS")
kfsSeasTri <- KFS(fitSeasTri$model)

# �ő�ΐ��ޓx
likSeasDummy0 <- kfsSeasDummy0$logLik - sum(kfsSeasDummy0$Finf>0) * log(2*pi)/2
likSeasDummy <- kfsSeasDummy$logLik - sum(kfsSeasDummy$Finf>0) * log(2*pi)/2
likSeasTri0 <- kfsSeasTri0$logLik - sum(kfsSeasTri0$Finf>0) * log(2*pi)/2
likSeasTri <- kfsSeasTri$logLik - sum(kfsSeasTri$Finf>0) * log(2*pi)/2

# 1 ����\���̕��ϓ��덷
mseSeasDummy0 <- sum(kfsSeasDummy0$v[14:144]^2) / 131
mseSeasDummy <- sum(kfsSeasDummy$v[14:144]^2) / 131
mseSeasTri0 <- sum(kfsSeasTri0$v[14:144]^2) / 131
mseSeasTri <- sum(kfsSeasTri$v[14:144]^2) / 131

# �U���ΐ��ޓx�̏C���i�G�ߕϓ����Œ肵�������ȃ��f���̖ޓx�𑵂���j
likSeasTri <- likSeasTri - (likSeasTri0 - likSeasDummy0)
likSeasTri0 <- likSeasDummy0

# AIC (�Ԓr���ʋK��)
aicSeasDummy0 <- -2*likSeasDummy0 + 2*(2+13)
aicSeasDummy <- -2*likSeasDummy + 2*(3+13)
aicSeasTri0 <- -2*likSeasTri0 + 2*(2+13)
aicSeasTri <- -2*likSeasTri + 2*(8+13)

# �}3.2�̕`��
mod1<-SSModel(c(sales$Fabric[1:120],rep(NA,24)) ~ SSMtrend(1,Q=NA) + SSMseasonal(12,Q=0 ), H=NA)
mod2<-SSModel(c(sales$Fabric[1:120],rep(NA,24)) ~ SSMtrend(1,Q=NA) + SSMseasonal(12,Q=NA), H=NA)
mod3<-SSModel(c(sales$Fabric[1:120],rep(NA,24)) ~ SSMtrend(2,Q=list(0,NA)) + SSMseasonal(12,Q=0 ), H=NA)
mod4<-SSModel(c(sales$Fabric[1:120],rep(NA,24)) ~ SSMtrend(2,Q=list(0,NA)) + SSMseasonal(12,Q=NA), H=NA)

fit1 <- fitSSM(mod1, numeric(2), method = "BFGS")
fit2 <- fitSSM(mod2, numeric(3), method = "BFGS")
fit3 <- fitSSM(mod3, numeric(2), method = "BFGS")
fit4 <- fitSSM(mod4, numeric(3), method = "BFGS")

kfs1 <- KFS(fit1$model)
kfs2 <- KFS(fit2$model)
kfs3 <- KFS(fit3$model)
kfs4 <- KFS(fit4$model) 

logLik1 <- kfs1$logLik - sum(kfs1$Finf>0) * log(2*pi)/2
logLik2 <- kfs2$logLik - sum(kfs2$Finf>0) * log(2*pi)/2
logLik3 <- kfs3$logLik - sum(kfs3$Finf>0) * log(2*pi)/2
logLik4 <- kfs4$logLik - sum(kfs4$Finf>0) * log(2*pi)/2

AIC1 <- -2*logLik1 + 2*( 2 + 12 )
AIC2 <- -2*logLik2 + 2*( 3 + 12 )
AIC3 <- -2*logLik3 + 2*( 2 + 13 )
AIC4 <- -2*logLik4 + 2*( 8 + 13 )

par(mfrow=c(3,1))
par(ps=16)
par(mar=c(3.5,3.5,1,1))
par(mgp=c(2,0.5,0))
plot(sales$Fabric, type="l", lty=1, ylab = "�̔��z�i10���~�j",xaxt="n",xaxs="i",col=1,xlab="(a) ��������")
axis(side=1,at=1+0:11*12,labels=c("02/1","03/1","04/1","05/1","06/1","07/1","08/1","09/1","10/1","11/1","12/1","13/1"))
lines(kfs2$alphahat[,"level"], col=3)
lines(kfs4$alphahat[,"level"], col=4)
abline(v=120.5,lty=3)

plot(kfs2$alphahat[,"sea_dummy1"],type="l", ylab = "�̔��z�i10���~�j",xaxt="n",xaxs="i",yaxs="i",col=3,xlab="(b) �G�ߐ���")
lines(kfs4$alphahat[,"sea_dummy1"],col=4)
axis(side=1,at=1+0:11*12,labels=c("02/1","03/1","04/1","05/1","06/1","07/1","08/1","09/1","10/1","11/1","12/1","13/1"))
abline(v=120.5,lty=3)

plot(sales$Fabric-kfs2$muhat,type="l", ylab="�̔��z�i10���~�j",�@xaxt="n",xaxs="i",yaxs="i",col=3,ylim=c(-150,150),xlab="(c) �������ϑ��l�h�����ƒ����\���덷")
lines(sales$Fabric-kfs4$muhat,col=4)
axis(side=1,at=1+0:11*12,labels=c("02/1","03/1","04/1","05/1","06/1","07/1","08/1","09/1","10/1","11/1","12/1","13/1"))
abline(v=120.5,lty=3)
abline(h=0,col=8)


### ��3.3 ###

# AR(1) �������f��
modAR1 <- SSModel(sales$Fabric ~ SSMtrend(2, Q = c(list(0), list(NA)))
  + SSMarima(ar = 0, Q = 0)
  + SSMseasonal(12, sea.type="dummy"), H = NA)
updatefn <- function(pars, model){
model <- SSModel(sales$Fabric ~
  SSMtrend(2, Q = c(list(0), list(exp(pars[1]))))
    + SSMarima(ar = artransform(pars[2]), Q = exp(pars[3]))
    + SSMseasonal(12, sea.type="dummy"), H = exp(pars[4]))
  return(model)
}
fitAR1 <- fitSSM(modAR1, c(-1,0,6,3), updatefn, method = "BFGS")
kfsAR1 <- KFS(fitAR1$model)

# AR(2) �������f��
modAR2 <- SSModel(sales$Fabric ~ SSMtrend(2, Q = c(list(0), list(NA)))
  + SSMarima(ar = c(0, 0), Q = 0)
  + SSMseasonal(12, sea.type="dummy"), H = NA)
updatefn <- function(pars, model){
  model <- SSModel(sales$Fabric ~ SSMtrend(2, Q = c(list(0), list(exp(pars[1]))))
    + SSMarima(ar = artransform(pars[2:3]), Q = exp(pars[4]))
    + SSMseasonal(12, sea.type="dummy"), H = exp(pars[5]))
  return(model)
}
fitAR2 <- fitSSM(modAR2, c(-1,0.1,0,6,3), updatefn, method = "BFGS")
kfsAR2 <- KFS(fitAR2$model)

# AR(12) �������f��(���O12 �ȊO�̎��ȉ�A�W���̓[���Ƃ���)
modAR12 <- SSModel(sales$Fabric ~ SSMtrend(2, Q = c(list(0), list(NA)))
  + SSMarima(ar = rep(0, 12), Q = 0)
  + SSMseasonal(12, sea.type="dummy"), H = NA)
updatefn <- function(pars, model){
  model <- SSModel(sales$Fabric ~ SSMtrend(2, Q = c(list(0), list(exp(pars[1]))))
    + SSMarima(ar = c(rep(0,11), artransform(pars[2])), Q = exp(pars[3]))
    + SSMseasonal(12, sea.type="dummy"), H = exp(pars[4]))
return(model)
}
fitAR12 <- fitSSM(modAR12, c(-1,0.4,6,0), updatefn, method = "BFGS")


### ��3.4 ###

Gasoline <- ts(scan("Gasoline.dat")) # �K�\�����P���f�[�^�̓ǂݍ���

# ��A�W�������ԕω�������ꍇ
modRegression <- SSModel(log(sales$Fuel) ~ SSMtrend(1, Q = NA)
  + SSMseasonal(12, sea.type="dummy")
  + SSMregression(~ log(Gasoline), Q = NA), H = NA)
fitRegression <- fitSSM(modRegression, numeric(3), method = "BFGS")
kfsRegression <- KFS(fitRegression$model)

# ��A�W�����Œ肷��ꍇ
modRegression0 <- SSModel(log(sales$Fuel) ~ SSMtrend(1, Q = NA)
  + SSMseasonal(12, sea.type="dummy")
  + log(Gasoline), H = NA)
fitRegression0 <- fitSSM(modRegression0, numeric(2), method = "BFGS")
kfsRegression0 <- KFS(fitRegression0$model)

### ��蒼�����R�[�h
# Reg0:��A�����Ȃ�
# Reg1:��A�W�����Œ肷��ꍇ
# Reg2:��A�W�������ԕω�������ꍇ

fuel <- sales[,4]                # �R�������Ƃ̌����̔��z
gas  <- ts(scan("Gasoline.dat")) # �K�\�����P���f�[�^�̓ǂݍ���

# ���f����`
modReg0 <- SSModel(log(fuel) ~ SSMtrend(1,Q=NA) + SSMseasonal(12,Q=0), H=NA)

modReg1 <- SSModel(log(fuel) ~ log(gas) + 
  SSMtrend(1,Q=NA) + SSMseasonal(12,Q=0), H=NA)

modReg2 <- SSModel(log(fuel) ~ SSMregression(~ log(gas), Q=NA) +
  SSMtrend(1,Q=NA) + SSMseasonal(12,Q=0), H=NA)

# ���m�p�����[�^�̐���
fitReg0 <- fitSSM(modReg0, rep(-8,2), method = "BFGS")
fitReg1 <- fitSSM(modReg1, rep(-8,2), method = "BFGS")
fitReg2 <- fitSSM(modReg2, rep(-8,3), method = "BFGS")

# �J���}���t�B���^�E�J���}���X���[�U�̎��s
kfsReg0 <- KFS(fitReg0$model)
kfsReg1 <- KFS(fitReg1$model)
kfsReg2 <- KFS(fitReg2$model)

# �ő�ΐ��ޓx
logLikReg0 <- kfsReg0$logLik - sum(kfsReg0$Finf>0) * log(2*pi)/2
logLikReg1 <- kfsReg1$logLik - sum(kfsReg1$Finf>0) * log(2*pi)/2
logLikReg2 <- kfsReg2$logLik - sum(kfsReg2$Finf>0) * log(2*pi)/2

# AIC (�Ԓr���ʋK��)
AICReg0 <- -2*logLikReg0 + 2*( 2 + 12 )
AICReg1 <- -2*logLikReg1 + 2*( 2 + 13 )
AICReg2 <- -2*logLikReg2 + 2*( 3 + 13 )

# ���������ꂽ���������C��A�����̕`��
par(mfrow=c(3,1))   # �`��̈���R����
par(mar=c(2,4,1,1)) # �`��̈�̗]���ݒ�
plot(log(fuel),type="o",lty=3,xaxt="n",xlab="",ylab="�̔��z�i�ΐ��j")
lines(kfsReg1$alphahat[,"level"]+kfsReg1$alphahat[,"log(gas)"]*log(gas),lwd=2)
axis(side=1,at=1+0:11*12,labels=c("02/1","03/1","04/1","05/1","06/1","07/1","08/1","09/1","10/1","11/1","12/1","13/1"))
plot(kfsReg1$alphahat[,"level"],xaxt="n",xlab="",ylab="��������")
axis(side=1,at=1+0:11*12,labels=c("02/1","03/1","04/1","05/1","06/1","07/1","08/1","09/1","10/1","11/1","12/1","13/1"))
plot(kfsReg1$alphahat[,"log(gas)"]*log(gas),ylim=c(1.4,1.6),xaxt="n",xlab="",ylab="��A����")
axis(side=1,at=1+0:11*12,labels=c("02/1","03/1","04/1","05/1","06/1","07/1","08/1","09/1","10/1","11/1","12/1","13/1"))


### ��3.5 ###

# �e���̗j���W�v
dates <- seq(as.Date("2002-01-01"), as.Date("2013-12-31"), by = 1)
weeks <- table(substr(dates,1,7), weekdays(dates, T))
sun <- weeks[,"��"]
mon <- weeks[,"��"]-sun; tue <- weeks[,"��"]-sun; wed <- weeks[,"��"]-sun
thu <- weeks[,"��"]-sun; fry <- weeks[,"��"]-sun; sat <- weeks[,"�y"]-sun
calendar <- cbind(mon, tue, wed, thu, fry, sat)

# ���邤�N�Q���̃_�~�[�ϐ�
leapyear <- rownames(weeks) %in% c("2004-02","2008-02","2012-02")

# �J�����_�[����(�j���E���邤�N) �̂��郂�f��
modCalender <- SSModel(sales$Fabric ~ SSMtrend(2, Q = c(list(0), list(NA)))
  + SSMseasonal(12, sea.type="dummy")
  + leapyear + calendar, H = NA)
fitCalender <- fitSSM(modCalender, numeric(2), method = "BFGS")
kfsCalender <- KFS(fitCalender$model)

# �}3.7�̕`��
plot(kfsCalender$muhat - kfsCalender$alphahat[,"level"],type="l",xaxs="i",xaxt="n",xlab="",ylab="�̔��z�i10���~�j")
axis(side=1,at=1+0:11*12,labels=c("02/1","03/1","04/1","05/1","06/1","07/1","08/1","09/1","10/1","11/1","12/1","13/1"))

### ��3.6 ###

# 2010�N11���C12���̃f�[�^�����O
salesNA <- sales$Machinery
salesNA[sales$month %in% c("2010�N11��","2010�N12��")] <- NA

# 2011�N8���ȍ~�̐����V�t�g���ϐ��̒�`
ShiftLevel <- (1:nrow(sales) >= which(sales$month=="2011�N8��"))

# �����V�t�g���ϐ������������f��
modShift <- SSModel(log(salesNA) ~ SSMtrend(1, Q = NA)
  + SSMseasonal(12, Q = NA, sea.type="dummy")
  + ShiftLevel, H = NA)
fitShift <- fitSSM(modShift, numeric(3))
kfsShift <- KFS(fitShift$model, smoothing=c("state","mean","disturbance"))

# �⏕�c���̃v���b�g��
plot(rstandard(kfsShift, "pearson")) # �ϑ��l�h����
plot(rstandard(kfsShift, "state")[,1]) # ��Ԋh����(��������)


# �o�[�W�����X�V�ŃG���[���o���ꍇ�͈ȉ������s��������
# �@�B�����Ɣ̔��z
machine <- sales[,3]

# �����V�t�g���ϐ������������f��
mod0 <- SSModel(log(machine) ~ SSMtrend(1, Q = NA)
  + SSMseasonal(12, Q = NA, sea.type="dummy"), H = NA)
fit0 <- fitSSM(mod0, numeric(3))
kfs0 <- KFS(fit0$model, smoothing=c("state","mean","disturbance"))

mod01 <- SSModel(log(machine) ~ SSMtrend(2, Q = list(NA,NA)), H = NA)
fit01 <- fitSSM(mod01, numeric(3))
kfs01 <- KFS(fit01$model, smoothing=c("state","mean","disturbance"))


# ��ԝ������̕W�����c���ikfs0$rstandard_state�j�̎Z�o
kfs0$V_eta_inv <- array(apply(kfs0$V_eta,3,solve),dim(kfs0$V_eta))
kfs0$V_eta_inv_chol <- array(apply(kfs0$V_eta_inv,3,chol),dim(kfs0$V_eta))
kfs0$rstandard_state <- kfs0$etahat*0
for(i in 1:nrow(kfs0$etahat)) kfs0$rstandard_state[i,] = kfs0$V_eta_inv_chol[,,i]%*%kfs0$etahat[i,]


# �}3.9�̕`��
par(mfrow=c(3,1))   # �`��̈���R����
par(mar=c(3,3,1,1)) # �`��̈�̗]���ݒ�
par(mgp=c(2,0.5,0))
plot(log(machine),type="l",xaxt="n",xlab="(a) �̔��z�̑ΐ��n��ƕ�������Ԃ̐�������",ylab="�̔��z�i�ΐ��j")
lines(kfs0$alphahat[,"level"],col=8)
axis(side=1,at=1+0:11*12,labels=c("02/1","03/1","04/1","05/1","06/1","07/1","08/1","09/1","10/1","11/1","12/1","13/1"))
plot(rstandard(kfs0, "pearson"),ylim=c(-6,6),xaxt="n",xlab="(b) �ϑ��l�������̕⏕�c��",ylab="�̔��z�i�ΐ��j")
abline(h=c(-1.96,1.96),lty=3)
axis(side=1,at=1+0:11*12,labels=c("02/1","03/1","04/1","05/1","06/1","07/1","08/1","09/1","10/1","11/1","12/1","13/1"))
plot(kfs0$rstandard_state[,1],ylim=c(-6,6),xaxt="n",xlab="(c) ��ԝ������i���������j�̕⏕�c��",ylab="�̔��z�i�ΐ��j")
abline(h=c(-1.96,1.96),lty=3)
axis(side=1,at=1+0:11*12,labels=c("02/1","03/1","04/1","05/1","06/1","07/1","08/1","09/1","10/1","11/1","12/1","13/1"))

# 2010 �N11 ���C12 ���̃f�[�^�����O
machineNA <- machine
machineNA[sales[,1] %in% c("2010�N11��","2010�N12��")] <- NA

# 2011 �N8 ���ȍ~�̐����V�t�g���ϐ��̒�`
ShiftLevel <- (1:nrow(sales) >= which(sales[,1]=="2011�N8��"))

# �����V�t�g���ϐ������������f��
modShift <- SSModel(log(machineNA) ~ SSMtrend(1, Q = NA)
  + SSMseasonal(12, Q = NA, sea.type="dummy")
  + ShiftLevel, H = NA)
fitShift <- fitSSM(modShift, numeric(3))
kfsShift <- KFS(fitShift$model, smoothing=c("state","mean","disturbance"))

# ��ԝ������̕W�����c���ikfsShift$rstandard_state�j�̎Z�o
kfsShift$V_eta_inv <- array(apply(kfsShift$V_eta,3,solve),dim(kfsShift$V_eta))
kfsShift$V_eta_inv_chol <- array(apply(kfsShift$V_eta_inv,3,chol),dim(kfsShift$V_eta))
kfsShift$rstandard_state <- kfsShift$etahat*0
for(i in 1:nrow(kfsShift$etahat)) kfsShift$rstandard_state[i,] = kfsShift$V_eta_inv_chol[,,i]%*%kfsShift$etahat[i,]

# �}3.10�̕`��
par(mfrow=c(3,1))   # �`��̈���R����
par(mar=c(3,3,1,1)) # �`��̈�̗]���ݒ�
par(mgp=c(2,0.5,0))
plot(log(machine),type="l",xaxt="n",xlab="(a) �̔��z�̑ΐ��n��ƕ�������Ԃ̐�������",ylab="�̔��z�i�ΐ��j")
lines(kfsShift$muhat - kfsShift$alphahat[,"sea_dummy1"],col=8)
axis(side=1,at=1+0:11*12,labels=c("02/1","03/1","04/1","05/1","06/1","07/1","08/1","09/1","10/1","11/1","12/1","13/1"))
plot(rstandard(kfsShift, "pearson"),ylim=c(-6,6),xaxt="n",xlab="(b) �ϑ��l�������̕⏕�c��",ylab="�̔��z�i�ΐ��j")
abline(h=c(-1.96,1.96),lty=3)
axis(side=1,at=1+0:11*12,labels=c("02/1","03/1","04/1","05/1","06/1","07/1","08/1","09/1","10/1","11/1","12/1","13/1"))
plot(kfsShift$rstandard_state[,1],ylim=c(-6,6),xaxt="n",xlab="(c) ��ԝ������i���������j�̕⏕�c��",ylab="�̔��z�i�ΐ��j")
abline(h=c(-1.96,1.96),lty=3)
axis(side=1,at=1+0:11*12,labels=c("02/1","03/1","04/1","05/1","06/1","07/1","08/1","09/1","10/1","11/1","12/1","13/1"))


### ��3.7 ###

Weight <- ts(scan("Weight.dat"))   # �̏d�f�[�^�̓ǂݍ���
Bodyfat <- ts(scan("Bodyfat.dat")) # �̎��b���f�[�^�̓ǂݍ���

WeightNA <- Weight
WeightNA[21:40] <- NA

modSUTSE <- SSModel(cbind(WeightNA, Bodyfat) ~
  SSMtrend(1, Q = matrix(NA,2,2), type = "distinct"), H = matrix(NA,2,2))
fitSUTSE <- fitSSM(modSUTSE, numeric(6), method="BFGS")
kfsSUTSE <- KFS(fitSUTSE$model)
