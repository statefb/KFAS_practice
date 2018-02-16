# install.packages("KFAS") # �~���[�T�[�o���擾���ăC���X�g�[��(1 �x����)
library(KFAS) # R �R���\�[���ւ̃p�b�P�[�W�̓ǂݍ���(R �̋N�����ƂɕK�v)

Weight <- ts(scan("Weight.dat")) # �̏d�f�[�^�̓ǂݍ���

# ���[�J�����x�����f���̒�`
mod <- SSModel(Weight ~ SSMtrend(1, Q = NA), H = NA)

# ���m�p�����[�^�̐���
fit <- fitSSM(mod, numeric(2), method = "BFGS")

# �J���}���t�B���^�ƕ������̎��s
kfs <- KFS(fit$model)

# �t�B���^������ʂƂ��̐M�����
afilt <- kfs$a[-1]
Pfilt <- kfs$P[,,-1] - fit$model$Q
afiltconf <- cbind(afilt+sqrt(Pfilt)*qnorm(0.025),afilt+sqrt(Pfilt)*qnorm(0.975))

# �}2.2�̕`��
plot(Weight,type="o",lty=3,xlab="�o�ߓ���",ylab="�̏d�ikg�j",ylim=c(83,87))
lines(afilt,lwd=2)
lines(afiltconf[,1])
lines(afiltconf[,2])

# ��������Ԃ̐M�����(kfs$alphahat+sqrt(kfs$V)*pnorm(0.025) �Ȃǂł��悢)
alphahatconf <- predict(fit$model, interval = "confidence", level = 0.95)

# �}2.3�̕`��
plot(Weight,type="o",lty=3,xlab="�o�ߓ���",ylab="�̏d�ikg�j",ylim=c(83,87))
lines(alphahatconf[,1],lwd=2)
lines(alphahatconf[,2])
lines(alphahatconf[,3])

# �����l�̕��
modNA <- SSModel(Weight[c(1:20,rep(NA,20),41:60)] ~ SSMtrend(1, Q = NA), H = NA)
fitNA <- fitSSM(modNA, numeric(2), method = "BFGS")
preNA <- predict(fitNA$model, interval = "prediction", level = 0.95)
confNA <- predict(fitNA$model, interval = "confidence", level = 0.95)

# �}2.4�̕`��
plot(Weight,type="o",lty=3,xlab="�o�ߓ���",ylab="�̏d�ikg�j",ylim=c(83,87))
lines(21:40,Weight[21:40],type="o",lty=3,col=8)
lines(confNA[,1],lwd=2)
lines(confNA[,2])
lines(confNA[,3])
lines(21:40,preNA[21:40,2],lty=2)
lines(21:40,preNA[21:40,3],lty=2)

# �����\��
mod50 <- SSModel(Weight[1:50] ~ SSMtrend(1, Q = NA), H = NA)
fit50 <- fitSSM(mod50, numeric(2), method = "BFGS")
pre50 <- predict(fit50$model, interval ="prediction", n.ahead = 10,
level = 0.95)
conf50 <- predict(fit50$model, interval ="confidence", n.ahead = 10)

# �}2.5�̕`��
plot(Weight,type="o",lty=3,xlab="�o�ߓ���",ylab="�̏d�ikg�j",ylim=c(83,87))
lines(51:60,Weight[51:60],type="o",lty=3,col=8)
lines(51:60,conf50[,1],lwd=2)
lines(51:60,conf50[,2])
lines(51:60,conf50[,3])
lines(51:60,pre50[,2],lty=2)
lines(51:60,pre50[,3],lty=2)

# �����\���i�����lNA�Ƃ��ė\����������j
mod50NA <- SSModel(Weight[c(1:50,rep(NA,10))] ~ SSMtrend(1, Q = NA), H = NA)
fit50NA <- fitSSM(mod50NA, numeric(2), method = "BFGS")
conf50NA <- predict(fit50NA$model, interval="confidence", level=0.95)
pre50NA  <- predict(fit50NA$model, interval="prediction", level=0.95)

