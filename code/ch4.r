library(KFAS); library(Nippon)

fire     <- ts(scan("fire.dat"))     # �΍Ќ����f�[�^�̓ǂݍ���
humidity <- ts(scan("humidity.dat")) # �����ώ��x�f�[�^�̓ǂݍ���
celsius  <- ts(scan("celsius.dat"))  # �����ϋC���f�[�^�̓ǂݍ���
rain     <- ts(scan("rain.dat"))     # ���~���ʃf�[�^�̓ǂݍ���

dates <- seq(as.Date("2005-01-01"), as.Date("2014-12-31"), by=1) # ���t��쐬
weekday <- weekdays(dates, T) # �j�����ʊ֐�("��"�`"��"��Ԃ�)
weekday[is.jholiday(dates)] <- "�j" # �p�b�P�[�WNippon �̏j�����ʊ֐�
weekday <- factor(weekday, c("��","��","��","��","��","�y","��","�j"))

# ���q��
modPois <- SSModel(fire ~ SSMtrend(2, Q = c(list(0), list(NA)))
  + SSMcycle(365.25) + SSMcycle(365.25/2) + SSMcycle(365.25/3)
  + SSMcycle(365.25/4) + SSMcycle(365.25/5) + SSMcycle(365.25/6)
  + weekday + humidity + celsius + rain,
distribution = "poisson")
diag(modPois$P1inf) <- 0 # �U������������肭�����Ȃ��������߁C
diag(modPois$P1) <- 100^2 # ����ɏ\���傫�ȏ������U��p����

# ��K�E�X���f���ł̓V�~�����[�V������nsim ��ݒ肷��
# �i���ԁE�������[�̐���ŏ�肭�����Ȃ��ꍇ�Ansim�����炵�Ď��s���Ă��������j
fitPois <- fitSSM(modPois, 0, method="Brent", lower=-40, upper=0, nsim=1000)
kfsPois <- KFS(fitPois$model, nsim=1000)

# �C���|�[�^���X�E�T���v�����O�Ɨ��p��(���x����(��8 ��Ԑ���) ��95%�M�����)
impPois <- importanceSSM(fitPois$model, type="states", nsim=4000)
emp <- cumsum(impPois$weights[order(impPois$samples[1,8,])])/
sum(impPois$weight)
conf <- sort(impPois$samples[1,8,])[rank(c(0.025, 0.975, emp))[1:2]]

# �ϑ��l��98���\�����
prePois <- predict(fitPois$model, interval="prediction", level=0.98, nsim=10000)
# ��͌��ʂ̕\��

kfsPois # ��A�W���̐���l��W���덷��\��

# ��A�W����95%�M����Ԃ�\��
sort(impPois$samples[1,1,])[rank(c(0.025, 0.975, emp))[1:2]] # �Ηj
sort(impPois$samples[1,2,])[rank(c(0.025, 0.975, emp))[1:2]] # ���j
sort(impPois$samples[1,3,])[rank(c(0.025, 0.975, emp))[1:2]] # �ؗj
sort(impPois$samples[1,4,])[rank(c(0.025, 0.975, emp))[1:2]] # ���j
sort(impPois$samples[1,5,])[rank(c(0.025, 0.975, emp))[1:2]] # �y�j
sort(impPois$samples[1,6,])[rank(c(0.025, 0.975, emp))[1:2]] # ���j
sort(impPois$samples[1,7,])[rank(c(0.025, 0.975, emp))[1:2]] # �j��
sort(impPois$samples[1,8,])[rank(c(0.025, 0.975, emp))[1:2]] # ���x
sort(impPois$samples[1,9,])[rank(c(0.025, 0.975, emp))[1:2]] # �C��
sort(impPois$samples[1,10,])[rank(c(0.025, 0.975, emp))[1:2]] # �~����

# ��������Ԃɂ��΍Ќ����̊��Ғl�A98%�\����ԂƊO��l
plot(prePois[,1], ylim=c(0,52), xaxs="i", yaxs="i", ylab="�΍Ќ���", xaxt="n")
lines(prePois[,2], col=8)
lines(prePois[,3], col=8)
points(which(fire>prePois[,3]),fire[fire>prePois[,3]])
points(which(fire<prePois[,2]),fire[fire<prePois[,2]])
axis(side=1,at=which(substr(dates,6,10)=="01-01"),
 labels=c("05/1/1","06/1/1","07/1/1","08/1/1","09/1/1","10/1/1","11/1/1","12/1/1","13/1/1","14/1/1"))

# ���̓񍀕��z���f���ɂ����

hum_cel = humidity * celsius # ���x�Ɖ��x�̌��ݍ�p
fireNA <- fire; fireNA[dates=="2011-03-11"] <- NA # �O��l�������l�ɑւ����O
modNegbin <- SSModel(fireNA ~ SSMtrend(2 , Q = c(list(0), list(NA)))
  + SSMcycle(365.25) + SSMcycle(365.25/2) + SSMcycle(365.25/3)
  + SSMcycle(365.25/4) + SSMcycle(365.25/5) + SSMcycle(365.25/6)
  + weekday + humidity + celsius + hum_cel,
  distribution = "negative binomial", u = 1)
diag(modNegbin$P1inf) <- 0; diag(modNegbin$P1) <- 100^2
updatefn <- function(pars, model){
  model$Q[2,2,] = exp(pars[1])
  model$u[,] = exp(pars[2])
  return(model)
}
fitNegbin <- fitSSM(modNegbin, c(-28,0), updatefn, method="BFGS", nsim=1000)
kfsNegbin <- KFS(fitNegbin$model, nsim=1000)

# �C���|�[�^���X�E�T���v�����O�Ɨ��p��(���x����(��8 ��Ԑ���) ��95%�M�����)
impNegbin <- importanceSSM(fitNegbin$model, type="states", nsim=1000)
emp <- cumsum(impNegbin$weights[order(impNegbin$samples[1,8,])])/
sum(impNegbin$weight)
conf <- sort(impNegbin$samples[1,8,])[rank(c(0.025, 0.975, emp))[1:2]]

# �ϑ��l��98���\�����
preNegbin <- predict(fitNegbin$model, interval="prediction", level=0.98, nsim=1000)


# ��͌��ʂ̕\��

kfsNegbin # ��A�W���̐���l��W���덷��\��

# ��A�W����95%�M����Ԃ�\��
sort(impNegbin$samples[1,1,])[rank(c(0.025, 0.975, emp))[1:2]] # �Ηj
sort(impNegbin$samples[1,2,])[rank(c(0.025, 0.975, emp))[1:2]] # ���j
sort(impNegbin$samples[1,3,])[rank(c(0.025, 0.975, emp))[1:2]] # �ؗj
sort(impNegbin$samples[1,4,])[rank(c(0.025, 0.975, emp))[1:2]] # ���j
sort(impNegbin$samples[1,5,])[rank(c(0.025, 0.975, emp))[1:2]] # �y�j
sort(impNegbin$samples[1,6,])[rank(c(0.025, 0.975, emp))[1:2]] # ���j
sort(impNegbin$samples[1,7,])[rank(c(0.025, 0.975, emp))[1:2]] # �j��
sort(impNegbin$samples[1,8,])[rank(c(0.025, 0.975, emp))[1:2]] # ���x
sort(impNegbin$samples[1,9,])[rank(c(0.025, 0.975, emp))[1:2]] # �C��
sort(impNegbin$samples[1,10,])[rank(c(0.025, 0.975, emp))[1:2]] # ���x�~�C��

# ��������Ԃɂ��΍Ќ����̊��Ғl�A98%�\����ԂƊO��l
plot(preNegbin[,1], ylim=c(0,52), xaxs="i", yaxs="i", ylab="�΍Ќ���", xaxt="n")
lines(preNegbin[,2], col=8)
lines(preNegbin[,3], col=8)
points(which(fire>preNegbin[,3]),fire[fire>preNegbin[,3]])
points(which(fire<preNegbin[,2]),fire[fire<preNegbin[,2]])
axis(side=1,at=which(substr(dates,6,10)=="01-01"),
 labels=c("05/1/1","06/1/1","07/1/1","08/1/1","09/1/1","10/1/1","11/1/1","12/1/1","13/1/1","14/1/1"))

