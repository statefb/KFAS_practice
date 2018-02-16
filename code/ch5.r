libor <- matrix(scan("libor.dat"),ncol=2,byrow=T) # �����f�[�^�̓ǂݍ���

# �T���v���T�C�YN, ���U�ߎ��̕�����J ����ѕ�����Ԓ�Dt=1/J
N <- 10^6; J <- 10; Dt <- 1/J

# libor ��1,�c,n �T�ڂɂ����閞��tau[1],�c,tau[p] �̋���������n �~ p �s��
n <- nrow(libor); p <- ncol(libor); tau <- c(3,6)/12 # ������3 ������6 ����

# �������z����̃T���v�����O
set.seed(1)
alp <- runif(N)*0.001; kap <- runif(N) ; the <- runif(N)*0.25
sig <- runif(N)*0.25 ; lam <- runif(N)*-1 ; sigeps <- runif(N)*0.0001

# CIR ���f���ɂ�闝�_���̌W��
gam <- sqrt((kap+lam)^2+2*sig^2)
Atau <- Btau <- matrix(0,N,p)
for(i in 1:p){
  Atau[,i] <- 2*kap*the/sig^2*log(2*gam*exp((gam+kap+lam)/2*tau[i])/
              ((gam+kap+lam)*(exp(gam*tau[i])-1)+2*gam))
  Btau[,i] <- 2*(exp(gam*tau[i])-1)/((gam+kap+lam)*(exp(gam*tau[i])-1)+2*gam)
}
alpPsi <- cbind(alp,kap,the,sig,lam,sigeps,Atau,Btau) # �g���ԃx�N�g��

signal <- matrix(0,p,n) # �`��p�Ɋe���_�̋������_�l�̃t�B���^������l��ۑ�

# ���q�t�B���^
for(s in 1:n){
  # �ޓx�ɂ��d��w �̎Z�o
  sigeps <- alpPsi[,6]; Atau <- alpPsi[,7:8]; Btau <- alpPsi[,9:10]
  w <- ifelse(alp>0, 1, 0) # �u���I�Z��������0 �ȉ��ł�����̂͏d��0 �Ƃ���
  for(i in 1:p) w <- w*dnorm(libor[s,i],(-Atau[,i]+Btau[,i]*alp)/tau[i], sigeps)
  # �t�B���^������l�i�������_�l�̏d�ݕt�����ρj�̎Z�o
  for(i in 1:p) signal[i,s] <- sum(w*(-Atau[,i]+Btau[,i]*alp)/tau[i])/sum(w)
  # �C���|�[�^���X�E���T���v�����O
  r <- rank(c((1:N-runif(1))/N,cumsum(w)/sum(w)),ties="random")[1:N]-1:N+1
  alpPsi <- alpPsi[r,]
  # 1 ����\���T���v�����O
  alp <- alpPsi[,1]; kap <- alpPsi[,2]; the <- alpPsi[,3];
  sig <- alpPsi[,4]
  for(j in 1:J) alp <- alp+kap*(the-alp)*Dt+sig*sqrt(alp)*rnorm(N,0,sqrt(Dt))
  alpPsi[,1] <- alp <- ifelse(is.nan(alp), 0, alp) # NaN ��0 �ɂ���������
}

par(mar=c(3,3,1,1))
par(mgp=c(2,.5,0))
plot(libor[,1]*100,type="l",ylim=c(0.07,0.16),xlab="Week",ylab="Interest (%)",xaxs="i",col=4)
legend(30,0.16,lty=1,col=c(2,4),legend=c("6m JPY LIBOR","3m JPY LIBOR"))
lines(signal[1,]*100,col=5,lwd=1)
lines(signal[2,]*100,col=6,lwd=1)
lines(libor[,2]*100,col=2)
lines(libor[,1]*100,col=4)

table(alpPsi[,2]) # kappa
table(alpPsi[,3]) # theta
table(alpPsi[,4]) # sigma
table(alpPsi[,5]) # lambda
table(alpPsi[,6]) # sigeps

