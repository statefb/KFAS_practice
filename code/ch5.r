libor <- matrix(scan("libor.dat"),ncol=2,byrow=T) # 金利データの読み込み

# サンプルサイズN, 離散近似の分割数J および分割区間長Dt=1/J
N <- 10^6; J <- 10; Dt <- 1/J

# libor は1,…,n 週目における満期tau[1],…,tau[p] の金利をもつn × p 行列
n <- nrow(libor); p <- ncol(libor); tau <- c(3,6)/12 # 満期は3 ヶ月と6 ヶ月

# 初期分布からのサンプリング
set.seed(1)
alp <- runif(N)*0.001; kap <- runif(N) ; the <- runif(N)*0.25
sig <- runif(N)*0.25 ; lam <- runif(N)*-1 ; sigeps <- runif(N)*0.0001

# CIR モデルによる理論式の係数
gam <- sqrt((kap+lam)^2+2*sig^2)
Atau <- Btau <- matrix(0,N,p)
for(i in 1:p){
  Atau[,i] <- 2*kap*the/sig^2*log(2*gam*exp((gam+kap+lam)/2*tau[i])/
              ((gam+kap+lam)*(exp(gam*tau[i])-1)+2*gam))
  Btau[,i] <- 2*(exp(gam*tau[i])-1)/((gam+kap+lam)*(exp(gam*tau[i])-1)+2*gam)
}
alpPsi <- cbind(alp,kap,the,sig,lam,sigeps,Atau,Btau) # 拡大状態ベクトル

signal <- matrix(0,p,n) # 描画用に各時点の金利理論値のフィルタ化推定値を保存

# 粒子フィルタ
for(s in 1:n){
  # 尤度による重みw の算出
  sigeps <- alpPsi[,6]; Atau <- alpPsi[,7:8]; Btau <- alpPsi[,9:10]
  w <- ifelse(alp>0, 1, 0) # 瞬時的短期金利が0 以下であるものは重み0 とする
  for(i in 1:p) w <- w*dnorm(libor[s,i],(-Atau[,i]+Btau[,i]*alp)/tau[i], sigeps)
  # フィルタ化推定値（金利理論値の重み付き平均）の算出
  for(i in 1:p) signal[i,s] <- sum(w*(-Atau[,i]+Btau[,i]*alp)/tau[i])/sum(w)
  # インポータンス・リサンプリング
  r <- rank(c((1:N-runif(1))/N,cumsum(w)/sum(w)),ties="random")[1:N]-1:N+1
  alpPsi <- alpPsi[r,]
  # 1 期先予測サンプリング
  alp <- alpPsi[,1]; kap <- alpPsi[,2]; the <- alpPsi[,3];
  sig <- alpPsi[,4]
  for(j in 1:J) alp <- alp+kap*(the-alp)*Dt+sig*sqrt(alp)*rnorm(N,0,sqrt(Dt))
  alpPsi[,1] <- alp <- ifelse(is.nan(alp), 0, alp) # NaN は0 におきかえる
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

