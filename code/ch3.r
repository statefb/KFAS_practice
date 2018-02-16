library(KFAS) # R コンソールへのパッケージの読み込み(R の起動ごとに必要)

### 例3.1 ###

Weight <- ts(scan("Weight.dat")) # 体重データの読み込み

# ローカルレベルモデル
modLocallevel <- SSModel(Weight ~ SSMtrend(1, Q = NA), H = NA)
fitLocallevel <- fitSSM(modLocallevel, numeric(2), method = "BFGS")
kfsLocallevel <- KFS(fitLocallevel$model)

# 2 次のトレンドモデル
modTrend <- SSModel(Weight ~ SSMtrend(2, Q = c(list(0),list(NA))), H = NA)
fitTrend <- fitSSM(modTrend, numeric(2), method = "BFGS")
kfsTrend <- KFS(fitTrend$model)

# ローカル線形トレンドモデル
modLocaltrend <- SSModel(Weight ~ SSMtrend(2, Q = c(list(NA),list(NA))), H = NA)
fitLocaltrend <- fitSSM(modLocaltrend, numeric(3), method = "BFGS")
kfsLocaltrend <- KFS(fitLocaltrend$model)

# 水準成分と傾き成分の平滑化状態の描画（図3.1）
plot(Weight, lty = 3, type = "o", ylab = "水準成分")
lines(kfsLocallevel$alphahat[,"level"], lwd = 2, col = 8)
lines(kfsTrend$alphahat[,"level"], lwd = 2)
lines(kfsLocaltrend$alphahat[,"level"], lwd = 2, lty = 2)
plot(kfsTrend$alphahat[,"slope"], lwd = 2, ylab="傾き成分")
lines(kfsLocaltrend$alphahat[,"slope"], lwd = 2, lty = 2)

# 最大対数尤度
likLocallevel <- kfsLocallevel$logLik - sum(kfsLocallevel$Finf>0) * log(2*pi)/2
likTrend <- kfsTrend$logLik - sum(kfsTrend$Finf>0) * log(2*pi)/2
likLocaltrend <- kfsLocaltrend$logLik - sum(kfsLocaltrend$Finf>0) * log(2*pi)/2

# AIC (赤池情報量規準)
aicLocallevel <- -2*likLocallevel + 2*(2+1)
aicTrend <- -2*likTrend + 2*(2+2)
aicLocaltrend <- -2*likLocaltrend + 2*(3+2)

# 1 期先予測の平均二乗誤差
mseLocallevel <- sum(kfsLocallevel$v[3:60]^2) / 58
mseTrend <- sum(kfsTrend$v[3:60]^2) / 58
mseLocaltrend <- sum(kfsLocaltrend$v[3:60]^2) / 58

# 図3.1の描画
par(mfrow=c(2,1))   # 描画領域を２分割
par(mar=c(2,4,1,1)) # 描画領域の余白設定
plot(Weight,type="o",lty=3,xlab="",ylab="水準成分")
lines(kfsLocallevel$alphahat[,"level"],lwd=2,col=8)
lines(kfsTrend$alphahat[,"level"],lwd=2)
lines(kfsLocaltrend$alphahat[,"level"],lwd=2,lty=5)

plot(kfsTrend$alphahat[,"slope"],type="l",lwd=2,xlab="",ylab="傾き成分")
lines(kfsLocaltrend$alphahat[,"slope"],lwd=2,lty=5)

### 例3.2 ###

sales <- read.csv("sales.csv")

#ダミー変数型(固定) の季節成分モデル(季節変動が固定の場合)
modSeasDummy0 <- SSModel(sales$Fabric ~ SSMtrend(2, Q = c(list(0), list(NA)))
  + SSMseasonal(12, sea.type="dummy"), H = NA)
fitSeasDummy0 <- fitSSM(modSeasDummy0, numeric(2), method = "BFGS")
kfsSeasDummy0 <- KFS(fitSeasDummy0$model)

#ダミー変数型(変化) の季節成分モデル(季節変動が変化する場合)
modSeasDummy <- SSModel(sales$Fabric ~ SSMtrend(2, Q = c(list(0), list(NA)))
  + SSMseasonal(12, sea.type="dummy", Q = NA), H = NA)
fitSeasDummy <- fitSSM(modSeasDummy, numeric(3), method = "BFGS")
kfsSeasDummy <- KFS(fitSeasDummy$model)

#三角関数型(固定) の季節成分モデル(季節変動が固定の場合)
modSeasTri0 <- SSModel(sales$Fabric ~ SSMtrend(2, Q = c(list(0), list(NA)))
  + SSMseasonal(12, sea.type="trigonometric"), H = NA)
fitSeasTri0 <- fitSSM(modSeasTri0, numeric(2), method = "BFGS")
kfsSeasTri0 <- KFS(fitSeasTri0$model)

#三角関数型(変化) の季節成分モデル(季節変動が変化する場合)
modSeasTri <- SSModel(sales$Fabric ~ SSMtrend(2, Q = c(list(0), list(NA)))
  + SSMseasonal(12, sea.type="trigonometric", Q = NA), H = NA)
updatefn <- function(pars, model){
model$H[] <- exp(pars[1])
diag(model$Q[,,1]) <- c(0,exp(pars[2]),rep(exp(pars[3:8]),c(rep(2,5),1)))
return(model)
}
fitSeasTri <- fitSSM(modSeasTri, c(6,0,1,2,0,0,0,0), updatefn, method="BFGS")
kfsSeasTri <- KFS(fitSeasTri$model)

# 最大対数尤度
likSeasDummy0 <- kfsSeasDummy0$logLik - sum(kfsSeasDummy0$Finf>0) * log(2*pi)/2
likSeasDummy <- kfsSeasDummy$logLik - sum(kfsSeasDummy$Finf>0) * log(2*pi)/2
likSeasTri0 <- kfsSeasTri0$logLik - sum(kfsSeasTri0$Finf>0) * log(2*pi)/2
likSeasTri <- kfsSeasTri$logLik - sum(kfsSeasTri$Finf>0) * log(2*pi)/2

# 1 期先予測の平均二乗誤差
mseSeasDummy0 <- sum(kfsSeasDummy0$v[14:144]^2) / 131
mseSeasDummy <- sum(kfsSeasDummy$v[14:144]^2) / 131
mseSeasTri0 <- sum(kfsSeasTri0$v[14:144]^2) / 131
mseSeasTri <- sum(kfsSeasTri$v[14:144]^2) / 131

# 散漫対数尤度の修正（季節変動を固定した同等なモデルの尤度を揃える）
likSeasTri <- likSeasTri - (likSeasTri0 - likSeasDummy0)
likSeasTri0 <- likSeasDummy0

# AIC (赤池情報量規準)
aicSeasDummy0 <- -2*likSeasDummy0 + 2*(2+13)
aicSeasDummy <- -2*likSeasDummy + 2*(3+13)
aicSeasTri0 <- -2*likSeasTri0 + 2*(2+13)
aicSeasTri <- -2*likSeasTri + 2*(8+13)

# 図3.2の描画
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
plot(sales$Fabric, type="l", lty=1, ylab = "販売額（10億円）",xaxt="n",xaxs="i",col=1,xlab="(a) 水準成分")
axis(side=1,at=1+0:11*12,labels=c("02/1","03/1","04/1","05/1","06/1","07/1","08/1","09/1","10/1","11/1","12/1","13/1"))
lines(kfs2$alphahat[,"level"], col=3)
lines(kfs4$alphahat[,"level"], col=4)
abline(v=120.5,lty=3)

plot(kfs2$alphahat[,"sea_dummy1"],type="l", ylab = "販売額（10億円）",xaxt="n",xaxs="i",yaxs="i",col=3,xlab="(b) 季節成分")
lines(kfs4$alphahat[,"sea_dummy1"],col=4)
axis(side=1,at=1+0:11*12,labels=c("02/1","03/1","04/1","05/1","06/1","07/1","08/1","09/1","10/1","11/1","12/1","13/1"))
abline(v=120.5,lty=3)

plot(sales$Fabric-kfs2$muhat,type="l", ylab="販売額（10億円）",　xaxt="n",xaxs="i",yaxs="i",col=3,ylim=c(-150,150),xlab="(c) 平滑化観測値撹乱項と長期予測誤差")
lines(sales$Fabric-kfs4$muhat,col=4)
axis(side=1,at=1+0:11*12,labels=c("02/1","03/1","04/1","05/1","06/1","07/1","08/1","09/1","10/1","11/1","12/1","13/1"))
abline(v=120.5,lty=3)
abline(h=0,col=8)


### 例3.3 ###

# AR(1) 成分モデル
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

# AR(2) 成分モデル
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

# AR(12) 成分モデル(ラグ12 以外の自己回帰係数はゼロとする)
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


### 例3.4 ###

Gasoline <- ts(scan("Gasoline.dat")) # ガソリン単価データの読み込み

# 回帰係数を時間変化させる場合
modRegression <- SSModel(log(sales$Fuel) ~ SSMtrend(1, Q = NA)
  + SSMseasonal(12, sea.type="dummy")
  + SSMregression(~ log(Gasoline), Q = NA), H = NA)
fitRegression <- fitSSM(modRegression, numeric(3), method = "BFGS")
kfsRegression <- KFS(fitRegression$model)

# 回帰係数を固定する場合
modRegression0 <- SSModel(log(sales$Fuel) ~ SSMtrend(1, Q = NA)
  + SSMseasonal(12, sea.type="dummy")
  + log(Gasoline), H = NA)
fitRegression0 <- fitSSM(modRegression0, numeric(2), method = "BFGS")
kfsRegression0 <- KFS(fitRegression0$model)

### 作り直したコード
# Reg0:回帰成分なし
# Reg1:回帰係数を固定する場合
# Reg2:回帰係数を時間変化させる場合

fuel <- sales[,4]                # 燃料小売業の月次販売額
gas  <- ts(scan("Gasoline.dat")) # ガソリン単価データの読み込み

# モデル定義
modReg0 <- SSModel(log(fuel) ~ SSMtrend(1,Q=NA) + SSMseasonal(12,Q=0), H=NA)

modReg1 <- SSModel(log(fuel) ~ log(gas) + 
  SSMtrend(1,Q=NA) + SSMseasonal(12,Q=0), H=NA)

modReg2 <- SSModel(log(fuel) ~ SSMregression(~ log(gas), Q=NA) +
  SSMtrend(1,Q=NA) + SSMseasonal(12,Q=0), H=NA)

# 未知パラメータの推定
fitReg0 <- fitSSM(modReg0, rep(-8,2), method = "BFGS")
fitReg1 <- fitSSM(modReg1, rep(-8,2), method = "BFGS")
fitReg2 <- fitSSM(modReg2, rep(-8,3), method = "BFGS")

# カルマンフィルタ・カルマンスムーザの実行
kfsReg0 <- KFS(fitReg0$model)
kfsReg1 <- KFS(fitReg1$model)
kfsReg2 <- KFS(fitReg2$model)

# 最大対数尤度
logLikReg0 <- kfsReg0$logLik - sum(kfsReg0$Finf>0) * log(2*pi)/2
logLikReg1 <- kfsReg1$logLik - sum(kfsReg1$Finf>0) * log(2*pi)/2
logLikReg2 <- kfsReg2$logLik - sum(kfsReg2$Finf>0) * log(2*pi)/2

# AIC (赤池情報量規準)
AICReg0 <- -2*logLikReg0 + 2*( 2 + 12 )
AICReg1 <- -2*logLikReg1 + 2*( 2 + 13 )
AICReg2 <- -2*logLikReg2 + 2*( 3 + 13 )

# 平滑化された水準成分，回帰成分の描画
par(mfrow=c(3,1))   # 描画領域を３分割
par(mar=c(2,4,1,1)) # 描画領域の余白設定
plot(log(fuel),type="o",lty=3,xaxt="n",xlab="",ylab="販売額（対数）")
lines(kfsReg1$alphahat[,"level"]+kfsReg1$alphahat[,"log(gas)"]*log(gas),lwd=2)
axis(side=1,at=1+0:11*12,labels=c("02/1","03/1","04/1","05/1","06/1","07/1","08/1","09/1","10/1","11/1","12/1","13/1"))
plot(kfsReg1$alphahat[,"level"],xaxt="n",xlab="",ylab="水準成分")
axis(side=1,at=1+0:11*12,labels=c("02/1","03/1","04/1","05/1","06/1","07/1","08/1","09/1","10/1","11/1","12/1","13/1"))
plot(kfsReg1$alphahat[,"log(gas)"]*log(gas),ylim=c(1.4,1.6),xaxt="n",xlab="",ylab="回帰成分")
axis(side=1,at=1+0:11*12,labels=c("02/1","03/1","04/1","05/1","06/1","07/1","08/1","09/1","10/1","11/1","12/1","13/1"))


### 例3.5 ###

# 各月の曜日集計
dates <- seq(as.Date("2002-01-01"), as.Date("2013-12-31"), by = 1)
weeks <- table(substr(dates,1,7), weekdays(dates, T))
sun <- weeks[,"日"]
mon <- weeks[,"月"]-sun; tue <- weeks[,"火"]-sun; wed <- weeks[,"水"]-sun
thu <- weeks[,"木"]-sun; fry <- weeks[,"金"]-sun; sat <- weeks[,"土"]-sun
calendar <- cbind(mon, tue, wed, thu, fry, sat)

# うるう年２月のダミー変数
leapyear <- rownames(weeks) %in% c("2004-02","2008-02","2012-02")

# カレンダー効果(曜日・うるう年) のあるモデル
modCalender <- SSModel(sales$Fabric ~ SSMtrend(2, Q = c(list(0), list(NA)))
  + SSMseasonal(12, sea.type="dummy")
  + leapyear + calendar, H = NA)
fitCalender <- fitSSM(modCalender, numeric(2), method = "BFGS")
kfsCalender <- KFS(fitCalender$model)

# 図3.7の描画
plot(kfsCalender$muhat - kfsCalender$alphahat[,"level"],type="l",xaxs="i",xaxt="n",xlab="",ylab="販売額（10億円）")
axis(side=1,at=1+0:11*12,labels=c("02/1","03/1","04/1","05/1","06/1","07/1","08/1","09/1","10/1","11/1","12/1","13/1"))

### 例3.6 ###

# 2010年11月，12月のデータを除外
salesNA <- sales$Machinery
salesNA[sales$month %in% c("2010年11月","2010年12月")] <- NA

# 2011年8月以降の水準シフト干渉変数の定義
ShiftLevel <- (1:nrow(sales) >= which(sales$month=="2011年8月"))

# 水準シフト干渉変数を加えたモデル
modShift <- SSModel(log(salesNA) ~ SSMtrend(1, Q = NA)
  + SSMseasonal(12, Q = NA, sea.type="dummy")
  + ShiftLevel, H = NA)
fitShift <- fitSSM(modShift, numeric(3))
kfsShift <- KFS(fitShift$model, smoothing=c("state","mean","disturbance"))

# 補助残差のプロット例
plot(rstandard(kfsShift, "pearson")) # 観測値撹乱項
plot(rstandard(kfsShift, "state")[,1]) # 状態撹乱項(水準成分)


# バージョン更新でエラーが出た場合は以下を実行ください
# 機械器具小売業販売額
machine <- sales[,3]

# 水準シフト干渉変数を加えたモデル
mod0 <- SSModel(log(machine) ~ SSMtrend(1, Q = NA)
  + SSMseasonal(12, Q = NA, sea.type="dummy"), H = NA)
fit0 <- fitSSM(mod0, numeric(3))
kfs0 <- KFS(fit0$model, smoothing=c("state","mean","disturbance"))

mod01 <- SSModel(log(machine) ~ SSMtrend(2, Q = list(NA,NA)), H = NA)
fit01 <- fitSSM(mod01, numeric(3))
kfs01 <- KFS(fit01$model, smoothing=c("state","mean","disturbance"))


# 状態攪乱項の標準化残差（kfs0$rstandard_state）の算出
kfs0$V_eta_inv <- array(apply(kfs0$V_eta,3,solve),dim(kfs0$V_eta))
kfs0$V_eta_inv_chol <- array(apply(kfs0$V_eta_inv,3,chol),dim(kfs0$V_eta))
kfs0$rstandard_state <- kfs0$etahat*0
for(i in 1:nrow(kfs0$etahat)) kfs0$rstandard_state[i,] = kfs0$V_eta_inv_chol[,,i]%*%kfs0$etahat[i,]


# 図3.9の描画
par(mfrow=c(3,1))   # 描画領域を３分割
par(mar=c(3,3,1,1)) # 描画領域の余白設定
par(mgp=c(2,0.5,0))
plot(log(machine),type="l",xaxt="n",xlab="(a) 販売額の対数系列と平滑化状態の水準成分",ylab="販売額（対数）")
lines(kfs0$alphahat[,"level"],col=8)
axis(side=1,at=1+0:11*12,labels=c("02/1","03/1","04/1","05/1","06/1","07/1","08/1","09/1","10/1","11/1","12/1","13/1"))
plot(rstandard(kfs0, "pearson"),ylim=c(-6,6),xaxt="n",xlab="(b) 観測値攪乱項の補助残差",ylab="販売額（対数）")
abline(h=c(-1.96,1.96),lty=3)
axis(side=1,at=1+0:11*12,labels=c("02/1","03/1","04/1","05/1","06/1","07/1","08/1","09/1","10/1","11/1","12/1","13/1"))
plot(kfs0$rstandard_state[,1],ylim=c(-6,6),xaxt="n",xlab="(c) 状態攪乱項（水準成分）の補助残差",ylab="販売額（対数）")
abline(h=c(-1.96,1.96),lty=3)
axis(side=1,at=1+0:11*12,labels=c("02/1","03/1","04/1","05/1","06/1","07/1","08/1","09/1","10/1","11/1","12/1","13/1"))

# 2010 年11 月，12 月のデータを除外
machineNA <- machine
machineNA[sales[,1] %in% c("2010年11月","2010年12月")] <- NA

# 2011 年8 月以降の水準シフト干渉変数の定義
ShiftLevel <- (1:nrow(sales) >= which(sales[,1]=="2011年8月"))

# 水準シフト干渉変数を加えたモデル
modShift <- SSModel(log(machineNA) ~ SSMtrend(1, Q = NA)
  + SSMseasonal(12, Q = NA, sea.type="dummy")
  + ShiftLevel, H = NA)
fitShift <- fitSSM(modShift, numeric(3))
kfsShift <- KFS(fitShift$model, smoothing=c("state","mean","disturbance"))

# 状態攪乱項の標準化残差（kfsShift$rstandard_state）の算出
kfsShift$V_eta_inv <- array(apply(kfsShift$V_eta,3,solve),dim(kfsShift$V_eta))
kfsShift$V_eta_inv_chol <- array(apply(kfsShift$V_eta_inv,3,chol),dim(kfsShift$V_eta))
kfsShift$rstandard_state <- kfsShift$etahat*0
for(i in 1:nrow(kfsShift$etahat)) kfsShift$rstandard_state[i,] = kfsShift$V_eta_inv_chol[,,i]%*%kfsShift$etahat[i,]

# 図3.10の描画
par(mfrow=c(3,1))   # 描画領域を３分割
par(mar=c(3,3,1,1)) # 描画領域の余白設定
par(mgp=c(2,0.5,0))
plot(log(machine),type="l",xaxt="n",xlab="(a) 販売額の対数系列と平滑化状態の水準成分",ylab="販売額（対数）")
lines(kfsShift$muhat - kfsShift$alphahat[,"sea_dummy1"],col=8)
axis(side=1,at=1+0:11*12,labels=c("02/1","03/1","04/1","05/1","06/1","07/1","08/1","09/1","10/1","11/1","12/1","13/1"))
plot(rstandard(kfsShift, "pearson"),ylim=c(-6,6),xaxt="n",xlab="(b) 観測値攪乱項の補助残差",ylab="販売額（対数）")
abline(h=c(-1.96,1.96),lty=3)
axis(side=1,at=1+0:11*12,labels=c("02/1","03/1","04/1","05/1","06/1","07/1","08/1","09/1","10/1","11/1","12/1","13/1"))
plot(kfsShift$rstandard_state[,1],ylim=c(-6,6),xaxt="n",xlab="(c) 状態攪乱項（水準成分）の補助残差",ylab="販売額（対数）")
abline(h=c(-1.96,1.96),lty=3)
axis(side=1,at=1+0:11*12,labels=c("02/1","03/1","04/1","05/1","06/1","07/1","08/1","09/1","10/1","11/1","12/1","13/1"))


### 例3.7 ###

Weight <- ts(scan("Weight.dat"))   # 体重データの読み込み
Bodyfat <- ts(scan("Bodyfat.dat")) # 体脂肪率データの読み込み

WeightNA <- Weight
WeightNA[21:40] <- NA

modSUTSE <- SSModel(cbind(WeightNA, Bodyfat) ~
  SSMtrend(1, Q = matrix(NA,2,2), type = "distinct"), H = matrix(NA,2,2))
fitSUTSE <- fitSSM(modSUTSE, numeric(6), method="BFGS")
kfsSUTSE <- KFS(fitSUTSE$model)
