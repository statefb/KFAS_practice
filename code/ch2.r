# install.packages("KFAS") # ミラーサーバより取得してインストール(1 度きり)
library(KFAS) # R コンソールへのパッケージの読み込み(R の起動ごとに必要)

Weight <- ts(scan("Weight.dat")) # 体重データの読み込み

# ローカルレベルモデルの定義
mod <- SSModel(Weight ~ SSMtrend(1, Q = NA), H = NA)

# 未知パラメータの推定
fit <- fitSSM(mod, numeric(2), method = "BFGS")

# カルマンフィルタと平滑化の実行
kfs <- KFS(fit$model)

# フィルタ化推定量とその信頼区間
afilt <- kfs$a[-1]
Pfilt <- kfs$P[,,-1] - fit$model$Q
afiltconf <- cbind(afilt+sqrt(Pfilt)*qnorm(0.025),afilt+sqrt(Pfilt)*qnorm(0.975))

# 図2.2の描画
plot(Weight,type="o",lty=3,xlab="経過日数",ylab="体重（kg）",ylim=c(83,87))
lines(afilt,lwd=2)
lines(afiltconf[,1])
lines(afiltconf[,2])

# 平滑化状態の信頼区間(kfs$alphahat+sqrt(kfs$V)*pnorm(0.025) などでもよい)
alphahatconf <- predict(fit$model, interval = "confidence", level = 0.95)

# 図2.3の描画
plot(Weight,type="o",lty=3,xlab="経過日数",ylab="体重（kg）",ylim=c(83,87))
lines(alphahatconf[,1],lwd=2)
lines(alphahatconf[,2])
lines(alphahatconf[,3])

# 欠測値の補間
modNA <- SSModel(Weight[c(1:20,rep(NA,20),41:60)] ~ SSMtrend(1, Q = NA), H = NA)
fitNA <- fitSSM(modNA, numeric(2), method = "BFGS")
preNA <- predict(fitNA$model, interval = "prediction", level = 0.95)
confNA <- predict(fitNA$model, interval = "confidence", level = 0.95)

# 図2.4の描画
plot(Weight,type="o",lty=3,xlab="経過日数",ylab="体重（kg）",ylim=c(83,87))
lines(21:40,Weight[21:40],type="o",lty=3,col=8)
lines(confNA[,1],lwd=2)
lines(confNA[,2])
lines(confNA[,3])
lines(21:40,preNA[21:40,2],lty=2)
lines(21:40,preNA[21:40,3],lty=2)

# 長期予測
mod50 <- SSModel(Weight[1:50] ~ SSMtrend(1, Q = NA), H = NA)
fit50 <- fitSSM(mod50, numeric(2), method = "BFGS")
pre50 <- predict(fit50$model, interval ="prediction", n.ahead = 10,
level = 0.95)
conf50 <- predict(fit50$model, interval ="confidence", n.ahead = 10)

# 図2.5の描画
plot(Weight,type="o",lty=3,xlab="経過日数",ylab="体重（kg）",ylim=c(83,87))
lines(51:60,Weight[51:60],type="o",lty=3,col=8)
lines(51:60,conf50[,1],lwd=2)
lines(51:60,conf50[,2])
lines(51:60,conf50[,3])
lines(51:60,pre50[,2],lty=2)
lines(51:60,pre50[,3],lty=2)

# 長期予測（欠測値NAとして予測するやり方）
mod50NA <- SSModel(Weight[c(1:50,rep(NA,10))] ~ SSMtrend(1, Q = NA), H = NA)
fit50NA <- fitSSM(mod50NA, numeric(2), method = "BFGS")
conf50NA <- predict(fit50NA$model, interval="confidence", level=0.95)
pre50NA  <- predict(fit50NA$model, interval="prediction", level=0.95)

