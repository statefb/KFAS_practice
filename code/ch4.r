library(KFAS); library(Nippon)

fire     <- ts(scan("fire.dat"))     # 火災件数データの読み込み
humidity <- ts(scan("humidity.dat")) # 日平均湿度データの読み込み
celsius  <- ts(scan("celsius.dat"))  # 日平均気温データの読み込み
rain     <- ts(scan("rain.dat"))     # 日降水量データの読み込み

dates <- seq(as.Date("2005-01-01"), as.Date("2014-12-31"), by=1) # 日付列作成
weekday <- weekdays(dates, T) # 曜日判別関数("月"～"日"を返す)
weekday[is.jholiday(dates)] <- "祝" # パッケージNippon の祝日判別関数
weekday <- factor(weekday, c("月","火","水","木","金","土","日","祝"))

# 因子化
modPois <- SSModel(fire ~ SSMtrend(2, Q = c(list(0), list(NA)))
  + SSMcycle(365.25) + SSMcycle(365.25/2) + SSMcycle(365.25/3)
  + SSMcycle(365.25/4) + SSMcycle(365.25/5) + SSMcycle(365.25/6)
  + weekday + humidity + celsius + rain,
distribution = "poisson")
diag(modPois$P1inf) <- 0 # 散漫初期化が上手くいかなかったため，
diag(modPois$P1) <- 100^2 # 代わりに十分大きな初期分散を用いる

# 非ガウスモデルではシミュレーション数nsim を設定する
# （時間・メモリーの制約で上手くいかない場合、nsimを減らして実行してください）
fitPois <- fitSSM(modPois, 0, method="Brent", lower=-40, upper=0, nsim=1000)
kfsPois <- KFS(fitPois$model, nsim=1000)

# インポータンス・サンプリングと利用例(湿度効果(第8 状態成分) の95%信頼区間)
impPois <- importanceSSM(fitPois$model, type="states", nsim=4000)
emp <- cumsum(impPois$weights[order(impPois$samples[1,8,])])/
sum(impPois$weight)
conf <- sort(impPois$samples[1,8,])[rank(c(0.025, 0.975, emp))[1:2]]

# 観測値の98％予測区間
prePois <- predict(fitPois$model, interval="prediction", level=0.98, nsim=10000)
# 解析結果の表示

kfsPois # 回帰係数の推定値を標準誤差を表示

# 回帰係数の95%信頼区間を表示
sort(impPois$samples[1,1,])[rank(c(0.025, 0.975, emp))[1:2]] # 火曜
sort(impPois$samples[1,2,])[rank(c(0.025, 0.975, emp))[1:2]] # 水曜
sort(impPois$samples[1,3,])[rank(c(0.025, 0.975, emp))[1:2]] # 木曜
sort(impPois$samples[1,4,])[rank(c(0.025, 0.975, emp))[1:2]] # 金曜
sort(impPois$samples[1,5,])[rank(c(0.025, 0.975, emp))[1:2]] # 土曜
sort(impPois$samples[1,6,])[rank(c(0.025, 0.975, emp))[1:2]] # 日曜
sort(impPois$samples[1,7,])[rank(c(0.025, 0.975, emp))[1:2]] # 祝日
sort(impPois$samples[1,8,])[rank(c(0.025, 0.975, emp))[1:2]] # 湿度
sort(impPois$samples[1,9,])[rank(c(0.025, 0.975, emp))[1:2]] # 気温
sort(impPois$samples[1,10,])[rank(c(0.025, 0.975, emp))[1:2]] # 降水量

# 平滑化状態による火災件数の期待値、98%予測区間と外れ値
plot(prePois[,1], ylim=c(0,52), xaxs="i", yaxs="i", ylab="火災件数", xaxt="n")
lines(prePois[,2], col=8)
lines(prePois[,3], col=8)
points(which(fire>prePois[,3]),fire[fire>prePois[,3]])
points(which(fire<prePois[,2]),fire[fire<prePois[,2]])
axis(side=1,at=which(substr(dates,6,10)=="01-01"),
 labels=c("05/1/1","06/1/1","07/1/1","08/1/1","09/1/1","10/1/1","11/1/1","12/1/1","13/1/1","14/1/1"))

# 負の二項分布モデルによる解析

hum_cel = humidity * celsius # 湿度と温度の交互作用
fireNA <- fire; fireNA[dates=="2011-03-11"] <- NA # 外れ値を欠測値に替え除外
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

# インポータンス・サンプリングと利用例(湿度効果(第8 状態成分) の95%信頼区間)
impNegbin <- importanceSSM(fitNegbin$model, type="states", nsim=1000)
emp <- cumsum(impNegbin$weights[order(impNegbin$samples[1,8,])])/
sum(impNegbin$weight)
conf <- sort(impNegbin$samples[1,8,])[rank(c(0.025, 0.975, emp))[1:2]]

# 観測値の98％予測区間
preNegbin <- predict(fitNegbin$model, interval="prediction", level=0.98, nsim=1000)


# 解析結果の表示

kfsNegbin # 回帰係数の推定値を標準誤差を表示

# 回帰係数の95%信頼区間を表示
sort(impNegbin$samples[1,1,])[rank(c(0.025, 0.975, emp))[1:2]] # 火曜
sort(impNegbin$samples[1,2,])[rank(c(0.025, 0.975, emp))[1:2]] # 水曜
sort(impNegbin$samples[1,3,])[rank(c(0.025, 0.975, emp))[1:2]] # 木曜
sort(impNegbin$samples[1,4,])[rank(c(0.025, 0.975, emp))[1:2]] # 金曜
sort(impNegbin$samples[1,5,])[rank(c(0.025, 0.975, emp))[1:2]] # 土曜
sort(impNegbin$samples[1,6,])[rank(c(0.025, 0.975, emp))[1:2]] # 日曜
sort(impNegbin$samples[1,7,])[rank(c(0.025, 0.975, emp))[1:2]] # 祝日
sort(impNegbin$samples[1,8,])[rank(c(0.025, 0.975, emp))[1:2]] # 湿度
sort(impNegbin$samples[1,9,])[rank(c(0.025, 0.975, emp))[1:2]] # 気温
sort(impNegbin$samples[1,10,])[rank(c(0.025, 0.975, emp))[1:2]] # 湿度×気温

# 平滑化状態による火災件数の期待値、98%予測区間と外れ値
plot(preNegbin[,1], ylim=c(0,52), xaxs="i", yaxs="i", ylab="火災件数", xaxt="n")
lines(preNegbin[,2], col=8)
lines(preNegbin[,3], col=8)
points(which(fire>preNegbin[,3]),fire[fire>preNegbin[,3]])
points(which(fire<preNegbin[,2]),fire[fire<preNegbin[,2]])
axis(side=1,at=which(substr(dates,6,10)=="01-01"),
 labels=c("05/1/1","06/1/1","07/1/1","08/1/1","09/1/1","10/1/1","11/1/1","12/1/1","13/1/1","14/1/1"))

