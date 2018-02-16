Weight <- ts(scan("Weight.dat")) # 体重データの読み込み

arima(Weight, c(1,1,0)) # ARIMA(1,1,0) モデルの当てはめと推定
arima(Weight, c(0,1,1)) # ARIMA(0,1,1) モデルの当てはめと推定
