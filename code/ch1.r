Weight <- ts(scan("Weight.dat")) # �̏d�f�[�^�̓ǂݍ���

arima(Weight, c(1,1,0)) # ARIMA(1,1,0) ���f���̓��Ă͂߂Ɛ���
arima(Weight, c(0,1,1)) # ARIMA(0,1,1) ���f���̓��Ă͂߂Ɛ���
