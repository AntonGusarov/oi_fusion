# Data assimilation using ensemble optimal interpolation 

This is the code used in the published paper 
[(link)](https://pdf.sciencedirectassets.com/280203/1-s2.0-S1877050917X00185/1-s2.0-S1877050917324018/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEOv%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJHMEUCIF3HyEj4mBd52IK3zURu9duVe8bGG307y0ilSWDBO2ddAiEA8YWMCaHkavyyyVaCcua6nmxsMZWq7JQsa%2Fw7AU2c6G8qtAMIcxACGgwwNTkwMDM1NDY4NjUiDDnfD9In0KNUouQGkiqRA36VqpAJZXvfNnC2ooX6B8PC%2FVcZhkS82EMgsfgRTYEGqrZGqMJ%2BUKA3%2FlkbeJsoixCE1cl%2FijXhEJxehplFadwbJylGZm5Nt9K%2BFxRUYT%2BCEjv5MYX8Wz0mPvB4WSFERzd4TnN0wsXztgdlOaMMeM8m2RhvwyJa70Fg5irDeClgONDuHV6bjOtSy0WIhp8bywcijpx%2B5mVOs81FIAggaDlGS1EG27dlzncYLED8TMQqCp6cdii%2FRe66e%2FevpKgutDkG3ebvk4AabCkDDRlBgwTJXvCeB%2FzJAY%2Bua8mXWKerqnlqCEGob%2FJChxHY9aC3ZvbW3nvA7MXc3rujoNmxdlkjjaylk88QzwRzBvrvSjaZTntCCEZDh9Of9%2Fpq7lr2cfCkqFTDGYmOZYMqfplfQLD8jrTzjuuugWrUJvy3ou01SOFJ8EmzDPRHEItmIXPyg9PPXheAhGTT7FRyw5DtrX01cnF3whlPsbH1kQhW1mKvxpfBAVnm9QvPH0FihpnkrbRgLHCPmXQQQCvjySsX1E%2FvMMTN%2B%2FAFOusBmmVohI%2FZF8R416EzBvOm6g0TeWqMgrDSIP%2F1iDG031qdPQ%2B58sEPO4ttDvQO7YXO8bLEQFPcxFvB2CFWcHKNYDLygugpeaNRQARoR2naqL61BNI8EPS7%2FM4XXPJ3pDV39xnVyUwILuOXv2cGPZRZj3nPDB8lKRet48FweeRROVHAPwT%2Bg7UNvYhsUtKzqWcFi2YITM1SpBbzDGphk%2BfiM3x1k4tPmUfQVKKdiCX%2BIjuQ9WEgVFWXfTnm5Q3UntVB1%2Br7uBeGP9%2FahDQx03M9UAxm9vOLTSHAhIpcF1tTjv%2BJx6LzINzcAxblCw%3D%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20200115T112738Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYSTT3PKHB%2F20200115%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=11104111acd2c5b5bbde5995b14ea96773c4a68f5d1420f133a2769de3d96376&hash=ef367e610d7f5354ba7e1b7e81a829c78afed3f7f5b528b52e75465ee62a5945&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S1877050917324018&tid=spdf-5612e98e-b063-4c91-a35b-d67c416c63ce&sid=c6eb608775fd5340483b0a5-04e8b0c018afgxrqb&type=client)
and implemented in MATLAB.

## 1. Introduction
 **Ensemble optimal interpolation (EnOI)** is an efficient, having a relatively low computational cost yet powerful
 method for correcting numerical simulation models in according to the field measurements. Basic EnOI is used in different areas 
 of geoscience including meteorology and oceanography. 
 
 The goal was to apply the method of EnOI to post-process the numerically simulated wind field in the way that it would maximally 
 fit the actual measurements at the weather stations. The basic principle behind this technique is to interpolate discrepancies 
 between measurements and simulation following the spatiotemporal consistency 
 of the modeled wind field. However, its statistical properties are in general non-stationary, spatially inhomogeneous 
 as well as the quality of measurements vary from station to station. 
 Thus, such a correction can introduce undesirable artifacts (outliers) into the output. We adressed this problem by 
 adding the ridge regression into the basic EnOI scheme.
 
 The algorithm was applied to assimilate observations to the wind field simulated for the Polar Arctic region.
 
![EnOI work example for Arctic](https://github.com/AntonGusarov/oi_fusion/blob/master/Fields.jpg)
 
 ## 2. Mathematical model
 
### Ensemble optimal interpolation

Resulting correction ![x_i^{a}](https://render.githubusercontent.com/render/math?math=x_i%5E%7Ba%7D) can be expressed as a weighted linear combination of the differences between observed and modeled time series 
![\Delta_n(\tau) = y_n^{o}(\tau)-Hx_n^{b}(\tau) ](https://render.githubusercontent.com/render/math?math=%5CDelta_n(%5Ctau)%20%3D%20y_n%5E%7Bo%7D(%5Ctau)-Hx_n%5E%7Bb%7D(%5Ctau)%20) added then to each grid point of the modeled (background) field ![x_i^{b}](https://render.githubusercontent.com/render/math?math=x_i%5E%7Bb%7D):

![x_i^{a}(t) = x_i^{b}(t) + \sum_{n=1}^{N} \sum_{\tau=t-T}^{t+T} w_n^{(i,t)}(\tau)\Delta_n(\tau)](https://render.githubusercontent.com/render/math?math=x_i%5E%7Ba%7D(t)%20%3D%20x_i%5E%7Bb%7D(t)%20%2B%20%5Csum_%7Bn%3D1%7D%5E%7BN%7D%20%5Csum_%7B%5Ctau%3Dt-T%7D%5E%7Bt%2BT%7D%20w_n%5E%7B(i%2Ct)%7D(%5Ctau)%5CDelta_n(%5Ctau))

where:\
![N](https://render.githubusercontent.com/render/math?math=N)- number of observation points,\
![t](https://render.githubusercontent.com/render/math?math=t) - time instance of the background model,\
![T](https://render.githubusercontent.com/render/math?math=T) - maximum time lag that depends on the variability of the studied process.

The conventional OI results in the Best Linear Unbiased Estimate (BLUE) ![\textbf{x}^{a}](https://render.githubusercontent.com/render/math?math=%5Ctextbf%7Bx%7D%5E%7Ba%7D) of the true system state column vector 
![$\textbf{x}^{t}$](https://render.githubusercontent.com/render/math?math=%24%5Ctextbf%7Bx%7D%5E%7Bt%7D%24) having available vector of observations ![\textbf{y}^{o}](https://render.githubusercontent.com/render/math?math=%5Ctextbf%7By%7D%5E%7Bo%7D), 
initial background field state vector ![\textbf{x}^{b}](https://render.githubusercontent.com/render/math?math=%5Ctextbf%7Bx%7D%5E%7Bb%7D), and estimated error covariances:

![\textbf{x}^{a}=\textbf{x}^{b}+\textbf{K}(\textbf{y}^{o}-\textbf{H}\textbf{x}^{b})](https://render.githubusercontent.com/render/math?math=%5Ctextbf%7Bx%7D%5E%7Ba%7D%3D%5Ctextbf%7Bx%7D%5E%7Bb%7D%2B%5Ctextbf%7BK%7D(%5Ctextbf%7By%7D%5E%7Bo%7D-%5Ctextbf%7BH%7D%5Ctextbf%7Bx%7D%5E%7Bb%7D))

![\textbf{K}=\textbf{P}^{b}\textbf{H}^{T}(\textbf{H}\textbf{P}^b\textbf{H}^{T}+\textbf{R})^{-1}](https://render.githubusercontent.com/render/math?math=%5Ctextbf%7BK%7D%3D%5Ctextbf%7BP%7D%5E%7Bb%7D%5Ctextbf%7BH%7D%5E%7BT%7D(%5Ctextbf%7BH%7D%5Ctextbf%7BP%7D%5Eb%5Ctextbf%7BH%7D%5E%7BT%7D%2B%5Ctextbf%7BR%7D)%5E%7B-1%7D)

where:\
![\textbf{P}^b = \langle (\textbf{x}^b - \textbf{x}^t){(\textbf{x}^b - \textbf{x}^t)}^T \rangle](https://render.githubusercontent.com/render/math?math=%5Ctextbf%7BP%7D%5Eb%20%3D%20%5Clangle%20(%5Ctextbf%7Bx%7D%5Eb%20-%20%5Ctextbf%7Bx%7D%5Et)%7B(%5Ctextbf%7Bx%7D%5Eb%20-%20%5Ctextbf%7Bx%7D%5Et)%7D%5ET%20%5Crangle)

and 

![\textbf{R} = \langle (\textbf{y}^o - \textbf{x}^t){(\textbf{y}^o - \textbf{x}^t)}^T \rangle](https://render.githubusercontent.com/render/math?math=%5Ctextbf%7BR%7D%20%3D%20%5Clangle%20(%5Ctextbf%7By%7D%5Eo%20-%20%5Ctextbf%7Bx%7D%5Et)%7B(%5Ctextbf%7By%7D%5Eo%20-%20%5Ctextbf%7Bx%7D%5Et)%7D%5ET%20%5Crangle) 

are the matrices of estimated covariances of the background field and observation errors, respectively;\
![\textbf{H}](https://render.githubusercontent.com/render/math?math=%5Ctextbf%7BH%7D) - the representation matrix of mapping observations onto the corresponding model grid;\
![\mathbf{K}](https://render.githubusercontent.com/render/math?math=%5Cmathbf%7BK%7D) - the weight matrix referred to as the Kalman gains, since the OI corresponds to the analysis stage of the Kalman filter. 
