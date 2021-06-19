import pandas as pd
import numpy as np
import math

df_all = pd.read_csv('retty_gourmet_open_q1.csv')

df_answer = df_all['answer']
df_latitude = df_all['latitude']
df_longitude = df_all['longitude']  # 経度
df_category = df_all['category']


for i in range(len(df_all)):
    if df_category[i] == '餃子':
        latitude = float(df_latitude[i])
        longitude = float(df_longitude[i])
        
        def calc_xy(phi_deg, lambda_deg, phi0_deg, lambda0_deg):
            # 緯度経度・平面直角座標系原点をラジアンに直す
            phi_rad = np.deg2rad(phi_deg)
            lambda_rad = np.deg2rad(lambda_deg)
            phi0_rad = np.deg2rad(phi0_deg)
            lambda0_rad = np.deg2rad(lambda0_deg)

            # 補助関数
            def A_array(n):
                A0 = 1 + (n ** 2) / 4. + (n ** 4) / 64.
                A1 = -(3. / 2) * (n - (n ** 3) / 8. - (n ** 5) / 64.) 
                A2 = (15. / 16) * (n**2 - (n ** 4) / 4.)
                A3 = -(35. / 48) * (n ** 3 - (5. / 16) * (n**5))
                A4 = (315. / 512) * (n**4)
                A5 = -(693. / 1280) * (n**5)
                return np.array([A0, A1, A2, A3, A4, A5])

            def alpha_array(n):
                a0 = np.nan  # dummy
                a1 = (1. / 2) * n - (2. / 3) * (n ** 2) + (5. / 16) * (n ** 3) + (41. / 180) * (n ** 4) - (127. / 288) * (n ** 5)
                a2 = (13. / 48) * (n ** 2) - (3. / 5) * (n ** 3) + (557. / 1440) * (n ** 4) + (281. / 630) * (n ** 5)
                a3 = (61. / 240) * (n ** 3) - (103. / 140) * (n ** 4) + (15061. / 26880) * (n ** 5)
                a4 = (49561. / 161280) * (n ** 4) - (179. / 168) * (n ** 5)
                a5 = (34729. / 80640) * (n ** 5)
                return np.array([a0, a1, a2, a3, a4, a5])

            # 定数 (a, F: 世界測地系-測地基準系1980（GRS80）楕円体)
            m0 = 0.9999
            a = 6378137.
            F = 298.257222101

            # (1) n, A_i, alpha_iの計算
            n = 1. / (2 * F - 1)
            A_array = A_array(n)
            alpha_array = alpha_array(n)

            # (2), S, Aの計算
            A_ = ((m0 * a) / (1. + n)) * A_array[0]  # [m]
            S_ = ((m0 * a) / (1. + n)) * (A_array[0] * phi0_rad + np.dot(A_array[1:], np.sin(2 * phi0_rad * np.arange(1, 6))))  # [m]

            # (3) lambda_c, lambda_sの計算
            lambda_c = np.cos(lambda_rad - lambda0_rad)
            lambda_s = np.sin(lambda_rad - lambda0_rad)

            # (4) t, t_の計算
            t = np.sinh(np.arctanh(np.sin(phi_rad)) - ((2 * np.sqrt(n)) / (1 + n)) * np.arctanh(((2 * np.sqrt(n)) / (1 + n)) * np.sin(phi_rad)))
            t_ = np.sqrt(1 + t * t)

            # (5) xi', eta'の計算
            xi2 = np.arctan(t / lambda_c)  # [rad]
            eta2 = np.arctanh(lambda_s / t_)

            # (6) x, yの計算
            x = A_ * (xi2 + np.sum(np.multiply(alpha_array[1:], np.multiply(np.sin(2 * xi2 * np.arange(1, 6)), np.cosh(2 * eta2 * np.arange(1, 6)))))) - S_  # [m]
            y = A_ * (eta2 + np.sum(np.multiply(alpha_array[1:], np.multiply(np.cos(2 * xi2 * np.arange(1, 6)), np.sinh(2 * eta2 * np.arange(1, 6))))))  # [m]
            # return
            return x, y  # [m]
            
        x, y = calc_xy(latitude, longitude, 35.6721903, 139.7363287)
        
        distance = math.sqrt(x ** 2 + y ** 2) / 1000
        if distance <= 1.5:
            print(df_answer[i])
