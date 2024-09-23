# result.csvを読んでグラフ描画

import csv
import matplotlib.pyplot as plt

x = []
y_true = []
y_orig = []
y_fined = []
error_orig = []
error_fined = []

xj_x = []
xj_y = []

# CSVからデータを読み出して配列に追加
with open('./result.csv') as file:
    reader = csv.reader(file)
    for row in reader:
        nums = [float(v) for v in row]
        x.append(nums[0])
        y_true.append(nums[1])
        y_orig.append(nums[2])
        y_fined.append(nums[3])
        error_orig.append(nums[4])
        error_fined.append(nums[5])

with open('./points.csv') as file:
    reader = csv.reader(file)
    for row in reader:
        nums = [float(v) for v in row]
        xj_x.append(nums[0])
        xj_y.append(nums[1])

# --- 描画 --- #
# 近似曲線
fig1 = plt.figure(figsize = (10, 7))  # 横, 縦
ax1 = fig1.add_subplot(211)
ax1.set_title('Approximate curve', fontsize=20)
ax1.step(x, y_true, label="true", linestyle="dashdot")
ax1.step(x, y_orig, label="origin", linestyle="--", color="green")
ax1.step(x, y_fined, label="fined", linestyle="-", color="orange")
ax1.set_xlim(x[0], x[-1])
ax1.tick_params(labelsize=13)  # 軸目盛の大きさ
ax1.legend(fontsize=15, loc="upper right")
ax1.set_ylabel("sinc(x)", fontsize=15)

# 誤差曲線
ax2 = fig1.add_subplot(212)
ax2.step(x, error_orig,  label="origin", linestyle="--", color="green")
ax2.step(x, error_fined, label="fined", linestyle="-", color="orange")
ax2.scatter(xj_x, xj_y, label="deviation point", marker='o', zorder=2) # 一番手前に表示
ax2.set_xlim(x[0], x[-1])
ax2.tick_params(labelsize=13)
ax2.legend(fontsize=15, loc="upper right")
ax2.set_ylabel("error",   fontsize=15)
ax2.set_xlabel("x [rad]", fontsize=15)

plt.show()