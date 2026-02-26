import numpy as np
import matplotlib.pyplot as plt

np.random.seed(42)

# Simulated predicted underdetermination scores
U_pred = np.sort(np.random.uniform(0, 1, 100))

# Simulated true ambiguity reduction correlated with U_pred
Delta_true = U_pred + 0.1*np.random.randn(100)

# Bin into deciles
bins = np.linspace(0,1,6)
digitized = np.digitize(U_pred, bins)

bin_means = [Delta_true[digitized == i].mean() for i in range(1,len(bins))]

plt.figure()
plt.plot(bins[:-1], bin_means, marker='o')
plt.xlabel("Predicted underdetermination $U_k$")
plt.ylabel("Observed ambiguity reduction $\Delta U_k$")
plt.title("Calibration under synthetic masking")
plt.savefig("Fig_calibration_curve.pdf")
plt.show()
