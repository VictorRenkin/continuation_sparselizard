import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

backward = pd.read_csv("data/backward.csv")
forward = pd.read_csv("data/forward.csv")
freq_backward = backward["u"]
amp_backward = backward["w"]

freq_forward = forward["u"]
amp_forward = forward["w"]

# Tracé du graphe
plt.figure(figsize=(10, 6))
plt.scatter(freq_forward/2/np.pi, amp_forward, color='blue', label='forward_sweep')
plt.scatter(freq_backward/2/np.pi, amp_backward, color='red', label='backward_sweep')

plt.xlabel("Frequency [Hz]")
plt.ylabel("Displacement u [m]")
plt.legend()

plt.savefig("fig/ClampedBeam_NLFR.pdf", bbox_inches='tight', dpi=300)