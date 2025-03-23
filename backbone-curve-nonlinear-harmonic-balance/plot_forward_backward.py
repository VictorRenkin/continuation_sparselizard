import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

backward = pd.read_csv("data/backward.csv")
forward = pd.read_csv("data/forward.csv")
freq_backward = backward["w"]
amp_backward = backward["u"]

freq_forward = forward["w"]
amp_forward = forward["u"]

# Tracé du graphe
plt.figure(figsize=(10, 6))
plt.scatter(freq_forward/2/np.pi, amp_forward, color='blue', label='forward_sweep')
plt.scatter(freq_backward/2/np.pi, amp_backward, color='red', label='backward_sweep')

plt.xlabel("Frequency [Hz]")
plt.ylabel("Displacement u [m]")
plt.legend()

plt.savefig("fig/ClampedBeam_NLFR.pdf", bbox_inches='tight', dpi=300)