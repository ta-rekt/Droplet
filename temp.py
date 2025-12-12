import numpy as np
import matplotlib.pyplot as plt

labels = ["tf", "rnap", "enhancer", "dna"]
intmat = np.array([[-1.0, -5.0, -2.0, -0.5], 
                   [-5.0, -1.0, -0.2, -0.2], 
                   [-2.0, -0.2,  0.0,  0.0], 
                   [-0.5, -0.2,  0.0,  0.0]])

fig, ax = plt.subplots(figsize=(3, 3))
ax.matshow(intmat, cmap="coolwarm", alpha=0.5)
ax.set_xticks(range(4)); ax.set_yticks(range(4))
ax.set_xticklabels(labels, rotation=45, ha="left")
ax.set_yticklabels(labels, rotation=45, va="bottom")

for (i, j), val in np.ndenumerate(intmat):
    ax.text(j, i, val, ha='center', va='center')
plt.show()
