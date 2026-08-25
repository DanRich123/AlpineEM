import numpy as np
import matplotlib.pyplot as plt

# 1. Generate continuous spiral parameters
a = 0.5
b = 1
theta = np.linspace(0, 8 * np.pi, 5000)  # High resolution prevents gaps

# 2. Compute coordinates for Arm 1
r1 = a + b * theta
x1 = r1 * np.cos(theta)
y1 = r1 * np.sin(theta)

# 3. Compute coordinates for Arm 2 (180-degree / pi radians phase shift)
x2 = r1 * np.cos(theta + np.pi)
y2 = r1 * np.sin(theta + np.pi)

# Concatenate both arms into single arrays
x = np.concatenate([x1, x2])
y = np.concatenate([y1, y2])

# 4. Define discretized grid parameters
grid_size = 140  # Resolution: 140x140 array
max_val = np.max(np.abs([x, y])) * 1.05  # Grid extent bound with a small margin

# 5. Map continuous (x, y) coordinates to discrete matrix indices [0, grid_size - 1]
col_indices = np.clip(np.floor((x + max_val) / (2 * max_val) * (grid_size - 1)).astype(int), 0, grid_size - 1)
row_indices = np.clip(np.floor((y + max_val) / (2 * max_val) * (grid_size - 1)).astype(int), 0, grid_size - 1)

# 6. Create the binary mask array (0s and 1s)
mask = np.zeros((grid_size, grid_size), dtype=int)
mask[row_indices, col_indices] = 1

# Inspect the result array
print("Mask shape:", mask.shape)
print("Total elements set to 1:", np.sum(mask))

# 7. Plot the 2-arm binary mask
plt.figure(figsize=(6, 6))
plt.imshow(mask, cmap='binary', origin='lower')
plt.title('Discretized 2-Arm Binary Mask')
plt.xlabel('Grid Column Index (X)')
plt.ylabel('Grid Row Index (Y)')
plt.colorbar(label='Mask Value (0 or 1)')
plt.show()
plt.savefig('spiral.png')

# 8. Ouput the file with the mask as a numpy array
np.save('mask', mask)