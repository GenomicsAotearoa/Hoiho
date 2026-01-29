import matplotlib.pyplot as plt
import re

# If you're reading from a file, change this to: open("reads.sam")
import sys
ds_values = []

pattern = re.compile(r'\bDS:f:([0-9.+-eE]+)')

for line in sys.stdin:
    match = pattern.search(line)
    if match:
        ds_values.append(float(match.group(1)))

# Plot
plt.hist(ds_values, bins=50, edgecolor='black')
plt.title("Histogram of DS:f: values")
plt.xlabel("DS value")
plt.ylabel("Frequency")
plt.grid(True)
plt.tight_layout()
plt.show()
