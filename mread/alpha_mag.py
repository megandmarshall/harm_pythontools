### This code finds all of the magnetic stresses \alpha that are present in the simulation data.
### Then, it plots the magnetic stresses over time
### Author: Brett Hajdaj '27

import matplotlib.pyplot as plt
import numpy as np
import os, sys
import time

# Specific directories for each user
MREAD_DIR = r"C:\Users\brhaj\harm_pythontools\mread"
DATA_DIR = r"D:\MADNT"

sys.path.insert(0, MREAD_DIR)
os.chdir(DATA_DIR)
import __init_student__ as ist

# Bin file start/end numbers
FILE_START = 7371   
FILE_END = 7753     
FILE_STEP = 1

# Main processing loop to read the data files and collect the results
results = []
for file_num in range(FILE_START, FILE_END + 1, FILE_STEP):
    result = ist.stressdecompvtime(file_num)
    results.append(result)
    print(f"Processed file {file_num} of {FILE_END}")
print(f"Collected {len(results)} results in total.")
records = np.array(results)

# Time = file # * 4
t = records[:, 0] * 4  

# Plotting configuration of the magnetic stresses over time
np.savetxt("alpha_mag_vs_time.csv", records, delimiter=",",
           header="fnumber,alpha_tot,alpha_mean,alpha_cross1,alpha_cross2,alpha_pert,"
                  "alpha_bub_tot,alpha_bub_cross1,alpha_bub_cross2,alpha_bub_pert,"
                  "alpha_disk_tot,alpha_disk_cross1,alpha_disk_cross2,alpha_disk_pert",
           comments="")
plt.style.use("dark_background")
plt.figure()
plt.plot(t, records[:, 1], label="total")
plt.plot(t, records[:, 2], label="mean field (MAD)")
plt.plot(t, records[:, 5], label="turbulent (MRI)")
plt.plot(t, records[:, 6], label="bubble/corona, total")
plt.plot(t, records[:, 10], label="disk, total")
plt.xlabel("time", fontsize=12)
plt.ylabel(r"$\alpha_{mag}$", fontsize=14)
plt.title("Magnetic stress vs. time", fontsize=14)
plt.legend(frameon=False)
plt.grid(alpha=0.2)
plt.tight_layout()
plt.savefig("alpha_mag_vs_time.png")
print("Saved alpha_mag_vs_time.csv and alpha_mag_vs_time.png")

