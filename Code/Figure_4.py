import pandas as pd
import matplotlib.pyplot as plt

# Path to the NMDS data file
nmds_data_path = "/Users/christopherbivins/Desktop/Gall_Project_Neater/Beta_Diversity/nmds_bray_df.csv"

# Load the NMDS data
nmds_df = pd.read_csv(nmds_data_path)

# Define the color and marker mapping
color_mapping = {
    'Leaf_galled_by_Andricus': ('#FDB462', '+'),
    'Andricus_Gall': ('#FDB462', 'o'),
    'Leaf_galled_by_Cynips': ('#FCCDE5', '+'),
    'Cynips_Gall': ('#FCCDE5', 'o'),
    'Ungalled_leaf': ('#80B1D3', '_')
}

# Create the NMDS plot
plt.figure(figsize=(10, 8))

for group, (hex_color, marker) in color_mapping.items():
    subset = nmds_df[nmds_df['Treatment'] == group]
    size = 300 if marker in ['+', '_'] else 100  # Larger size for plus and minus
    plt.scatter(
        subset['NMDS1'], subset['NMDS2'],
        label=group.replace('_', ' '), c=hex_color, marker=marker, s=size, linewidth=2
    )

# Customize the plot
plt.title("NMDS Ordination Plot - Bray-Curtis Distance Matrix")
plt.xlabel("NMDS1")
plt.ylabel("NMDS2")
plt.legend(loc='upper right', title='Sample Groups')

# Remove all gridlines and central axes
ax = plt.gca()
ax.grid(False)
ax.axhline(0, color='none', linewidth=0)
ax.axvline(0, color='none', linewidth=0)
ax.set_frame_on(False)

# Save the plot as an SVG file
plt.tight_layout()
plt.savefig("nmds_ordination_plot_final.svg", format="svg")
plt.close()

print("Plot saved as 'nmds_ordination_plot_final.svg' in the current directory.")

