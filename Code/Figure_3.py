import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Verify current working directory
print("Current Working Directory:", os.getcwd())

# Path to the CSV file
file_path = "/Users/christopherbivins/Desktop/Gall_Project_Neater/Alpha_Diversity/richness_shannon_diversity.csv"

# Load the dataset
data = pd.read_csv(file_path)

# Define the color palette and label mapping
palette = {
    'Leaf_galled_by_Andricus': '#FDB462',  # Pastel orange
    'Andricus_Gall': '#FDB462',           # Same pastel orange
    'Leaf_galled_by_Cynips': '#FCCDE5',   # Same pastel pink as Cynips_Gall
    'Cynips_Gall': '#FCCDE5',             # Pastel pink
    'Ungalled_leaf': '#80B1D3'            # Pastel blue
}

# Shortened labels for the x-axis
label_mapping = {
    'Leaf_galled_by_Andricus': 'FGL',
    'Andricus_Gall': 'FGG',
    'Leaf_galled_by_Cynips': 'CQL',
    'Cynips_Gall': 'CQG',
    'Ungalled_leaf': 'UGL'
}

# Function to create and save a jittered strip plot
def create_strip_plot(data, y_column, title, output_filename):
    plt.figure(figsize=(12, 6))
    sns.stripplot(
        x='Treatment',
        y=y_column,
        data=data,
        jitter=True,
        size=8,
        palette=palette,
        marker='o',  # Use circles as markers
        edgecolor='gray',
        linewidth=0.5
    )
    # Customize the plot
    plt.title(title, fontsize=16)
    plt.xlabel('Treatment', fontsize=14)
    plt.ylabel(y_column, fontsize=14)
    plt.xticks(
        ticks=range(len(label_mapping)),
        labels=[label_mapping[t] for t in data['Treatment'].unique()],
        rotation=0,  # No rotation needed for short labels
        fontsize=12
    )
    plt.ylim(0, None)  # Ensure the y-axis starts at 0
    
    # Remove the top and right spines
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.grid(False)  # Remove gridlines
    plt.tight_layout()
    plt.savefig(output_filename, format='svg')
    plt.close()

# Create and save both plots
create_strip_plot(
    data=data,
    y_column='Shannon',
    title='Jittered Strip Plot of Shannon Diversity by Treatment',
    output_filename='jittered_strip_plot_shannon.svg'
)

create_strip_plot(
    data=data,
    y_column='Richness',
    title='Jittered Strip Plot of Richness by Treatment',
    output_filename='jittered_strip_plot_richness.svg'
)

print("Plots saved as 'jittered_strip_plot_shannon.svg' and 'jittered_strip_plot_richness.svg'")

