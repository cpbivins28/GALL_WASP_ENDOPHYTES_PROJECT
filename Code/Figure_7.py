import pandas as pd
import matplotlib.pyplot as plt

def generate_otu4_plot():
    # File paths
    presence_absence_path = '/Users/christopherbivins/Desktop/Gall_Project_Neater/Presence_Absence/OTU_presence_absence_table.csv'
    metadata_path = '/Users/christopherbivins/Desktop/Gall_Project_Neater/Relative_Abundance/Gall_project_metadata.csv'

    # Load the data
    presence_absence = pd.read_csv(presence_absence_path)
    metadata = pd.read_csv(metadata_path)

    # Reshape the presence-absence data to long format
    presence_absence_long = presence_absence.melt(
        id_vars='OTU_ID', 
        var_name='SampleID', 
        value_name='Presence'
    )

    # Filter for OTU4 data
    otu4_data_long = presence_absence_long[presence_absence_long['OTU_ID'] == 'OTU4']

    # Merge with metadata to get treatment information
    otu4_with_metadata = pd.merge(otu4_data_long, metadata, left_on='SampleID', right_on='Sample_ID', how='inner')

    # Filter for samples where OTU4 is present (assuming presence is marked as 1)
    otu4_present = otu4_with_metadata[otu4_with_metadata['Presence'] == 1]

    # Count occurrences of OTU4 in each Treatment group
    otu4_counts_by_treatment = otu4_present['Treatment'].value_counts()

    # Calculate the total number of samples for each treatment group
    total_samples_by_treatment = metadata['Treatment'].value_counts()

    # Calculate the proportion of samples in each treatment group with OTU4
    ordered_treatments = [
        'Andricus_Gall',
        'Leaf_galled_by_Andricus',
        'Cynips_Gall',
        'Leaf_galled_by_Cynips',
        'Ungalled_leaf'
    ]
    otu4_proportions = (otu4_counts_by_treatment / total_samples_by_treatment).reindex(ordered_treatments)

    # Define custom colors for each treatment group
    custom_colors = {
        'Ungalled_leaf': '#80B1D3',
        'Leaf_galled_by_Cynips': '#FCCDE5',
        'Leaf_galled_by_Andricus': '#FDB462',
        'Andricus_Gall': '#FDB462',
        'Cynips_Gall': '#FCCDE5'
    }

    # Create the horizontal bar plot
    plt.figure(figsize=(8, 6))
    otu4_proportions.plot(
        kind='barh', 
        color=[custom_colors[treatment] for treatment in ordered_treatments]
    )
    plt.xlabel('Proportion of Samples')
    plt.ylabel('Treatment Group')
    plt.title('Proportion of Samples with OTU4 in Each Treatment Group (Reordered)')
    plt.grid(False)  # Remove gridlines
    plt.tight_layout()

    # Save the plot as an SVG file
    plt.savefig('OTU4_Proportions_Plot.svg', format='svg')
    print("Plot saved as 'OTU4_Proportions_Plot.svg'")

if __name__ == "__main__":
    generate_otu4_plot()

