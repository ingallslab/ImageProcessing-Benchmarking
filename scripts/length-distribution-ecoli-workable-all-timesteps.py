import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob


def draw_dist():
    """
    Generate and save histograms to compare the distribution of lengths across different data analysis tools.

    This function creates histograms for datasets corresponding to different tools, comparing the distribution of
    bacterial lengths. It calculates a common set of bin edges based on the data range and plots the distributions
    for each tool.
    """

    # Determine the common bin edges for all subplots based on the data range
    bin_edges = np.linspace(min(min(CP_Omnipose_length), min(DeLTA_length), min(FAST_length), min(SuperSegger_length)),
                            max(max(CP_Omnipose_length), max(DeLTA_length), max(FAST_length), max(SuperSegger_length)),
                            21)

    # Create the subplots
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))

    # Plot histograms for each tool
    axs[0, 0].hist(CP_Omnipose_length, bins=bin_edges, color='red', edgecolor='white', label="CP-Omnipose")
    axs[0, 1].hist(SuperSegger_length, bins=bin_edges, color='black', edgecolor='white', label="SS-Omnipose")
    axs[1, 0].hist(DeLTA_length, bins=bin_edges, color='blue', edgecolor='white', label='DeLTA')
    axs[1, 1].hist(FAST_length, bins=bin_edges, color='#3cb44b', edgecolor='white', label='FAST')

    # Set common labels for the subplots
    axs[0, 0].set_ylabel('Frequency', fontfamily='serif', fontsize=14, labelpad=10)
    axs[1, 0].set_ylabel('Frequency', fontfamily='serif', fontsize=14, labelpad=10)

    axs[1, 0].set_xlabel('length (um)', fontfamily='serif', fontsize=14)
    axs[1, 1].set_xlabel('length (um)', fontfamily='serif', fontsize=14)

    # Calculate the maximum and minimum y-axis limits from all subplots
    ymax_list = []
    ymin_list = []

    for ax in axs.flat:
        ymin, ymax = ax.get_ylim()
        ymax_list.append(ymax)
        ymin_list.append(ymin)

    max_ymax = int(max(ymax_list)) + 1
    min_ymin = int(min(ymin_list))

    # Standardize y-axis limits and x-tick labels across all subplots
    for ax in axs.flat:
        # Adjust the x-ticks and labels
        tick_positions = [(bin_edges[i] + bin_edges[i + 1]) / 2 for i in range(len(bin_edges) - 1)]
        ax.set_xticks(tick_positions)

        bin_labels = [round((bin_edges[i] + bin_edges[i + 1]) / 2) for i in range(len(bin_edges) - 1)]
        # Slice tick_positions and bin_labels to show every alternate bin
        tick_positions = tick_positions[::4] + [tick_positions[-1]]
        bin_labels = bin_labels[::4] + [bin_labels[-1]]

        ax.set_xticks(tick_positions)
        ax.set_xticklabels(bin_labels, rotation=0, fontsize=12, fontfamily='serif')
        ax.set_xlim([min(tick_positions) - 0.4, max(tick_positions) + 0.4])

        # Create a list of even y-ticks within the range
        ax.set_ylim([min_ymin, max_ymax])
        yticks = [i for i in range(min_ymin, max_ymax + 1) if i % 250 == 0]
        # Set the y-ticks to the even values
        ax.set_yticks(yticks, fontsize=12, fontfamily='serif')

    # Increase margin around the plot
    axs[0, 0].legend(prop={'size': 13, 'family': 'serif'})
    axs[0, 1].legend(prop={'size': 13, 'family': 'serif'})
    axs[1, 0].legend(prop={'size': 13, 'family': 'serif'})
    axs[1, 1].legend(prop={'size': 13, 'family': 'serif'})

    plt.tight_layout(pad=2.0)

    # Show the plot
    # plt.show()
    plt.savefig('distribution of length - ecoli babies.pdf')


if __name__ == '__main__':
    babies_folder = 'G:/ecoli babies/done'
    tools = ['CP_Omnipose', 'DeLTA', 'FAST', 'SuperSegger']

    CP_Omnipose_length = []
    DeLTA_length = []
    FAST_length = []
    SuperSegger_length = []

    for tool in tools:
        for folder in glob.glob(babies_folder + '/*'):
            if tool == 'CP_Omnipose':
                bacteria_properties_path = folder + '/' + tool + '/processing/CP_Omnipose_bacteria_feature_analysis.csv'
                bacteria_properties_df = pd.read_csv(bacteria_properties_path)
                CP_Omnipose_length.extend([v * 0.144 for v in bacteria_properties_df['Major_axis'].values.tolist()])
            elif tool == 'DeLTA':
                bacteria_properties_path = folder + '/' + tool + '/results/DeLTA_bacteria_feature_analysis.csv'
                bacteria_properties_df = pd.read_csv(bacteria_properties_path)
                DeLTA_length.extend([v * 0.144 for v in bacteria_properties_df['Major_axis'].values.tolist()])
            elif tool == 'FAST':
                bacteria_properties_path = folder + '/' + tool + '/results/FAST_bacteria_feature_analysis.csv'
                bacteria_properties_df = pd.read_csv(bacteria_properties_path)
                FAST_length.extend([v * 0.144 for v in bacteria_properties_df['Major_axis'].values.tolist()])
            elif tool == 'SuperSegger':
                bacteria_properties_path = folder + '/' + tool + '/results/SuperSegger_bacteria_feature_analysis.csv'
                bacteria_properties_df = pd.read_csv(bacteria_properties_path)
                SuperSegger_length.extend([v * 0.144 for v in bacteria_properties_df['Major_axis'].values.tolist()])

    # filtration
    CP_Omnipose_length = [v for v in CP_Omnipose_length if 0.68 < v < 13]
    DeLTA_length = [v for v in DeLTA_length if 0.68 < v < 13]
    FAST_length = [v for v in FAST_length if 0.68 < v < 13]
    SuperSegger_length = [v for v in SuperSegger_length if 0.68 < v < 13]

    df = pd.DataFrame.from_dict({
        'CP_Omnipose': CP_Omnipose_length,
        'DeLTA': DeLTA_length,
        'FAST': FAST_length,
        'SuperSegger': SuperSegger_length
    }, orient='index').transpose()

    # Compute bins and histograms for each column
    bin_edges = np.linspace(df.min().min(), df.max().max(), num=21)

    hist_data = {
        'bin_start': bin_edges[:-1],  # Start of each bin
        'bin_end': bin_edges[1:]  # End of each bin
    }

    for col in df.columns:
        freq, _ = np.histogram(df[col].dropna(), bins=bin_edges)
        hist_data[col] = freq

    # Convert histogram data to dataframe
    hist_df = pd.DataFrame(hist_data)

    print(df.mean())
    print(df.std())
    draw_dist()
