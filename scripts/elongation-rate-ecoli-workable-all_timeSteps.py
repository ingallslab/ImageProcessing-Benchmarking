import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob


# Transformation function
def transform(y):
    return (scale_y_upper_than + 1) + np.sqrt(y - scale_y_upper_than) if y > scale_y_upper_than else y


# Inverse transformation for y-tick labels
def inverse_transform(y):
    return (y - (scale_y_upper_than + 1)) ** 2 + scale_y_upper_than if y > scale_y_upper_than else y


def draw_dist():
    """
    Generate and save histograms to compare the distribution of elongation rates across different data analysis tools.

    This function creates histograms for each tool's dataset to compare the distribution of bacterial elongation
    rates. It calculates a common set of bin edges based on the data range and plots the distributions for each
    tool.
    """

    # Determine the common bin edges for all subplots based on the data range
    bin_edges = np.linspace(min(min(CP_Omnipose_growth_rate), min(DeLTA_growth_rate), min(FAST_growth_rate),
                                min(SuperSegger_growth_rate)),
                            max(max(CP_Omnipose_growth_rate), max(DeLTA_growth_rate), max(FAST_growth_rate),
                                max(SuperSegger_growth_rate)), 31)

    # Create the subplots
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))

    # Setup data and parameters for each tool
    axes = [axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1]]
    data = [CP_Omnipose_growth_rate, SuperSegger_growth_rate, DeLTA_growth_rate, FAST_growth_rate]
    tool_names = ['CP-Omnipose', 'SS-Omnipose', 'DeLTA', 'FAST']
    colors = ['red', 'black', 'blue', '#3cb44b']

    # Variable to store the maximum height after transformation
    max_transformed_height = 0

    # Plot histograms for each dataset
    for ax, d, color, tool_name in zip(axes, data, colors, tool_names):
        counts, _, patches = ax.hist(d, bins=bin_edges, color=color, edgecolor='white', label=tool_name)

        # Adjust the heights of the bars
        for patch in patches:
            height = patch.get_height()
            new_height = transform(height)
            patch.set_height(new_height)
            max_transformed_height = max(max_transformed_height, new_height)

    # Set common labels and titles
    axs[0, 0].set_ylabel('Frequency', fontfamily='serif', fontsize=14, labelpad=10)
    axs[1, 0].set_ylabel('Frequency', fontfamily='serif', fontsize=14, labelpad=10)

    axs[1, 0].set_xlabel('elongation rate (um/min)', fontfamily='serif', fontsize=14,
                         labelpad=10)  # Increase space between x label and x axis
    axs[1, 1].set_xlabel('elongation rate (um/min)', fontfamily='serif', fontsize=14, labelpad=10)

    # Standardize axes labels and limits
    for ax in axs.flat:
        # Set x-ticks and labels
        tick_positions = [(bin_edges[i] + bin_edges[i + 1]) / 2 for i in range(len(bin_edges) - 1)]
        bin_labels = [round((bin_edges[i] + bin_edges[i + 1]) / 2, 3) for i in range(len(bin_edges) - 1)]
        # Slice tick_positions and bin_labels to show every alternate bin
        tick_positions = tick_positions[::6] + [tick_positions[-1]]
        bin_labels = bin_labels[::6] + [bin_labels[-1]]

        ax.set_xticks(tick_positions)
        ax.set_xticklabels(bin_labels, rotation=0, fontsize=12, fontfamily='serif')
        ax.set_xlim([min(tick_positions) - 0.08, max(tick_positions) + 0.08])

        # Set y-ticks and limits based on transformed data
        transformed_yticks = np.linspace(0, (np.round(max_transformed_height // 4.95) + 1) * 4.95,
                                         (int(np.round(max_transformed_height // 12)) + 1))

        original_yticks = [inverse_transform(y) for y in transformed_yticks]

        ax.set_yticks(transformed_yticks)
        ax.set_yticklabels([f"{int(y)}" for y in original_yticks], fontsize=12, fontfamily='serif')

        # Get the current y-axis limits
        ymin, ymax = ax.get_ylim()

        ax.set_ylim([ymin, (np.round(max(transformed_yticks) // 4) + 1) * 4])

    # Increase margin around the plot
    axs[0, 0].legend(prop={'size': 13, 'family': 'serif'})
    axs[0, 1].legend(prop={'size': 13, 'family': 'serif'})
    axs[1, 0].legend(prop={'size': 13, 'family': 'serif'})
    axs[1, 1].legend(prop={'size': 13, 'family': 'serif'})

    # Increase margin around the plot
    plt.tight_layout(pad=2.0)

    # Show the plot
    # plt.show()
    plt.savefig('distribution of elongation rate - ecoli babies.pdf')


if __name__ == '__main__':
    babies_folder = 'G:/ecoli babies/done'
    tools = ['CP_Omnipose', 'DeLTA', 'FAST', 'SuperSegger']
    num_bin = 30
    scale_y_upper_than = 50

    CP_Omnipose_growth_rate = []
    DeLTA_growth_rate = []
    FAST_growth_rate = []
    SuperSegger_growth_rate = []

    for tool in tools:
        for folder in glob.glob(babies_folder + '/*'):
            if tool == 'CP_Omnipose':
                bacteria_properties_path = folder + '/' + tool + '/processing/CP_Omnipose_LifeHistory_based_Analysis.csv'
                bacteria_properties_df = pd.read_csv(bacteria_properties_path)
                CP_Omnipose_growth_rate.extend(bacteria_properties_df['GrowthRate'].values.tolist())
            elif tool == 'DeLTA':
                bacteria_properties_path = folder + '/' + tool + '/results/DeLTA_LifeHistory_based_Analysis.csv'
                bacteria_properties_df = pd.read_csv(bacteria_properties_path)
                DeLTA_growth_rate.extend(bacteria_properties_df['GrowthRate'].values.tolist())
            elif tool == 'FAST':
                bacteria_properties_path = folder + '/' + tool + '/results/FAST_LifeHistory_based_Analysis.csv'
                bacteria_properties_df = pd.read_csv(bacteria_properties_path)
                FAST_growth_rate.extend(bacteria_properties_df['GrowthRate'].values.tolist())
            elif tool == 'SuperSegger':
                bacteria_properties_path = folder + '/' + tool + '/results/SuperSegger_LifeHistory_based_Analysis.csv'
                bacteria_properties_df = pd.read_csv(bacteria_properties_path)
                SuperSegger_growth_rate.extend(bacteria_properties_df['GrowthRate'].values.tolist())

    # filtration (remove zeros and nan)
    CP_Omnipose_growth_rate = [v for v in CP_Omnipose_growth_rate if v != 0 and v != 0.0 and str(v) != 'nan'
                               and 2.551672 >= v >= -1.042678]
    DeLTA_growth_rate = [v for v in DeLTA_growth_rate if v != 0 and v != 0.0 and str(v) != 'nan'
                         and 2.551672 >= v >= -1.042678]
    FAST_growth_rate = [v for v in FAST_growth_rate if v != 0 and v != 0.0 and str(v) != 'nan'
                        and 2.551672 >= v >= -1.042678]
    SuperSegger_growth_rate = [v for v in SuperSegger_growth_rate if v != 0 and v != 0.0 and str(v) != 'nan'
                               and 2.551672 >= v >= -1.042678]

    df = pd.DataFrame.from_dict({
        'CP_Omnipose': CP_Omnipose_growth_rate,
        'DeLTA': DeLTA_growth_rate,
        'FAST': FAST_growth_rate,
        'SuperSegger': SuperSegger_growth_rate
    }, orient='index').transpose()

    # Compute bins and histograms for each column
    bin_edges = np.linspace(df.min().min(), df.max().max(), num=31)

    hist_data = {
        'bin_start': bin_edges[:-1],  # Start of each bin
        'bin_end': bin_edges[1:]  # End of each bin
    }

    for col in df.columns:
        freq, _ = np.histogram(df[col].dropna(), bins=bin_edges)
        hist_data[col] = freq

    # Convert histogram data to dataframe
    hist_df = pd.DataFrame(hist_data)

    print(hist_df)
    print(df.mean())
    print(df.std())

    draw_dist()
