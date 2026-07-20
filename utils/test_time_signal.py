import argparse
import random
import numpy as np
import pandas as pd
from ete3 import Tree
from scipy.stats import linregress
import matplotlib.pyplot as plt


# Function to extract root-to-tip distances
def get_root_to_tip_distances(tree, exclude=[]):
    distances = {}
    for leaf in tree.iter_leaves():
        if leaf.name in exclude:
            continue
        distances[leaf.name] = tree.get_distance(leaf)
    return distances


# Function to compute R² correlation
def compute_r2(distances, dates):
    common_samples = list(set(distances.keys()) & set(dates.keys()))
    if len(common_samples) < 2:
        return None  # Not enough data
    
    x = np.array([dates[sample] for sample in common_samples])  # Decimal years
    y = np.array([distances[sample] for sample in common_samples])
    
    slope, intercept, r_value, _, _ = linregress(x, y)
    return r_value**2


# Function to perform randomization test
def randomization_test(tree, meta, fname_plot=None, num_shuffles=100_000, exclude=[]):
    distances = get_root_to_tip_distances(tree, exclude=exclude)
    dates = {
        seqid: meta.loc[seqid, 'Date']
        for seqid in distances.keys()
    }

    if len(dates) < 2:
        print("Not enough valid dates for analysis.")
        return
    
    observed_r2 = compute_r2(distances, dates)
    
    if observed_r2 is None:
        print("Could not compute observed R² (insufficient data).")
        return

    randomized_r2_values = []
    sample_names = list(dates.keys())

    for _ in range(num_shuffles):
        random.shuffle(sample_names)
        shuffled_dates = dict(zip(sample_names, dates.values()))
        r2 = compute_r2(distances, shuffled_dates)
        if r2 is not None:
            randomized_r2_values.append(r2)

    # Compute p-value
    p_value = sum(r >= observed_r2 for r in randomized_r2_values) / num_shuffles

    print(f"Observed R²: {observed_r2:.4f}")
    print(f"P-value: {p_value:.4e} (Fraction of randomized R² values ≥ observed R²)")

    # Plot histogram of randomized R² values
    plt.hist(randomized_r2_values, bins=30, alpha=0.7, label="Randomized R²")
    plt.axvline(observed_r2, color="red", linestyle="dashed", linewidth=2, label="Observed R²")
    plt.xlabel("R² value")
    plt.ylabel("Frequency")
    plt.legend()
    plt.title("Randomization Test for Temporal Signal")

    if fname_plot:
        plt.savefig(fname_plot)
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tree', required=True, help='ML tree.')
    parser.add_argument('--meta', required=True, help='CSV file with metadata include .')
    parser.add_argument('--plot_out', default=None, help='Filename to write the plot to.')
    parser.add_argument('--num_shuffles', default=100_000, help='Number of samples drawn.')
    parser.add_argument('--exclude', default=None, help='Filename with accessions to exclude. One accession per line.')
    args = parser.parse_args()

    # Load tree and run the test
    tree = Tree(args.tree, format=1)  # Format 1 preserves branch lengths
    meta = pd.read_csv(args.meta, sep='\t').set_index('Accession')

    if args.exclude:
        with open(args.exclude) as f:
            exlude_lst = [line.strip() for line in f]
    else:
      exlude_lst = []
      
    randomization_test(tree, meta, args.plot_out, args.num_shuffles, exclude=exlude_lst)


if __name__ == '__main__':
    main()
