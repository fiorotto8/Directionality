import ROOT
import argparse

def split_tree(input_filename, tree_name, n_subsets=4):
    # Open the original ROOT file
    original_file = ROOT.TFile.Open(input_filename, "READ")
    if not original_file or original_file.IsZombie():
        print(f"Error: Cannot open file {input_filename}")
        return

    original_tree = original_file.Get(tree_name)
    if not original_tree:
        print(f"Error: Cannot find tree {tree_name} in file {input_filename}")
        return
    
    # Determine the number of entries for each subset
    total_entries = original_tree.GetEntries()
    entries_per_subset = total_entries // n_subsets
    
    for i in range(n_subsets):
        # Create a new ROOT file for each subset
        subset_filename = f"subset_{i+1}.root"
        subset_file = ROOT.TFile(subset_filename, "RECREATE")
        
        # Clone the original tree structure
        subset_tree = original_tree.CloneTree(0)
        
        # Define the entry range for the current subset
        start_entry = i * entries_per_subset
        end_entry = start_entry + entries_per_subset if i < n_subsets - 1 else total_entries
        
        # Copy the entries for the current subset
        for j in range(start_entry, end_entry):
            original_tree.GetEntry(j)
            subset_tree.Fill()
        
        # Write the subset tree to the new file
        subset_tree.AutoSave()
        subset_file.Close()
        print(f"Created {subset_filename} with entries from {start_entry} to {end_entry-1}")

    original_file.Close()


def main():
    parser = argparse.ArgumentParser(description="Split a ROOT TTree into multiple subsets.")
    parser.add_argument("input_file", help="The input ROOT file containing the TTree.")
    parser.add_argument("--tree_name", help="The name of the TTree to split.", default="Events")
    parser.add_argument("--n_subsets", type=int, default=4, help="Number of subsets to split the TTree into. Default is 4.")
    
    args = parser.parse_args()
    
    split_tree(args.input_file, args.tree_name, args.n_subsets)

if __name__ == "__main__":
    main()
