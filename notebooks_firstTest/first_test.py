import uproot



# Define the relative or absolute path to the ROOT file
# Example: '../data/myfile.root' or '/home/user/data/myfile.root'
root_file_path = "/data/cgarcia_2002/WCTE/data/WCTE_offline_R2385S0P0.root"

# Open the ROOT file
with uproot.open(root_file_path) as file:
    print(file)
    print(file.classnames())

    tree_v2 = file["WCTEReadoutWindows;3"]
    tree_v3 = file["WCTEReadoutWindows;2"]

    print("tree 2", tree_v2.num_entries)  # Or tree_v3.keys()
    print("tree 3", tree_v3.num_entries)  # Or tree_v3.keys()
    print(tree_v2.keys())

    print(tree_v2["window_time"].array(library="np"))

    #print("Number of events in WCTEReadoutWindows;2:", tree_v2.num_entries)
    #print("Number of events in WCTEReadoutWindows;3:", tree_v3.num_entries)
"""
    print("\nDetalles de branches:")
    for branch_name, branch_obj in tree_v2.items():
        print(f"Branch: {branch_name}")
        print(f"  Tipo: {branch_obj.classname}")


"""