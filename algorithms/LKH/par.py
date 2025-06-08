import os
import glob

def write_par_file(tsp_file, output_folder):
    base = os.path.splitext(os.path.basename(tsp_file))[0]
    par_file = os.path.join(output_folder, base + ".par")
    with open(par_file, 'w') as f:
        f.write(f"PROBLEM_FILE = {tsp_file}\n")
        f.write(f"OUTPUT_TOUR_FILE = {os.path.join(output_folder, base + '.tour')}\n")
        f.write("MOVE_TYPE = 5\n")
        f.write("PATCHING_C = 3\n")
        f.write("PATCHING_A = 2\n")
        f.write("RUNS = 1\n")
        f.write("SEED = 42\n")
    print(f"Generated {par_file}")

def main(tsp_folder, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    tsp_files = glob.glob(os.path.join(tsp_folder, "*.tsp"))
    for tsp_file in tsp_files:
        write_par_file(tsp_file, output_folder)

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python generate_par_files.py <tsp_folder> <par_output_folder>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
