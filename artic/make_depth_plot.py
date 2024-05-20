import matplotlib.pyplot as plt
import argparse


def create_depth_plot(depth_file, output_file):
    # Read in the depth file
    depths = []
    with open(depth_file, "r") as infile:
        for line in infile:
            if line:
                depths.append(line.strip().split("\t"))

    # Plot the depth
    fig, ax = plt.subplots()
    ax.plot([int(x[1]) for x in depths], [int(x[2]) for x in depths])
    ax.set_yscale("log")
    # Add mindepth line
    ax.axhline(y=20, color="r", linestyle="-")

    # Labels
    ax.set_ylabel("Read depth")
    ax.set_xlabel("Position")
    ax.set_ylim(1, max([int(x[2]) for x in depths]) * 1.1)

    ax.figure.set_size_inches(16 / 2, 9 / 2)

    # Save plot
    plt.savefig(output_file + ".png")


def main():
    parser = argparse.ArgumentParser(
        description="Plot the depth of coverage from a depth file"
    )
    parser.add_argument("--depth", help="Depth file", required=True, type=str)
    parser.add_argument("--min-depth", help="Depth file", required=True, type=int)
    parser.add_argument("--output", help="Output file", required=True, type=str)

    args = parser.parse_args()
    create_depth_plot(args.depth, args.output)


if __name__ == "__main__":
    main()
