import numpy as np
from itertools import takewhile, dropwhile


def read_slurm_file(file):
    """Reads slurm output file and skips initial newlines and the info text at the bottom."""
    with open(file) as f:
        lines = dropwhile(lambda line: line == '\n', f)
        lines = takewhile(lambda line: line != '\n', lines)
        lines = (line.strip().split(' ') for line in lines)
        lines = list(lines)
        return np.array(lines, dtype=float)


def read(file1, file2):
    data1 = read_slurm_file(file1)
    data2 = read_slurm_file(file2)

    return data1, data2


def draw(data, ax, filename):
    x, y = data[:, 0], data[:, 1]
    ax.scatter(x, y)
    ax.set_yscale('log')
    # ax.set_xscale('log')
    ax.set_xlabel('exp(log_k(i)/n_power(i))')
    ax.set_ylabel('power(i)/n_power(i)')
    ax.title.set_text(filename)



if __name__ == "__main__":
    from sys import argv
    file1, file2 = argv[1:3]
    data1, data2 = read(file1, file2)
    print(file1)
    print(data1)
    print()
    print(file2)
    print(data2)
    print()
    print(f'{(data1==data2).all() = }')
    if not (data1==data2).all():
        print(f'{(data1==data2) = }')
        print()
        print(file1 + ' / ' + file2)
        print(data1/data2)

    try:
        if argv[3] == 'd':
            import matplotlib.pyplot as plt
            fig, (ax1, ax2) = plt.subplots(1, 2)
            draw(data1, ax1, file1)
            draw(data2, ax2, file2)
            plt.show()
    except IndexError:
        print("Not drawing.")
    except ModuleNotFoundError:
        print("Couldn't draw, matplotlib not available.")
