import numpy as np
from itertools import takewhile, dropwhile


def read_file(file):
    with open(file) as f:
        lines = dropwhile(lambda line: line == '\n', f)
        lines = takewhile(lambda line: line != '\n', lines)
        lines = (line.strip().split(' ') for line in lines)
        lines = list(lines)
        return np.array(list(lines), dtype=float)


def read(file1, file2):
    data1 = read_file(file1)
    data2 = read_file(file2)

    return data1, data2


def draw(data, ax):
    x, y = data[:, 0], data[:, 1]
    ax.scatter(x, y)
    ax.set_yscale('log')
    # ax.set_xscale('log')
    ax.set_xlabel('exp(log_k(i)/n_power(i))')
    ax.set_ylabel('power(i)/n_power(i)')


def draw_relative(data1, data2, ax):
    assert((data1[:, 0] == data2[:, 0]).all())
    x = data1[:, 0]
    y1 = data1[:, 1]
    y2 = data2[:, 1]

    ax.scatter(x, y1/y2)
    ax.set_xlabel('exp(log_k(i)/n_power(i))')
    ax.set_ylabel('relative power')


if __name__ == "__main__":
    from sys import argv
    file1, file2 = argv[1:3]
    data1, data2 = read(file1, file2)
    print(data1/data2)

    try:
        if argv[3] == 'd':
            import matplotlib.pyplot as plt
            fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
            draw(data1, ax1)
            draw(data2, ax2)
            draw_relative(data1, data2, ax3)
            plt.show()
    except IndexError:
        print("Not drawing.")
    except ModuleNotFoundError:
        print("Couldn't draw, matplotlib not available.")