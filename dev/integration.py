import matplotlib.pyplot as plt
import numpy as np


def quad():

    n = 100
    x, y = np.linspace(0, 1, num=n), np.linspace(0, 1, num=n)

    X, Y = np.meshgrid(x, y)

    f1 = (1 - X) * (1 - Y)
    f2 = X * (1 - Y)
    f3 = (1 - X) * Y
    f4 = X * Y

    f_sum = f1 + f2 + f3 + f4

    plt.figure()
    plt.title('f1 [0, 0]')
    plt.imshow(f1)
    plt.colorbar()
    plt.axis('equal')

    plt.figure()
    plt.title('f2 [1, 0]')
    plt.imshow(f2)
    plt.colorbar()
    plt.axis('equal')

    plt.figure()
    plt.title('f3 [0, 1]')
    plt.imshow(f3)
    plt.colorbar()
    plt.axis('equal')

    plt.figure()
    plt.title('f4 [1, 1]')
    plt.imshow(f4)
    plt.colorbar()
    plt.axis('equal')

    plt.figure()
    plt.title('f_sum')
    plt.imshow(f_sum)
    plt.colorbar()
    plt.axis('equal')

    plt.show()

    return

if __name__ == '__main__':

    n = 100
    x, y = np.linspace(0, 1, num=n), np.linspace(0, 1, num=n)

    X, Y = np.meshgrid(x, y)

    f1 = (1 - X) * (1 - Y)
    f2 = X * (1 - np.abs(0.5 - Y))
    f3 = (1 - X) * Y

    f_sum = f1 + f2 + f3

    f1 = f1 / f_sum
    f2 = f2 / f_sum
    f3 = f3 / f_sum

    f_sum = f1 + f2 + f3

    plt.figure()
    plt.title('f1 [0, 0]')
    plt.contourf(f1)
    plt.colorbar()
    plt.axis('equal')

    plt.figure()
    plt.title('f2 [1, 0]')
    plt.contourf(f2)
    plt.colorbar()
    plt.axis('equal')

    plt.figure()
    plt.title('f3 [1, 1]')
    plt.contourf(f3)
    plt.colorbar()
    plt.axis('equal')

    plt.figure()
    plt.title('f_sum')
    plt.contourf(f_sum)
    plt.colorbar()
    plt.axis('equal')

    plt.show()