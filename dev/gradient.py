import typing as tp
import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(formatter={'float': '{: 0.5e}'.format})


def func(x: tp.Iterable) -> float:
    return (x[0] + 1.5) * (x[0] + 1.5) + (x[0] + 1.5) * x[1] + x[1] * x[1]

def show(f: tp.Callable[[tp.Iterable], float], ps: np.ndarray, x0: np.ndarray) -> None:
    
    n = 100
    l = 1.0
    
    xs = np.linspace(x0[0] - l, x0[0] + l, num=n)
    ys = np.linspace(x0[1] - l, x0[1] + l, num=n)

    X, Y = np.meshgrid(xs, ys)
    F = np.zeros_like(X)

    for i in range(n):
        for j in range(n):
            F[i, j] = f([X[i, j], Y[i, j]])
    
    x = [p[0] for p in ps]; x.append(ps[0, 0])
    y = [p[1] for p in ps]; y.append(ps[0, 1])

    plt.figure()
    plt.contourf(X, Y, F)
    plt.colorbar()
    plt.scatter(x0[0], x0[1], color='black')
    plt.plot(x, y, 'k')
    plt.axis('equal')
    plt.show()

    return

def num_grad(f: tp.Callable[[np.ndarray], float], x0: np.ndarray) -> np.ndarray:

    dx = 1e-12

    grad = np.array([
        (f(x0 + np.array([dx, 0.0])) - f(x0 - np.array([dx, 0.0]))) / (2 * dx),
        (f(x0 + np.array([0.0, dx])) - f(x0 - np.array([0.0, dx]))) / (2 * dx)
    ])

    return grad

def least_square_grad(f: tp.Callable[[tp.Iterable], float], xs: np.ndarray, x0: tp.Iterable) -> np.ndarray:

    n = xs.shape[0]
    lhs = np.empty((n, 2), dtype=np.double)
    rhs = np.empty(n, dtype=np.double)

    for i in range(n):
        lhs[i, :] = xs[i, :] - x0[:]
        rhs[i] = f(xs[i, :]) -f(x0[:])

    sol = np.linalg.lstsq(lhs, -rhs, rcond=None)

    return sol[0]

def subpanels_grad(f: tp.Callable[[tp.Iterable], float], xs: np.ndarray, x0: tp.Iterable) -> np.ndarray:
    
    f0 = f(x0)
    fs = np.asarray([f(x) for x in xs])

    n = np.zeros(3, dtype=np.double)

    imax = xs.shape[0]
    for i in range(imax):
        if i == imax - 1:
            v1 = np.array([xs[i, 0] - x0[0], xs[i, 1] - x0[1], fs[i] - f0])
            v2 = np.array([xs[0, 0] - x0[0], xs[0, 1] - x0[1], fs[0] - f0])
            n[:] = n[:] + np.cross(v1, v2)
        else:
            v1 = np.array([xs[i, 0] - x0[0], xs[i, 1] - x0[1], fs[i] - f0])
            v2 = np.array([xs[i + 1, 0] - x0[0], xs[i + 1, 1] - x0[1], fs[i + 1] - f0])
            n[:] = n[:] + np.cross(v1, v2)
    
    return np.array([n[0] / n[2], n[1] / n[2]])

def quad_panel(side) -> tp.List[np.ndarray]:
    
    factor = 0.5 * side

    p1 = np.array([ side * 1.01 + factor * np.random.random(1), side * 1.01 + factor * np.random.random(1)])
    p2 = np.array([-side * 1.01 - factor * np.random.random(1), side * 1.01 + factor * np.random.random(1)])
    p3 = np.array([-side * 1.01 - factor * np.random.random(1),-side * 1.01 - factor * np.random.random(1)])
    p4 = np.array([ side * 1.01 + factor * np.random.random(1),-side * 1.01 - factor * np.random.random(1)])

    xs = np.asarray([p1, p2, p3, p4]).reshape((4, 2))
    x0 = np.mean(xs, axis=0).reshape(2)

    x = [p[0] for p in xs]; x.append(xs[0, 0])
    y = [p[1] for p in xs]; y.append(xs[0, 1])

    return [x0, xs]

def tri_panel(side) -> tp.List[np.ndarray]:
    
    factor = 0.5 * side

    p1 = np.array([ side * 1.01 + factor * np.random.random(1), side * 1.01 + factor * np.random.random(1)])
    p2 = np.array([-side * 1.01 - factor * np.random.random(1), side * 1.01 + factor * np.random.random(1)])
    p3 = np.array([-side * 1.01 - factor * np.random.random(1),-side * 1.01 - factor * np.random.random(1)])

    xs = np.asarray([p1, p2, p3]).reshape((3, 2))
    x0 = np.mean(xs, axis=0).reshape(2)

    return [x0, xs]

if __name__ == '__main__':

    x0, xs = tri_panel(0.01)

    grad_1 = num_grad(func, x0)
    grad_2 = least_square_grad(func, xs, x0)
    grad_3 = subpanels_grad(func, xs, x0)

    print(grad_1)
    print(grad_2)
    print(grad_3)

    show(func, xs, x0)