import numpy as np

def generate_nonlinear_system(n, seed, mode):
    np.random.seed(seed)

    n = n + 100

    if mode == "l-uni":
        # n = 20
        X = np.zeros((n, 1))
        Y = np.zeros((n, 1))
        X[0] = np.random.rand()
        Y[0] = np.random.rand()
        gamma_yx = 0.32
        r_x = 3.7
        r_y = 3.8
        for i in range(1, n):
            X[i] = X[i - 1] * (r_x - r_x * X[i - 1])
            Y[i] = Y[i - 1] * (r_y - r_y * Y[i - 1] - gamma_yx * X[i - 1])

    elif mode == "l-bi":
        # n = 20
        X = np.zeros((n, 1))
        Y = np.zeros((n, 1))
        X[0] = np.random.rand()
        Y[0] = np.random.rand()
        gamma_xy = 0.02
        gamma_yx = 0.1
        r_x = 3.7
        r_y = 3.8
        for i in range(1, n):
            X[i] = X[i - 1] * (r_x - r_x * X[i - 1] - gamma_xy * Y[i - 1])
            Y[i] = Y[i - 1] * (r_y - r_y * Y[i - 1] - gamma_yx * X[i - 1])

    elif mode == "RL":
        # n = 50
        alpha = 6
        C = 2
        dt = 0.01
        X = np.zeros((n, 3))
        Y = np.zeros((n, 3))
        X[0, :] = (np.random.rand(1, 3) * 2 - 1) * 10
        Y[0, :] = (np.random.rand(1, 3) * 2 - 1) * 10

        cause_start = 1
        for i in range(1, n):
            X[i, 0] = X[i-1, 0] - alpha * (X[i-1, 1] + X[i-1, 2]) * dt
            X[i, 1] = X[i-1, 1] + alpha * (X[i-1, 0] + 0.2 * X[i-1, 1]) * dt
            X[i, 2] = X[i-1, 2] + alpha * (0.2 + X[i-1, 2] * (X[i-1, 0] - 5.7)) * dt
            Y[i, 0] = Y[i-1, 0] + 10 * (-Y[i-1, 0] + Y[i-1, 1]) * dt
            if i > cause_start:
                Y[i, 1] = Y[i-1, 1] + (28 * Y[i-1, 0] - Y[i-1, 1] - Y[i-1, 0] * Y[i-1, 2] + C * X[i-cause_start, 1]**2) * dt
            else:
                Y[i, 1] = Y[i-1, 1] + (28 * Y[i-1, 0] - Y[i-1, 1] - Y[i-1, 0] * Y[i-1, 2]) * dt
            Y[i, 2] = Y[i-1, 2] + (Y[i-1, 0] * Y[i-1, 1] - 8/3 * Y[i-1, 2]) * dt

            X[i, :] = X[i, :] + (np.random.rand(1, 3) - 0.5) * 0.1
            Y[i, :] = Y[i, :] + (np.random.rand(1, 3) - 0.5) * 0.1

        X = X.T
        Y = Y.T

    X = X[100:, :]
    Y = Y[100:, :]

    return X, Y
