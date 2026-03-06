import numpy as np

def lu_no_pivoting(A):
    n = A.shape[0]
    L = np.eye(n)
    U = A.copy().astype(float)
    try:
        for i in range(n):
            for j in range(i + 1, n):
                if U[i, i] == 0:
                    raise ValueError("Zero pivot encountered. No pivoting fails.")
                factor = U[j, i] / U[i, i]
                U[j, i:] -= factor * U[i, i:]
                L[j, i] = factor
    except ValueError as e:
        print("LU won't exist:", e)
        return None, None
    return L, U

def lu_partial_pivoting(A):
    n = A.shape[0]
    P = np.eye(n)
    L = np.zeros((n, n))
    U = A.copy().astype(float)
    for i in range(n):
        max_row = np.argmax(np.abs(U[i:, i])) + i
        U[[i, max_row], :] = U[[max_row, i], :]
        P[[i, max_row], :] = P[[max_row, i], :]
        for j in range(i + 1, n):
            factor = U[j, i] / U[i, i]
            U[j, i:] -= factor * U[i, i:]
            L[j, i] = factor
    np.fill_diagonal(L, 1)
    return P, L, U

def lu_complete_pivoting(A):
    n = A.shape[0]
    P = np.eye(n)
    Q = np.eye(n)
    L = np.zeros((n, n))
    U = A.copy().astype(float)
    for i in range(n):
        max_index = np.unravel_index(np.argmax(np.abs(U[i:, i:])), U[i:, i:].shape)
        max_row, max_col = max_index[0] + i, max_index[1] + i
        U[[i, max_row], :] = U[[max_row, i], :]
        P[[i, max_row], :] = P[[max_row, i], :]
        U[:, [i, max_col]] = U[:, [max_col, i]]
        Q[:, [i, max_col]] = Q[:, [max_col, i]]
        for j in range(i + 1, n):
            factor = U[j, i] / U[i, i]
            U[j, i:] -= factor * U[i, i:]
            L[j, i] = factor
    np.fill_diagonal(L, 1)
    return P, Q, L, U

def generate_spd_matrix(n):
    A = np.random.rand(n, n)
    return np.dot(A, A.T) + n * np.eye(n)

def cholesky_from_lu(A):
    L, U = lu_no_pivoting(A)
    if L is None or U is None:
        print("LU factorization failed. No Cholesky factor can be derived.")
        return
    D = np.diag(np.diag(U))
    S = np.sqrt(D)
    L0 = L @ S
    print("\nCholesky Factor L0 from LU:\n", L0)
    print("\nVerification (L0 @ L0.T):\n", L0 @ L0.T)
    print("\nOriginal Matrix A:\n", A)

def perform_lu(A):
    print("\nOriginal Matrix A:\n", A)
    L, U = lu_no_pivoting(A)
    if L is not None and U is not None:
        print("\nLU without Pivoting:\nL:\n", L, "\nU:\n", U)
    P, L, U = lu_partial_pivoting(A)
    print("\nLU with Partial Pivoting:\nP:\n", P, "\nL:\n", L, "\nU:\n", U)
    P, Q, L, U = lu_complete_pivoting(A)
    print("\nLU with Complete Pivoting:\nP:\n", P, "\nQ:\n", Q, "\nL:\n", L, "\nU:\n", U)

def main():
    A_diag_increasing = np.diag([1, 2, 3, 4, 5])
    print("\nDiagonal Matrix (Increasing Values):")
    perform_lu(A_diag_increasing)

    A_antidiag_decreasing = np.fliplr(np.diag([5, 4, 3, 2, 1]))
    print("\nAntidiagonal Matrix (Decreasing Values):")
    perform_lu(A_antidiag_decreasing)

    A_sum = A_diag_increasing + A_antidiag_decreasing
    print("\nSum of Diagonal and Antidiagonal Matrices:")
    perform_lu(A_sum)

    A_spd = generate_spd_matrix(5)
    print("\nSPD Matrix for Cholesky Example:")
    cholesky_from_lu(A_spd)

if __name__ == "__main__":
    main()
