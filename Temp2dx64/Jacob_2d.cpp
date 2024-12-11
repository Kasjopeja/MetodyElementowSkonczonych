#include "Jacob_2d.h"

void Jacob_2d(std::vector<std::vector<double>>& J_,
              std::vector<std::vector<double>>& J_inv,
              int P, int N_p, int NBN,
              const std::vector<std::vector<double>>& N1,
              const std::vector<std::vector<double>>& N2,
              const std::vector<double>& X,
              const std::vector<double>& Y,
              double& DetJ) {
    // Initialize Jacobian and its inverse to zero
    J_ = {{0.0, 0.0}, {0.0, 0.0}};
    J_inv = {{0.0, 0.0}, {0.0, 0.0}};

    // Compute Jacobian matrix
    for (int i = 0; i < NBN; ++i) {
        J_[0][0] += N1[i][P] * X[i];
        J_[0][1] += N1[i][P] * Y[i];
        J_[1][0] += N2[i][P] * X[i];
        J_[1][1] += N2[i][P] * Y[i];
    }

    // Compute determinant of the Jacobian
    DetJ = J_[0][0] * J_[1][1] - J_[0][1] * J_[1][0];

    // Compute inverse of the Jacobian
    Inv_MAT(2, J_, J_inv);
}

void Inv_MAT(int N, const std::vector<std::vector<double>>& Mat, std::vector<std::vector<double>>& Inv) {
    // Augmented matrix for Gaussian elimination
    std::vector<std::vector<double>> A(N, std::vector<double>(2 * N, 0.0));

    // Initialize the augmented matrix
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A[i][j] = Mat[i][j];
        }
        A[i][N + i] = 1.0;
    }

    // Gaussian elimination to compute the inverse
    for (int PivRow = 0; PivRow < N; ++PivRow) {
        double PivElt = A[PivRow][PivRow];

        if (std::abs(PivElt) < 1e-12) {
            // Find a non-zero pivot and swap rows
            int K = PivRow + 1;
            while (K < N && std::abs(A[K][PivRow]) < 1e-12) {
                ++K;
            }

            if (K == N) {
                std::cerr << "Error: Singular matrix, cannot compute inverse." << std::endl;
                return;
            }

            std::swap(A[PivRow], A[K]);
            PivElt = A[PivRow][PivRow];
        }

        // Normalize the pivot row
        for (int j = 0; j < 2 * N; ++j) {
            A[PivRow][j] /= PivElt;
        }

        // Eliminate other rows
        for (int TarRow = 0; TarRow < N; ++TarRow) {
            if (TarRow != PivRow) {
                double TarElt = A[TarRow][PivRow];
                for (int j = 0; j < 2 * N; ++j) {
                    A[TarRow][j] -= A[PivRow][j] * TarElt;
                }
            }
        }
    }

    // Extract the inverse matrix
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            Inv[i][j] = A[i][N + j];
        }
    }
}