/* REQUIRED */
#include "nelder_mead.h"
using namespace std;
/* REQUIRED */

/* CLASS: Matrix */
Matrix::Matrix() {}
Matrix::Matrix(float C0, float C1) {
    coordinates.push_back(C0);
    coordinates.push_back(C1);
}
Matrix::Matrix(float C0, float C1, float C2) {
    coordinates.push_back(C0);
    coordinates.push_back(C1);
    coordinates.push_back(C2);
}
Matrix::Matrix(float C0, float C1, float C2, float C3) {
    coordinates.push_back(C0);
    coordinates.push_back(C1);
    coordinates.push_back(C2);
    coordinates.push_back(C3);
}
Matrix::Matrix(float C0, float C1, float C2, float C3, float C4) {
    coordinates.push_back(C0);
    coordinates.push_back(C1);
    coordinates.push_back(C2);
    coordinates.push_back(C3);
    coordinates.push_back(C4);
}
float& Matrix::operator[](int i) {
    return coordinates[i];
}
float Matrix::at(int i) const {
    return coordinates[i];
}
int Matrix::dimension() const {
    return coordinates.size();
}
void Matrix::prepare(int size) {
    for (int i(0); i < size; ++i) {
        coordinates.push_back(0);
    }
}
Matrix Matrix::operator+(Matrix operand) {
    Matrix output;
    output.prepare(dimension());
    for (int i(0); i < dimension(); ++i) {
        output[i] = coordinates[i] + operand[i];
    }
    return output;
}
void Matrix::operator+=(Matrix operand) {
    for (int i(0); i < dimension(); ++i) {
        coordinates[i] += operand[i];
    }
}
Matrix Matrix::operator-(Matrix operand) {
    Matrix output;
    output.prepare(dimension());
    for (int i(0); i < dimension(); ++i) {
        output[i] = coordinates[i] - operand[i];
    }
    return output;
}
bool Matrix::operator==(Matrix operand) {
    if (dimension() != operand.dimension()) {
        return false;
    }
    for (int i(0); i < dimension(); ++i) {
        if (coordinates[i] != operand[i]) {
            return false;
        }
    }
    return true;
}
Matrix Matrix::operator*(float factor) {
    Matrix output;
    output.prepare(dimension());
    for (int i(0); i < dimension(); ++i) {
        output[i] = coordinates[i] * factor;
    }
    return output;
}
Matrix Matrix::operator/(float denominator) {
    Matrix output;
    output.prepare(dimension());
    for (int i(0); i < dimension(); ++i) {
        output[i] = coordinates[i] / denominator;
    }
    return output;
}
void Matrix::operator/=(float denominator) {
    for (int i(0); i < dimension(); ++i) {
        coordinates[i] /= denominator;
    }
}
bool Matrix::operator<(const Matrix operand) const {
    for (int i(0); i < dimension(); ++i) {
        if (at(i) < operand.at(i)) {
            return false;
        }
        else if (at(i) > operand.at(i)) {
            return true;
        }
    }
}
float Matrix::length() {
    float absolute(0);
    for (int i(0); i < dimension(); ++i) {
        absolute += coordinates[i] * coordinates[i];
    }
    return pow(absolute, 0.5f);
}
/* CLASS: Matrix */

/* CLASS: Memo */
Memo::Memo() {}
float Memo::lookup(Matrix A) {
    if (!contains(A)) {
        throw A;
    }
    else {
        return values[A];
    }
}
void Memo::insert(Matrix A, float value) {
    values[A] = value;
}
bool Memo::contains(Matrix A) {
    map<Matrix, float>::iterator it = values.find(A);
    return it != values.end();
}
/* CLASS: Memo */

/* CLASS: NelderMead */
NelderMead::NelderMead(int dimension, float termination_distance = 0.001) {
    this->dimension = dimension;
    srand(time(NULL));
    alpha = 1;
    gamma = 2;
    rho = -0.5;
    sigma = 0.5;
    this->termination_distance = termination_distance;
}
bool NelderMead::operator()(const Matrix& A, const Matrix& B) {
    return dp.lookup(A) < dp.lookup(B);
}
bool NelderMead::done() {
    if (matrices.size() < dimension) {
        return false;
    }
    for (int i(0); i < dimension + 1; ++i) {
        for (int j(0); j < dimension + 1; ++j) {
            if (i == j) continue;
            if ((matrices[i] - matrices[j]).length() > termination_distance) {
                return false;
            }
        }
    }
    return true;
}
void NelderMead::insert(Matrix A) {
    if (matrices.size() < dimension + 1) {
        matrices.push_back(A);
    }
}
Matrix NelderMead::step(Matrix A, float score) {
    dp.insert(A, score);
    try {
        if (matrices.size() < dimension + 1) {
            matrices.push_back(A);
        }
        if (matrices.size() == dimension + 1) {
            while(!done()) {
                sort(matrices.begin(), matrices.end(), *this);
                /* Center of gravity */
                Matrix cog;
                cog.prepare(dimension);
                for (int i(1); i <= dimension; ++i) {
                    cog += matrices[i];
                }
                cog /= dimension;
                Matrix best = matrices[dimension];
                Matrix worst = matrices[0];
                Matrix snd_worst = matrices[1];
                /* Reflect */
                Matrix reflected = cog + (cog - worst) * alpha;
                if (f(reflected) > f(snd_worst) && f(reflected) < f(best)) {
                    matrices[0] = reflected;
                }
                else if (f(reflected) > f(best)) {
                    /* Expand */
                    Matrix expanded = cog + (cog - worst) * gamma;
                    if (f(expanded) > f(reflected)) {
                        matrices[0] = expanded;
                    }
                    else {
                        matrices[0] = reflected;
                    }
                }
                else {
                    /* Contract */
                    Matrix contracted = cog + (cog - worst) * rho;
                    if (f(contracted) > f(worst)) {
                        matrices[0] = contracted;
                    }
                    else {
                        for (int i(0); i < dimension; ++i) {
                            matrices[i] = best + (matrices[i] - best) * sigma;
                        }
                    }
                }
            }
                /* Output cog */
                Matrix cog;
                for (int i(0); i <= dimension; ++i) {
                    cog += matrices[i];
                }
                return cog / (dimension + 1);
            }
            else {
                /* output random vectors */
                Matrix output;
                output.prepare(dimension);
                for (int i(0); i < dimension; ++i) {
                    output[i] = 0.001 * (rand() % 1000);
                }
                return output;
            }
            
        }
        catch (Matrix A) {
            return A;
    }
}
float NelderMead::f(Matrix A) {
    return dp.lookup(A);
}
/* CLASS: NelderMead */
