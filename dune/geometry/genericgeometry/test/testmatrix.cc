// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <iomanip>
#include <limits>

#include <dune/common/gmpfield.hh>
#include <dune/common/parallel/mpihelper.hh>

#include "../matrixhelper.hh"
#include "../geometrytraits.hh"

using namespace Dune;

template <class Matrix>
void printAsSingleLine(Matrix const &m, std::string name)
{
  std::cout << name << " = [ ";
  for (size_t i = 0; i<m.N() - 1; ++i)
    std::cout << m[i] << ", ";
  std::cout << m[m.N() - 1] << " ]" << std::endl;
}

template <class Matrix>
typename FieldTraits<typename DenseMatVecTraits<Matrix>::value_type>::real_type
matrixNormDifference(Matrix const &m1, Matrix const &m2)
{
  Matrix tmp = m1;
  tmp -= m2;
  return tmp.infinity_norm();
}

template <typename Field, int rows, int cols>
void examineMatrix(std::string name, Dune::FieldMatrix<Field, rows, cols> const &A) {
  using MatrixHelper =
      GenericGeometry::MatrixHelper<GenericGeometry::DuneCoordTraits<Field>>;
  std::cout << "Examining matrix " << name << std::endl;

  Field sqrtDetAAT = MatrixHelper::template sqrtDetAAT< rows, cols >( A );
  std::cout << "sqrtDetAAT  = " << sqrtDetAAT << std::endl;

  FieldMatrix< Field, cols, rows > invA;
  Field detA = MatrixHelper::template rightInvA< rows, cols >( A, invA );
  std::cout << "detA        = " << detA << std::endl;

  // Check the identity A*invA*A = A
  const FieldMatrix<Field, rows, cols> AinvAA = A.rightmultiplyany(invA).rightmultiplyany(A);
  printAsSingleLine(A,      "A          ");
  printAsSingleLine(AinvAA, "A*invA*A   ");
  std::cout << "Identity    A*invA*A = A    satisfied up to error: "
            << matrixNormDifference(AinvAA, A) << std::endl;

  // Check the identity invA*A*invA = A
  const FieldMatrix<Field, cols, rows> invAAinvA = invA.rightmultiplyany(A).rightmultiplyany(invA);
  printAsSingleLine(invA,      "invA       ");
  printAsSingleLine(invAAinvA, "invA*A*invA");
  std::cout << "Identity invA*A*invA = invA satisfied up to error: "
            << matrixNormDifference(invAAinvA, invA) << std::endl;
  std::cout << std::endl;
}

template <typename Field>
void runTest() {
  // Test whether I can compute the square root of the determinant of A A^T of a nearly singular matrix.
  // This particular matrix makes the sqrtDetAAT method abort in dune-grid revision 6631.
  FieldMatrix< Field, 2, 2 > A;
  A[0][0] =  0.099999999999999867;
  A[0][1] = -0.010000000000002118;
  A[1][0] =  0.099999999999999867;
  A[1][1] = -0.0099999999999998979;
  examineMatrix("#1", A);

  // Lets do the same for a non-square matrix.
  FieldMatrix< Field, 2, 3 > B;
  B[0][0] = A[0][0];
  B[0][1] = A[0][1];
  B[0][2] = 0;
  B[1][0] = A[1][0];
  B[1][1] = A[1][1];
  B[1][2] = 0;
  examineMatrix("#2", B);
}

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);
  std::cout << "Running tests for type: double" << std::endl;
  runTest<double>();
  std::cout << std::endl;

  std::cout << "Running tests for type: long double" << std::endl;
  runTest<long double>();
  std::cout << std::endl;

#if HAVE_GMP
  std::cout << "Running tests for type: GMPField<72>" << std::endl;
  runTest<GMPField<72>>();
  std::cout << std::endl;

  std::cout << "Running tests for type: GMPField<160>" << std::endl;
  runTest<GMPField<160>>();
  std::cout << std::endl;
#endif
}
