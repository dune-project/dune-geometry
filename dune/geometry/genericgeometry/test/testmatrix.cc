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

template <typename Field>
void runTest() {

  typedef GenericGeometry::MatrixHelper< GenericGeometry::DuneCoordTraits< Field > > MatrixHelper;

  // Test whether I can compute the square root of the determinant of A A^T of a nearly singular matrix.
  // This particular matrix makes the sqrtDetAAT method abort in dune-grid revision 6631.
  FieldMatrix< Field, 2, 2 > A;
  A[0][0] =  0.099999999999999867;
  A[0][1] = -0.010000000000002118;
  A[1][0] =  0.099999999999999867;
  A[1][1] = -0.0099999999999998979;

  Field sqrtDetAAT = MatrixHelper::template sqrtDetAAT< 2, 2 >( A );
  std::cout << "sqrtDetAAT = " << sqrtDetAAT << std::endl;

  FieldMatrix< Field, 2, 2 > invA;
  Field detA = MatrixHelper::template rightInvA< 2, 2 >( A, invA );
  std::cout << "detA       = " << detA << std::endl;

  // Check the identity invA*A*invA = A
  const FieldMatrix<Field, 2, 2> invAAinvA = invA.rightmultiplyany(A).rightmultiplyany(invA);
  printAsSingleLine(invA,      "invA       ");
  printAsSingleLine(invAAinvA, "invA*A*invA");
  std::cout << "Identity invA*A*invA = invA satisfied up to error: "
            << matrixNormDifference(invAAinvA, invA) << std::endl;

  // Check the identity A*invA*A = A
  const FieldMatrix<Field, 2, 2> AinvAA = A.rightmultiplyany(invA).rightmultiplyany(A);
  printAsSingleLine(A,      "A       ");
  printAsSingleLine(AinvAA, "A*invA*A");
  std::cout << "Identity    A*invA*A = A    satisfied up to error: "
            << matrixNormDifference(AinvAA, A) << std::endl;

  // Lets do the same for a non-square matrix.
  FieldMatrix< Field, 2, 3 > B;
  B[0][0] = A[0][0];
  B[0][1] = A[0][1];
  B[0][2] = 0;
  B[1][0] = A[1][0];
  B[1][1] = A[1][1];
  B[1][2] = 0;

  Field sqrtDetBBT = MatrixHelper::template sqrtDetAAT< 2, 3 >( B );
  std::cout << "sqrtDetBBT = " << sqrtDetBBT << std::endl;

  FieldMatrix< Field, 3, 2 > invB;
  Field detB = MatrixHelper::template rightInvA< 2, 3 >( B, invB );
  std::cout << "detB       = " << detB << std::endl;

  // Check the identity invB*B*invB = B
  const FieldMatrix<Field, 3, 2> invBBinvB = invB.rightmultiplyany(B).rightmultiplyany(invB);
  printAsSingleLine(invB,      "invB       ");
  printAsSingleLine(invBBinvB, "invB*B*invB");
  std::cout << "Identity invB*B*invB = invB satisfied up to error: "
            << matrixNormDifference(invBBinvB, invB) << std::endl;

  // Check the identity B*invB*B = B
  const FieldMatrix<Field, 2, 3> BinvBB = B.rightmultiplyany(invB).rightmultiplyany(B);
  printAsSingleLine(B,      "B       ");
  printAsSingleLine(BinvBB, "B*invB*B");
  std::cout << "Identity    B*invB*B = B    satisfied up to error: "
            << matrixNormDifference(BinvBB, B) << std::endl;
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
