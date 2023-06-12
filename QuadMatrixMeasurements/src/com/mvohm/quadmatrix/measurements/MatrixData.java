/*

 Copyright 2021-2023 M.Vokhmentsev

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

*/

package com.mvohm.quadmatrix.measurements;

import static com.mvohm.quadmatrix.measurements.AuxMethods.*;
import static com.mvohm.quadmatrix.measurements.MatrixDataGenerators.*;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.Random;

import com.mvohm.quadmatrix.BigDecimalMatrix;
import com.mvohm.quadmatrix.DoubleMatrix;
import com.mvohm.quadmatrix.Matrix;
import com.mvohm.quadmatrix.QuadrupleMatrix;
import com.mvohm.quadmatrix.measurements.AuxMethods.ErrorSet;

import com.mvohm.quadruple.Quadruple;

import Jama.CholeskyDecomposition;

/**
 * A helper class to generate random matrix data and examine accuracy of the solutions of the corresponding matrix
 */
public class MatrixData {

  protected static final double RANGE_BOTTOM          = -1.0;
  protected static final double RANGE_TOP             = 1.0;

  private int size;
  protected double[][] matrixData;
  protected double[] vector;
  protected double[] solution;

  private Quadruple[] quadrupleVector;
  private Quadruple[] quadrupleSolution;

  private BigDecimal[] bigDecimalVector;
  private BigDecimal[] bigDecimalSolution;

  protected double[][] matrixB;
  protected double[][] matrixX;

  private Quadruple[][] quadrupleMatrixB;
  private Quadruple[][] quadrupleMatrixX;

  private BigDecimal[][] bigDecimalMatrixB;
  private BigDecimal[][] bigDecimalMatrixX;

  private static double[][] unityMatrix = null;
  private static Quadruple[][] quadrupleUnityMatrix = null;
  private static BigDecimal[][] bigDecimalUnityMatrix = null;

  protected long time;

  protected enum Purpose { VECTOR_SOLUTION, MATRIX_SOLUTION, INVERSION };

  protected Purpose purpose;

  private String method_1_name = null;
  private String method_2_name = null;

  private final MathContext mc = new MathContext( BigDecimalMatrix.getDefaultPrecision(),
                                                  RoundingMode.HALF_EVEN);

  /* ***************************************************************************
   *********** Creating datasets ***********************************************
   *****************************************************************************/

  /**
   * fills matrixData with random values providing that it is SPD,
   * fills the solution with random values in the range defined by corresponding constants,
   * and finds the vector such that solving the equation system with this vector would give the found solution
   * @param size the size of the matrix to create
   * @param random contains a random seed and ensures reproducibility
   * @return and instance of RandomMatrixData with the generated data
   */
  public static MatrixData makeDataSetForSPDSolutions(int size, Random random) {
    // TODO 2023-06-05 18:18:54 Don't forget to remove them all -- added for debugging CollectStatistics
    // printMethodName();
    final MatrixData data = new MatrixData();
    data.purpose = Purpose.VECTOR_SOLUTION;
    MatrixDataGenerators.setRandomSeed(random);
    data.setMatrixData(randomSpdMatrix(size, RANGE_BOTTOM, RANGE_TOP, 1.0));
    data.solution = randomVector(size, RANGE_BOTTOM, RANGE_TOP);
    data.vector = multiply(data.matrixData, data.solution);
    return data;
  }

  /**
   * fills matrixData and solution with random values in the range defined by corresponding constants
   * and finds the vector such that solving the system with this vector would give the found solution
   * @param size the size of the matrix to create
   * @param random contains a random seed and ensures reproducibility
   * @return and instance of RandomMatrixData with the generated data
   */
  public static MatrixData makeDataSetForVectorSolutions(int size, Random random) {
    // printMethodName();
    final MatrixData data = new MatrixData();
    data.purpose = Purpose.VECTOR_SOLUTION;
    MatrixDataGenerators.setRandomSeed(random);
    data.setMatrixData(randomMatrix(size, RANGE_BOTTOM, RANGE_TOP));
    data.solution = randomVector(size, RANGE_BOTTOM, RANGE_TOP);
    data.vector = multiply(data.matrixData, data.solution);
    return data;
  }

  /**
   * Generates matrix data whose rows differ significantly in magnitude from each other.
   * fills matrixData A and a target solution x with random values in the range -1.0 .. 1.0,
   * then scales the rows of the matrix by a scale gradually growing from 1/sqrt(range) to sqrt(range),
   * and finds the vector b such that solving the system with the scaled matrix
   * and this vector would give the required solution.
   * @param size the size of the matrix to create
   * @param random contains a random seed and ensures reproducibility
   * @return and instance of RandomMatrixData with the generated data
   * @param size the size of the matrix to create
   * @param random contains a random seed and ensures reproducibility
   * @param scaleRange
   * @return
   */
  public static MatrixData makeLargeRangeDataSetForVectorSolutions(int size, Random random, double scaleRange) {
    // printMethodName();
    final MatrixData data = new MatrixData();
    data.purpose = Purpose.VECTOR_SOLUTION;
    MatrixDataGenerators.setRandomSeed(random);
    data.setMatrixData(randomMatrix(size, RANGE_BOTTOM, RANGE_TOP));

    // Multiply rows by an increasing scale
    final double scaleIncrease = Math.pow(scaleRange, 1.0/(size));  //
    double scale = 1.0 / Math.sqrt(scaleRange);
    for (int i = 0; i < size; i++) {
      data.matrixData[i] = multiply(data.matrixData[i], scale);
      scale *= scaleIncrease;
    }

    data.solution = randomVector(size, RANGE_BOTTOM, RANGE_TOP);
    data.vector = multiply(data.matrixData, data.solution);
    return data;
  }

  /**
   * fills matrixData and solution with non-uniform random values in the range defined by corresponding constants
   * and finds the vector such that solving the system with this vector would give the found solution
   * @param size the size of the matrix to create
   * @param random contains a random seed and ensures reproducibility
   * @return and instance of RandomMatrixData with the generated data
   */
  public static MatrixData makeNonUniformDataSet(int size, Random random, double density, double power, double slope) {
    // printMethodName();
    final MatrixData data = new MatrixData();
    data.purpose = Purpose.VECTOR_SOLUTION;
    MatrixDataGenerators.setRandomSeed(random);
    data.setMatrixData(randomSparsePowPlusLinearMatrix(size, density, power, slope));
    data.solution = randomVector(size, RANGE_BOTTOM, RANGE_TOP);
    data.vector = multiply(data.matrixData, data.solution);
    return data;
  }

  public static MatrixData makeDataSetForMatrixSolutions(int size, Random random) {
    // printMethodName();
    final MatrixData data = new MatrixData();
    data.purpose = Purpose.MATRIX_SOLUTION;
    MatrixDataGenerators.setRandomSeed(random);
    data.setMatrixData(randomMatrix(size, RANGE_BOTTOM, RANGE_TOP));
    data.matrixX = randomMatrix(size, RANGE_BOTTOM, RANGE_TOP);
    data.matrixB = multiply(data.matrixData, data.matrixX);
    return data;
  }

  public static MatrixData makeDataSetForInversions(int size, Random random) {
    // printMethodName();
    final MatrixData data = new MatrixData();
    data.purpose = Purpose.INVERSION;
    MatrixDataGenerators.setRandomSeed(random);
    data.setMatrixData(randomMatrix(size, RANGE_BOTTOM, RANGE_TOP));
    if (unityMatrix == null || unityMatrix.length != size) {
      unityMatrix = unityMatrix(size);
    }
    return data;
  }

  /* ************************************************************************
   ****** Perform operations and return errors ******************************
   **************************************************************************/

  // LU-decomposition-based simple vector solutions

  public ErrorSet doubleLuSolutionWithScalingErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    setTestedMethodName("(scaled) DoubleMatrix.solve(double[])");
    final DoubleMatrix matrix = new DoubleMatrix(matrixData, true);
    time = -System.nanoTime();
    matrix.solve(vector);
    time += System.nanoTime();
    final double[] actualSolution = matrix.getDoubleSolution();
    return findErrors(this.solution, actualSolution).setTime(time);
  }

  public ErrorSet jamaLuSolutionErrors() {
    // printMethodName();
    setTestedMethodName("Jama Matrix.solve(double[])");
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    final Jama.Matrix matrixA = new Jama.Matrix(matrixData);
    final Jama.Matrix matrixB = new Jama.Matrix(vector, size);
    time = -System.nanoTime();
    final Jama.Matrix matrixX = matrixA.solve(matrixB);
    time += System.nanoTime();
    final double[] actualSolution = matrixX.getColumnPackedCopy();
    return findErrors(this.solution, actualSolution).setTime(time);
  }

  public ErrorSet quadrupleLuSolutionWithScalingErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    setTestedMethodName("(scaled) QuadrupleMatrix.solve(Quadruple[])");
    final QuadrupleMatrix matrix = new QuadrupleMatrix(matrixData, true);
    quadrupleSolution = makeQuadrupleSolution();
    quadrupleVector = multiply(matrix.getQuadrupleData(), quadrupleSolution);
    time = -System.nanoTime();
    matrix.solve(quadrupleVector);
    time += System.nanoTime();
    final Quadruple[] actualSolution = matrix.getQuadrupleSolution();
    return findErrors(quadrupleSolution, actualSolution, false).setTime(time);
  }

  public ErrorSet bigDecimalLuSolutionWithScalingErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    setTestedMethodName("(scaled) BigDecimalMatrix.solve(BigDecimal[])");
    final BigDecimalMatrix matrix = new BigDecimalMatrix(matrixData, true);
    bigDecimalSolution = makeBigDecimalSolution();
    bigDecimalVector = multiply(matrix.getBigDecimalData(), bigDecimalSolution );
    time = -System.nanoTime();
    matrix.solve(bigDecimalVector);
    time += System.nanoTime();
    final BigDecimal[] actualSolution = matrix.getBigDecimalSolution();
    return findErrors(bigDecimalSolution, actualSolution, false).setTime(time);
  }


  public ErrorSet doubleLuSolutionWithoutScalingErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    setTestedMethodName("(unscaled) DoubleMatrix.solve(double[])");
    final DoubleMatrix matrix = new DoubleMatrix(matrixData, false);
    time = -System.nanoTime();
    matrix.solve(vector);
    time += System.nanoTime();
    final double[] actualSolution = matrix.getDoubleSolution();
    return findErrors(this.solution, actualSolution).setTime(time);
  }

  public ErrorSet quadrupleLuSolutionWithoutScalingErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    setTestedMethodName("(unscaled) QuadrupleMatrix.solve(Quadruple[])");
    final QuadrupleMatrix matrix = new QuadrupleMatrix(matrixData, false);
    quadrupleSolution = makeQuadrupleSolution();
    quadrupleVector = multiply(matrix.getQuadrupleData(), quadrupleSolution);
    time = -System.nanoTime();
    matrix.solve(quadrupleVector);
    time += System.nanoTime();
    final Quadruple[] actualSolution = matrix.getQuadrupleSolution();
    return findErrors(quadrupleSolution, actualSolution).setTime(time);
  }

  public ErrorSet bigDecimalLuSolutionWithoutScalingErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    setTestedMethodName("(unscaled) BigDecimalMatrix.solve(BigDecimal[])");
    final BigDecimalMatrix matrix = new BigDecimalMatrix(matrixData, false);
    bigDecimalSolution = makeBigDecimalSolution();
    bigDecimalVector = multiply(matrix.getBigDecimalData(), bigDecimalSolution );
    time = -System.nanoTime();
    matrix.solve(bigDecimalVector);
    time += System.nanoTime();
    final BigDecimal[] actualSolution = matrix.getBigDecimalSolution();
    return findErrors(bigDecimalSolution, actualSolution, false).setTime(time);
  }

  // Cholesky-decomposition-based accurate vector solutions

  public ErrorSet doubleAccurateLUSolutionWithScalingErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    setTestedMethodName("DoubleMatrix.solveAccurately(double[])");
    final DoubleMatrix matrix = new DoubleMatrix(matrixData);
    time = -System.nanoTime();
    matrix.solveAccurately(vector);
    time += System.nanoTime();
    final double[] actualSolution = matrix.getDoubleSolution();
    return findErrors(this.solution, actualSolution).setTime(time);
  }

  public ErrorSet quadrupleAccurateLUSolutionWithScalingErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    setTestedMethodName("QuadrupleMatrix.solveAccurately(Quadruple[])");
    final QuadrupleMatrix matrix = new QuadrupleMatrix(matrixData, true);
    if (quadrupleSolution == null) {
      quadrupleSolution = makeQuadrupleSolution();
      quadrupleVector = multiply(matrix.getQuadrupleData(), quadrupleSolution);
    }
    time = -System.nanoTime();
    matrix.solveAccurately(quadrupleVector);
    time += System.nanoTime();
    final Quadruple[] actualSolution = matrix.getQuadrupleSolution();
    return findErrors(quadrupleSolution, actualSolution).setTime(time);
  }

  public ErrorSet bigDecimalAccurateLUSolutionWithScalingErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    setTestedMethodName("BigDecimalMatrix.solveAccurately(BigDecimal[])");
    final BigDecimalMatrix matrix = new BigDecimalMatrix(matrixData, true);
    if (bigDecimalSolution == null) {
      bigDecimalSolution = makeBigDecimalSolution();
      bigDecimalVector = multiply(matrix.getBigDecimalData(), bigDecimalSolution );
    }
    time = -System.nanoTime();
    matrix.solveAccurately(bigDecimalVector);
    time += System.nanoTime();
    final BigDecimal[] actualSolution = matrix.getBigDecimalSolution();
    return findErrors(bigDecimalSolution, actualSolution, false).setTime(time);
  }

  // Cholesky-decomposition-based simple vector solutions

  public ErrorSet doubleSpdSolutionErrors() {
    // printMethodName();
    setTestedMethodName("DoubleMatrix.solveSPD(double[])");
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    final DoubleMatrix matrix = new DoubleMatrix(matrixData);
    time = -System.nanoTime();
    matrix.solveSPD(vector);
    time += System.nanoTime();
    final double[] actualSolution = matrix.getDoubleSolution();
    return findErrors(this.solution, actualSolution).setTime(time);
  }

  public ErrorSet jamaSpdSolutionErrors() {
    // printMethodName();
    setTestedMethodName("Jama Matrix.solveSPD(double[])");
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    final Jama.Matrix matrixA = new Jama.Matrix(matrixData);
    final Jama.Matrix matrixB = new Jama.Matrix(vector, size);
    time = -System.nanoTime();
    final CholeskyDecomposition d = new CholeskyDecomposition(matrixA);
    final Jama.Matrix matrixX = d.solve(matrixB);
    time += System.nanoTime();
    final double[] actualSolution = matrixX.getColumnPackedCopy();
    return findErrors(this.solution, actualSolution).setTime(time);
  }

  public ErrorSet quadrupleSPDSolutionErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    setTestedMethodName("QuadrupleMatrix.solveSPD(Quadruple[])");
    final QuadrupleMatrix matrix = new QuadrupleMatrix(matrixData, true);
    quadrupleSolution = makeQuadrupleSolution();
    quadrupleVector = multiply(matrix.getQuadrupleData(), quadrupleSolution);
    time = -System.nanoTime();
    matrix.solveSPD(quadrupleVector);
    time += System.nanoTime();
    final Quadruple[] actualSolution = matrix.getQuadrupleSolution();
    return findErrors(quadrupleSolution, actualSolution, false).setTime(time);
  }

  public ErrorSet bigDecimalSPDSolutionErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    setTestedMethodName("BigDecimalMatrix.solveSPD(BigDecimal[])");
    final BigDecimalMatrix matrix = new BigDecimalMatrix(matrixData, true);
    bigDecimalSolution = makeBigDecimalSolution();
    bigDecimalVector = multiply(matrix.getBigDecimalData(), bigDecimalSolution);
    time = -System.nanoTime();
    matrix.solveSPD(bigDecimalVector);
    time += System.nanoTime();
    final BigDecimal[] actualSolution = matrix.getBigDecimalSolution();
    return findErrors(bigDecimalSolution, actualSolution).setTime(time);
  }

  // Cholesky-decomposition-based precise vector solutions

  public ErrorSet doubleAccurateSPDSolutionErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    setTestedMethodName("DoubleMatrix.solveSPDAccurately(double[])");
    final DoubleMatrix matrix = new DoubleMatrix(matrixData);
    time = -System.nanoTime();
    matrix.solveSPDAccurately(vector);
    time += System.nanoTime();
    final double[] actualSolution = matrix.getDoubleSolution();
    return findErrors(this.solution, actualSolution).setTime(time);
  }

  public ErrorSet quadrupleAccurateSPDSolutionErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    setTestedMethodName("QuadrupleMatrix.solveSPDAccurately(Quadruple[])");
    final QuadrupleMatrix matrix = new QuadrupleMatrix(matrixData, true);
    if (quadrupleSolution == null) {
      quadrupleSolution = makeQuadrupleSolution();
      quadrupleVector = multiply(matrix.getQuadrupleData(), quadrupleSolution);
    }
    time = -System.nanoTime();
    matrix.solveSPDAccurately(quadrupleVector);
    time += System.nanoTime();
    final Quadruple[] actualSolution = matrix.getQuadrupleSolution();
    return findErrors(quadrupleSolution, actualSolution).setTime(time);
  }

  public ErrorSet bigDecimalAccurateSPDSolutionErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.VECTOR_SOLUTION);
    setTestedMethodName("BigDecimalMatrix.solveSPDAccurately(BigDecimal[])");
    final BigDecimalMatrix matrix = new BigDecimalMatrix(matrixData, true);
    if (bigDecimalSolution == null) {
      bigDecimalSolution = makeBigDecimalSolution();
      bigDecimalVector = multiply(matrix.getBigDecimalData(), bigDecimalSolution);
    }
    time = -System.nanoTime();
    matrix.solveSPDAccurately(bigDecimalVector);
    time += System.nanoTime();
    final BigDecimal[] actualSolution = matrix.getBigDecimalSolution();
    return findErrors(bigDecimalSolution, actualSolution).setTime(time);
  }

  // LU-decomposition-based simple matrix solutions

  public ErrorSet doubleMatrixSolutionErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.MATRIX_SOLUTION);
    setTestedMethodName("DoubleMatrix.solve(double[][])");
    final DoubleMatrix matrix = new DoubleMatrix(matrixData, true);
    time = -System.nanoTime();
    final double[][] actualSolution = matrix.solve(matrixB).getDoubleData();
    time += System.nanoTime();
    return findErrors(this.matrixX, actualSolution).setTime(time);
  }

  public ErrorSet jamaMatrixSolutionErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.MATRIX_SOLUTION);
    setTestedMethodName("Jama Matrix.solve(double[][])");
    final Jama.Matrix matrixA = new Jama.Matrix(matrixData);
    final Jama.Matrix matrixB = new Jama.Matrix(this.matrixB);
    time = -System.nanoTime();
    final Jama.Matrix matrixX = matrixA.solve(matrixB);
    time += System.nanoTime();
    final double[][] actualSolution = matrixX.getArray();
    return findErrors(this.matrixX, actualSolution).setTime(time);
  }

  public ErrorSet quadrupleMatrixSolutionErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.MATRIX_SOLUTION);
    setTestedMethodName("QuadrupleMatrix.solve(Quadruple[][])");
    final QuadrupleMatrix matrix = new QuadrupleMatrix(matrixData, true);

    quadrupleMatrixX = makeQuadrupleMatrixX();
    quadrupleMatrixB = multiply(matrix.getQuadrupleData(), quadrupleMatrixX);

    time = -System.nanoTime();
    final Quadruple[][] actualSolution = matrix.solve(quadrupleMatrixB).getQuadrupleData();
    time += System.nanoTime();
    return findErrors(quadrupleMatrixX, actualSolution).setTime(time);
  }

  public ErrorSet bigDecimalMatrixSolutionErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.MATRIX_SOLUTION);
    setTestedMethodName("BigDecimalMatrix.solve(BigDecimal[][])");
    final BigDecimalMatrix matrix = new BigDecimalMatrix(matrixData, true);

    bigDecimalMatrixX = makeBigDecimalMatrixX();
    bigDecimalMatrixB = multiply(matrix.getBigDecimalData(), bigDecimalMatrixX);

    time = -System.nanoTime();
    final BigDecimal[][] actualSolution = matrix.solve(bigDecimalMatrixB).getBigDecimalData();
    time += System.nanoTime();
    return findErrors(bigDecimalMatrixX, actualSolution).setTime(time);
  }

  // LU-decomposition-based accurate matrix solutions

  public ErrorSet doubleAccurateMatrixSolutionErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.MATRIX_SOLUTION);
    setTestedMethodName("DoubleMatrix.solveAccurately(double[][])");
    final DoubleMatrix matrix = new DoubleMatrix(matrixData, true);
    time = -System.nanoTime();
    final double[][] actualSolution = matrix.solveAccurately(matrixB).getDoubleData();
    time += System.nanoTime();
    return findErrors(this.matrixX, actualSolution).setTime(time);
  }

  public ErrorSet quadrupleAccurateMatrixSolutionErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.MATRIX_SOLUTION);
    setTestedMethodName("QuadrupleMatrix.solveAccurately(Quadruple[][])");
    final QuadrupleMatrix matrix = new QuadrupleMatrix(matrixData, true);
    if (quadrupleMatrixX == null) {
      quadrupleMatrixX = makeQuadrupleMatrixX();
      quadrupleMatrixB = multiply(matrix.getQuadrupleData(), quadrupleMatrixX);
    }
    time = -System.nanoTime();
    final Quadruple[][] actualSolution = matrix.solveAccurately(quadrupleMatrixB).getQuadrupleData();
    time += System.nanoTime();
    return findErrors(quadrupleMatrixX, actualSolution).setTime(time);
  }

  public ErrorSet bigDecimalAccurateMatrixSolutionErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.MATRIX_SOLUTION);
    setTestedMethodName("BigDecimalMatrix.solveAccurately(BigDecimal[][])");
    final BigDecimalMatrix matrix = new BigDecimalMatrix(matrixData, true);
    if (bigDecimalMatrixX == null) {
      bigDecimalMatrixX = makeBigDecimalMatrixX();
      bigDecimalMatrixB = multiply(matrix.getBigDecimalData(), bigDecimalMatrixX);
    }
    time = -System.nanoTime();
    final BigDecimal[][] actualSolution = matrix.solveAccurately(bigDecimalMatrixB).getBigDecimalData();
    time += System.nanoTime();
    return findErrors(bigDecimalMatrixX, actualSolution).setTime(time);
  }

  // Matrix inversion

  public ErrorSet doubleMatrixInversionErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.INVERSION);
    setTestedMethodName("DoubleMatrix.inverse()");
    final DoubleMatrix matrix = new DoubleMatrix(matrixData, true);
    time = -System.nanoTime();
    final Matrix inverse = matrix.inverse();
    time += System.nanoTime();
    final double[][] product = multiply(matrixData, inverse.getDoubleData());
    return findErrors(unityMatrix, product, 1).setTime(time);
  }

  public ErrorSet jamaMatrixInversionErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.INVERSION);
    setTestedMethodName("Jama Matrix.inverse()");
    final Jama.Matrix matrixA = new Jama.Matrix(matrixData);
    time = -System.nanoTime();
    final Jama.Matrix inverse = matrixA.inverse();
    time += System.nanoTime();
    final double[][] product = multiply(matrixData, inverse.getArray());
    return findErrors(unityMatrix, product, 1).setTime(time);
  }

  public ErrorSet quadrupleMatrixInversionErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.INVERSION);
    setTestedMethodName("QuadrupleMatrix.inverse()");
    quadrupleUnityMatrix = convertToQuadruples(unityMatrix);
    final QuadrupleMatrix matrix = new QuadrupleMatrix(matrixData, true);
    time = -System.nanoTime();
    final Matrix inverse = matrix.inverse();
    time += System.nanoTime();
    final Quadruple[][] product = multiply(matrix.getQuadrupleData(), inverse.getQuadrupleData());
    return findErrors(quadrupleUnityMatrix, product).setTime(time);
  }

  public ErrorSet bigDecimalMatrixInversionErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.INVERSION);
    setTestedMethodName("BigDecimalMatrix.inverse()");
    bigDecimalUnityMatrix = convertToBigDecimals(unityMatrix);
    final BigDecimalMatrix matrix = new BigDecimalMatrix(matrixData, true);
    time = -System.nanoTime();
    final Matrix inverse = matrix.inverse();
    time += System.nanoTime();
    final BigDecimal[][] product = multiply(matrix.getBigDecimalData(), inverse.getBigDecimalData());
    return findErrors(bigDecimalUnityMatrix, product).setTime(time);
  }

  // Accurate matrix inversion

  public ErrorSet doubleAccurateMatrixInversionErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.INVERSION);
    setTestedMethodName("DoubleMatrix.inverseAccurately()");
    final DoubleMatrix matrix = new DoubleMatrix(matrixData, true);
    time = -System.nanoTime();
    final Matrix inverse = matrix.inverseAccurately();
    time += System.nanoTime();
    final double[][] product = multiply(matrixData, inverse.getDoubleData());
    return findErrors(unityMatrix, product, 1).setTime(time);
  }

  public ErrorSet quadrupleAccurateMatrixInversionErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.INVERSION);
    if (quadrupleUnityMatrix == null) {
      quadrupleUnityMatrix = convertToQuadruples(unityMatrix);
    }
    setTestedMethodName("QuadrupleMatrix.inverseAccurately()");
    final QuadrupleMatrix matrix = new QuadrupleMatrix(matrixData, true);
    time = -System.nanoTime();
    final Matrix inverse = matrix.inverseAccurately();
    time += System.nanoTime();
    final Quadruple[][] product = multiply(matrix.getQuadrupleData(), inverse.getQuadrupleData());
    return findErrors(quadrupleUnityMatrix, product).setTime(time);
  }

  public ErrorSet bigDecimalAccurateMatrixInversionErrors() {
    // printMethodName();
    checkPurpose(performerMethodName(), Purpose.INVERSION);
    setTestedMethodName("BigDecimalMatrix.inverseAccurately()");
    if (bigDecimalUnityMatrix == null)
      bigDecimalUnityMatrix = convertToBigDecimals(unityMatrix);
    final BigDecimalMatrix matrix = new BigDecimalMatrix(matrixData, true);
    time = -System.nanoTime();
    final Matrix inverse = matrix.inverseAccurately();
    time += System.nanoTime();
    final BigDecimal[][] product = multiply(matrix.getBigDecimalData(), inverse.getBigDecimalData());
    return findErrors(bigDecimalUnityMatrix, product).setTime(time);
  }

  // matrix multiplications

  public ErrorSet doubleMultiplicationErrors() {
    checkPurpose(performerMethodName(), Purpose.MATRIX_SOLUTION);
    setTestedMethodName("DoubleMatrix.inverse()");
    final Matrix matrix = new DoubleMatrix(matrixData, true);
    time = -System.nanoTime();
    final Matrix product = matrix.multiply(matrixX);
    time += System.nanoTime();
    return findErrors(matrixB, product.getDoubleData()).setTime(time);
  }

  public ErrorSet jamaMultiplicationErrors() {
    checkPurpose(performerMethodName(), Purpose.MATRIX_SOLUTION);
    setTestedMethodName("Jama Matrix.inverse()");
    final Jama.Matrix matrixA = new Jama.Matrix(matrixData);
    final Jama.Matrix matrixX = new Jama.Matrix(this.matrixX);
    time = -System.nanoTime();
    final Jama.Matrix matrixB = matrixA.times(matrixX);
    time += System.nanoTime();
    return findErrors(this.matrixB, matrixB.getArray()).setTime(time);
  }

  public ErrorSet quadrupleMultiplicationErrors() {
    checkPurpose(performerMethodName(), Purpose.MATRIX_SOLUTION);
    setTestedMethodName("DoubleMatrix.inverse()");
    final QuadrupleMatrix matrix = new QuadrupleMatrix(matrixData, true);

    quadrupleMatrixX = makeQuadrupleMatrixX();
    final Quadruple[][] expectedProduct = multiply(matrix.getQuadrupleData(), quadrupleMatrixX);

    time = -System.nanoTime();
    final Matrix actualProduct = matrix.multiply(quadrupleMatrixX);
    time += System.nanoTime();

    return findErrors(expectedProduct, actualProduct.getQuadrupleData()).setTime(time);
  }

  public ErrorSet bigDecimalMultiplicationErrors() {
    checkPurpose(performerMethodName(), Purpose.MATRIX_SOLUTION);
    setTestedMethodName("DoubleMatrix.inverse()");
    final BigDecimalMatrix matrix = new BigDecimalMatrix(convertToQuadruples(matrixData), true);

    // If it were converted from doubles directly, the precision would be too low
    bigDecimalMatrixX = convertToBigDecimals(convertToQuadruples(matrixX));
    final BigDecimal[][] expectedProduct = multiply(matrix.getBigDecimalData(), bigDecimalMatrixX);

    time = -System.nanoTime();
    final Matrix actualProduct = matrix.multiply(bigDecimalMatrixX);
    time += System.nanoTime();

    return findErrors(expectedProduct, actualProduct.getBigDecimalData()).setTime(time);
  }



  /* *************************************************************************
  ******** Private methods ***************************************************
  ***************************************************************************/

  protected void setMatrixData(double[][] data) {
    matrixData = data;
    size = data.length;
  }

  int getSize() {
    return size;
 }

  protected void setTestedMethodName(String methodName) {
    if (method_1_name == null) {
      method_1_name = methodName;
    } else if (method_2_name == null) {
      method_2_name = methodName;
    }
  }

  String testedMethod1name() {
    return method_1_name;
  }

  String testedMethod2name() {
    return method_2_name;
  }

  protected static String performerMethodName() {
    final String s = new Exception().getStackTrace()[1].toString();
    return s.replaceFirst("\\(.*\\)", "()").replaceFirst(".*\\.", "");
  }

  protected void checkPurpose(String methodName, Purpose impliedPurpose) {
    if (this.purpose != impliedPurpose) {
      throw new IllegalArgumentException(String.format(
                  "Calling %s on MatrixData implies %s, while this instance was creted for %s",
                  methodName, impliedPurpose, this.purpose));
    }
  }

  public long  getTime() {
    return time;
  }

  private static double[] multiply(double[] vector, double factor) {
    final double[] result = new double[vector.length];
    for (int i = 0; i < vector.length; i++) {
      result[i] = vector[i] * factor;
    }
    return result;
  }

  private Quadruple[] makeQuadrupleSolution() {
    return convertToQuadruples(solution);
  }

  private BigDecimal[] makeBigDecimalSolution() {
    return convertToBigDecimals(solution);
  }

  private Quadruple[][] makeQuadrupleMatrixX() {
    return convertToQuadruples(matrixX);
  }

  private BigDecimal[][] makeBigDecimalMatrixX() {
    return convertToBigDecimals(matrixX);
  }

  private static double[][] multiply(double[][] matrixA, double[][] matrixB) {
    final int size = matrixA.length;
    final double[][] result = new double[size][size];
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        result[i][j] = multiplyRowByColumn(matrixA, i, matrixB, j);
      }
    }
    return result;
  }

  private static double multiplyRowByColumn(double[][] matrixA, int rowIndex, double[][] matrixB, int colIndex) {
    final int size = matrixA.length;
    final Quadruple productElement = new Quadruple();
    for (int k = 0; k < size; k++) {
      productElement.add(new Quadruple(matrixA[rowIndex][k]).multiply(matrixB[k][colIndex]));
    }
    return productElement.doubleValue();
  }

  protected static double[] multiply(double[][] matrix, double[] vector) {
    final int length = matrix.length;
    final double[] result = new double[length];
    for (int i = 0; i < length; i++) {
      final Quadruple productElement = new Quadruple();
      for (int j = 0; j < length; j++) {
        productElement.add(new Quadruple(matrix[i][j]).multiply(vector[j]));
      }
      result[i] = productElement.doubleValue();
    }
    return result;
  }

  private static Quadruple[] multiply(Quadruple[][] matrix, Quadruple[] vector) {
    final int length = matrix.length;
    final Quadruple[] result = new Quadruple[length];
    for (int i = 0; i < length; i++) {
      final Quadruple productElement = new Quadruple();
      for (int j = 0; j < length; j++) {
        productElement.add(new Quadruple(matrix[i][j]).multiply(vector[j]));
      }
      result[i] = productElement;
    }
    return result;
  }

  private BigDecimal[] multiply(BigDecimal[][] matrix, BigDecimal[] vector) {
    final int length = matrix.length;
    final BigDecimal[] result = new BigDecimal[length];
    for (int i = 0; i < length; i++) {
      BigDecimal productElement = BigDecimal.ZERO;
      for (int j = 0; j < length; j++) {
        productElement = productElement.add(matrix[i][j].multiply(vector[j], mc), mc);
      }
      result[i] = productElement;
    }
    return result;
  }

  private static Quadruple[][] multiply(Quadruple[][] matrixA, Quadruple[][] matrixB) {
    final int length = matrixA.length;
    final Quadruple[][] result = new Quadruple[length][length];
    for (int i = 0; i < length; i++) {
      for (int j = 0; j < length; j++) {
        result[i][j] = multiplyRowByColumn(matrixA, i, matrixB, j);
      }
    }
    return result;
  }

  private static BigDecimal[][] multiply(BigDecimal[][] matrixA, BigDecimal[][] matrixB) {
    final int length = matrixA.length;
    final BigDecimal[][] result = new BigDecimal[length][length];
    for (int i = 0; i < length; i++) {
      for (int j = 0; j < length; j++) {
        result[i][j] = multiplyRowByColumn(matrixA, i, matrixB, j);
      }
    }
    return result;
  }

  private static Quadruple multiplyRowByColumn(Quadruple[][] matrixA, int rowIndex, Quadruple[][] matrixB, int colIndex) {
    final int size = matrixA.length;
    final Quadruple productElement = new Quadruple();
    for (int k = 0; k < size; k++) {
      productElement.add(Quadruple.multiply(matrixA[rowIndex][k], matrixB[k][colIndex]));
    }
    return productElement;
  }

  private static BigDecimal multiplyRowByColumn(BigDecimal[][] matrixA, int rowIndex, BigDecimal[][] matrixB, int colIndex) {
    final int size = matrixA.length;
    BigDecimal productElement = BigDecimal.ZERO;
    for (int k = 0; k < size; k++) {
      productElement = productElement.add(matrixA[rowIndex][k].multiply(matrixB[k][colIndex]));
    }
    return productElement;
  }

} // class MatrixData {