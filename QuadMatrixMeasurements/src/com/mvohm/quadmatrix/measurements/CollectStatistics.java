/*

 Copyright 2021-2025 M.Vokhmentsev

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

import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.Locale;
import java.util.Random;

import com.mvohm.quadmatrix.BigDecimalMatrix;

/**
 * Estimates accuracy and execution times of the most common operations 
 * performed on matrices of different types. 
 * The matrices being tested are BigDecimalMatrix, DoubleMatrix, and QuadrupleMatrix 
 * from project QuadMartix (https://github.com/m-vokhm/QuadMatrix),
 * and JAMA from https://math.nist.gov/javanumerics/jama/ (used as reference library). 
 * 
 * Complete execution may take up to 48 hours or even more,
 * depending on the machine performance.   
 */

public class CollectStatistics {

  private static final int RAND_SEED = 123;
  static Random random = new Random(RAND_SEED);

  /** An object of a class implementing this interface is encapsulated in an OperationTester object 
   * and is responsible for generating a data sample of a certain size to test an operation of a certain type. */
  public interface DataGenerator {
    MatrixData generate();
  }

  /** An object of a class implementing this interface is encapsulated in an OperationTester object 
   * and is responsible for executing the method under test and collecting statistics 
   * on the performance and accuracy of the method under test */
  public interface OperationPerformer {
    ErrorSet perform(MatrixData matrixData);
  }

  /** Types of operations performed on matrices. */
  enum Operations {
    SIMPLE_VECTOR_SOLUTION,       // Solves A * x = b, where x and b are vectors
    ACCURATE_VECTOR_SOLUTION,     // with iterative refinement
    SIMPLE_SPD_SOLUTION,          // solve symmetric  positive-definite matrix using Cholesky decomposition 
    ACCURATE_SPD_SOLUTION,        // with iterative refinement
    SIMPLE_MATRIX_SOLUTION,       // Solves A * X = B, where X and B are matrices
    ACCURATE_MATRIX_SOLUTION,     // with iterative refinement  
    SIMPLE_INVERSION,
    ACCURATE_INVERSION,           // with iterative refinement
    MULTIPLICATION,
  };

  /** Types of matrices being tested */
  enum MatrixTypes {
    JAMA,
    DOUBLE_MATRIX,                // com.mvohm.quadmatrix.DoubleMatrix                                
    QUADRUPLE_MATRIX,             // com.mvohm.quadmatrix.QuadrupleMatrix
    BIGDECIMAL_MATRIX_40,         // com.mvohm.quadmatrix.BigDecimalMatrix with the matrix precision set to 40 decimal digits 
    BIGDECIMAL_MATRIX_80,         // com.mvohm.quadmatrix.BigDecimalMatrix with the matrix precision set to 80 decimal digits
  };

  /** Sizes of matrices to test */
  int [] sizes = new int[] {
     50,
    100,
    200,
    400,
  };

  static final long WARMUP_TIME =     3_000; // Start timing after warmup of 3 seconds after the start of the execution
  static final int WARMUP_COUNT =     1_000; // Start timing after warmup of 1000 test iterations, if they take less than WARMUP_TIME 
  static final int MAX_ITERATIONS = 100_000; // For fast methods, do max 100,000 iterations
  static final int MIN_ITERATIONS =      20; // Do at least 20 iterations
  static final long MAXTIME_MS =     30_000; // max 30 seconds per every type + operation;

  interface TesterMaker {
    OperationTester make(int size);
  }

  interface GeneratorMaker {
    DataGenerator make(int size);
  }

  /**
   * For each operation type and size, stores the corresponding method for dataset generation.
   * Remove the types of operations you do not want to perform.
   */
  HashMap<Operations, GeneratorMaker> generatorMakers = new HashMap<>()
  {{
    put(Operations.SIMPLE_VECTOR_SOLUTION,    size -> () -> MatrixData.makeDataSetForVectorSolutions(size, random));
    put(Operations.ACCURATE_VECTOR_SOLUTION,  size -> () -> MatrixData.makeDataSetForVectorSolutions(size, random));
    put(Operations.SIMPLE_SPD_SOLUTION,       size -> () -> MatrixData.makeDataSetForSPDSolutions(size, random));
    put(Operations.ACCURATE_SPD_SOLUTION,     size -> () -> MatrixData.makeDataSetForSPDSolutions(size, random));
    put(Operations.SIMPLE_MATRIX_SOLUTION,    size -> () -> MatrixData.makeDataSetForMatrixSolutions(size, random));
    put(Operations.ACCURATE_MATRIX_SOLUTION,  size -> () -> MatrixData.makeDataSetForMatrixSolutions(size, random));
    put(Operations.SIMPLE_INVERSION,          size -> () -> MatrixData.makeDataSetForInversions(size, random));
    put(Operations.ACCURATE_INVERSION,        size -> () -> MatrixData.makeDataSetForInversions(size, random));
    put(Operations.MULTIPLICATION,            size -> () -> MatrixData.makeDataSetForMatrixSolutions(size, random));
  }};

  /**
   * For each operation type and matrix type, stores the corresponding method that 
   * performs the operation and collects statistics regarding execution time and errors.
   * Remove the types of operations and matrices you do not want to test.
   */
  HashMap<Operations, HashMap<MatrixTypes, OperationPerformer>> performers = new HashMap<>() {{
    put(Operations.SIMPLE_VECTOR_SOLUTION, new HashMap<>() {{
      put(MatrixTypes.JAMA,                 MatrixData::jamaLuSolutionErrors);
      put(MatrixTypes.DOUBLE_MATRIX,        MatrixData::doubleLuSolutionWithScalingErrors);
      put(MatrixTypes.QUADRUPLE_MATRIX,     MatrixData::quadrupleLuSolutionWithScalingErrors);
      put(MatrixTypes.BIGDECIMAL_MATRIX_40, MatrixData::bigDecimalLuSolutionWithScalingErrors);
      put(MatrixTypes.BIGDECIMAL_MATRIX_80, MatrixData::bigDecimalLuSolutionWithScalingErrors);
    }});
    put(Operations.ACCURATE_VECTOR_SOLUTION, new HashMap<>() {{
      put(MatrixTypes.DOUBLE_MATRIX,        MatrixData::doubleAccurateLUSolutionWithScalingErrors);
      put(MatrixTypes.QUADRUPLE_MATRIX,     MatrixData::quadrupleAccurateLUSolutionWithScalingErrors);
      put(MatrixTypes.BIGDECIMAL_MATRIX_40, MatrixData::bigDecimalAccurateLUSolutionWithScalingErrors);
      put(MatrixTypes.BIGDECIMAL_MATRIX_80, MatrixData::bigDecimalAccurateLUSolutionWithScalingErrors);
    }});
    put(Operations.SIMPLE_SPD_SOLUTION, new HashMap<>() {{
      put(MatrixTypes.JAMA,                 MatrixData::jamaSpdSolutionErrors);
      put(MatrixTypes.DOUBLE_MATRIX,        MatrixData::doubleSpdSolutionErrors);
      put(MatrixTypes.QUADRUPLE_MATRIX,     MatrixData::quadrupleSPDSolutionErrors);
      put(MatrixTypes.BIGDECIMAL_MATRIX_40, MatrixData::bigDecimalSPDSolutionErrors);
      put(MatrixTypes.BIGDECIMAL_MATRIX_80, MatrixData::bigDecimalSPDSolutionErrors);
    }});
    put(Operations.ACCURATE_SPD_SOLUTION, new HashMap<>() {{
      put(MatrixTypes.DOUBLE_MATRIX,        MatrixData::doubleAccurateSPDSolutionErrors);
      put(MatrixTypes.QUADRUPLE_MATRIX,     MatrixData::quadrupleAccurateSPDSolutionErrors);
      put(MatrixTypes.BIGDECIMAL_MATRIX_40, MatrixData::bigDecimalAccurateSPDSolutionErrors);
      put(MatrixTypes.BIGDECIMAL_MATRIX_80, MatrixData::bigDecimalAccurateSPDSolutionErrors);
    }});
    put(Operations.SIMPLE_MATRIX_SOLUTION, new HashMap<>() {{
      put(MatrixTypes.JAMA,                 MatrixData::jamaMatrixSolutionErrors);
      put(MatrixTypes.DOUBLE_MATRIX,        MatrixData::doubleMatrixSolutionErrors);
      put(MatrixTypes.QUADRUPLE_MATRIX,     MatrixData::quadrupleMatrixSolutionErrors);
      put(MatrixTypes.BIGDECIMAL_MATRIX_40, MatrixData::bigDecimalMatrixSolutionErrors);
      put(MatrixTypes.BIGDECIMAL_MATRIX_80, MatrixData::bigDecimalMatrixSolutionErrors);
    }});
    put(Operations.ACCURATE_MATRIX_SOLUTION, new HashMap<>() {{
      put(MatrixTypes.DOUBLE_MATRIX,        MatrixData::doubleAccurateMatrixSolutionErrors);
      put(MatrixTypes.QUADRUPLE_MATRIX,     MatrixData::quadrupleAccurateMatrixSolutionErrors);
      put(MatrixTypes.BIGDECIMAL_MATRIX_40, MatrixData::bigDecimalAccurateMatrixSolutionErrors);
      put(MatrixTypes.BIGDECIMAL_MATRIX_80, MatrixData::bigDecimalAccurateMatrixSolutionErrors);
    }});
    put(Operations.SIMPLE_INVERSION, new HashMap<>() {{
      put(MatrixTypes.JAMA,                 MatrixData::jamaMatrixInversionErrors);
      put(MatrixTypes.DOUBLE_MATRIX,        MatrixData::doubleMatrixInversionErrors);
      put(MatrixTypes.QUADRUPLE_MATRIX,     MatrixData::quadrupleMatrixInversionErrors);
      put(MatrixTypes.BIGDECIMAL_MATRIX_40, MatrixData::bigDecimalMatrixInversionErrors);
      put(MatrixTypes.BIGDECIMAL_MATRIX_80, MatrixData::bigDecimalMatrixInversionErrors);
    }});
    put(Operations.ACCURATE_INVERSION, new HashMap<>() {{
      put(MatrixTypes.DOUBLE_MATRIX,        MatrixData::doubleAccurateMatrixInversionErrors);
      put(MatrixTypes.QUADRUPLE_MATRIX,     MatrixData::quadrupleAccurateMatrixInversionErrors);
      put(MatrixTypes.BIGDECIMAL_MATRIX_40, MatrixData::bigDecimalAccurateMatrixInversionErrors);
      put(MatrixTypes.BIGDECIMAL_MATRIX_80, MatrixData::bigDecimalAccurateMatrixInversionErrors);
    }});
    put(Operations.MULTIPLICATION, new HashMap<>() {{
      put(MatrixTypes.JAMA,                 MatrixData::jamaMultiplicationErrors);
      put(MatrixTypes.DOUBLE_MATRIX,        MatrixData::doubleMultiplicationErrors);
      put(MatrixTypes.QUADRUPLE_MATRIX,     MatrixData::quadrupleMultiplicationErrors);
      put(MatrixTypes.BIGDECIMAL_MATRIX_40, MatrixData::bigDecimalMultiplicationErrors);
      put(MatrixTypes.BIGDECIMAL_MATRIX_80, MatrixData::bigDecimalMultiplicationErrors);
    }});
  }};

  private PrintStream output = null;

  public static void main(String[] args) throws IOException {
    new CollectStatistics().run();
  }

  /** Traverses through operation types, matrix types and sizes, */ 
  private void run() throws IOException {
    Locale.setDefault(Locale.US);
    output = openOutput();
    for (final Operations operation: Operations.values()) {
      testOperation(operation);
    }
    output.close();
    say("Done!");
  }

  /**
   * Perform the specified operation on all types of matrices of all sizes
   * @param operation
   * @param output
   */
  private void testOperation(Operations operation) {
    for (final MatrixTypes matrixType: MatrixTypes.values()) {
      testOperationOnType(operation, matrixType);
    }
    write();
  }

  /**
   * Perform the specified operation on the specified type of matrix of all sizes
   * @param operation
   * @param output
   */
  private void testOperationOnType(Operations operation, MatrixTypes matrixType) {
    setBigDecimalMatrixPrecision(matrixType);

    final ErrorSet[] results = collectStatsOnSizes(operation, matrixType);

    if (results != null) {
      writeResults(results, operation, matrixType);
    }
    write();
  }

  private static void setBigDecimalMatrixPrecision(MatrixTypes matrixType) {
    if (matrixType == MatrixTypes.BIGDECIMAL_MATRIX_40) {
      BigDecimalMatrix.setDefaultPrecision(40);
    } else if (matrixType == MatrixTypes.BIGDECIMAL_MATRIX_80) {
      BigDecimalMatrix.setDefaultPrecision(80);
    }
  }

  /**
   * Collect statistics in respect of errors and times for different sizes 
   * of matrices of the specified type, performing the specified type of operation 
   */
  private ErrorSet[] collectStatsOnSizes(Operations operation, MatrixTypes matrixType) {
    final ErrorSet[] results = new ErrorSet[sizes.length];

    for (int i = 0; i < sizes.length; i++) {
      results[i] = testOperationOnTypeOfSize(operation, matrixType, sizes[i]);
      if (results[i] == null) {
        write(operation + " for " + matrixType + " not implemented");
        return null;
      } else {
        showResults(results[i]);
      }
    }
    return results;
  }

  /** Creates a tester object for the given combination of the operation type, matrix type, and size,
   * and runs it.
   */
  private ErrorSet testOperationOnTypeOfSize(Operations operation, MatrixTypes matrixType, int size) {
    final OperationTester tester = makeTester(operation, matrixType, size);
    if (tester == null) {
      return null;
    }
    showSectionHeader(operation, matrixType, size);
    runTester(tester);
    return tester.getStatistics();
  }

  private static void showSectionHeader(Operations operation, MatrixTypes matrixType, int size) {
    say("Operation:     " + operation);
    say("  Matrix type: " + matrixType);
    say("    Matrix size = %4s", size);
  }


  private OperationTester makeTester(Operations operation, MatrixTypes matrixType, int size) {
    final GeneratorMaker generatorMaker = generatorMakers.get(operation);
    final HashMap<MatrixTypes, OperationPerformer> performerTable = performers .get(operation);
    final OperationPerformer performer = performerTable == null? null: performerTable.get(matrixType);
    if (generatorMaker == null || performer == null)
      return null;

    final DataGenerator generator = generatorMaker.make(size);
    return new OperationTester(generator, performer);
  }

  private void runTester(final OperationTester tester) {
    resetTime();
    long lastTime = 0;

    // run it not more than ITERATIONS times 
    // and not longer than MAXTIME_MS milliseconds,
    // but not less than MIN_ITERATIONS times 
    for (int i = 1; i <= MAX_ITERATIONS; i++) {
      tester.perform();
      final long currentTime = System.currentTimeMillis();
      if (currentTime - lastTime > 2000) {  // Show progress every 2 sec
        showProgress(tester, i);
        lastTime = currentTime;
      }
      if (elapsedTime() > MAXTIME_MS && i >= MIN_ITERATIONS)
        break;
    }
  }

  private static void showProgress(final OperationTester tester, int i) {
    say("      %,6d: Max error: %7.3e, MSE: %7.3e, time: %,11.3f ms",
        i, tester.getLastMaxError(), tester.getLastMse(), tester.getLastTimeMs());
  }

  private static void showResults(ErrorSet result) {
    say("========================");
    say("     Average: Max error: %7.3e, MSE: %7.3e, time: %,11.3f ms\n",
        result.maxError(), result.mse(), result.getTime() * 1e-6);
  }

  private void writeResults(ErrorSet[] results, Operations operation, MatrixTypes matrixType) {
    write("Statistics for %s on %s", operation, matrixType);
    write_("    Size:    ");
    for (int i = 0; i < results.length; i++) {
      write_("\t%12s", sizes[i]);
    }
    write();
    write_("    Errors:  ");
    for (int i = 0; i < results.length; i++) {
      write_("\t%12.3e", results[i].mse());
    }
    write();
    write_("    Time, ms:");
    for (int i = 0; i < results.length; i++) {
      write_("\t%12.3f", results[i].getTime() * 1e-6);
    }
    write();
  }

  private long testStartTime;

  private void resetTime() {
    testStartTime = System.currentTimeMillis();
  }

  /** Returns the amount of time passed since the beginning of the current test, in milliseconds */
  private long elapsedTime() {
    return System.currentTimeMillis() - testStartTime;
  }

  private static PrintStream openOutput() throws IOException {
    final String dataPath = System.getProperty("user.dir") + "\\Results\\";
    // Bug fix 2024-12-30 18:53:43: Create the folder if it does not exist
    final Path path = Paths.get(dataPath, new String[] {});
    if (!Files.exists(path)) {
      Files.createDirectory(path);
    }
    final SimpleDateFormat dateFormat = new SimpleDateFormat("yyMMdd_HHmm");
    final Date now = new Date();
    final String dateStr = dateFormat.format(now);
    final String fileName = dataPath + "stats_" + dateStr + ".txt";
    final PrintStream out = new PrintStream(fileName);
    return out;
  };

  @SuppressWarnings("unused")
  private void write(String format, Object ... args) {
    if (output != null) {
      output.printf(format + "\n", args);
    }
    say(format, args);
  }

  private void write(Object data) {
    if (output != null) {
      output.println(data);
    }
    say(data);
  }

  private void write_(String format, Object ... args) {
    if (output != null) {
      output.printf(format, args);
    }
    say_(format, args);
  }

  private void write_(Object data) {
    if (output != null) {
      output.print(data);
    }
    say_(data);
  }

  private void write() {
    if (output != null) {
      output.println();
    }
    say();
  }

  /* ***************************************************************************
   ***** OperationTester *******************************************************
   *****************************************************************************/
  
  /**
   * An instance of this class encapsulates objects for generating test data 
   * and performing the method under test on that data, 
   * an object of type MatrixData that stores the test data, 
   * and an object of type ErrorSet that keeps the results of the most recent 
   * performed operation (execution time and error values) 
   * 
   * Method perform() generates a new random data set and performs the tested operation that data set.    
   */
  private class OperationTester  {

    private DataGenerator generator;
    private OperationPerformer performer;

    private MatrixData matrixData;
    private ErrorSet errorSet;

    private int trialCount;
    private int timedTrialCount;

    private double accumulatedMse;
    private double accumulatedMeanErr;
    private double accumulatedMaxErr;
    private long accumulatedTime;

    public OperationTester(DataGenerator generator, OperationPerformer performer) {
      this.generator = generator;
      this.performer = performer;
    }

    public double getLastTimeMs() {
      return errorSet.getTime() * 1e-6;
    }

    public double getLastMse() {
      return errorSet.mse();
    }

    public double getLastMaxError() {
      return errorSet.maxError();
    }

    public void perform() {
      trialCount++;
      matrixData = generator.generate();

      errorSet = performer.perform(matrixData);
      accumulatedMse += errorSet.mse();
      accumulatedMaxErr = Math.max(accumulatedMaxErr, errorSet.maxError());

      if (    trialCount > WARMUP_COUNT
          || (elapsedTime() > WARMUP_TIME && trialCount > 1)) {
        accumulatedTime += errorSet.getTime();
        timedTrialCount++;
      }
    }

    public ErrorSet getStatistics() {
      final double avrMse = accumulatedMse / trialCount;
      final double avrMaxErr = accumulatedMaxErr;
      final double avrMeanErr = accumulatedMeanErr / trialCount;
      final long avrTime = Math.round((double)accumulatedTime / timedTrialCount);
      return new ErrorSet(avrMse, avrMeanErr, avrMaxErr).setTime(avrTime);
    }

  } // private class OperationTester

}
