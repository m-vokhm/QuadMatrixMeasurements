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

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.Locale;
import java.util.Random;

import javax.sound.sampled.AudioFormat.Encoding;

import com.mvohm.quadmatrix.BigDecimalMatrix;

public class CollectStatistics {

  static final Random random = new Random();

  public interface DataGenerator {
    MatrixData generate();
  }

  public interface OperationPerformer {
    ErrorSet perform(MatrixData matrixData);
  }



  enum Operations {
    SIMPLE_VECTOR_SOLUTION,
    ACCURATE_VECTOR_SOLUTION,
    SIMPLE_SPD_SOLUTION,
    ACCURATE_SPD_SOLUTION,
    SIMPLE_MATRIX_SOLUTION,
    ACCURATE_MATRIX_SOLUTION,
    SIMPLE_INVERSION,
    ACCURATE_INVERSION,
    MULTIPLICATION,
  };

  enum MatrixTypes {
    JAMA,
    DOUBLE_MATRIX,
    QUADRUPLE_MATRIX,
    BIGDECIMAL_MATRIX_40,
    BIGDECIMAL_MATRIX_80,
  };

  int [] sizes = new int[] {
     50,
    100,
//    200,
//    400,
  };

  static final long WARMUP_TIME =   3_000; // Start timing after warmup of 3 seconds after the start of the execution
  static final int WARMUP_COUNT =   1_000; // Start timing after warmup of 1000 test iterations, if they take less than 3 a
  static final int ITERATIONS =    10_000; // For fast methods, do at least 10,000 iterations
  static final int MIN_ITERATIONS =    10; // Do at least 10 iterations
  static final long MAXTIME_MS =   20_000; // max 20 seconds per every type + operation;

  interface TesterMaker {
    OperationTester make(int size);
  }

  interface GeneratorMaker {
    DataGenerator make(int size);
  }

  /**
   * For each operation type, stores the corresponding method for data generation
   */
  HashMap<Operations, GeneratorMaker> generatorMakers = new HashMap<>()
  {{
    put(Operations.SIMPLE_VECTOR_SOLUTION, size -> () -> MatrixData.makeDataSetForVectorSolutions(size, random));
    put(Operations.MULTIPLICATION, size -> () -> MatrixData.makeDataSetForMatrixSolutions(size, random));
  }};

  /**
   * For each operation type and matrix type, stores the corresponding method performing the operation to be tested
   */
  HashMap<Operations, HashMap<MatrixTypes, OperationPerformer>> performers = new HashMap<>() {{
//    put(Operations.SIMPLE_VECTOR_SOLUTION, new HashMap<>() {{
//      put(MatrixTypes.DOUBLE_MATRIX,        MatrixData::luSolutionWithScalingErrors);
//      put(MatrixTypes.QUADRUPLE_MATRIX,     MatrixData::quadrupleLuSolutionWithScalingErrors);
//      put(MatrixTypes.BIGDECIMAL_MATRIX_40, MatrixData::bigDecimalLuSolutionWithScalingErrors);
//      put(MatrixTypes.BIGDECIMAL_MATRIX_80, MatrixData::bigDecimalLuSolutionWithScalingErrors);
//    }});
//    put(Operations.ACCURATE_VECTOR_SOLUTION, new HashMap<>() {{
//      put(MatrixTypes.DOUBLE_MATRIX,        MatrixData::accurateLUSolutionWithScalingErrors);
//      put(MatrixTypes.QUADRUPLE_MATRIX,     MatrixData::quadrupleAccurateLUSolutionWithScalingErrors);
//      put(MatrixTypes.BIGDECIMAL_MATRIX_40, MatrixData::bigDecimalAccurateLUSolutionWithScalingErrors);
//      put(MatrixTypes.BIGDECIMAL_MATRIX_80, MatrixData::bigDecimalAccurateLUSolutionWithScalingErrors);
//    }});
//    put(Operations.SIMPLE_SPD_SOLUTION, new HashMap<>() {{
//      put(MatrixTypes.DOUBLE_MATRIX,        MatrixData::spdSolutionErrors);
//      put(MatrixTypes.QUADRUPLE_MATRIX,     MatrixData::quadrupleSPDSolutionErrors);
//      put(MatrixTypes.BIGDECIMAL_MATRIX_40, MatrixData::bigDecimalSPDSolutionErrors);
//      put(MatrixTypes.BIGDECIMAL_MATRIX_80, MatrixData::bigDecimalSPDSolutionErrors);
//    }});
//    put(Operations.ACCURATE_SPD_SOLUTION, new HashMap<>() {{
//      put(MatrixTypes.DOUBLE_MATRIX,        MatrixData::accurateSPDSolutionErrors);
//      put(MatrixTypes.QUADRUPLE_MATRIX,     MatrixData::quadrupleAccurateSPDSolutionErrors);
//      put(MatrixTypes.BIGDECIMAL_MATRIX_40, MatrixData::bigDecimalAccurateSPDSolutionErrors);
//      put(MatrixTypes.BIGDECIMAL_MATRIX_80, MatrixData::bigDecimalAccurateSPDSolutionErrors);
//    }});
//    put(Operations.SIMPLE_MATRIX_SOLUTION, new HashMap<>() {{
//      put(MatrixTypes.DOUBLE_MATRIX,        MatrixData::matrixSolutionErrors);
//      put(MatrixTypes.QUADRUPLE_MATRIX,     MatrixData::quadrupleMatrixSolutionErrors);
//      put(MatrixTypes.BIGDECIMAL_MATRIX_40, MatrixData::bigDecimalMatrixSolutionErrors);
//      put(MatrixTypes.BIGDECIMAL_MATRIX_80, MatrixData::bigDecimalMatrixSolutionErrors);
//    }});
//    put(Operations.ACCURATE_MATRIX_SOLUTION, new HashMap<>() {{
//      put(MatrixTypes.DOUBLE_MATRIX,        MatrixData::accurateMatrixSolutionErrors);
//      put(MatrixTypes.QUADRUPLE_MATRIX,     MatrixData::quadrupleAccurateMatrixSolutionErrors);
//      put(MatrixTypes.BIGDECIMAL_MATRIX_40, MatrixData::bigDecimalAccurateMatrixSolutionErrors);
//      put(MatrixTypes.BIGDECIMAL_MATRIX_80, MatrixData::bigDecimalAccurateMatrixSolutionErrors);
//    }});
//    put(Operations.SIMPLE_INVERSION, new HashMap<>() {{
//      put(MatrixTypes.DOUBLE_MATRIX,        MatrixData::matrixInversionErrors);
//      put(MatrixTypes.QUADRUPLE_MATRIX,     MatrixData::quadrupleMatrixInversionErrors);
//      put(MatrixTypes.BIGDECIMAL_MATRIX_40, MatrixData::bigDecimalMatrixInversionErrors);
//      put(MatrixTypes.BIGDECIMAL_MATRIX_80, MatrixData::bigDecimalMatrixInversionErrors);
//    }});
//    put(Operations.ACCURATE_INVERSION, new HashMap<>() {{
//      put(MatrixTypes.DOUBLE_MATRIX,        MatrixData::accurateMatrixInversionErrors);
//      put(MatrixTypes.QUADRUPLE_MATRIX,     MatrixData::quadrupleAccurateMatrixInversionErrors);
//      put(MatrixTypes.BIGDECIMAL_MATRIX_40, MatrixData::bigDecimalAccurateMatrixInversionErrors);
//      put(MatrixTypes.BIGDECIMAL_MATRIX_80, MatrixData::bigDecimalAccurateMatrixInversionErrors);
//    }});
    put(Operations.MULTIPLICATION, new HashMap<>() {{
//      put(MatrixTypes.DOUBLE_MATRIX,        MatrixData::multiplicationErrors);
//      put(MatrixTypes.QUADRUPLE_MATRIX,     MatrixData::quadrupleMultiplicationErrors);
      put(MatrixTypes.BIGDECIMAL_MATRIX_40, MatrixData::bigDecimalMultiplicationErrors);
      put(MatrixTypes.BIGDECIMAL_MATRIX_80, MatrixData::bigDecimalMultiplicationErrors);
    }});
  }};


  private PrintStream output = null;

  public static void main(String[] args) throws FileNotFoundException {
    new CollectStatistics().run();
  }

  private void run() throws FileNotFoundException {
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
    write("Operation: " + operation);
    write("  Matrix type: " + matrixType);

    setBigDecimalMatrixPrecision(matrixType);

    final ErrorSet[] results = collectStatsOnSizes(operation, matrixType);

    if (results != null) {
      writeResults(results);
    }
    write();
  }

  private ErrorSet[] collectStatsOnSizes(Operations operation, MatrixTypes matrixType) {
    final ErrorSet[] results = new ErrorSet[sizes.length];

    // Collect errors and times for different sizes
    for (int i = 0; i < sizes.length; i++) {
      results[i] = testOperationOnTypeOfSize(operation, matrixType, sizes[i]);
      if (results[i] == null) {
        write("  Not implemented");
        return null;
      } else {
        showResults(results[i]);
      }
    }
    return results;
  }

  private static void setBigDecimalMatrixPrecision(MatrixTypes matrixType) {
    if (matrixType == MatrixTypes.BIGDECIMAL_MATRIX_40) {
      BigDecimalMatrix.setDefaultPrecision(40);
    } else if (matrixType == MatrixTypes.BIGDECIMAL_MATRIX_80) {
      BigDecimalMatrix.setDefaultPrecision(80);
    }
  }

  private void writeResults(ErrorSet[] results) {
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

  private ErrorSet testOperationOnTypeOfSize(Operations operation, MatrixTypes matrixType, int size) {
    final OperationTester tester = makeTester(operation, matrixType, size);
    if (tester == null) {
      return null;
    }
    say("    Matrix size = %4s", size);
    runTester(tester);
    return tester.getStatistics();
  }


  TesterMaker petMaker =
      size -> new OperationTester (
                () -> MatrixData.makeDataSetForVectorSolutions(size, random),
                MatrixData::quadrupleLuSolutionWithScalingErrors);


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

    // run it ITERATIONS times or MAXTIME milliseconds, which comes first
    for (int i = 1; i <= ITERATIONS; i++) {
      tester.perform();
      final long currentTime = System.currentTimeMillis();
      if (currentTime - lastTime > 1000) {
        showProgress(tester, i);
        lastTime = currentTime;
      }
      if (elapsedTime() > MAXTIME_MS && i >= MIN_ITERATIONS)
        break;
    }
  }

  private static void showProgress(final OperationTester tester, int i) {
    say("        %4d: Max error: %7.3e, MSE: %7.3e, time: %,9.3f ms",
        i, tester.getLastMaxError(), tester.getLastMse(), tester.getLastTimeMs());
  }

  private static void showResults(ErrorSet result) {
    say("========================");
    say("     Average: Max error: %7.3e, MSE: %7.3e, time: %,9.3f ms\n",
        result.maxError(), result.mse(), result.getTime() * 1e-6);
  }

  private long testStartTime;

  private void resetTime() {
    testStartTime = System.currentTimeMillis();
  }

  /** Return amount of time passed since the beginning of the current test, in milliseconds */
  private long elapsedTime() {
    return System.currentTimeMillis() - testStartTime;
  }

  private static PrintStream openOutput() throws FileNotFoundException {
    final String dataPath = System.getProperty("user.dir") + "\\Results\\";
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

  private class OperationTester  {

    DataGenerator generator;
    OperationPerformer performer;

    MatrixData matrixData;
    ErrorSet errorSet;

    int trialCount;
    int timedTrialCount;

    double accumulatedMse;
    double accumulatedMeanErr;
    double accumulatedMaxErr;
    long accumulatedTime;

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
      final double avrMaxErr = accumulatedMaxErr / trialCount;
      final double avrMeanErr = accumulatedMeanErr / trialCount;
      final long avrTime = Math.round((double)accumulatedTime / timedTrialCount);
      return new ErrorSet(avrMse, avrMeanErr, avrMaxErr).setTime(avrTime);
    }

  }

}
