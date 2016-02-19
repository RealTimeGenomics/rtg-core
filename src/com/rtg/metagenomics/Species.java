/*
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.metagenomics;


import com.rtg.launcher.GlobalFlags;
import com.rtg.metagenomics.matrix.Matrix;
import com.rtg.metagenomics.matrix.MatrixSymmetric;
import com.rtg.metagenomics.matrix.MatrixUtils;
import com.rtg.metagenomics.matrix.Vector;
import com.rtg.metagenomics.matrix.VectorSimple;
import com.rtg.util.Pair;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

import Jama.EigenvalueDecomposition;

/**
 * Generates a list of estimates of the frequency of each species.
 */
public class Species extends IntegralAbstract {
  /** Used to detect when L increases */
  private static final double INCREASE_THRESHOLD = -0.1;  //TODO something sensible about this arbitrary constant

  private static final double L_TERMINATION = (Double) GlobalFlags.getFlag(GlobalFlags.SPECIES_LTERMINATION_FLAG).getValue();
  private static final Double L_TERMINATION_TARGET = (Double) GlobalFlags.getFlag(GlobalFlags.SPECIES_TERMINATION_TARGET_FLAG).getValue();

  private static final int[] EMPTY = new int[0];

  private final BlockInfo mBlockInfo;

  private final int[][] mMembersOf;

  private double mLastL;
  private Vector mLastR;
  private final SimpleTerminator mSimpleTerminator;

  /**
   * Primarily for testing - creates a default taxonomy membership assuming no taxonomy structure.
   * @param blockInfo parameters for the block currently being processed
   */
  public Species(BlockInfo blockInfo) {
    this(makeFlatMembership(blockInfo.getGlobalSpeciesMap() == null ? blockInfo.getSpeciesMap().size() : blockInfo.getGlobalSpeciesMap().size()), blockInfo);
  }

  /**
   *
   * @param membersOf list of local genome ids members for each global id species
   * @param blockInfo parameters for the block currently being processed
   */
  public Species(int[][] membersOf, BlockInfo blockInfo) {
    if (membersOf == null) {
      throw new NullPointerException();
    }
    mBlockInfo = blockInfo;
    mMembersOf = membersOf;
    mSimpleTerminator = new SimpleTerminator(ltermination());
  }

  private double ltermination() {
    return L_TERMINATION;
  }

  /**
   * @param minIter minimum number of iterations
   * @return result of species calculation
   */
  public SubBlockResult solve(final int minIter) {
    final int numSpecies = mBlockInfo.getN();
    final Vector m = new VectorSimple(numSpecies);
    for (final Frag frag : mBlockInfo.getFrags()) {
      frag.init(m);
    }
    final Vector initialR = new VectorSimple(numSpecies);
    //initialize r
    for (int i = 0; i < numSpecies; i++) {
      final long l = mBlockInfo.getGenomeLength(i);
      if (l == 0) {
        Diagnostic.warning("Could not determine length of taxonId \"" + mBlockInfo.getTaxonId(i) + "\", its either not present or the length is 0");
        continue;
      }
      final double rv = m.get(i) / l;
      //rv should be valid so long as m.get(i) is finite and positive (see comments in Frag.init).
      initialR.set(i, rv);
    }
    solve(initialR, EMPTY, minIter, mSimpleTerminator);
    final Matrix hessian = hessian(mLastR);
    //System.err.println("Hessian:");
    //Extract eigenvectors and eigenvalues.
    final long startTime = System.currentTimeMillis();
    final EigenvalueDecomposition ed = hessian.toJama().eig();
    final long endTime = System.currentTimeMillis();
    Diagnostic.developerLog("B:" + mBlockInfo.id() + " EigenValue Decomp took: " + ((endTime - startTime) / 1000) + "s" + " Hessian dimensions: " + hessian.dimension());
    final int totalGenomes =  mBlockInfo.getGlobalSpeciesMap() == null ? mBlockInfo.getSpeciesMap().size() : mBlockInfo.getGlobalSpeciesMap().size();
    return new SubBlockResult(mLastR, variance(mMembersOf, mBlockInfo, mLastR, ed), new VectorSimple(totalGenomes), mLastL);

  }

  private boolean checkRemainedFixed(final Vector v, final int[] fixedIds) {
    for (final int i : fixedIds) {
      final long l = mBlockInfo.getGenomeLength(i);
      if (l == 0) {
        continue;
      }
      final double rv = 1.0 / l;
      if (Math.abs(v.get(i) - rv) > 1e-12) {
        return false;
      }
    }
    return true;
  }

  /**
   * Perform objective function solving with some genomes fixed.
   * @param previousR vector containing existing solution.
   * @param fixedIds list of local genome ids to be fixed
   * @param minIter minimum number of iterations
   * @param targetL target L value
   * @return the final L value
   */
  public double solveFixed(Vector previousR, int[] fixedIds, final int minIter, final double targetL) {
    final Vector initialR = new VectorSimple(previousR);
    //initialize r
    for (final int i : fixedIds) {
      final long l = mBlockInfo.getGenomeLength(i);
      if (l == 0) {
        Diagnostic.warning("Could not determine length of taxonId \"" + mBlockInfo.getTaxonId(i) + "\", its either not present or the length is 0");
        continue;
      }
      final double rv = 1.0 / l;
      //rv should be valid so long as m.get(i) is finite and positive (see comments in Frag.init).
      initialR.set(i, rv);
    }
    final Terminator terminator = L_TERMINATION_TARGET == null ? mSimpleTerminator : new TargetTerminator(ltermination(), L_TERMINATION_TARGET, targetL);
    solve(initialR, fixedIds, minIter, terminator);
    assert checkRemainedFixed(mLastR, fixedIds);
    return mLastL;
  }

  private void solve(Vector initialR, int[] fixedIds, final int minIter, final Terminator terminator) {
    if (!MatrixUtils.isFinitePositive(initialR)) {
      message("Invalid initial estimate:", initialR);
      throw new RuntimeException("Invalid initial estimate.");
    }
    mLastR = initialR;
    final int n = mBlockInfo.getN();
    assert n >= 1;
    int iter = 0;
    //compute first derivative (in log space)
    //TODO if calculation of jacobian fails then return the current x.
    //TODO if calculation of delta or x fails then return the last x
    final Pair<Vector, Double> pair0 = jacobian(mLastR, fixedIds);
    Vector jacobian = pair0.getA();
    mLastL = pair0.getB();
    if (n == 1) {
      return; //avoid warning messages
    }

    if (!MatrixUtils.isFinite(jacobian)) {
      message("Invalid initial derivatives:", jacobian);
    } else {
      Vector lastJacobian = jacobian;
      Vector delta = MatrixUtils.negative(jacobian);
      double smooth = 0.0;
      int floundered = 0;
      while (true) {
        final String msg = "B:" + mBlockInfo.id() + " Iteration:" + iter;
        if (BlockInfo.VERY_VERBOSE) {
          Diagnostic.developerLog(msg + " R: " + mLastR.toString());
        }
        final Pair<Vector, Double> solved = solveLine(mLastR, delta);
        if (solved == null) {
          Diagnostic.developerLog(msg + " Solve line unable to make an estimate");
          break;
        }
        final Vector rnew = solved.getA();
        if (!MatrixUtils.isFinitePositive(rnew)) {
          message(msg + " Invalid estimate:", rnew);
          break;
        }
        final Pair<Vector, Double> pair = jacobian(rnew, fixedIds);
        final Double ll = pair.getB();
        final double ldelta;
        final double ldelta0 = mLastL - ll;
        if (ldelta0 < INCREASE_THRESHOLD) {
          //System.err.println("L increased.");
          Diagnostic.developerLog(msg + " L increased." + " LastL:" + mLastL + " L:" + ll);
          bug(mLastR, delta, solved.getB());
          //try using more robust (but probably slower) solver
          final Pair<Vector, Double> solvedR = solveLineRobust(mLastR, delta);
          final Vector rnewR = solvedR.getA();
          if (!MatrixUtils.isFinitePositive(rnewR)) {
            message(msg + " Invalid estimate:", rnewR);
            break;
          }
          final Pair<Vector, Double> pairR = jacobian(rnewR, fixedIds);
          final Double llR = pairR.getB();
          final double ldelta0R = mLastL - llR;
          if (ldelta0R < INCREASE_THRESHOLD) {
            if (floundered > 2) {
              break;
            }
            floundered++;
            //System.err.println("L increased again.");
            Diagnostic.developerLog(msg + " L increased again." + " LastL:" + mLastL + " L:" + llR);
            bug(mLastR, delta, solvedR.getB());
            //leave polakRibiere to determine a new direction - ?? not sure if this works
            ldelta = 0.0;
          } else {
            ldelta = ldelta0R < 0.0 ? 0.0 : ldelta0R;
            lastJacobian = jacobian;
            jacobian = pairR.getA();
            mLastR = rnewR; //must be set before the breaks below
            mLastL = llR;
          }
        } else {
          floundered = 0;
          ldelta = ldelta0 < 0.0 ? 0.0 : ldelta0;
          lastJacobian = jacobian;
          jacobian = pair.getA();
          mLastR = rnew; //must be set before the breaks below
          mLastL = ll;
        }

        final double s = 1.0 / Math.sqrt(iter + 1);
        smooth = ldelta * s + (1 - s) * smooth;
        final double est = smooth * (iter + 1);
        Diagnostic.developerLog(msg + " L: " + Utils.realFormat(mLastL, 8) + " ldelta: " + ldelta + " smooth: " + Utils.realFormat(smooth, 8) + " est: " + Utils.realFormat(est, 8));

        if (iter >= minIter && terminator.terminate(est, mLastL)) {
          Diagnostic.developerLog("B:" + mBlockInfo.id() + " Termination minIter: " + minIter);
          if (BlockInfo.VERY_VERBOSE) {
            Diagnostic.developerLog("R: " + rnew.toString());
          }
          break;
        }

        delta = polakRibiere(lastJacobian, jacobian, delta);
        if (!MatrixUtils.isFinite(delta)) {
          message("B:" + mBlockInfo.id() + " Invalid delta:", delta);
          break;
        }
        if (BlockInfo.VERY_VERBOSE) {
          Diagnostic.developerLog("B:" + mBlockInfo.id() + " D: " + delta.toString());
        }
        iter++;
      }
    }
  }

  private void bug(final Vector r, final Vector deltaV, final double d) {
    final int intDelta = 10;
    final double delta = d / intDelta;
    final Line line = new LLine(r, deltaV, mBlockInfo);
    for (int i = -intDelta; i <= 2 * intDelta; i++) {
      final double b = i * delta;
      final double[] vs = line.values(b);
      Diagnostic.developerLog(b + " " + vs[0] + " " + vs[1]);
    }
  }

  static Vector variance(int[][] membersOf, BlockInfo info, Vector r, EigenvalueDecomposition ed) {
    final Jama.Matrix eig = ed.getV();
    final double[] realEigenvalues = ed.getRealEigenvalues();
    final boolean isGlobal = info.getGlobalSpeciesMap() == null;
    final int blockSize = realEigenvalues.length;
    final int totalGenomes =  isGlobal ? blockSize : info.getGlobalSpeciesMap().size();
    //System.err.println("Eigen matrix:");
    //System.err.println(IntegralAbstract.toString(eig.getArray()));
    //System.err.println("Eigen values:");
    //System.err.println(IntegralAbstract.toString(realEigenvalues));

    final Vector v = new VectorSimple(totalGenomes);
    for (int j = 0; j < totalGenomes; j++) { // This guy should iterate over all members of the taxonomy (i.e. including clades), using global ids
      double sum = 0.0;
      for (int i = 0; i < blockSize; i++) {
        //final double viold = eig.get(j, i);
        final double la = realEigenvalues[i];
        if (la < 0.0) { //protect against errors in eigenvector code
          sum = Double.POSITIVE_INFINITY;
          continue;
        }
        double vi = 0;
        for (int k = 0; k < membersOf[j].length; k++) { // This guy iterates over the species in the block, using local ids
          final int j2 = membersOf[j][k];
          vi += eig.get(j2, i) * r.get(j2);
        }

        final double vi2 = Math.pow(vi, 2);
        if (vi2 > 0.0) { //ignore cases where dot product is zero - because it is a square cant be negative
          final double t = vi2 / la;

          if (Double.isInfinite(t) || Double.isNaN(t)) {
            Diagnostic.developerLog("i=" + i + " j=" + j + " vi=" + vi + " vi^2=" + vi2 + " lamb_i=" + la + " t=" + t);
          }
          //If sum ever ends up as +infinity or NaN then the lower and higher bounds are set to 0.0 and 1.0
          //System.err.println("t=" + t + " vi2=" + vi2 + " la=" + la);
          sum += t;
        }
      }
      //System.err.println("sum=" + sum);
      v.set(j, sum);
    }
    return v;
  }

  /* Make a default membership matrix corresponding to a flat taxonomy */
  static int[][] makeFlatMembership(final int totalGenomes) {
    final int[][] membersOf = new int[totalGenomes][];
    for (int j = 0; j < totalGenomes; j++) {
      membersOf[j] = new int[] {j};
    }
    return membersOf;
  }

  private void message(final String msg, final Vector v) {
    //final String str = msg + (mBlockInfo.isVerbose() && v != null ? LS + v.toString() : "");
    final String str = msg + (v != null ? LS + v.toString() : "");
    Diagnostic.developerLog(str);
  }

  private static Vector polakRibiere(Vector lastJacobian, Vector jacobian, Vector delta) {
    final Vector diff = MatrixUtils.subtract(jacobian, lastJacobian);
    final double a = MatrixUtils.multiply(jacobian, diff);
    final double b = MatrixUtils.multiply(lastJacobian, lastJacobian);
    final double beta = a / b;
    if (beta < 0.0) {
      return MatrixUtils.negative(jacobian);
    }

    final Vector bd = MatrixUtils.multiply(beta, delta);
    return MatrixUtils.subtract(bd, jacobian);
  }

  /**
   * Compute Jacobian in log space.
   * @param r current position (in frequency space).
   * @param fixedIds Set of ids which have been forced to have a derivative of 0.0 and a value close to 0.0.
   * @return the Jacobian in log space.
   */
  private Pair<Vector, Double>  jacobian(Vector r, final int[] fixedIds) {
    final Pair<Vector, Double> pair = jacobianR(mBlockInfo, r, fixedIds);
    final Vector jr = pair.getA();
    final Vector jacobian = MatrixUtils.pointProduct(r, jr);
    if (BlockInfo.VERY_VERBOSE) {
      Diagnostic.developerLog("Jr: " + jr.toString());
      Diagnostic.developerLog("Js: " + jacobian.toString());
    }
    return new Pair<>(jacobian, pair.getB());
  }

  /**
   * Compute Jacobian in frequency space.
   * @param blockInfo information about the genomes etc.
   * @param r current position (in frequency space).
   * @param fixedIds Set of ids which have been forced to have a derivative of 0.0 and a value close to 0.0.
   *
   * @return the Jacobian in frequency space.
   */
  static Pair<Vector, Double> jacobianR(BlockInfo blockInfo, final Vector r, final int[] fixedIds) {
    final Vector jacobian = new VectorSimple(blockInfo.getN());
    double ll = 0.0;
    for (int i = 0; i < blockInfo.getN(); i++) {
      final long length = blockInfo.getGenomeLength(i);
      final double lr = length * r.get(i);
      ll += lr;
      jacobian.set(i, length);
    }
    for (final Frag frag : blockInfo.getFrags()) {
      ll += frag.increment(r, jacobian);
    }
    for (final int id : fixedIds) {
      jacobian.set(id, 0);
    }
    return new Pair<>(jacobian, ll);
  }

  /**
   * Compute L in frequency space.
   *
   * @param r current position (in frequency space).
   * @param blockInfo current block genomes and fragment counts
   * @return L.
   */
  static double ll(final Vector r, final BlockInfo blockInfo) {
    double ll = 0.0;
    for (int i = 0; i < blockInfo.getN(); i++) {
      final long length = blockInfo.getGenomeLength(i);
      final double lr = length * r.get(i);
      //System.err.println("lr=" + lr);
      ll += lr;
    }
    for (final Frag frag : blockInfo.getFrags()) {
      final double lf = frag.l(r);
      //System.err.println("lf=" + lf);
      ll += lf;
    }
    return ll;
  }

  /**
   * Compute Hessian in frequency space.
   * @param r current position (in frequency space).
   * @return the Jacobian in frequency space.
   */
  Matrix hessianR(final Vector r) {
    final Matrix hessian = new MatrixSymmetric(mBlockInfo.getN());
    for (final Frag frag : mBlockInfo.getFrags()) {
      frag.incrementR(r, hessian);
    }
    return hessian;
  }

  /**
   * Compute Jacobian and Hessian in frequency space.
   * @param r current position (in frequency space).
   * @return the Jacobian in frequency space.
   */
  Matrix hessian(final Vector r) {
    final Matrix hessian = new MatrixSymmetric(mBlockInfo.getN());
    final Vector jacobian = new VectorSimple(mBlockInfo.getN());
    for (int i = 0; i < mBlockInfo.getN(); i++) {
      final long length = mBlockInfo.getGenomeLength(i);
      final double lr = length * r.get(i);
      jacobian.set(i, lr);
      hessian.set(i, i, lr);
    }

    for (final Frag frag : mBlockInfo.getFrags()) {
      frag.increment(r, jacobian, hessian);
    }
    return hessian;
  }


  private Pair<Vector, Double> solveLine(final Vector r, final Vector delta) {
    final Line line = new SpeciesLine(r, delta, mBlockInfo);
    //final Line line = new SpeciesLineLinearDeriv(r, delta, mBlockInfo);
    final LineSolver solver = new LinearInterpolationSolver();
    //final LineSolver solver = new NewtonRaphsonSolver(false);
    final LineSolver ls = new LoggingLineSolver(solver, mBlockInfo.isVerbose());
    final double d = ls.solveLine(line, LoggingLineSolver.RELATIVE_THRESHOLD);
    //final double d = ls.solveLine(line, L_TERMINATION_TARGET);
    if (d == 0.0) {
      return null;
    }
    final Vector soln = incrS(mBlockInfo, r, delta, d);
    return new Pair<>(soln, d);
  }

  private Pair<Vector, Double> solveLineRobust(final Vector r, final Vector delta) {
    final Line line = new LLine(r, delta, mBlockInfo);
    final Minimizer ls = new Minimizer();
    final double d = ls.solveLine(line, LoggingLineSolver.RELATIVE_THRESHOLD);
    if (d == 0.0) {
      return null;
    }
    final Vector soln = incrS(mBlockInfo, r, delta, d);
    return new Pair<>(soln, d);
  }

  static Vector incrS(BlockInfo blockInfo, final Vector r, final Vector delta, final double d) {
    final Vector soln = new VectorSimple(blockInfo.getN());
    for (int i = 0; i < blockInfo.getN(); i++) {
      final double e = Math.exp(delta.get(i) * d);
      final double rv = r.get(i);
      soln.set(i, rv * e);
    }
    return soln;
  }


  @Override
  public void toString(final StringBuilder sb) {
    sb.append("genomes:").append(LS);
    for (int i = 0; i < mBlockInfo.getN(); i++) {
      sb.append("[").append(i).append("]").append(mBlockInfo.getTaxonId(i)).append(LS);
    }

    sb.append("fragments:").append(LS);
    for (int i = 0; i < mBlockInfo.getFrags().length; i++) {
      sb.append("[").append(i).append("]").append(mBlockInfo.getFrags()[i].toString()).append(LS);
    }

  }

  @Override
  public boolean integrity() {
    if (BlockInfo.VERY_VERBOSE) {
      Exam.assertTrue(mBlockInfo.isVerbose());
    }
    Exam.assertTrue(L_TERMINATION > 0.0);
    return true;
  }

  @Override
  public boolean globalIntegrity() {
    for (final Frag frag : mBlockInfo.getFrags()) {
      frag.integrity();
    }
    return true;
  }

}
