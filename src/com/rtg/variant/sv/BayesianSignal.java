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

package com.rtg.variant.sv;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.License;
import com.rtg.util.Pair;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * A factory for bayesian signals.
 *
 */
@TestClass(value = {"com.rtg.variant.sv.BayesianSignalTest", "com.rtg.variant.sv.DeleteBoundaryBayesianSignalTest"})
public abstract class BayesianSignal {

  /**
   * <code>DistrSam</code> bundles up a distribution along with the
   * window over which the distribution will be placed.
   */
  public static final class DistrSam extends Pair<Distribution, SamCounts> {
    DistrSam(Distribution a, SamCounts b) {
      super(a, b);
    }
  }

  static double offsetLeft(ReadGroupStats stats, double offset, boolean reverse) {
    if (reverse) {
      final double length = stats.meanLength();
      return offset - length - 1;
    } else {
      return -offset;
    }
  }

  static double offsetRight(ReadGroupStats stats, double offset, boolean reverse) {
    if (reverse) {
      final double length = stats.meanLength();
      return -offset - length - 1;
    } else {
      return offset;
    }
  }

  protected boolean mDebug;

  BayesianSignal(boolean debug) {
    mDebug = debug;
  }

  BayesianSignal() {
    mDebug = false;
  }

  abstract Distribution leftArmProper(ReadGroupStats stats, boolean reverse);

  abstract Distribution leftArmDiscordant(ReadGroupStats stats, boolean reverse);

  abstract Distribution leftArmUnmated(ReadGroupStats stats, boolean reverse);

  abstract Distribution rightArmProper(ReadGroupStats stats, boolean reverse);

  abstract Distribution rightArmDiscordant(ReadGroupStats stats, boolean reverse);

  abstract Distribution rightArmUnmated(ReadGroupStats stats, boolean reverse);

  /**
   * Create pairings of distributions with SamArray objects for the
   * supplied read group. This implementation uses the distributions
   * returned by the above methods.
   * @param state the read-group state
   * @param reverse if true, the signal should reverse the window inputs.
   * @param columnLabel a string prefix to use in debug output file naming.
   * @return pairings suitable for using as inputs to a bayesian signal.
   */
  protected List<DistrSam> makeInputs(ReadGroupState state, boolean reverse, String columnLabel) {
    final ReadGroupStats stats = state.stats();
    Diagnostic.userLog("Creating model for " + columnLabel + " - " + state.stats().id());
    final String dir = reverse ? "-rev" : "";
    if (mDebug && License.isDeveloper()) {
      DistributionUtils.plot(new File("debug-dist-" + columnLabel + "-leftProper" + dir + ".txt"), "leftProper" + dir, leftArmProper(stats, reverse));
      DistributionUtils.plot(new File("debug-dist-" + columnLabel + "-leftDiscordant" + dir + ".txt"), "leftDiscordant" + dir, leftArmDiscordant(stats, reverse));
      DistributionUtils.plot(new File("debug-dist-" + columnLabel + "-leftUnmated" + dir + ".txt"), "leftUnmated" + dir, leftArmUnmated(stats, reverse));
      DistributionUtils.plot(new File("debug-dist-" + columnLabel + "-rightProper" + dir + ".txt"), "rightProper" + dir, rightArmProper(stats, reverse));
      DistributionUtils.plot(new File("debug-dist-" + columnLabel + "-rightDiscordant" + dir + ".txt"), "rightDiscordant" + dir, rightArmDiscordant(stats, reverse));
      DistributionUtils.plot(new File("debug-dist-" + columnLabel + "-rightUnmated" + dir + ".txt"), "rightUnmated" + dir, rightArmUnmated(stats, reverse));
    }

    return Arrays.asList(// left arm
      new DistrSam(compactDistribution(leftArmProper(stats, reverse), "left-proper" + dir), state.properLeftArm(reverse)),
      new DistrSam(compactDistribution(leftArmDiscordant(stats, reverse), "left-discordant" + dir), state.discordantLeftArm(reverse)),
      new DistrSam(compactDistribution(leftArmUnmated(stats, reverse), "left-unmated" + dir), state.unmatedLeftArm(reverse)),

      // right arm
      new DistrSam(compactDistribution(rightArmProper(stats, reverse), "right-proper" + dir), state.properRightArm(reverse)),
      new DistrSam(compactDistribution(rightArmDiscordant(stats, reverse), "right-discordant" + dir), state.discordantRightArm(reverse)),
      new DistrSam(compactDistribution(rightArmUnmated(stats, reverse), "right-unmated" + dir), state.unmatedRightArm(reverse)));
  }

  static Distribution compactDistribution(Distribution a, String label) {
    if (a instanceof DistributionConstant) {
      return a;
    } else if (a instanceof DistributionStep) {
      return a;
    }
    final int lo = a.lo();
    final int hi = a.hi();
    final double precision = 0.0000000001;
    final double constantval = a.get(lo);
    boolean constant = true;
    int steppos = 0;
    double stepval = 0;
    boolean step = false;
    for (int i = lo + 1; i < hi; ++i) {
      final double current = a.get(i);
      if (constant && (Math.abs(current - constantval) > precision)) {
        constant = false;
        step = true;
        steppos = i - 1;
        stepval = current;
      } else if (step && (Math.abs(current - stepval) > precision)) {
        step = false;
        break;
      }
    }
    if (constant) {
      Diagnostic.developerLog("Upgraded " + label + " from "  + a.getClass().getSimpleName() + " to DistributionConstant");
      return new DistributionConstant(lo, hi, constantval);
    } else if (step) {
      Diagnostic.developerLog("Upgraded " + label + " from " + a.getClass().getSimpleName() + " to DistributionStep");
      return new DistributionStep(lo, hi, steppos, constantval, stepval);
    } else {
      return a;
    }
  }

  /**
   * Creates a bayesian signal across the input read group states
   *
   * @param groups a separate <code>ReadGroupState</code> object for each read group
   * @param reverse if true, the model will be reversed on the
   * template. Typically this will swap a left breakpoint signal to
   * become a right breakpoint signal.
   * @param columnLabel the column label for the output file.
   * @return a <code>Signal</code> value
   */
  public final Signal makeSignal(ReadGroupState[] groups, boolean reverse, String columnLabel) {
    final ArrayList<DistrSam> inputs = new ArrayList<>();
    for (final ReadGroupState state : groups) {
      inputs.addAll(makeInputs(state, reverse, columnLabel));
    }
    return makeSignal(inputs, columnLabel);
  }

  private Signal makeSignal(List<DistrSam> pairs, String columnLabel) {
    final Signal[] signals = new Signal[pairs.size()];
    int i = 0;
    for (final DistrSam ds : pairs) {
      final Distribution di = ds.getA();
      signals[i] = di.getSignalLn(ds.getB(), columnLabel);
      ++i;
    }
    if (mDebug) {
      try {
        return new SignalSum(new PrintStream(new FileOutputStream(new File("debug-" + columnLabel + ".txt"))), columnLabel, signals);
        //return new SignalSum(System.err, columnLabel, signals);
      } catch (final IOException e) {
        throw new RuntimeException("eek", e);
      }
    } else {
      return new SignalSum(columnLabel, signals);
    }
  }
}
