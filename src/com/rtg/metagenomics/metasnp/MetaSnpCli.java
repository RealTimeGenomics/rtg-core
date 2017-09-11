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
package com.rtg.metagenomics.metasnp;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.LoggedCli;
import com.rtg.util.PosteriorUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.Flag;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LogStream;
import com.rtg.variant.util.arithmetic.LogPossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;
import com.rtg.vcf.VcfWriter;
import com.rtg.vcf.VcfWriterFactory;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfHeader;
import com.rtg.vcf.header.VcfNumber;

/**
 * Main entry point for <code>metasnp</code>.
 */
public class MetaSnpCli extends LoggedCli {

  /** Flag name for low frequency filter */
  static final String MIN_FREQ_FLAG = "min-freq";
  /** Flag name for low coverage filter */
  static final String MIN_TOTAL_COVERAGE = "min-total-coverage";
  /** Flag name for high frequency filter */
  static final String MAX_FREQ_FLAG = "max-freq";
  /** Flag name for high coverage filter */
  static final String MAX_TOTAL_COVERAGE = "max-total-coverage";
  private static final String BETA = "beta";
  private static final String MODULE_NAME = "metasnp";
  private static final String STRAINS = "strains";
  private static final String VISUALISATION = "visualisation";
  private static final String VISUALISATION_PREFIX = "visual";
  private static final String ITERATIONS = "iterations";
  private static final String ERROR_RATE = "error-rate";
  private static final String XI_PRIORS = "Xxi";

  private static final String VCF_OUTPUT = "strains.vcf";
  private static final String XI_FILE = "xi.txt";
  private static final String LIKE = "LIKE";
  private static final String SYNDROME = "SYNDROME";

  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  @Override
  public String description() {
    return "estimate fractions and variants of species in multiple samples";
  }

  @Override
  protected void initFlags() {
    CommonFlags.initOutputDirFlag(mFlags);
    mFlags.registerRequired(File.class, CommonFlags.FILE, "allele stats file");
    mFlags.registerRequired('s', STRAINS, Integer.class, CommonFlags.INT, "number of strains");
    mFlags.registerOptional('e', ERROR_RATE, Double.class, CommonFlags.FLOAT, "read/mapping error rate", 0.01);
    mFlags.registerRequired('m', MIN_FREQ_FLAG, Integer.class, CommonFlags.INT, "minimum allele frequency");
    mFlags.registerOptional('M', MAX_FREQ_FLAG, Integer.class, CommonFlags.INT, "maximum allele frequency");
    mFlags.registerOptional(MAX_TOTAL_COVERAGE, Integer.class, CommonFlags.INT, "maximum coverage threshold");
    mFlags.registerOptional(MIN_TOTAL_COVERAGE, Integer.class, CommonFlags.INT, "minimum coverage threshold");
    mFlags.registerOptional('v', VISUALISATION, "produce visualisation files");
    mFlags.registerOptional('i', ITERATIONS, Integer.class, CommonFlags.INT, "number of iterations to attempt convergence", 10);
    mFlags.registerOptional(XI_PRIORS, String.class, "FLOAT...", "initial values for the xi matrix");
    final Flag<String> betaType = mFlags.registerOptional(BETA, String.class, CommonFlags.STRING, "hypothesis probability method", "reestimate");
    betaType.setParameterRange(betaFlagValues());
    mFlags.setValidator(new Validator());
  }

  static String[] betaFlagValues() {
    final EmIterate.BetaType[] values = EmIterate.BetaType.values();
    final String[] str = new String[values.length];
    for (int i = 0; i < values.length; ++i) {
      str[i] = values[i].toString().toLowerCase(Locale.getDefault());
    }
    return str;
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(CommonFlags.OUTPUT_FLAG);
  }

  private MetaSnpReader getReader(final File file) throws IOException {
    return VcfUtils.isVcfExtension(file) ? new VcfMetaSnpReader(file) : new AlleleStatReader(file);
  }

  private static double[][] initXi(final int samples, final int strains, PossibilityArithmetic arith, final String xiPriors) {
    final String[] xis = StringUtils.split(xiPriors, ',');
    if (xis.length != samples * strains) {
      throw new NoTalkbackSlimException("Initial values for xi must have length " + samples + "*" + strains + " = " + samples * strains);
    }
    final double[][] xi = new double[samples][strains];
    try {
      for (int k = 0, i = 0; k < samples; ++k) {
        for (int j = 0; j < strains; ++j, ++i) {
          final double v = Double.parseDouble(xis[i]);
          if (v <= 0 || Double.isInfinite(v)) {
            throw new NoTalkbackSlimException("Invalid xi value: " + v);
          }
          xi[k][j] = arith.prob2Poss(v);
        }
      }
    } catch (final NumberFormatException e) {
      throw new NoTalkbackSlimException("Invalid xi values: " + xiPriors);
    }
    return xi;
  }

  static double[][] initXi(int samples, int strains, PossibilityArithmetic arith) {
    final double[][] xi = new double[samples][strains];
    for (int i = 0; i < xi.length; ++i) {
      xi[i] = RandomWalkXiFinder.uniformDistribution(strains, arith);
    }
    return xi;
  }

  @Override
  protected int mainExec(OutputStream out, LogStream err) throws IOException {

    final File outputDirectory = outputDirectory();
    final File f = (File) mFlags.getAnonymousValue(0);
    final int strains = (Integer) mFlags.getValue(STRAINS);
    final int minFreq = (Integer) mFlags.getValue(MIN_FREQ_FLAG);
    final int maxFreq = mFlags.isSet(MAX_FREQ_FLAG) ? (Integer) mFlags.getValue(MAX_FREQ_FLAG) : Integer.MAX_VALUE;

    final int minCov = mFlags.isSet(MIN_TOTAL_COVERAGE) ? (Integer) mFlags.getValue(MIN_TOTAL_COVERAGE) : 0;
    final int maxCov = mFlags.isSet(MAX_TOTAL_COVERAGE) ? (Integer) mFlags.getValue(MAX_TOTAL_COVERAGE) : Integer.MAX_VALUE;
    final EmIterate.BetaType updateBeta = EmIterate.BetaType.valueOf(((String) mFlags.getValue(BETA)).toUpperCase(Locale.getDefault()));

    final double error = (Double) mFlags.getValue(ERROR_RATE);
    int approxLength = 0; // number of lines of input approximates length of genome
    try (final MetaSnpReader reader = getReader(f)) {
      final List<Integer> ref = new ArrayList<>();
      final List<double[][]> evidence = new ArrayList<>();
      final List<MetaSnpLine> lines = new ArrayList<>();
      MetaSnpLine line;
      while ((line = reader.nextLine()) != null) {
        ++approxLength;
        final int refAllele = line.getReferenceIndex();
        int nonRefCount = 0;
        int total = 0;
        final double[][] evidenceArray = new double[line.mCounts[0].length][line.mCounts.length];
        for (int i = 0; i < evidenceArray.length; ++i) {
          for (int j = 0; j < evidenceArray[i].length; ++j) {
            evidenceArray[i][j] = line.mCounts[j][i];
            if (j != refAllele) {
              nonRefCount += line.mCounts[j][i];
            }
            total += line.mCounts[j][i];
          }
        }
        if (refAllele >= 0 && nonRefCount >= minFreq && nonRefCount < maxFreq && total >= minCov && total < maxCov) {
          ref.add(refAllele);
          evidence.add(evidenceArray);
          lines.add(line);
        }
      }
      final PossibilityArithmetic arith = LogPossibility.SINGLETON;
      final int samples = reader.samples().size();
      final double[][] xiPriors = mFlags.isSet(XI_PRIORS) ? initXi(samples, strains, arith, (String) mFlags.getValue(XI_PRIORS)) : initXi(samples, strains, arith);
      Diagnostic.info(evidence.size() + "/" + approxLength + " positions passed initial thresholding");
      final List<EmIterate.EmResult> iterations = EmIterate.iterate(ref, evidence, strains, approxLength, arith, new EmIterate.FixedIterations((Integer) mFlags.getValue(ITERATIONS)), updateBeta, error, xiPriors);
      final EmIterate.EmResult result = iterations.get(iterations.size() - 1);
      final double[][] xi = result.mXi;
      try (PrintStream xiOut = new PrintStream(FileUtils.createOutputStream(new File(outputDirectory, XI_FILE)))) {
        writeXi(xi, arith, xiOut, reader.samples());
      }
      try (final ByteArrayOutputStream xiBytes = new ByteArrayOutputStream()) {
        try (PrintStream xiString = new PrintStream(xiBytes)) {
          writeXi(xi, arith, xiString, reader.samples());
        }
        Diagnostic.info("Estimated strain proportions: " + StringUtils.LS + xiBytes);
      }
      writeVcf(ref, lines, iterations.get(iterations.size() - 1), new File(outputDirectory, VCF_OUTPUT), arith);
      if (mFlags.isSet(VISUALISATION)) {
        final File visual = new File(outputDirectory, VISUALISATION_PREFIX);
        try (FileOutputStream visualStream = new FileOutputStream(visual)) {
          outputVisualisation(ref, lines, evidence, result, visualStream);
        }
        for (int i = 0; i < iterations.size(); ++i) {
          final File visualIt = new File(visual.getPath() + i);
          try (OutputStream visualStream = FileUtils.createOutputStream(visualIt)) {
            outputVisualisation(ref, lines, evidence, iterations.get(i), visualStream);
          }
        }
      }
    }
    return 0;
  }

  static void writeXi(double[][] xi, PossibilityArithmetic arith, PrintStream xiOut, List<String> samples) {
    if (samples == null) {
      for (int j = 0; j < xi.length; ++j) {
        xiOut.print("\tSample" + j);
      }
    } else {
      for (String sample : samples) {
        xiOut.print("\t" + sample);
      }
    }
    xiOut.println();
    for (int i = 0; i < xi[0].length; ++i) {
      xiOut.print("Strain" + i);
      for (final double[] sample : xi) {
        final double value = sample[i];
        xiOut.print("\t" + Utils.realFormat(arith.poss2Prob(value), 3));
      }
      xiOut.println();
    }
  }

  static void outputVisualisation(List<Integer> refBytes, List<MetaSnpLine> lines, List<double[][]> evidence, EmIterate.EmResult result, OutputStream out) throws IOException {
    for (int i = 0; i < refBytes.size(); ++i) {
      final double[][] currentEvidence = evidence.get(i);
      final int[] currentAssignments = result.mAssignments.get(i).mCalls;
      final int[] totals = new int[currentEvidence.length];
      final int numAlleles = currentEvidence[0].length;
      for (int sample = 0; sample < currentEvidence.length; ++sample) {
        for (byte allele = 0; allele < numAlleles; ++allele) {
          totals[sample] += currentEvidence[sample][allele];
        }
      }
      for (byte allele = 0; allele < numAlleles; ++allele) {
        if (allele == refBytes.get(i)) {
          continue;
        }
        boolean unassigned = true;
        for (int currentAssignment : currentAssignments) {
          if (currentAssignment == allele) {
            unassigned = false;
            break;
          }
        }
        if (unassigned) {
          continue;
        }
        int color = 0;
        for (int currentAssignment : currentAssignments) {
          color <<= 1;
          if (currentAssignment == allele) {
            ++color;
          }
        }
        final double[] coordinates = new double[currentEvidence.length];
        for (int sample = 0; sample < currentEvidence.length; ++sample) {
          coordinates[sample] = currentEvidence[sample][allele] / (double) totals[sample];
        }
        final StringBuilder output = new StringBuilder();
        output.append(lines.get(i).getSequence()).append('\t');
        output.append(lines.get(i).getPosition() + 1).append('\t');
        output.append(lines.get(i).getReferenceAllele()).append('\t');
        output.append(lines.get(i).mAlleles[allele]).append('\t');
        output.append(color);
        for (double coordinate : coordinates) {
          output.append("\t").append(Utils.realFormat(coordinate));
        }
        output.append(StringUtils.LS);
        out.write(output.toString().getBytes());
      }
    }
  }

  static void writeVcf(List<Integer> refBytes, List<MetaSnpLine> lines, EmIterate.EmResult res, File out, PossibilityArithmetic arith) throws IOException {
    final VcfHeader header = new VcfHeader();
    header.addCommonHeader();
    for (int i = 0; i < res.mAssignments.get(0).mCalls.length; ++i) {
      header.addSampleName("Strain" + i);
    }
    header.addInfoField(LIKE, MetaType.FLOAT, VcfNumber.DOT, "phred scaled likelihood of genotype assignments");
    header.addInfoField(SYNDROME, MetaType.STRING, VcfNumber.DOT, "packed representation of strain assignment");
    final List<Integer> alts = new ArrayList<>();
    try (final VcfWriter writer = new VcfWriterFactory().zip(false).make(header, out)) {
      for (int i = 0; i < lines.size(); ++i) {
        final MetaSnpLine line = lines.get(i);
        final int[] assignments = res.mAssignments.get(i).mCalls;
        final int ref = refBytes.get(i);
        alts.clear();
        alts.add(ref);
        for (int assignment1 : assignments) {
          if (!alts.contains(assignment1)) {
            alts.add(assignment1);
          }
        }
        final VcfRecord record = new VcfRecord(line.getSequence(), line.getPosition(), line.getReferenceAllele());
        final double phred = PosteriorUtils.phredIfy(arith.poss2Ln(res.mAssignments.get(i).mLikelihood));
        record.setInfo("LIKE", Utils.realFormat(phred, 3));
        record.setNumberOfSamples(assignments.length);
        for (int alt = 1; alt < alts.size(); ++alt) {
          record.addAltCall(line.mAlleles[alts.get(alt)]);
        }
        final StringBuilder syndrome = new StringBuilder();
        for (int assignment : assignments) {
          record.addFormatAndSample(VcfUtils.FORMAT_GENOTYPE, String.valueOf(alts.indexOf(assignment)));
          syndrome.append(alts.indexOf(assignment) == 0 ? '0' : '1');
        }
        record.setInfo("SYNDROME", syndrome.toString());
        writer.write(record);
      }
    }
  }

  private static class Validator implements com.rtg.util.cli.Validator {
    @Override
    public boolean isValid(CFlags flags) {
      if (!CommonFlags.validateOutputDirectory(flags)) {
        return false;
      }
      final int strains = (Integer) flags.getValue(STRAINS);
      if (strains < 2) {
        flags.setParseMessage("It makes no sense to run this with less than 2 strains");
        return false;
      }
      final File f = (File) flags.getAnonymousValue(0);
      if (!f.exists()) {
        flags.setParseMessage("The file '" + f + "' doesn't exist");
        return false;
      }
      final double error = (Double) flags.getValue(ERROR_RATE);
      if (error <= 0 || error >= 1.0) {
        flags.setParseMessage("--" + ERROR_RATE + " should be between 0.0 and 1.0 (exclusive)");
      }
      return true;
    }
  }

}
