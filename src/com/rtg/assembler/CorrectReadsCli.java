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

package com.rtg.assembler;

import static com.rtg.launcher.CommonFlags.OUTPUT_FLAG;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.LoggedCli;
import com.rtg.util.cli.CFlags;
import com.rtg.util.io.LogStream;

/**
 */
public class CorrectReadsCli extends LoggedCli {
  private static final String INPUT = "input";
  private static final String THRESHOLD = "threshold";

  @Override
  public String moduleName() {
    return "correctreads";

  }

  @Override
  public String description() {
    return "correct read errors using kmer analysis";
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(OUTPUT_FLAG);
  }

  @Override
  protected int mainExec(OutputStream out, LogStream log) throws IOException {
    final File input = (File) mFlags.getValue(INPUT);
    final File output = (File) mFlags.getValue(OUTPUT_FLAG);
    final int kmerSize = (Integer) mFlags.getValue(DeBruijnAssemblerCli.KMER_SIZE);
    final int threshold = mFlags.isSet(THRESHOLD) ? (Integer) mFlags.getValue(THRESHOLD) : -1;
    CorrectReads.correct(input, output, kmerSize, threshold);
    return 0;
  }

  @Override
  protected void initFlags() {
    initFlagsLocal(mFlags);
  }
  static void initFlagsLocal(CFlags flags) {
    flags.setDescription("Attempt to correct read errors by locating low frequency kmers and mutating them into high frequency ones.");
    CommonFlags.initOutputDirFlag(flags);
    flags.registerRequired('i', INPUT, File.class, CommonFlags.SDF, "read SDF to correct");
    flags.registerRequired('k', DeBruijnAssemblerCli.KMER_SIZE, Integer.class, CommonFlags.INT, "size of kmer to use in correction");
    flags.registerOptional('c', THRESHOLD, Integer.class, CommonFlags.INT, "override the calculated frequency threshold");
  }

  /**
   * main
   * @param args command line arguments
   */
  public static void main(String[] args) {
      new CorrectReadsCli().mainExit(args);
  }

}
