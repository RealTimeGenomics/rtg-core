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
package com.rtg.zooma;

import java.io.File;
import java.io.FileDescriptor;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.ngs.MapFlags;
import com.rtg.util.Environment;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

/**
 * Provide integration with C Zooma read mapping
 */
public class ZoomaNativeMapReadsCli extends AbstractCli {

  private static final String INDEX_FLAG = "index";
  private static final String LEFT_FLAG = "left";
  private static final String RIGHT_FLAG = "right";
  private static final String OUT_PREFIX_FLAG = "out";
  private static final String NO_QUICK_FLAG = "Xno-quick-cache";
  private static final String OUTPUT_CHUNK_FLAG = "output-chunk";
  private static final String INPUT_CHUNK_FLAG = "input-chunk";
  private static final String SAM_RG_FLAG = "sam-rg";
  private static final String E_SCORE = "e-score";
  private static final String STEP = "step";

  @Override
  public String moduleName() {
    return "zmap";
  }

  @Override
  public String description() {
    return "even faster read mapping";
  }

  @Override
  protected void initFlags() {
    mFlags.registerOptional('t', INDEX_FLAG, File.class, "FILE", "reference index file", new File("zooma.index.bin"));
    mFlags.registerRequired('l', LEFT_FLAG, File.class, "FILE", "left input file for FASTQ paired end data");
    mFlags.registerRequired('r', RIGHT_FLAG, File.class, "FILE", "right input file for FASTQ paired end data");
    mFlags.registerOptional('o', OUT_PREFIX_FLAG, String.class, "STRING", "prefix for output files", "outfile_zooma");
    mFlags.registerOptional('Q', NO_QUICK_FLAG, "do not build and use the \"quick cache\"");
    mFlags.registerOptional('k', OUTPUT_CHUNK_FLAG, Integer.class, "INT", "number of reads per chunk", 100000);
    mFlags.registerOptional('e', E_SCORE, Integer.class, CommonFlags.INT, "alignment score threshold per arm (2*e for mates)", 10);
    mFlags.registerOptional('s', STEP, Integer.class, CommonFlags.INT, "set the search step size (The w value from inside the index)");
    mFlags.registerOptional(SAM_RG_FLAG, String.class, "STRING", "string in the form \"@RG\\tID:RG1\\tSM:SAMPLE_NAME\\tPL:ILLUMINA\"");
    MapFlags.initFragmentSizeFlags(mFlags);
    CommonFlags.initThreadsFlag(mFlags);
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    //Diagnostic.setLogStream(err);
    final File indexFile = (File) mFlags.getValue(INDEX_FLAG);
    if (!indexFile.exists()) {
      throw new NoTalkbackSlimException("No such file: " + indexFile);
    }
    final File leftFile = (File) mFlags.getValue(LEFT_FLAG);
    if (!leftFile.exists()) {
      throw new NoTalkbackSlimException("No such file: " + leftFile);
    }
    final File rightFile = (File) mFlags.getValue(RIGHT_FLAG);
    if (!rightFile.exists()) {
      throw new NoTalkbackSlimException("No such file: " + rightFile);
    }

    if (!NativeZooma.isEnabled()) {
      err.println("Native library libZooma could not be loaded.");
      return 1;
    }
    final NativeZooma zooma = new NativeZooma();
    final Integer chunkSize = (Integer) mFlags.getValue(OUTPUT_CHUNK_FLAG);
    final Integer threads = CommonFlags.parseThreads((Integer) mFlags.getValue(CommonFlags.THREADS_FLAG));
    final int step = mFlags.isSet(STEP) ? (Integer) mFlags.getValue(STEP) : Integer.MAX_VALUE;
    double tsm = Environment.getTotalMemory() / 1024.0 / 1024.0 / 1024.0 - 128;
      // subtract chunksize * threads * average read size
    final double memforthreads = (chunkSize * 150L * 6 * threads) / 1024.0 / 1024.0 / 1024.0;
    tsm -= memforthreads;
    if (tsm < 100) {
      Diagnostic.warning("Not enough ram for default settings");
    }
    int icacheGB = (int) (tsm / 3);
    if (icacheGB > 20) {
      icacheGB = 20;
    }
    final int ocacheGB = (int) (tsm - icacheGB); // subtract input cache, the rest for output cache
    return zooma.mapReads(
      indexFile.getPath(), leftFile.getPath(), rightFile.getPath(), (String) mFlags.getValue(OUT_PREFIX_FLAG), threads,
      chunkSize, (Integer) mFlags.getValue(E_SCORE), (Integer) mFlags.getValue(CommonFlags.MIN_FRAGMENT_SIZE), (Integer) mFlags.getValue(CommonFlags.MAX_FRAGMENT_SIZE),
      step, (String) mFlags.getValue(SAM_RG_FLAG), true, true, false,
      false, false, 0, null, FileDescriptor.err, icacheGB, ocacheGB);
  }

}
