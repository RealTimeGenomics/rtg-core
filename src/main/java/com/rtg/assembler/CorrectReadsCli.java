/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
