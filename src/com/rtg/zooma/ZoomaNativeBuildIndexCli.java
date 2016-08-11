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
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

/**
 * Provide integration with C Zooma index building
 */
public class ZoomaNativeBuildIndexCli extends AbstractCli {

  private static final String REF_FLAG = "reference";
  private static final String OUTPUT_FLAG = "output";
  private static final String INCLUDE_FLAG = "include";
  private static final String EXCLUDE_FLAG = "exclude";
  private static final String WORD = "word";
  private static final String STEP = "step";

  @Override
  public String moduleName() {
    return "zbuild";
  }

  @Override
  public String description() {
    return "build an index for use with zmap";
  }

  @Override
  protected void initFlags() {
    mFlags.registerRequired('i', REF_FLAG, File.class, "FILE", "FASTA file for reference");
    mFlags.registerOptional('o', OUTPUT_FLAG, File.class, "FILE", "name of output index file", new File("zooma.index.bin"));
    mFlags.registerOptional('c', INCLUDE_FLAG, String.class, "STR", "include chromosomes with this string in the name");
    mFlags.registerOptional('e', EXCLUDE_FLAG, String.class, "STR", "exclude chromosomes with this string in the name");
    mFlags.registerOptional('w', WORD, Integer.class, CommonFlags.INT, "hash width (<= 21)", 18);
    mFlags.registerOptional('s', STEP, Integer.class, CommonFlags.INT, "step size", 1);
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    //Diagnostic.setLogStream(err);
    final File fastaFile = (File) mFlags.getValue(REF_FLAG);
    if (!fastaFile.exists()) {
      throw new NoTalkbackSlimException("No such file: " + fastaFile);
    }
    final File indexFile = (File) mFlags.getValue(OUTPUT_FLAG);
    if (indexFile.exists()) {
      throw new NoTalkbackSlimException("Output file " + indexFile + " already exists");
    }
    final String include = mFlags.isSet(INCLUDE_FLAG) ? (String) mFlags.getValue(INCLUDE_FLAG) : null;
    final String exclude = mFlags.isSet(EXCLUDE_FLAG) ? (String) mFlags.getValue(EXCLUDE_FLAG) : null;

    if (!NativeZooma.isEnabled()) {
      err.println("Native library libZooma could not be loaded.");
      return 1;
    }
    final NativeZooma zooma = new NativeZooma();
    return zooma.buildIndex(indexFile.getPath(), fastaFile.getPath(), include, exclude, (Integer) mFlags.getValue(WORD), (Integer) mFlags.getValue(STEP));
  }

}
