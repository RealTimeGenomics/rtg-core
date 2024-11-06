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
package com.rtg.zooma;

import static com.rtg.launcher.CommonFlags.FILE;
import static com.rtg.launcher.CommonFlags.INT;
import static com.rtg.launcher.CommonFlags.STRING;

import java.io.File;
import java.io.OutputStream;
import java.io.PrintStream;

import com.rtg.launcher.AbstractCli;
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
    mFlags.registerRequired('i', REF_FLAG, File.class, FILE, "FASTA file for reference");
    mFlags.registerOptional('o', OUTPUT_FLAG, File.class, FILE, "name of output index file", new File("zooma.index.bin"));
    mFlags.registerOptional('c', INCLUDE_FLAG, String.class, STRING, "include chromosomes with this string in the name");
    mFlags.registerOptional('e', EXCLUDE_FLAG, String.class, STRING, "exclude chromosomes with this string in the name");
    mFlags.registerOptional('w', WORD, Integer.class, INT, "hash width (<= 21)", 18);
    mFlags.registerOptional('s', STEP, Integer.class, INT, "step size", 1);
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) {
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
