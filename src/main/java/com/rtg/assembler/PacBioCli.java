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

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Collection;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.ParamsCli;
import com.rtg.util.IORunnable;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;

/**
 */
public class PacBioCli extends ParamsCli<PacBioParams> {

  static final String MODULE_NAME = "addpacbio";
  static final String TRIM = "trim";

  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  @Override
  public String description() {
    return "add Pacific Biosciences reads to an assembly";
  }

  @Override
  protected IORunnable task(PacBioParams params, OutputStream out) {
    return new PacBio(params, out);
  }

  @Override
  protected PacBioParams makeParams() throws IOException {
    return makeParamsLocal(mFlags);
  }

  protected static PacBioParams makeParamsLocal(CFlags flags) throws IOException {
    final PacBioParams.Builder builder = PacBioParams.builder();
    builder.directory((File) flags.getValue(CommonFlags.OUTPUT_FLAG))
        .reads(CommonFlags.getFileList(flags, CommonFlags.INPUT_LIST_FLAG, null, true))
        .graph((File) flags.getValue(GraphMapCli.GRAPH_FLAG))
        .trimGraph(flags.isSet(TRIM));

    return builder.create();
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(CommonFlags.OUTPUT_FLAG);
  }

  @Override
  protected void initFlags() {
    initLocalFlags(mFlags);
  }
  protected static void initLocalFlags(CFlags flags) {
    CommonFlagCategories.setCategories(flags);
    CommonFlags.initOutputDirFlag(flags);
    flags.registerRequired('g', GraphMapCli.GRAPH_FLAG, File.class, CommonFlags.DIR, "graph of the assembly to map against").setCategory(INPUT_OUTPUT);
    flags.registerOptional(TRIM, "before mapping remove any short disconnected sequences from the graph").setCategory(INPUT_OUTPUT);
    final Flag<File> inFlag = flags.registerRequired(File.class, CommonFlags.SDF, "SDF directories containing reads to map");
    inFlag.setCategory(INPUT_OUTPUT);
    inFlag.setMinCount(0);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    final Flag<File> listFlag = flags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, CommonFlags.FILE, "file containing a list of SDF directories (1 per line) containing sequences to assemble").setCategory(INPUT_OUTPUT);
    flags.addRequiredSet(inFlag);
    flags.addRequiredSet(listFlag);
    flags.setValidator(new PacBioValidator());
  }

  private static class PacBioValidator implements Validator {
    /**
     * Check the file list and anonymous file input flags.
     * @param flags the flags to check
     * @return <code>true</code> if all okay <code>false</code> otherwise
     */
    public static boolean checkSdfFileList(CFlags flags) {
      final Collection<File> files;
      try {
        files = CommonFlags.getFileList(flags, CommonFlags.INPUT_LIST_FLAG, null, true);
      } catch (final IOException e) {
        flags.setParseMessage("An error occurred reading " + flags.getValue(CommonFlags.INPUT_LIST_FLAG));
        return false;
      }
      if (files.isEmpty()) {
        flags.setParseMessage("No input files specified.");
        return false;
      }
      return true;
    }

    @Override
    public boolean isValid(CFlags flags) {
      if (!CommonFlags.validateOutputDirectory(flags)) {
        return false;
      }
      if (!checkSdfFileList(flags)) {
        return false;
      }
      return CommonFlags.validateThreads(flags);
    }

  }

  /**
   * Noddy main.
   * @param args see usage
   */
  public static void main(final String[] args) {
    new PacBioCli().mainInit(args, System.out, System.err);
  }
}
