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
package com.rtg.simulation.variants;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.LoggedCli;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LogStream;

/**
 * Generates a genome SDF corresponding to the sample genotype described in a VCF file.
 */
public class SampleReplayerCli extends LoggedCli {

  private static final String MODULE_NAME = "samplereplay";

  private static final String SAMPLE_VCF = "input";
  private static final String SAMPLE_NAME = "sample";
  private static final String OUTPUT_SDF = "output";
  private static final String REFERENCE_SDF = "reference";


  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  @Override
  protected void initFlags() {
    mFlags.setDescription("Generates the genome corresponding to a sample genotype.");
    mFlags.registerExtendedHelp();
    CommonFlagCategories.setCategories(mFlags);
    mFlags.registerRequired('t', REFERENCE_SDF, File.class, "SDF", "SDF containing reference genome").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerRequired('o', OUTPUT_SDF, File.class, "SDF", "name for output SDF").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerRequired('i', SAMPLE_VCF, File.class, "FILE", "input VCF containing the sample genotype").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerRequired('s', SAMPLE_NAME, String.class, "STRING", "name of the sample to select from the VCF").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.setValidator(new Validator() {
      @Override
      public boolean isValid(CFlags flags) {
        return CommonFlags.validateOutputDirectory(outputDirectory());
      }
    });

  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(OUTPUT_SDF);
  }

  @Override
  protected int mainExec(OutputStream out, LogStream log) throws IOException {
    final CFlags flags = mFlags;
    final File input = (File) flags.getValue(REFERENCE_SDF);
    final File sampleVcf = FileUtils.getZippedFileName(true, (File) flags.getValue(SAMPLE_VCF));
    final String sample = (String) flags.getValue(SAMPLE_NAME);
    try (SequencesReader dsr = SequencesReaderFactory.createMemorySequencesReader(input, true, LongRange.NONE)) {
      final SampleReplayer vr = new SampleReplayer(dsr);
      vr.replaySample(sampleVcf, outputDirectory(), sample);
      return 0;
    }
  }

  /**
   * @param args arguments
   */
  public static void main(String[] args) {
    new SampleReplayerCli().mainExit(args);
  }
}
