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
package com.rtg.metagenomics;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.OutputParams;
import com.rtg.launcher.ParamsCli;
import com.rtg.metagenomics.MetaPipelineParams.MetaPipelineParamsBuilder;
import com.rtg.mode.SequenceType;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.Environment;
import com.rtg.util.IORunnable;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.util.io.LogStream;

/**
 * Prototype species pipeline wrapper
 * Filtering + Species + Protein Search Pipeline
 */
public class MetagenomicsWrapperCli extends ParamsCli<MetaPipelineParams> {

  private static final String ENVIRONMENT_REFERENCES_DIR = "references.dir";
  private static final String FILTER_REFERENCE_DEFAULT = "filter";
  private static final String SPECIES_REFERENCE_DEFAULT = "species";
  private static final String PROTEIN_REFERENCE_DEFAULT = "protein";

  static final String PROTEIN = "protein";
  static final String SPECIES = "species";
  static final String INPUT = "input";
  static final String INPUT_LEFT = "input-left";
  static final String INPUT_RIGHT = "input-right";
  static final String PLATFORM = "platform";
  static final String FILTER = "filter";

  LogStream mLogStream;
  PrintStream mErr;
  MetagenomicsWrapperTask mTask = null;

  protected enum Platform {
    ILLUMINA,
    IONTORRENT
  }

  @Override
  public String moduleName() {
    return "composition-functional-meta-pipeline";
  }

  @Override
  public String description() {
    return "run metagenomic composition and functional pipelines";
  }

  protected static class MetaWrapperValidator implements Validator {
    @Override
    public boolean isValid(CFlags flags) {
      int inputCount = 0;
      boolean valid = true;
      if (flags.isSet(INPUT_LEFT) || flags.isSet(INPUT_RIGHT)) {
        ++inputCount;
        if (!flags.isSet(INPUT_LEFT) || !flags.isSet(INPUT_RIGHT)) {
          flags.setParseMessage("You must provide both --" + INPUT_LEFT + " and --" + INPUT_RIGHT + " or use --" + INPUT);
          valid = false;
        } else {
          final File left = (File) flags.getValue(INPUT_LEFT);
          final File right = (File) flags.getValue(INPUT_RIGHT);
          if (!left.isFile()) {
            flags.setParseMessage("--" + INPUT_LEFT + " should be a file");
            valid = false;
          }
          if (!right.isFile()) {
            flags.setParseMessage("--" + INPUT_RIGHT + " should be a file");
            valid = false;
          }
        }
      }
      if (flags.isSet(INPUT)) {
        ++inputCount;
        final File illumina = (File) flags.getValue(INPUT);
        if (!illumina.exists() || (!illumina.isDirectory() && !illumina.isFile())) {
          flags.setParseMessage("--" + INPUT + " file doesn't exist: " + illumina.getPath());
          valid = false;
        }
      }
      if (inputCount < 1) {
        flags.setParseMessage("Some read data is required (use --" + INPUT + "/--" + INPUT_LEFT + "/--" + INPUT_RIGHT + ")");
        valid = false;
      }
      if (inputCount > 1) {
        flags.setParseMessage("Too many read datasets supplied. Provide one of formatted SDF, a single end fastq or left and right paired fastq files");
        valid = false;
      }
      final File filterSdf = (File) flags.getValue(FILTER);
      if (filterSdf != null) {
        if (!CommonFlags.validateSDF(filterSdf)) {
          valid = false;
        }
      } else {
        flags.setParseMessage("Filter reference required (use --" + FILTER + ")");
        valid = false;
      }
      if (flags.getFlag(SPECIES) != null) {
        final File speciesSdf = (File) flags.getValue(SPECIES);
        if (speciesSdf != null) {
          if (!CommonFlags.validateSDF(speciesSdf)) {
            valid = false;
          }
        } else {
          flags.setParseMessage("Species reference required (use --" + SPECIES + ")");
          valid = false;
        }
      }
      if (flags.getFlag(PROTEIN) != null) {
        final File proteinSdf = (File) flags.getValue(PROTEIN);
        if (proteinSdf != null) {
          if (!CommonFlags.validateSDF(proteinSdf)) {
            valid = false;
          } else {
            try {
              try (SequencesReader r = SequencesReaderFactory.createDefaultSequencesReader(proteinSdf)) {
                if (r.type() != SequenceType.PROTEIN) {
                  flags.setParseMessage("Protein database should be an SDF containing formatted protein data");
                  valid = false;
                }
              }
            } catch (IOException e) {
              flags.setParseMessage("Invalid protein database.");
              valid = false;
            }
          }
        } else {
          flags.setParseMessage("Protein reference required (use --" + PROTEIN + ")");
          valid = false;
        }
      }
      return valid;
    }
  }

  boolean hasSpecies() {
    return true;
  }

  boolean hasProtein() {
    return true;
  }

  @Override
  protected void initFlags() {
    initFlags(mFlags, hasSpecies(), hasProtein());
  }

  protected static void initFlags(CFlags flags, boolean species, boolean protein) {
    CommonFlagCategories.setCategories(flags);
    flags.setDescription("Runs the metagenomic composition and functional pipelines. The pipelines consist of read filtering, read alignment then species composition, and protein searching.");
    CommonFlags.initOutputDirFlag(flags);
    final Flag<File> input = flags.registerOptional(INPUT, File.class, CommonFlags.SDF_OR_FILE, "input SDF or single end fastq input").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    final Flag<File> left = flags.registerOptional(INPUT_LEFT, File.class, CommonFlags.FILE, "left arm of paired end FASTQ input").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    final Flag<File> right = flags.registerOptional(INPUT_RIGHT, File.class, CommonFlags.FILE, "right arm of paired end FASTQ input").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerOptional(PLATFORM, Platform.class, CommonFlags.STRING, "platform of the input data", Platform.ILLUMINA).setCategory(CommonFlagCategories.INPUT_OUTPUT);
    final String referencesDir = Environment.getEnvironmentMap().get(ENVIRONMENT_REFERENCES_DIR);
    final Flag<File> filterFlag = flags.registerOptional(FILTER, File.class, CommonFlags.SDF, "filter reference sequence").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    if (referencesDir != null) {
      filterFlag.setParameterDefault(new File(referencesDir, FILTER_REFERENCE_DEFAULT));
    }
    if (species) {
      final Flag<File> speciesFlag = flags.registerOptional(SPECIES, File.class, CommonFlags.SDF, "SDF containing reference sequences").setCategory(CommonFlagCategories.INPUT_OUTPUT);
      if (referencesDir != null) {
        speciesFlag.setParameterDefault(new File(referencesDir, SPECIES_REFERENCE_DEFAULT));
      }
    }
    if (protein) {
      final Flag<File> proteinFlag = flags.registerOptional(PROTEIN, File.class, CommonFlags.SDF, "SDF containing protein database").setCategory(CommonFlagCategories.INPUT_OUTPUT);
      if (referencesDir != null) {
        proteinFlag.setParameterDefault(new File(referencesDir, PROTEIN_REFERENCE_DEFAULT));
      }
    }
    flags.addRequiredSet(input);
    flags.addRequiredSet(left, right);
    flags.setValidator(new MetaWrapperValidator());
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(CommonFlags.OUTPUT_FLAG);
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    mErr = err;
    return super.mainExec(out, err);
  }

  @Override
  protected int mainExec(OutputStream out, LogStream log) throws IOException {
    mLogStream = log;
    final int ret = super.mainExec(out, log);
    if (ret != 0) {
      return ret;
    }
    if (mTask != null) {
      return mTask.returnCode();
    } else {
      return 1;
    }
  }

  @Override
  protected IORunnable task(MetaPipelineParams params, OutputStream out) {
    mTask = new MetagenomicsWrapperTask(params, out, mUsageMetric, mLogStream, mErr);
    return mTask;
  }

  @Override
  protected MetaPipelineParams makeParams() {
    final MetaPipelineParamsBuilder builder = MetaPipelineParams.builder().name(mFlags.getName());
    final File output = (File) mFlags.getValue(CommonFlags.OUTPUT_FLAG);
    final OutputParams outParams = new OutputParams(output, false);
    builder.outputParams(outParams);
    if (mFlags.isSet(INPUT)) {
      builder.inputFile((File) mFlags.getValue(INPUT));
    } else {
      builder.inputLeft((File) mFlags.getValue(INPUT_LEFT))
             .inputRight((File) mFlags.getValue(INPUT_RIGHT));
    }
    builder.inputPlatform((Platform) mFlags.getValue(PLATFORM));
    builder.filterSdf((File) mFlags.getValue(FILTER));
    if (mFlags.getFlag(SPECIES) != null) {
      builder.speciesSdf((File) mFlags.getValue(SPECIES));
    }
    if (mFlags.getFlag(PROTEIN) != null) {
      builder.proteinSdf((File) mFlags.getValue(PROTEIN));
    }
    return builder.create();
  }
}
