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
package com.rtg.variant.cnv;

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.UTILITY;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.io.Reader;
import java.io.Writer;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;


/**
 * Compare detected CNVs with the generated CNVs
 *
 */
public class CnvStatsCli extends AbstractCli {
  //private static final String VALIDATE = "validate";
  static final String GENERATED = "generated-cnvs";
  static final String DETECTED = "detected-cnvs";
  static final String OUTPUT_DIR = "outdir";
  static final String BUCKET_SIZE = "bucket-size";
  static final String COLLAPSE_CNV = "collapse-cnv";
  //private static final String PRIORS_FLAG = "priors";


  /**
   * Compare detected CNVs with the generated CNVs
   *
   * @param args arguments
   */
  public static void main(final String[] args) {
    new CnvStatsCli().mainExit(args);
  }

  private static class FlagsValidator implements Validator {
    @Override
    public boolean isValid(final CFlags cflags) {
      final File generated = (File) cflags.getValue(GENERATED);
      if (!generated.exists() || generated.isDirectory()) {
        cflags.setParseMessage("CNV generation file doesn't exist");
        return false;
      }
      final File detected = (File) cflags.getValue(DETECTED);
      if (!detected.exists() || detected.isDirectory()) {
        cflags.setParseMessage("CNV detection file doesn't exist");
        return false;
      }
      if (cflags.isSet(OUTPUT_DIR)) {
        if (!CommonFlags.validateOutputDirectory((File) cflags.getValue(OUTPUT_DIR))) {
          return false;
        }
      }
       return true;
    }
  }

  @Override
  protected void initFlags() {
    initFlags(mFlags);
  }

  /**
   * initialize a flags object
   * @param flags the flags object to initialize
   */
  public void initFlags(final CFlags flags) {
    flags.registerExtendedHelp();
    CommonFlagCategories.setCategories(flags);
    flags.registerRequired('g', GENERATED, File.class, "file", "CNV generation file").setCategory(INPUT_OUTPUT);
    flags.registerRequired('s', DETECTED, File.class, "file", "CNV detection file").setCategory(INPUT_OUTPUT);
    flags.registerOptional('o', OUTPUT_DIR, File.class, "dir", "output directory for results").setCategory(INPUT_OUTPUT);
    flags.registerOptional('b', BUCKET_SIZE, Integer.class, "int", "size of buckets used for CNV", 100).setCategory(UTILITY);
    flags.registerOptional(COLLAPSE_CNV, "collapse adjacent generated regions with identical CNV value").setCategory(UTILITY);
    flags.setValidator(new FlagsValidator());
  }


  /**
   * Compare detected CNVs with the generated CNVs
   *
   * @param outStr stream to deliver output to
   * @param errStr stream to print errors to
   * @throws IOException should an IO error occur
   * @return 0 on success 1 on failure
   */
  @Override
  protected int mainExec(final OutputStream outStr, final PrintStream errStr) throws IOException {
    final CFlags flags = mFlags;
    final CnvStats cs = new CnvStats();
    final Writer output = new OutputStreamWriter(outStr);
    try {
      final File generatedSnpFile = (File) flags.getValue(GENERATED);
      final File detectedSnpFile = (File) flags.getValue(DETECTED);

      try (Reader generatedReader = getReader(generatedSnpFile)) {
        try (Reader detectedReader = getReader(detectedSnpFile)) {

          if (flags.isSet(OUTPUT_DIR)) {
            final File outputDirectory = (File) flags.getValue(OUTPUT_DIR);
            if (!outputDirectory.exists() && !outputDirectory.mkdirs()) {
              throw new IOException("Couldn't create output directory");
            }
            cs.setOutput(outputDirectory);
          }
          final int bucketSize = (Integer) flags.getValue(BUCKET_SIZE);
          cs.setThreshold(bucketSize);
          cs.setCollapse(flags.isSet(COLLAPSE_CNV));
          try {
            cs.getStats(generatedReader, detectedReader, output);
          } finally {
            cs.closeOutput();
          }
        }
      }
    } catch (final NumberFormatException e) {
      Diagnostic.error("Problem parsing numbers from file: " + e.getMessage());
      return 1;
    } finally {
      output.flush();
    }
    return 0;
  }


  @Override
  public String moduleName() {
     return "cnvsimeval";
  }

  @Override
  public String description() {
    return "evaluate accuracy of calling CNVs on simulated data";
  }

  private BufferedReader getReader(final File f) throws IOException {
    return new BufferedReader(new InputStreamReader(FileUtils.createInputStream(f, true)));
  }

}
