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
package com.rtg.vcf;


import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.REPORTING;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.util.StringUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.vcf.header.VcfHeader;

/**
 * Calculate and output statistics for a <code>VCF</code> file.
 */
public class VcfStatsCli extends AbstractCli {

  private static final String SAMPLE = "sample";
  private static final String LENGTHS = "allele-lengths";
  private static final String COUNTS = "Xallele-counts";
  private static final String KNOWN = "known";
  private static final String NOVEL = "novel";

  @Override
  protected void initFlags() {
    CommonFlagCategories.setCategories(mFlags);
    mFlags.registerExtendedHelp();
    mFlags.setDescription("Display statistics from a set of VCF files.");
    mFlags.registerOptional("Xvariant", "calculate statistics via Variant API rather than VcfRecord API")
      .setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(SAMPLE, String.class, "STRING", "only calculate statistics for the specified sample (Default is to include all samples)")
      .setMaxCount(Integer.MAX_VALUE)
      .setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(KNOWN, "only calculate statistics for known variants (Default is to ignore known/novel status)").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(NOVEL, "only calculate statistics for novel variants (Default is to ignore known/novel status)").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(LENGTHS, "output variant length histogram").setCategory(REPORTING);
    mFlags.registerOptional(COUNTS, "output alleles per variant site histogram").setCategory(REPORTING);
    mFlags.registerRequired(File.class, "FILE", "input VCF files from which to derive statistics or '-' to read from standard input")
      .setMinCount(1)
      .setMaxCount(Integer.MAX_VALUE)
      .setCategory(INPUT_OUTPUT);

    mFlags.setValidator(new FlagValidator());
  }

  private static final class FlagValidator implements Validator {
    @Override
    public boolean isValid(CFlags flags) {
      final List<Object> files = flags.getAnonymousValues(0);
      if (files.size() > 1) {
        for (Object file : files) {
          if (CommonFlags.isStdio((File) file)) {
            flags.setParseMessage("Reading from standard in, not expecting other FILE arguments.");
            return false;
          }
        }
      }
      return flags.checkNand(KNOWN, NOVEL);
    }
  }

  @Override
  public String moduleName() {
    return "vcfstats";
  }

  @Override
  protected int mainExec(final OutputStream out, final PrintStream err) throws IOException {
    final List<File> inputs;
    final List<Object> files = mFlags.getAnonymousValues(0);
    // Kludge to let '-' as a single filename through the file checking to make std in possible
    if (files.size() == 1 && CommonFlags.isStdio((File) files.get(0))) {
      inputs = new ArrayList<>();
      inputs.add((File) files.get(0));
    } else {
      inputs = CommonFlags.getFileList(mFlags, null, null, false);
    }
    final boolean showHistograms = mFlags.isSet(LENGTHS);
    final boolean showAlleleCountHistograms = mFlags.isSet(COUNTS);
    String[] samples = null;
    if (mFlags.isSet(SAMPLE)) {
      samples = new String[mFlags.getValues(SAMPLE).size()];
      int i = 0;
      for (Object o : mFlags.getValues(SAMPLE)) {
        samples[i++] = (String) o;
      }
    }
    for (final File vcffile : inputs) {
      try (VcfReader vr = VcfReader.openVcfReader(vcffile)) {
        final VariantStatistics stats = new VariantStatistics(null);
        if (mFlags.isSet(KNOWN)) {
          stats.onlyKnown(true);
        } else if (mFlags.isSet(NOVEL)) {
          stats.onlyKnown(false);
        }
        stats.showLengthHistograms(showHistograms);
        stats.showAlleleCountHistograms(showAlleleCountHistograms);
        final VcfHeader header = vr.getHeader();
        if (samples != null) {
          stats.onlySamples(samples);
          for (String sample : samples) {
            if (!header.getSampleNames().contains(sample)) {
              Diagnostic.warning("Specified sample '" + sample + "' is not contained in file: " + vcffile);
              continue;
            }
          }
        }
        while (vr.hasNext()) {
          //System.err.println(vr.getCurrent());
          stats.tallyVariant(header, vr.next());
        }
        out.write(("Location                     : " + vcffile + StringUtils.LS).getBytes());
        stats.printStatistics(out);
      }
    }
    return 0;
  }

  /**
   * Test from the command line - just dump statistics for an existing VCF file
   * @param args command line arguments.
   */
  public static void main(String... args) {
    new VcfStatsCli().mainExit(args);
  }
}
