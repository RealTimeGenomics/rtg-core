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

import static com.rtg.launcher.CommonFlags.NO_GZIP;
import static com.rtg.util.cli.CommonFlagCategories.FILTERING;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.REPORTING;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

import com.rtg.bed.BedUtils;
import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.sam.SamRangeUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.util.io.FileUtils;
import com.rtg.vcf.VcfFilterStatistics.Stat;

/**
 * Class inputs the variant results and filters them based on various criteria.
 *
 */
public final class VcfFilterCli extends AbstractCli {

  // flags
  private static final String MODULE_NAME = "vcffilter";
  private static final String INPUT = "input";
  private static final String OUTPUT = "output";
  private static final String RESTRICTION_FLAG = "region";

  // filter flags
  private static final String REMOVE_INFO = "remove-info";
  private static final String KEEP_INFO = "keep-info";

  private static final String REMOVE_FILTER = "remove-filter";
  private static final String KEEP_FILTER = "keep-filter";
  private static final String INCLUDE_BED = "include-bed";
  private static final String EXCLUDE_BED = "exclude-bed";
  private static final String INCLUDE_VCF = "include-vcf";
  private static final String EXCLUDE_VCF = "exclude-vcf";

  private static final String REMOVE_SAME_AS_REF = "remove-same-as-ref";
  private static final String REMOVE_ALL_SAME_AS_REF = "remove-all-same-as-ref";

  private static final String SNPS_ONLY = "snps-only";
  private static final String NON_SNPS_ONLY = "non-snps-only";

  private static final String MIN_DEPTH = "min-read-depth";
  private static final String MAX_DEPTH = "max-read-depth";
  private static final String MIN_COMBINED_DEPTH = "min-combined-read-depth";
  private static final String MAX_COMBINED_DEPTH = "max-combined-read-depth";
  private static final String MIN_GENOTYPE_QUALITY = "min-genotype-quality";
  private static final String MAX_GENOTYPE_QUALITY = "max-genotype-quality";
  private static final String MIN_QUALITY = "min-quality";
  private static final String MAX_QUALITY = "max-quality";
  private static final String MAX_AMBIGUITY_RATIO = "max-ambiguity-ratio";
  private static final String MIN_AVR_SCORE = "min-avr-score";
  private static final String MAX_AVR_SCORE = "max-avr-score";
  private static final String MIN_DENOVO_SCORE = "min-denovo-score";
  private static final String MAX_DENOVO_SCORE = "max-denovo-score";

  private static final String FAIL_FLAG = "fail";
  private static final String CLEAR_FAILED_SAMPLES = "clear-failed-samples";

  private static final String SAMPLE = "sample";
  private static final String ALL_SAMPLES = "all-samples";

  private static final String DENSITY_WINDOW = "density-window";
  private static final String REMOVE_OVERLAPPING = "remove-overlapping";

  // old filter flags
  private static final String MIN_POSTERIOR_SCORE = "Xmin-posterior-score";
  private static final String MAX_POSTERIOR_SCORE = "Xmax-posterior-score";


  private final VcfFilterTask mVcfFilterTask = new VcfFilterTask();


  @Override
  protected void initFlags() {
    mFlags.registerExtendedHelp();
    mFlags.setDescription("Filters VCF records based on various criteria. When filtering on multiple samples, if any of the specified samples fail the criteria, the record will be filtered.");
    CommonFlagCategories.setCategories(mFlags);

    mFlags.registerRequired('i', INPUT, File.class, "file", "VCF file containing variants to be filtered. Use '-' to read from standard input").setCategory(INPUT_OUTPUT);
    mFlags.registerRequired('o', OUTPUT, File.class, "file", "output VCF file. Use '-' to write to standard output").setCategory(INPUT_OUTPUT);
    CommonFlags.initNoGzip(mFlags);
    CommonFlags.initIndexFlags(mFlags);

    mFlags.registerOptional(RESTRICTION_FLAG, String.class, "string", "if set, only read VCF records within the specified range. The format is one of <template_name>, <template_name>:start-end or <template_name>:start+length").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(CommonFlags.BED_REGIONS_FLAG, File.class, "File", "if set, only read VCF records that overlap the ranges contained in the specified BED file").setCategory(INPUT_OUTPUT);

    // What to apply to, what to do with the results of filtering
    mFlags.registerOptional(FAIL_FLAG, String.class, "STRING", "instead of removing failed records set their filter field to the provided value").setCategory(REPORTING);
    mFlags.registerOptional(CLEAR_FAILED_SAMPLES, "instead of removing failed records set the sample GT fields to missing").setCategory(REPORTING);

    // Variant position
    mFlags.registerOptional(INCLUDE_BED, File.class, "File", "only keep variants within the regions in this BED file").setCategory(FILTERING);
    mFlags.registerOptional(EXCLUDE_BED, File.class, "File", "discard all variants within the regions in this BED file").setCategory(FILTERING);
    mFlags.registerOptional(INCLUDE_VCF, File.class, "File", "only keep variants that overlap with the ones in this file").setCategory(FILTERING);
    mFlags.registerOptional(EXCLUDE_VCF, File.class, "File", "discard all variants that overlap with the ones in this file").setCategory(FILTERING);
    mFlags.registerOptional('w', DENSITY_WINDOW, Integer.class, "INT", "window within which multiple variants are discarded").setCategory(FILTERING);
    mFlags.registerOptional(REMOVE_OVERLAPPING, "remove records that overlap with previous records").setCategory(FILTERING);

    // REF/ALT contents
    mFlags.registerOptional(SNPS_ONLY, "if set, will output simple SNPs only").setCategory(FILTERING);
    mFlags.registerOptional(NON_SNPS_ONLY, "if set, will output MNPs and INDELs only").setCategory(FILTERING);

    // Contents of FILTER
    mFlags.registerOptional('r', REMOVE_FILTER, String.class, "STRING", "remove variants with this FILTER tag").setCategory(FILTERING).setMinCount(0).setMaxCount(Integer.MAX_VALUE);
    mFlags.registerOptional('k', KEEP_FILTER, String.class, "STRING", "only keep variants with this FILTER tag").setCategory(FILTERING).setMinCount(0).setMaxCount(Integer.MAX_VALUE);

    // Contents of INFO
    mFlags.registerOptional('R', REMOVE_INFO, String.class, "STRING", "remove variants with this INFO tag").setCategory(FILTERING).setMinCount(0).setMaxCount(Integer.MAX_VALUE);
    mFlags.registerOptional('K', KEEP_INFO, String.class, "STRING", "only keep variants with this INFO tag").setCategory(FILTERING).setMinCount(0).setMaxCount(Integer.MAX_VALUE);

    // Which FORMAT columns do we look at
    mFlags.registerOptional(SAMPLE, String.class, "STRING", "apply sample-specific criteria to the named sample contained in the input VCF").setCategory(INPUT_OUTPUT).setMaxCount(Integer.MAX_VALUE);
    mFlags.registerOptional(ALL_SAMPLES, "apply sample-specific criteria to all samples contained in the input VCF").setCategory(INPUT_OUTPUT);

    // FORMAT/GT
    mFlags.registerOptional(REMOVE_SAME_AS_REF, "remove where sample is same as reference").setCategory(FILTERING);
    mFlags.registerOptional(REMOVE_ALL_SAME_AS_REF, "remove where all samples are same as reference").setCategory(FILTERING);

    // Other INFO fields
    mFlags.registerOptional('c', MIN_COMBINED_DEPTH, Integer.class, "INT", "minimum allowed combined read depth").setCategory(FILTERING);
    mFlags.registerOptional('C', MAX_COMBINED_DEPTH, Integer.class, "INT", "maximum allowed combined read depth").setCategory(FILTERING);

    // Other FORMAT fields
    mFlags.registerOptional('d', MIN_DEPTH, Integer.class, "INT", "minimum allowed sample read depth").setCategory(FILTERING);
    mFlags.registerOptional('D', MAX_DEPTH, Integer.class, "INT", "maximum allowed sample read depth").setCategory(FILTERING);
    mFlags.registerOptional('g', MIN_GENOTYPE_QUALITY, Double.class, "float", "minimum allowed genotype quality").setCategory(FILTERING);
    mFlags.registerOptional('G', MAX_GENOTYPE_QUALITY, Double.class, "float", "maximum allowed genotype quality").setCategory(FILTERING);
    mFlags.registerOptional('q', MIN_QUALITY, Double.class, "float", "minimum allowed quality").setCategory(FILTERING);
    mFlags.registerOptional('Q', MAX_QUALITY, Double.class, "float", "maximum allowed quality").setCategory(FILTERING);
    mFlags.registerOptional('A', MAX_AMBIGUITY_RATIO, Double.class, "float", "maximum allowed ambiguity ratio").setCategory(FILTERING);
    mFlags.registerOptional(MIN_AVR_SCORE, Double.class, "float", "minimum allowed AVR score").setCategory(FILTERING);
    mFlags.registerOptional(MAX_AVR_SCORE, Double.class, "float", "maximum allowed AVR score").setCategory(FILTERING);
    mFlags.registerOptional(MIN_DENOVO_SCORE, Double.class, "float", "minimum de novo score threshold").setCategory(FILTERING);
    mFlags.registerOptional(MAX_DENOVO_SCORE, Double.class, "float", "maximum de novo score threshold").setCategory(FILTERING);

    // Xflags
    mFlags.registerOptional('p', MIN_POSTERIOR_SCORE, Double.class, "float", "minimum allowed posterior score").setCategory(FILTERING);
    mFlags.registerOptional('P', MAX_POSTERIOR_SCORE, Double.class, "float", "maximum allowed posterior score").setCategory(FILTERING);

    mFlags.setValidator(new VcfFilterValidator());
  }

  private static class VcfFilterValidator implements Validator {
    @Override
    public boolean isValid(final CFlags flags) {
      final File input = (File) flags.getValue(INPUT);
      if (!CommonFlags.isStdio(input)) {
        if (!input.exists()) {
          flags.setParseMessage("Given file \"" + input.getPath() + "\" does not exist.");
          return false;
        }
        if (input.isDirectory()) {
          flags.setParseMessage("Given file \"" + input.getPath() + "\" is a directory.");
          return false;
        }
      }
      final File o = (File) flags.getValue(OUTPUT);
      if (!CommonFlags.isStdio(o)) {
        final File output = FileUtils.getZippedFileName(!flags.isSet(NO_GZIP), o);
        if (output.exists()) {
          flags.setParseMessage("The file \"" + output + "\" already exists. Please remove it first or choose a different file");
          return false;
        }
      }
      if ((flags.isSet(MIN_GENOTYPE_QUALITY) || flags.isSet(MAX_GENOTYPE_QUALITY)) && (flags.isSet(MIN_POSTERIOR_SCORE) || flags.isSet(MAX_POSTERIOR_SCORE))) {
        //Only possible if someone is monkeying with X-flags
        flags.setParseMessage("Use genotype-quality or posterior filters, not both.");
        return false;
      }
      if (!isMinMaxValid(flags, MIN_GENOTYPE_QUALITY, MAX_GENOTYPE_QUALITY, 0.0, Double.MAX_VALUE, true)) {
        return false;
      }
      if (!isMinMaxValid(flags, MIN_POSTERIOR_SCORE, MAX_POSTERIOR_SCORE, 0.0, Double.MAX_VALUE, true)) {
        return false;
      }
      if (!isMinMaxValid(flags, MIN_QUALITY, MAX_QUALITY, 0.0, Double.MAX_VALUE, true)) {
        return false;
      }
      if (!isMinMaxValid(flags, null, MAX_AMBIGUITY_RATIO, 0.0, 1.0, true)) {
        return false;
      }
      if (!isMinMaxValid(flags, MIN_AVR_SCORE, MAX_AVR_SCORE, 0.0, Double.MAX_VALUE, true)) {
        return false;
      }
      if (!isMinMaxValid(flags, MIN_DENOVO_SCORE, MAX_DENOVO_SCORE, 0.0, Double.MAX_VALUE, true)) {
        return false;
      }
      if (flags.isSet(MIN_DENOVO_SCORE) || flags.isSet(MAX_DENOVO_SCORE)) {
        if (flags.isSet(ALL_SAMPLES) || flags.getValues(SAMPLE).size() != 1) {
          flags.error("De Novo filtering requires a single sample to be specified");
          return false;
        }
      }
      if (flags.isSet(DENSITY_WINDOW) && !CommonFlags.validateFlagBetweenValues(flags, DENSITY_WINDOW, 1, Integer.MAX_VALUE)) {
        return false;
      }
      if (!validateDepthFlags(flags, MIN_DEPTH, MAX_DEPTH) || !validateDepthFlags(flags, MIN_COMBINED_DEPTH, MAX_COMBINED_DEPTH)) {
        return false;
      }
      for (final String flag : new String[] {EXCLUDE_BED, EXCLUDE_VCF, INCLUDE_BED, INCLUDE_VCF}) {
        if (validateFile(flags, flag)) {
          return false;
        }
      }
      if (!CommonFlags.validateRegions(flags)) {
        return false;
      }
      if (!flags.checkNand(EXCLUDE_BED, EXCLUDE_VCF)) {
        return false;
      }
      if (!flags.checkNand(INCLUDE_BED, INCLUDE_VCF)) {
        return false;
      }
      if (!flags.checkNand(SNPS_ONLY, NON_SNPS_ONLY)) {
        return false;
      }
      if (!flags.checkNand(SAMPLE, ALL_SAMPLES)) {
        return false;
      }
      if (!flags.checkNand(FAIL_FLAG, CLEAR_FAILED_SAMPLES)) {
        return false;
      }
      return true;
    }

    static boolean validateDepthFlags(CFlags flags, String minFlag, String maxFlag) {
      int c = 0;
      if (flags.isSet(minFlag)) {
        c = (Integer) flags.getValue(minFlag);
        if (c < 0) {
          flags.setParseMessage("The specified flag \"--" + minFlag + "\" has invalid value \"" + c + "\". It should be greater than or equal to 0.");
          return false;
        }
      }
      if (flags.isSet(maxFlag)) {
        final int maxc = (Integer) flags.getValue(maxFlag);
        if (maxc < c) {
          flags.setParseMessage("The specified flag \"--" + maxFlag + "\" has invalid value \"" + maxc + "\". It should be greater than or equal to " + c + ".");
          return false;
        }
      }
      return true;
    }

    static boolean validateFile(CFlags flags, String flag) {
      if (flags.isSet(flag)) {
        final File file = (File) flags.getValue(flag);
        if (!file.exists()) {
          flags.setParseMessage("The \"--" + flag + "\" file: \"" + file + "\" doesn't exist.");
          return true;
        }
        if (file.isDirectory()) {
          flags.setParseMessage("The \"--" + flag + "\" file: \"" + file + "\" is a directory.");
          return true;
        }
      }
      return false;
    }

    static boolean isMinMaxValid(final CFlags flags, final String minFlag, final String maxFlag, double min, double max, boolean maxInclusive) {
      double minp = min;
      if (minFlag != null && flags.isSet(minFlag)) {
        minp = (Double) flags.getValue(minFlag);
        if (minp < min || minp > max || (!maxInclusive && minp >= max)) {
          flags.setParseMessage(getMinMaxMessage(minFlag, minp, min, max, maxInclusive));
          return false;
        }
      }
      if (flags.isSet(maxFlag)) {
        final double maxp = (Double) flags.getValue(maxFlag);
        if (maxp < minp || maxp > max || (!maxInclusive && (maxp >= max || minp >= maxp))) {
          flags.setParseMessage(getMinMaxMessage(maxFlag, maxp, minp, max, maxInclusive));
          return false;
        }
      }
      return true;
    }

    static String getMinMaxMessage(final String flag, final double val, final double min, final double max, boolean maxInclusive) {
      return "The specified flag \"--" + flag + "\" has invalid value \"" + val + "\". It should be greater than or equal to " + min + " and less than " + (maxInclusive ? "or equal to " : "") + getDoubleString(max) + ".";
    }

    static String getDoubleString(final double val) {
      if (val >= Double.MAX_VALUE) {
        return "Infinity";
      }
      return String.valueOf(val);
    }
  }

  /**
   */
  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  @Override
  protected int mainExec(final OutputStream out, final PrintStream err) throws IOException {
    try {
      process(out);
      return 0;
    } finally {
      out.flush();
    }
  }

  private void process(final OutputStream output) throws IOException {

    if (mFlags.isSet(SAMPLE)) {
      for (Object o : mFlags.getValues(SAMPLE)) {
        mVcfFilterTask.mSampleNames.add((String) o);
      }
    }

    mVcfFilterTask.mResetFailedSampleGts =  mFlags.isSet(CLEAR_FAILED_SAMPLES);

    // These involve checking the specified sample field
    mVcfFilterTask.mRemoveSameAsRef = mFlags.isSet(REMOVE_SAME_AS_REF);
    mVcfFilterTask.mSnpsOnly = mFlags.isSet(SNPS_ONLY);
    mVcfFilterTask.mNonSnpsOnly = mFlags.isSet(NON_SNPS_ONLY);

    if (mFlags.isSet(MIN_QUALITY) || mFlags.isSet(MAX_QUALITY)) {
      final double minQuality = mFlags.isSet(MIN_QUALITY) ? (Double) mFlags.getValue(MIN_QUALITY) : 0.0;
      final double maxQuality = mFlags.isSet(MAX_QUALITY) ? (Double) mFlags.getValue(MAX_QUALITY) : Double.MAX_VALUE;
      mVcfFilterTask.mFilters.add(new VcfFilter.QualFilter(mVcfFilterTask.mVcfFilterStatistics, minQuality, maxQuality));
    }
    if (mFlags.isSet(MIN_GENOTYPE_QUALITY) || mFlags.isSet(MAX_GENOTYPE_QUALITY)
        || mFlags.isSet(MIN_POSTERIOR_SCORE) || mFlags.isSet(MAX_POSTERIOR_SCORE)) {
      final double minGq = mFlags.isSet(MIN_GENOTYPE_QUALITY) ? (Double) mFlags.getValue(MIN_GENOTYPE_QUALITY) : mFlags.isSet(MIN_POSTERIOR_SCORE) ? (Double) mFlags.getValue(MIN_POSTERIOR_SCORE) : 0.0;
      final double maxGq = mFlags.isSet(MAX_GENOTYPE_QUALITY) ? (Double) mFlags.getValue(MAX_GENOTYPE_QUALITY) : mFlags.isSet(MAX_POSTERIOR_SCORE) ? (Double) mFlags.getValue(MAX_POSTERIOR_SCORE) : Double.MAX_VALUE;
      final boolean posteriorFiltering = mFlags.isSet(MIN_POSTERIOR_SCORE) || mFlags.isSet(MAX_POSTERIOR_SCORE);
      mVcfFilterTask.mVcfFilterStatistics.setPosteriorFiltering(posteriorFiltering);
      mVcfFilterTask.mFilters.add(new VcfSampleFilter.GqFilter(mVcfFilterTask.mVcfFilterStatistics, minGq, maxGq, posteriorFiltering));
    }
    if (mFlags.isSet(MIN_DEPTH) || mFlags.isSet(MAX_DEPTH)) {
      final int minReadDepth = mFlags.isSet(MIN_DEPTH) ? (Integer) mFlags.getValue(MIN_DEPTH) : 0;
      final int maxReadDepth = mFlags.isSet(MAX_DEPTH) ? (Integer) mFlags.getValue(MAX_DEPTH) : Integer.MAX_VALUE;
      mVcfFilterTask.mFilters.add(new VcfSampleFilter.MinMaxIntFilter(mVcfFilterTask.mVcfFilterStatistics, Stat.READ_DEPTH_FILTERED_COUNT, minReadDepth, maxReadDepth, VcfUtils.FORMAT_SAMPLE_DEPTH));
    }
    if (mFlags.isSet(MIN_COMBINED_DEPTH) || mFlags.isSet(MAX_COMBINED_DEPTH)) {
      final int minCombinedReadDepth = mFlags.isSet(MIN_COMBINED_DEPTH) ? (Integer) mFlags.getValue(MIN_COMBINED_DEPTH) : 0;
      final int maxCombinedReadDepth = mFlags.isSet(MAX_COMBINED_DEPTH) ? (Integer) mFlags.getValue(MAX_COMBINED_DEPTH) : Integer.MAX_VALUE;
      mVcfFilterTask.mFilters.add(new VcfInfoFilter.MinMaxIntFilter(mVcfFilterTask.mVcfFilterStatistics, Stat.COMBINED_READ_DEPTH_FILTERED_COUNT, minCombinedReadDepth, maxCombinedReadDepth, VcfUtils.INFO_COMBINED_DEPTH));
    }
    if (mFlags.isSet(MIN_AVR_SCORE) || mFlags.isSet(MAX_AVR_SCORE)) {
      final double minAvrScore = mFlags.isSet(MIN_AVR_SCORE) ? (Double) mFlags.getValue(MIN_AVR_SCORE) : Double.NEGATIVE_INFINITY;
      final double maxAvrScore = mFlags.isSet(MAX_AVR_SCORE) ? (Double) mFlags.getValue(MAX_AVR_SCORE) : Double.POSITIVE_INFINITY;
      mVcfFilterTask.mFilters.add(new VcfSampleFilter.MinMaxDoubleFilter(mVcfFilterTask.mVcfFilterStatistics, Stat.AVR_SCORE_FILTERED_COUNT, minAvrScore, maxAvrScore, VcfUtils.FORMAT_AVR));
    }
    if (mFlags.isSet(MAX_AMBIGUITY_RATIO)) {
      final double maxAmbiguityRatio = mFlags.isSet(MAX_AMBIGUITY_RATIO) ? (Double) mFlags.getValue(MAX_AMBIGUITY_RATIO) : Double.POSITIVE_INFINITY;
      mVcfFilterTask.mFilters.add(new VcfSampleFilter.MinMaxDoubleFilter(mVcfFilterTask.mVcfFilterStatistics, Stat.AMBIGOUS_FILTERED_COUNT, Double.NEGATIVE_INFINITY, maxAmbiguityRatio, VcfUtils.FORMAT_AMBIGUITY_RATIO));
    }
    if (mFlags.isSet(MAX_DENOVO_SCORE) || mFlags.isSet(MIN_DENOVO_SCORE)) {
      mVcfFilterTask.mFilters.add(new VcfSampleFilter.DenovoFilter(mVcfFilterTask.mVcfFilterStatistics,
          mFlags.isSet(MIN_DENOVO_SCORE) ? (Double) mFlags.getValue(MIN_DENOVO_SCORE) : 0,
          mFlags.isSet(MAX_DENOVO_SCORE) ? (Double) mFlags.getValue(MAX_DENOVO_SCORE) : Double.MAX_VALUE));
    }

    final String[] sampleFlags = {
        REMOVE_SAME_AS_REF,
        SNPS_ONLY, NON_SNPS_ONLY,
        MIN_DEPTH, MAX_DEPTH,
        MIN_GENOTYPE_QUALITY, MAX_GENOTYPE_QUALITY,
        MIN_POSTERIOR_SCORE, MAX_POSTERIOR_SCORE,
        MAX_AMBIGUITY_RATIO,
        MIN_AVR_SCORE, MAX_AVR_SCORE,
        MIN_DENOVO_SCORE, MAX_DENOVO_SCORE,
    };
    mVcfFilterTask.mCheckingSample = false;
    for (final String flag : sampleFlags) {
      mVcfFilterTask.mCheckingSample |= mFlags.isSet(flag);
    }

    mVcfFilterTask.mRemoveOverlapping = mFlags.isSet(REMOVE_OVERLAPPING);
    mVcfFilterTask.mRemoveAllSameAsRef = mFlags.isSet(REMOVE_ALL_SAME_AS_REF);
    mVcfFilterTask.mDensityWindow = mFlags.isSet(DENSITY_WINDOW) ? (Integer) mFlags.getValue(DENSITY_WINDOW) : null;

    if (mFlags.isSet(KEEP_FILTER)) {
      for (final Object tag : mFlags.getValues(KEEP_FILTER)) {
        mVcfFilterTask.mKeepFilters.add((String) tag);
      }
    }
    if (mFlags.isSet(REMOVE_FILTER)) {
      for (final Object tag : mFlags.getValues(REMOVE_FILTER)) {
        mVcfFilterTask.mRemoveFilters.add((String) tag);
      }
    }
    if (mFlags.isSet(KEEP_INFO)) {
      for (final Object tag : mFlags.getValues(KEEP_INFO)) {
        mVcfFilterTask.mKeepInfos.add((String) tag);
      }
    }
    if (mFlags.isSet(REMOVE_INFO)) {
      for (final Object tag : mFlags.getValues(REMOVE_INFO)) {
        mVcfFilterTask.mRemoveInfos.add((String) tag);
      }
    }
    if (mFlags.isSet(INCLUDE_BED)) {
      mVcfFilterTask.mIncludeBed = BedUtils.regions((File) mFlags.getValue(INCLUDE_BED));
    } else if (mFlags.isSet(INCLUDE_VCF)) {
      mVcfFilterTask.mIncludeBed = VcfUtils.regionsVcf((File) mFlags.getValue(INCLUDE_VCF));
    }
    if (mFlags.isSet(EXCLUDE_BED)) {
      mVcfFilterTask.mExcludeBed = BedUtils.regions((File) mFlags.getValue(EXCLUDE_BED));
    } else if (mFlags.isSet(EXCLUDE_VCF)) {
      mVcfFilterTask.mExcludeBed = VcfUtils.regionsVcf((File) mFlags.getValue(EXCLUDE_VCF));
    }
    final RegionRestriction region;
    if (mFlags.isSet(CommonFlags.RESTRICTION_FLAG)) {
      region = new RegionRestriction((String) mFlags.getValue(RESTRICTION_FLAG));
    } else {
      region = null;
    }
    final ReferenceRanges ranges;
    if (mFlags.isSet(CommonFlags.BED_REGIONS_FLAG)) {
      Diagnostic.developerLog("Loading BED regions");
      ranges = SamRangeUtils.createBedReferenceRanges((File) mFlags.getValue(CommonFlags.BED_REGIONS_FLAG));
    } else {
      ranges = null;
    }
    mVcfFilterTask.mAllSamples = mFlags.isSet(ALL_SAMPLES);
    mVcfFilterTask.mFailFilterName = mFlags.isSet(FAIL_FLAG) ? (String) mFlags.getValue(FAIL_FLAG) : null;

    final File in = (File) mFlags.getValue(INPUT);
    final File out = (File) mFlags.getValue(OUTPUT);
    final boolean gzip = !mFlags.isSet(NO_GZIP);
    final boolean index = !mFlags.isSet(CommonFlags.NO_INDEX);
    final boolean stdout = CommonFlags.isStdio(out);
    Diagnostic.developerLog("Starting filter");
    try (VcfReader r = ranges != null ? VcfReader.openVcfReader(in, ranges) : VcfReader.openVcfReader(in, region)) {
      final File vcfFile = stdout ? null : FileUtils.getZippedFileName(gzip, out);
      try (VcfWriter w = new VcfWriter(r.getHeader(), vcfFile, output, gzip, index)) {
        mVcfFilterTask.filterVcf(r, w);
      }
      if (!stdout) {
        mVcfFilterTask.printStatistics(output);
      }
    }
  }

}
