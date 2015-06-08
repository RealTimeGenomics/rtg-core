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
package com.rtg.variant.eval;

import static com.rtg.launcher.BuildCommon.RESOURCE;
import static com.rtg.util.cli.CommonFlagCategories.FILTERING;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.REPORTING;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.Locale;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.NoStatistics;
import com.rtg.launcher.OutputParams;
import com.rtg.launcher.ParamsCli;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.IORunnable;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.variant.eval.VcfEvalParams.VcfEvalParamsBuilder;
import com.rtg.vcf.VcfUtils;
import com.rtg.vcf.annotation.DerivedAnnotations;

/**
 * Compare detected SNPs with the generated SNPs
 *
 */
public class VcfEvalCli extends ParamsCli<VcfEvalParams> {

  private static final String MODULE_NAME = "vcfeval";
  private static final String BASELINE = "baseline";
  private static final String CALLS = "calls";
  private static final String SORT_ORDER = "sort-order";
  private static final String ALL_RECORDS = "all-records";
  private static final String SAMPLE = "sample";
  private static final String SORT_FIELD = "vcf-score-field";
  private static final String SQUASH_PLOIDY = "squash-ploidy";
  private static final String BASELINE_TP = "baseline-tp";
  private static final String SLOPE_FILES = "slope-files";

  private static final String MAX_LENGTH = "Xmax-length";
  private static final String RTG_STATS = "Xrtg-stats";

  /**
   * Compare detected SNPs with the generated SNPs
   *
   * @param args arguments
   */
  public static void main(final String[] args) {
    new VcfEvalCli().mainExit(args);
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(CommonFlags.OUTPUT_FLAG);
  }

  @Override
  protected void initFlags() {
    initFlags(mFlags);
    mFlags.setName(applicationName() + " " + moduleName());
  }

  /**
   * initialize a flags object
   * @param flags the flags object to initialize
   */
  public static void initFlags(final CFlags flags) {
    CommonFlagCategories.setCategories(flags);
    flags.registerExtendedHelp();
    flags.setDescription("Evaluates called variants for genotype agreement with a baseline variant set irrespective of representational differences. Outputs a weighted ROC file which can be viewed with rtg rocplot and separate VCF files containing false positives (called variants not matched in the baseline), false negatives (baseline variants not matched in the call set), and true positives (called variants matched in the baseline).");
    flags.registerRequired('o', CommonFlags.OUTPUT_FLAG, File.class, "DIR", RESOURCE.getString("OUTPUT_DESC")).setCategory(INPUT_OUTPUT);
    flags.registerRequired('b', BASELINE, File.class, "file", "VCF file containing baseline variants").setCategory(INPUT_OUTPUT);
    flags.registerRequired('c', CALLS, File.class, "file", "VCF file containing called variants").setCategory(INPUT_OUTPUT);
    flags.registerRequired('t', CommonFlags.TEMPLATE_FLAG, File.class, "SDF", "SDF of the reference genome the variants are called against").setCategory(INPUT_OUTPUT);

    flags.registerOptional(CommonFlags.RESTRICTION_FLAG, String.class, "string", "if set, only read VCF records within the specified range. The format is one of <template_name>, <template_name>:start-end or <template_name>:start+length").setCategory(INPUT_OUTPUT);
    flags.registerOptional(CommonFlags.BED_REGIONS_FLAG, File.class, "File", "if set, only read VCF records that overlap the ranges contained in the specified BED file").setCategory(INPUT_OUTPUT);

    flags.registerOptional(SAMPLE, String.class, "STRING", "the name of the sample to select (required when using multi-sample VCF files)").setCategory(FILTERING);
    flags.registerOptional(ALL_RECORDS, "use all records regardless of FILTER status. Default is to only process records where FILTER is \".\" or \"PASS\"").setCategory(FILTERING);
    flags.registerOptional(SQUASH_PLOIDY, "treat heterozygous genotypes as homozygous ALT in both baseline and calls").setCategory(FILTERING);

    flags.registerOptional('f', SORT_FIELD, String.class, "STRING", "the name of the VCF FORMAT field to use as the ROC score. Also valid are \"QUAL\" or \"INFO=<name>\" to select the named VCF INFO field", VcfUtils.FORMAT_GENOTYPE_QUALITY).setCategory(REPORTING);
    flags.registerOptional('O', SORT_ORDER, RocSortOrder.class, "STRING", "the order in which to sort the ROC scores so that \"good\" scores come before \"bad\" scores", RocSortOrder.DESCENDING).setCategory(REPORTING);
    flags.registerOptional(MAX_LENGTH, Integer.class, "INT", "don't attempt to evaluate variant alternatives longer than this", 1000).setCategory(FILTERING);
    flags.registerOptional(RTG_STATS, "output RTG specific files and statistics").setCategory(REPORTING);
    flags.registerOptional(BASELINE_TP, "output an additional file containing the baseline version of true positive variants").setCategory(REPORTING);
    flags.registerOptional(SLOPE_FILES, "output files for ROC slope analysis").setCategory(REPORTING);

    CommonFlags.initThreadsFlag(flags);
    CommonFlags.initNoGzip(flags);
    flags.setValidator(new VcfEvalFlagsValidator());
  }

  /**
   * @return current name of the module
   */
  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  static class VcfEvalFlagsValidator implements Validator {

    @Override
    public boolean isValid(final CFlags flags) {
      if (!CommonFlags.validateOutputDirectory(flags)) {
        return false;
      }

      final File generated = (File) flags.getValue(BASELINE);
      if (!(generated.exists() && generated.isFile())) {
        flags.setParseMessage("baseline VCF file doesn't exist");
        return false;
      }
      final File detected = (File) flags.getValue(CALLS);
      if (!(detected.exists() && detected.isFile())) {
        flags.setParseMessage("calls VCF file doesn't exist");
        return false;
      }
      if (flags.isSet(SORT_FIELD)) {
        final String field = (String) flags.getValue(SORT_FIELD);
        final int pIndex = field.indexOf('=');
        if (pIndex != -1) {
          final String fieldTypeName = field.substring(0, pIndex).toUpperCase(Locale.getDefault());
          try {
            final RocScoreField f = RocScoreField.valueOf(fieldTypeName);
            if (f == RocScoreField.DERIVED) {
              try {
                final DerivedAnnotations ann = DerivedAnnotations.valueOf(field.substring(pIndex + 1).toUpperCase(Locale.getDefault()));
                if (!DerivedAnnotations.singleValueNumericAnnotations().contains(ann)) {
                  throw new IllegalArgumentException("Non single value numeric annotation");
                }
              } catch (IllegalArgumentException e) {
                flags.setParseMessage("Unrecognized derived annotation \"" + field + "\", must be one of " + Arrays.toString(DerivedAnnotations.singleValueNumericAnnotations().toArray()));
                return false;
              }
            }
          } catch (IllegalArgumentException e) {
            flags.setParseMessage("Unrecognized field type \"" + fieldTypeName + "\", must be one of " + Arrays.toString(RocScoreField.values()));
            return false;
          }
        }
      }
      if (!CommonFlags.validateThreads(flags)) {
        return false;
      }
      if (!CommonFlags.validateTemplate(flags)) {
        return false;
      }
      if (!CommonFlags.validateRegions(flags)) {
        return false;
      }
      return true;
    }
  }


  @Override
  protected IORunnable task(VcfEvalParams params, OutputStream out) {
    return new VcfEvalTask(params, out, new NoStatistics());
  }

  private void checkTabix(File vcfFile) throws IOException {
    final File index = TabixIndexer.indexFileName(vcfFile);
    if (!TabixIndexer.isBlockCompressed(vcfFile)) {
      throw new NoTalkbackSlimException(vcfFile + " is not in bgzip format");
    } else if (!index.exists()) {
      throw new NoTalkbackSlimException("Index not found for file: " + index.getPath() + " expected index called: " + index.getPath());
    }
  }

  @Override
  protected VcfEvalParams makeParams() throws InvalidParamsException, IOException {
    final VcfEvalParamsBuilder builder = VcfEvalParams.builder();
    builder.name(mFlags.getName());
    builder.outputParams(new OutputParams(outputDirectory(), false, !mFlags.isSet(CommonFlags.NO_GZIP)));
    builder.templateFile((File) mFlags.getValue(CommonFlags.TEMPLATE_FLAG));
    final File baseLine = (File) mFlags.getValue(BASELINE);
    final File calls = (File) mFlags.getValue(CALLS);
    checkTabix(baseLine);
    checkTabix(calls);
    builder.baseLineFile(baseLine).callsFile(calls);
    builder.sortOrder((RocSortOrder) mFlags.getValue(SORT_ORDER));
    builder.scoreField((String) mFlags.getValue(SORT_FIELD));
    builder.maxLength((Integer) mFlags.getValue(MAX_LENGTH));
    if (mFlags.isSet(CommonFlags.RESTRICTION_FLAG)) {
      builder.restriction(new RegionRestriction((String) mFlags.getValue(CommonFlags.RESTRICTION_FLAG)));
    }
    if (mFlags.isSet(CommonFlags.BED_REGIONS_FLAG)) {
      builder.bedRegionsFile((File) mFlags.getValue(CommonFlags.BED_REGIONS_FLAG));
    }
    if (mFlags.isSet(SAMPLE)) {
      builder.sampleName((String) mFlags.getValue(SAMPLE));
    }
    builder.useAllRecords(mFlags.isSet(ALL_RECORDS));
    builder.squashPloidy(mFlags.isSet(SQUASH_PLOIDY));
    builder.rtgStats(mFlags.isSet(RTG_STATS));
    builder.outputBaselineTp(mFlags.isSet(BASELINE_TP));
    builder.outputSlopeFiles(mFlags.isSet(SLOPE_FILES));
    builder.numberThreads(CommonFlags.parseThreads((Integer) mFlags.getValue(CommonFlags.THREADS_FLAG)));
    return builder.create();
  }
}
