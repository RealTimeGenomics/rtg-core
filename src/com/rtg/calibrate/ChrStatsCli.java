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
package com.rtg.calibrate;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Collection;
import java.util.Map;
import java.util.Set;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommandLineFiles;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.SequenceMode;
import com.rtg.ngs.MapFlags;
import com.rtg.reader.SdfUtils;
import com.rtg.reference.Sex;
import com.rtg.relation.GenomeRelationships;
import com.rtg.sam.SamUtils;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.variant.bayes.multisample.MultisampleCli;

import htsjdk.samtools.SAMFileHeader;

/**
 * Main program for <code>chrstats</code> module.
 */
public class ChrStatsCli extends AbstractCli {

  private static final String PEDIGREE_FLAG = "pedigree";
  private static final String SAMPLE_FLAG = "sample";

  private static final String SEX_Z_THRESHOLD_FLAG = "sex-z-threshold";
  private static final String Z_THRESHOLD_FLAG = "z-threshold";

  @Override
  public String moduleName() {
    return "chrstats";
  }

  @Override
  protected void initFlags() {
    mFlags.registerExtendedHelp();
    mFlags.setDescription("Check expected chromosome coverage levels from mapping calibration files.");
    CommonFlagCategories.setCategories(mFlags);
    MultisampleCli.addSamFileFlags(mFlags);
    mFlags.registerRequired('t', CommonFlags.TEMPLATE_FLAG, File.class, CommonFlags.SDF, "SDF containing reference genome").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerOptional('p', PEDIGREE_FLAG, File.class, "file", "genome relationships PED file").setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    mFlags.registerOptional(MapFlags.SEX_FLAG, Sex.class, "sex", "sex setting that the individual was mapped as (when not using pedigree)", Sex.EITHER).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    mFlags.registerOptional('s', SAMPLE_FLAG, String.class, "STRING", "the name of the sample to check (required when multiple samples present)").setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    mFlags.registerOptional(SEX_Z_THRESHOLD_FLAG, Double.class, "FLOAT", "the z-score deviation threshold for sex chromosome consistency", ChrStats.DEFAULT_MIN_SEX_DEVIATIONS).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    mFlags.registerOptional(Z_THRESHOLD_FLAG, Double.class, "FLOAT", "the z-score deviation threshold for chromosome consistency", ChrStats.DEFAULT_MIN_DEVIATIONS).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
//    final Flag anon = new AnonymousFlag("files to process", File.class, "file");
//    anon.setMaxCount(Integer.MAX_VALUE);
//    anon.setMinCount(1);
//    anon.setCategory(CommonFlagCategories.INPUT_OUTPUT);
//    mFlags.register(anon);
    mFlags.setValidator(new CcValidator());
  }

  private static class CcValidator implements Validator {
    @Override
    public boolean isValid(CFlags flags) {
      if (!flags.checkNand(MapFlags.SEX_FLAG, PEDIGREE_FLAG)) {
        return false;
      }
      if ((Double) flags.getValue(Z_THRESHOLD_FLAG) <= 0) {
        flags.setParseMessage("Z-score threshold must be positive");
        return false;
      }
      if ((Double) flags.getValue(SEX_Z_THRESHOLD_FLAG) <= 0) {
        flags.setParseMessage("Sex-z-score threshold must be positive");
        return false;
      }
      return true;
    }
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    final File genomeFile = (File) mFlags.getValue(CommonFlags.TEMPLATE_FLAG);
    SdfUtils.validateHasNames(genomeFile);
    final SequenceParams genomeParams = SequenceParams.builder().directory(genomeFile).useMemReader(false).mode(SequenceMode.UNIDIRECTIONAL).create();
    final Collection<File> inputFiles = new CommandLineFiles(CommonFlags.INPUT_LIST_FLAG, null, CommandLineFiles.EXISTS).getFileList(mFlags);
    final SamCalibrationInputs inputs = new SamCalibrationInputs(inputFiles, true);
    final Collection<File> samFiles = inputs.getSamFiles();
    final Collection<File> calibrationFiles = inputs.getCalibrationFiles();
    if (samFiles.size() == 0) {
      throw new InvalidParamsException("No SAM files provided for input.");
    }
    if (calibrationFiles.size() != samFiles.size()) {
      throw new InvalidParamsException("Number of calibration files (" + calibrationFiles.size() + ") does not match number of SAM files (" + samFiles.size() + ").");
    }
    Diagnostic.developerLog("Loading calibration information");
    final Calibrator c = Calibrator.initCalibrator(calibrationFiles);
    final SAMFileHeader uberHeader = SamUtils.getUberHeader(samFiles, false, null);
    if (c == null) {
      throw new RuntimeException("Could not load calibration information"); // Runtime as above checks have ensured there are calibration files
    } else {
      final ChrStats cc = new ChrStats(genomeParams.reader(), (Double) mFlags.getValue(SEX_Z_THRESHOLD_FLAG), (Double) mFlags.getValue(Z_THRESHOLD_FLAG));
      if (!cc.referenceOk()) {
        throw new NoTalkbackSlimException("The supplied reference does not contain genome chromosome information. For more information, see the user manual.");
      }

      final GenomeRelationships pedigree = mFlags.isSet(PEDIGREE_FLAG) ? GenomeRelationships.loadGenomeRelationships((File) mFlags.getValue(PEDIGREE_FLAG)) : null;

      final Map<String, String> readGroupToSampleId = SamUtils.getReadGroupToSampleId(uberHeader);
      final Map<String, Integer> sequenceLengthMap = c.hasLengths() ? c.getSequenceLengths() : Calibrator.getSequenceLengthMap(genomeParams.reader(), (RegionRestriction) null);
      final CalibratedPerSequenceExpectedCoverage expectedCoverages = new CalibratedPerSequenceExpectedCoverage(c, sequenceLengthMap, readGroupToSampleId, null);
      final Set<String> samples = expectedCoverages.samples();
      if (samples.size() == 0) {
        throw new NoTalkbackSlimException("No sample information contained in mapping headers");
      } else if (samples.size() == 1 || mFlags.isSet(SAMPLE_FLAG))  {
        // Single sample mode
        final String sample = mFlags.isSet(SAMPLE_FLAG) ? (String) mFlags.getValue(SAMPLE_FLAG) : samples.iterator().next();
        final Sex sex;
        if (pedigree != null) {
          sex = pedigree.getSex(sample);
        } else {
          sex = (Sex) mFlags.getValue(MapFlags.SEX_FLAG);
        }
        cc.chrStatsCheckAndReport(expectedCoverages, sample, sex);
      } else {
        // Multiple samples
        if (pedigree == null) {
          throw new NoTalkbackSlimException("Multiple samples present. Please specify sample sex information via --" + PEDIGREE_FLAG);
        }
        cc.chrStatsCheck(expectedCoverages, samples, GenomeRelationships.loadGenomeRelationships((File) mFlags.getValue(PEDIGREE_FLAG)));
      }
    }
    return 0;
  }
}
