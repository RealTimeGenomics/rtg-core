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

import static com.rtg.launcher.CommonFlags.PEDIGREE_FLAG;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
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
import com.rtg.relation.PedFileParser;
import com.rtg.sam.SamUtils;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LineWriter;
import com.rtg.variant.bayes.multisample.AbstractMultisampleCli;

import htsjdk.samtools.SAMFileHeader;

/**
 * Main program for <code>chrstats</code> module.
 */
public class ChrStatsCli extends AbstractCli {

  private static final String OUTPUT_PEDIGREE_FLAG = "output-pedigree";
  private static final String SAMPLE_FLAG = "sample";

  private static final String SEX_Z_THRESHOLD_FLAG = "sex-z-threshold";
  private static final String Z_THRESHOLD_FLAG = "z-threshold";
  static final String DESCRIPTION = "check expected chromosome coverage levels from mapping calibration files";

  @Override
  public String moduleName() {
    return "chrstats";
  }

  @Override
  public String description() {
    return DESCRIPTION;
  }

  @Override
  protected void initFlags() {
    mFlags.setDescription("Check expected chromosome coverage levels from mapping calibration files.");
    CommonFlagCategories.setCategories(mFlags);
    AbstractMultisampleCli.addSamFileFlags(mFlags);
    CommonFlags.initReferenceTemplate(mFlags, true);
    mFlags.registerOptional('p', PEDIGREE_FLAG, File.class, CommonFlags.FILE, "genome relationships PED file").setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    mFlags.registerOptional(OUTPUT_PEDIGREE_FLAG, File.class, CommonFlags.FILE, "output best guest of per-sample sex information to PED file").setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    mFlags.registerOptional(MapFlags.SEX_FLAG, Sex.class, "sex", "sex setting that the individual was mapped as (when not using pedigree)", Sex.EITHER).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    mFlags.registerOptional('s', SAMPLE_FLAG, String.class, CommonFlags.STRING, "the name of the sample to check (required when multiple samples present)").setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    mFlags.registerOptional(SEX_Z_THRESHOLD_FLAG, Double.class, CommonFlags.FLOAT, "the z-score deviation threshold for sex chromosome consistency", ChrStats.DEFAULT_MIN_SEX_DEVIATIONS).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    mFlags.registerOptional(Z_THRESHOLD_FLAG, Double.class, CommonFlags.FLOAT, "the z-score deviation threshold for chromosome consistency", ChrStats.DEFAULT_MIN_DEVIATIONS).setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
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
    final Collection<File> inputFiles = new CommandLineFiles(CommonFlags.INPUT_LIST_FLAG, null, CommandLineFiles.EXISTS).getFileList(mFlags);
    final SamCalibrationInputs inputs = new SamCalibrationInputs(inputFiles, true);
    final Collection<File> samFiles = inputs.getSamFiles();
    final Collection<File> calibrationFiles = inputs.getCalibrationFiles();
    if (samFiles.isEmpty()) {
      throw new InvalidParamsException("No SAM files provided for input.");
    }
    if (calibrationFiles.size() != samFiles.size()) {
      throw new InvalidParamsException("Number of calibration files (" + calibrationFiles.size() + ") does not match number of SAM files (" + samFiles.size() + ").");
    }
    Diagnostic.developerLog("Loading calibration information");
    final Calibrator c = Calibrator.initCalibrator(calibrationFiles);
    if (c == null) {
      throw new NoTalkbackSlimException("Could not load calibration information - ensure SAM files have been calibrated"); // Runtime as above checks have ensured there are calibration files
    } else {
      final SequenceParams.SequenceParamsBuilder builder = SequenceParams.builder().directory(genomeFile).useMemReader(false).mode(SequenceMode.UNIDIRECTIONAL);
      try (SequenceParams genomeParams = builder.create()) {
        final ChrStats cc = new ChrStats(genomeParams.reader(), (Double) mFlags.getValue(SEX_Z_THRESHOLD_FLAG), (Double) mFlags.getValue(Z_THRESHOLD_FLAG));
        if (!cc.referenceOk()) {
          throw new NoTalkbackSlimException("The supplied reference does not contain sufficient genome chromosome information. For more information, see the user manual.");
        }

        final SAMFileHeader uberHeader = SamUtils.getUberHeader(genomeParams.reader(), samFiles);
        final Map<String, String> readGroupToSampleId = SamUtils.getReadGroupToSampleId(uberHeader);
        final Map<String, Integer> sequenceLengthMap = c.hasLengths() ? c.getSequenceLengths() : Calibrator.getNonNSequenceLengthMap(genomeParams.reader(), null);
        final CalibratedPerSequenceExpectedCoverage expectedCoverages = new CalibratedPerSequenceExpectedCoverage(c, sequenceLengthMap, readGroupToSampleId, null);
        final Set<String> samples = expectedCoverages.samples();
        if (samples.isEmpty()) {
          throw new NoTalkbackSlimException("No sample information contained in mapping headers");
        }
        if (mFlags.isSet(OUTPUT_PEDIGREE_FLAG)) {
          outputPedigree(cc, expectedCoverages, samples);
        } else if (samples.size() == 1 || mFlags.isSet(SAMPLE_FLAG)) {
          singleSample(cc, expectedCoverages, samples);
        } else {
          multiSample(cc, expectedCoverages, samples);
        }
      }
    }
    return 0;
  }

  private void multiSample(ChrStats cc, CalibratedPerSequenceExpectedCoverage expectedCoverages, Set<String> samples) throws IOException {
    final GenomeRelationships pedigree = mFlags.isSet(PEDIGREE_FLAG) ? GenomeRelationships.loadGenomeRelationships((File) mFlags.getValue(PEDIGREE_FLAG)) : null;
    // Multiple samples
    if (pedigree == null) {
      throw new NoTalkbackSlimException("Multiple samples present. Please specify sample and sex information via --" + PEDIGREE_FLAG);
    }
    cc.chrStatsCheck(expectedCoverages, samples, pedigree, false);
  }

  private void singleSample(ChrStats cc, CalibratedPerSequenceExpectedCoverage expectedCoverages, Set<String> samples) throws IOException {
    // Single sample mode
    final GenomeRelationships pedigree = mFlags.isSet(PEDIGREE_FLAG) ? GenomeRelationships.loadGenomeRelationships((File) mFlags.getValue(PEDIGREE_FLAG)) : null;
    final String sample = mFlags.isSet(SAMPLE_FLAG) ? (String) mFlags.getValue(SAMPLE_FLAG) : samples.iterator().next();
    final Sex sex;
    if (pedigree != null) {
      sex = pedigree.getSex(sample);
    } else {
      sex = (Sex) mFlags.getValue(MapFlags.SEX_FLAG);
    }
    cc.chrStatsCheckAndReport(expectedCoverages, sample, sex);
  }

  private void outputPedigree(ChrStats cc, CalibratedPerSequenceExpectedCoverage expectedCoverages, Set<String> samples) throws IOException {
    final GenomeRelationships pedigree = (mFlags.isSet(PEDIGREE_FLAG)) ? GenomeRelationships.loadGenomeRelationships((File) mFlags.getValue(PEDIGREE_FLAG)) : new GenomeRelationships();
    samples.forEach(pedigree::addGenome);
    cc.chrStatsCheck(expectedCoverages, samples, pedigree, true);
    final File outFile = (File) mFlags.getValue(OUTPUT_PEDIGREE_FLAG);
    try (LineWriter w = new LineWriter(new OutputStreamWriter(FileUtils.createOutputStream(outFile)))) {
      w.write(PedFileParser.toString(pedigree, "source: " + CommandLine.getCommandLine()));
    }
    Diagnostic.info("PED file written to " + outFile + ", please check and add relationships where known.");
  }
}
