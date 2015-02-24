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
package com.rtg.sam;

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Collection;
import java.util.regex.Pattern;

import com.rtg.calibrate.Calibrator;
import com.rtg.calibrate.Recalibrate;
import com.rtg.calibrate.SamCalibrationInputs;
import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommandLineFiles;
import com.rtg.launcher.CommonFlags;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.util.SingletonPopulatorFactory;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.variant.SamRecordPopulator;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;

/**
 * Convert SAM to BAM
 */
public class Sam2Bam extends AbstractCli {

  private static final String MODULE_NAME = "sam2bam";
  //private static final String INDEX_FLAG = "index";
  private static final String OUTPUT_FLAG = "output";

  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  @Override
  protected void initFlags() {
    initFlags(mFlags);
  }

  // Separated out from above method for ease of testing
  static void initFlags(CFlags flags) {
    flags.registerExtendedHelp();
    CommonFlagCategories.setCategories(flags);
    flags.setDescription("Produces an indexed BAM file from coordinate-sorted SAM/BAM files.");
    final Flag inFlag = flags.registerRequired(File.class, "file", "SAM/BAM format files containing mapped reads");
    inFlag.setCategory(INPUT_OUTPUT);
    inFlag.setMinCount(1);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    flags.registerRequired('o', OUTPUT_FLAG, File.class, "file", "name for output BAM file. Use '-' to write to standard output").setCategory(INPUT_OUTPUT);
    //flags.registerOptional('i', INDEX_FLAG, File.class, "file", "Output filename for index").setCategory(INPUT_OUTPUT);

    flags.setValidator(VALIDATOR);
  }

  private static final Validator VALIDATOR = new Validator() {
    @Override
    public boolean isValid(final CFlags flags) {
      if (!CommonFlags.checkFileList(flags, null, null, Integer.MAX_VALUE)) {
        return false;
      }
      final File outFile = (File) flags.getValue(OUTPUT_FLAG);
      if (outFile.isDirectory()) {
        flags.setParseMessage(outFile.getPath() + " is a directory, must be a file");
        return false;
      }
      for (final Object o : flags.getAnonymousValues(0)) {
        final File f = (File) o;
        if (!f.exists()) {
          flags.setParseMessage("\"" + f.getPath() + "\" does not exist");
          return false;
        }
        if (f.isDirectory()) {
          flags.setParseMessage("\"" + f.getPath() + "\" is a directory");
          return false;
        }
      }
      return true;
    }
  };

  /**
   * Convert and merge <code>files</code> into one BAM file <code>outFile</code>
   * @param outFile output BAM file
   * @param indexFile output BAI file
   * @param files input files
   * @throws IOException if an IO error occurs
   */
  public static void convertSamToBam(final File outFile, final File indexFile, File... files) throws IOException {
    convertSamToBam(outFile, indexFile, Arrays.asList(files));
  }

  /**
   * Convert and merge <code>files</code> into one BAM file <code>outFile</code>
   * @param outFile output BAM file
   * @param indexFile output BAI file
   * @param files input files
   * @throws IOException if an IO error occurs
   */
  public static void convertSamToBam(final File outFile, final File indexFile, Collection<File> files) throws IOException {
    final SamCalibrationInputs inputs = new SamCalibrationInputs(files, true);
    final Collection<File> samFiles = inputs.getSamFiles();
    final Collection<File> calibrationFiles = inputs.getCalibrationFiles();
    if (calibrationFiles.size() > 0 && calibrationFiles.size() != samFiles.size()) {
      Diagnostic.warning("Number of calibration files does not match number of SAM files, will not merge calibration files.");
    }
    try (RecordIterator<SAMRecord> sam = new ThreadedMultifileIterator<>(samFiles, new SingletonPopulatorFactory<>(new SamRecordPopulator()))) {
      final SAMFileHeader header = sam.header().clone();
      header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
      SamUtils.addProgramRecord(header);
      SamUtils.updateRunId(header);

      final SAMFileWriterFactory fact = new SAMFileWriterFactory();
      try (SAMFileWriter writer = fact.makeBAMWriter(header, true, outFile)) {
        while (sam.hasNext()) {
          final SAMRecord record = sam.next();
          try {
            writer.addAlignment(record);
          } catch (final IllegalArgumentException iae) {
            throw new NoTalkbackSlimException(iae.getMessage().replaceAll(Pattern.quote("SAMFileWriterImpl.addAlignment for "), ""));
          }
        }
      }
    }
    try {
      BamIndexer.saveBamIndex(outFile, indexFile);
    } catch (final UnindexableDataException e) {
      Diagnostic.warning("Cannot create BAM index: " + e.getMessage());
    }
    if (calibrationFiles.size() > 0 && calibrationFiles.size() == samFiles.size()) {
      final Calibrator c = new Calibrator(Calibrator.getCovariateSet(calibrationFiles.iterator().next()), null);
      for (final File f : calibrationFiles) {
        c.accumulate(f);
      }
      c.writeToFile(new File(outFile.getParent(), outFile.getName() + Recalibrate.EXTENSION));
    }
  }

  static File getBamOutputFile(File file) {
    if (file.getName().endsWith(".bam")) {
      return file;
    }
    return new File(file.getPath() + ".bam");
  }

  /**
   * Generate genomes with specified coverage.
   *
   * @param outStr stream to deliver output to
   * @param errStr stream to print errors to
   * @return 0 on success 1 on failure
   */
  @Override
  protected int mainExec(final OutputStream outStr, final PrintStream errStr) throws IOException {
    final Collection<File> inputFiles = new CommandLineFiles(null, null, CommandLineFiles.EXISTS, CommandLineFiles.NOT_DIRECTORY).getFileList(mFlags);
    final File outFile = getBamOutputFile((File) mFlags.getValue(OUTPUT_FLAG));
    convertSamToBam(outFile, BamIndexer.indexFileName(outFile), inputFiles);
    return 0;
  }
}

