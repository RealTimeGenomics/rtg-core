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
package com.rtg.sam;

import static com.rtg.launcher.CommonFlags.FILE;
import static com.rtg.launcher.CommonFlags.OUTPUT_FLAG;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Collection;
import java.util.regex.Pattern;

import com.rtg.calibrate.Calibrator;
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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;

/**
 * Convert SAM to BAM
 */
public class Sam2Bam extends AbstractCli {

  @Override
  public String moduleName() {
    return "sam2bam";
  }

  @Override
  public String description() {
    return "convert SAM file to BAM file and create index";
  }

  @Override
  protected void initFlags() {
    initFlags(mFlags);
  }

  // Separated out from above method for ease of testing
  static void initFlags(CFlags flags) {
    CommonFlagCategories.setCategories(flags);
    flags.setDescription("Produces an indexed BAM file from coordinate-sorted SAM/BAM files.");
    CommonFlags.initForce(flags);
    final Flag<File> inFlag = flags.registerRequired(File.class, FILE, "SAM/BAM format files containing mapped reads");
    inFlag.setCategory(INPUT_OUTPUT);
    inFlag.setMinCount(1);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    //flags.registerOptional('t', TEMPLATE_FLAG, File.class, CommonFlags.SDF, "SDF of the reference genome the reads have been mapped against (required for CRAM input)").setCategory(INPUT_OUTPUT);
    flags.registerRequired('o', OUTPUT_FLAG, File.class, FILE, "name for output BAM file").setCategory(INPUT_OUTPUT);

    flags.setValidator(VALIDATOR);
  }

  private static final Validator VALIDATOR = new Validator() {
    @Override
    public boolean isValid(final CFlags flags) {
      if (!CommonFlags.checkFileList(flags, null, null, Integer.MAX_VALUE)) {
        return false;
      }
      if (!CommonFlags.validateOutputFile(flags, getBamOutputFile((File) flags.getValue(OUTPUT_FLAG)))) {
        return false;
      }
      for (final Object o : flags.getAnonymousValues(0)) {
        if (!CommonFlags.validateInputFile(flags, (File) o)) {
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
    convertSamToBamSimple(outFile, indexFile, samFiles);
    if (calibrationFiles.size() > 0 && calibrationFiles.size() == samFiles.size()) {
      final Calibrator c = new Calibrator(Calibrator.getCovariateSet(calibrationFiles.iterator().next()), null);
      for (final File f : calibrationFiles) {
        c.accumulate(f);
      }
      c.writeToFile(new File(outFile.getParent(), outFile.getName() + CommonFlags.RECALIBRATE_EXTENSION));
    }
  }

  /**
   * Convert and merge <code>samFiles</code> into one BAM file <code>outFile</code>, no calibration file handling
   * @param outFile output BAM file
   * @param indexFile output BAI file
   * @param samFiles input files
   * @throws IOException if an IO error occurs
   */
  public static void convertSamToBamSimple(File outFile, File indexFile, Collection<File> samFiles) throws IOException {
    try (ThreadedMultifileIterator<SAMRecord> sam = new ThreadedMultifileIterator<>(samFiles, new SingletonPopulatorFactory<>(new SamRecordPopulator()))) {
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
  }

  static File getBamOutputFile(File file) {
    if (file.getName().endsWith(".bam")) {
      return file;
    }
    return new File(file.getPath() + ".bam");
  }

  /**
   * Generate a BAM output file.
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

