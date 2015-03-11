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

import java.io.Closeable;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.List;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.reader.SdfId;
import com.rtg.reader.SequencesReader;
import com.rtg.sam.ReferenceSequenceRequiredException;
import com.rtg.sam.SamFilterParams;
import com.rtg.sam.SamMerger;
import com.rtg.sam.SamReadingContext;
import com.rtg.sam.SamUtils;
import com.rtg.sam.ThreadedMultifileIterator;
import com.rtg.tabix.IndexingStreamCreator;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.intervals.ReferenceRegions;
import com.rtg.util.io.AsynchInputStream;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

/**
 * Does the stand alone re-calibration.
 */
@TestClass("com.rtg.calibrate.RecalibrateCliTest")
public class Recalibrate implements Closeable {

  /**
   * Extension for calibration files
   */
  public static final String EXTENSION = ".calibration";

  private final SequencesReader mTemplate;
  private final SdfId mTemplateSdfId;
  private final ReferenceRegions mRegions;

  /**
   * Constructor
   * @param template reference SDF
   * @param regions bed file regions to restrict calibration
   * @throws IOException if an IO error occurs
   */
  public Recalibrate(SequencesReader template, ReferenceRegions regions) throws IOException {
    mRegions = regions;
    mTemplate = template;
    mTemplateSdfId = mTemplate.getSdfId();
  }

  void doRecalibrate(List<File> samFiles, List<CovariateEnum> covs, boolean force) throws IOException {
    try {
      for (final File f : samFiles) {
        doRecalibrate(f, covs, force);
      }
    } catch (final ReferenceSequenceRequiredException e) {
      throw new NoTalkbackSlimException("Template SDF must be supplied when using legacy cigars");
    }
  }

  private void doRecalibrate(File samFile, List<CovariateEnum> covs, boolean force) throws IOException {
    try (SamReader reader = SamUtils.makeSamReader(new AsynchInputStream(new FileInputStream(samFile)), mTemplate)) {
      SamUtils.checkReferenceGuid(reader.getFileHeader(), mTemplateSdfId);
      final Calibrator c = doRecalibrate(reader, CovariateEnum.getCovariates(covs, reader.getFileHeader()));
      final File calibrationFile = new File(samFile.getParent(), samFile.getName() + EXTENSION);
      if (!force && calibrationFile.exists()) {
        throw new IOException("File: " + calibrationFile + " already exists");
      }
      c.writeToFile(calibrationFile);
    }
  }

  private Calibrator doRecalibrate(SamReader reader, Covariate[] covs) throws IOException {
    final Calibrator c = new Calibrator(covs, mRegions);
    if (mRegions != null) {
      c.setSequenceLengths(Calibrator.getSequenceLengthMap(mTemplate, mRegions));
    }
    final CalibratingSamRecordPopulator p = new CalibratingSamRecordPopulator(c, mTemplate, false);
    try (SAMRecordIterator it = reader.iterator()) {
      while (it.hasNext()) {
        p.populate(it.next());
      }
    }
    return c;
  }

  void doMergeRecalibrate(File output, List<File> samFiles, List<CovariateEnum> covs, int threads, boolean force, boolean compress) throws IOException {
    final SAMFileHeader uberHeader = SamUtils.getUberHeader(mTemplate, samFiles);
    final Covariate[] covariates = CovariateEnum.getCovariates(covs, uberHeader);
    final CalibratingPopulatorFactory rpf = new CalibratingPopulatorFactory(covariates, mRegions, mTemplate);
    final boolean outputBam = output.getName().endsWith(SamUtils.BAM_SUFFIX);
    final File alignmentOutputFile = outputBam ? output : SamUtils.getZippedSamFileName(compress, output);
    final File calibrationFile = new File(alignmentOutputFile.getParent(), alignmentOutputFile.getName() + EXTENSION);
    if (!force && calibrationFile.exists()) {
      throw new IOException("File: " + calibrationFile + " already exists");
    }
    final TabixIndexer.IndexerFactory indexerFactory = outputBam ? null : new TabixIndexer.SamIndexerFactory();
    final SamReadingContext context = new SamReadingContext(samFiles, threads, SamFilterParams.builder().create(), SamUtils.getUberHeader(mTemplate, samFiles), mTemplate);
    try (final ThreadedMultifileIterator<SAMRecord> it = new ThreadedMultifileIterator<>(context, rpf)) {
      final SAMFileHeader header = it.header().clone();
      header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
      SamUtils.addProgramRecord(header);
      SamUtils.updateRunId(header);
      try (final IndexingStreamCreator streamHandler = new IndexingStreamCreator(alignmentOutputFile, System.out, compress, indexerFactory, true /* create index */)) {
        try (final SAMFileWriter writer = SamMerger.getSAMFileWriter(streamHandler, header, outputBam, true /* writeHeader */, true /* terminateBlockedGzip */)) {
          while (it.hasNext()) {
            writer.addAlignment(it.next());
          }
        }
      }
    }

    final Calibrator c = rpf.mergedCalibrator();
    c.writeToFile(calibrationFile);
  }

  @Override
  public void close() throws IOException {
    mTemplate.close();
  }
}
