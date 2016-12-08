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
package com.rtg.ngs;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.calibrate.CalibratingSamRecordPopulator;
import com.rtg.calibrate.Calibrator;
import com.rtg.calibrate.CovariateEnum;
import com.rtg.mode.DnaUtils;
import com.rtg.ngs.blocking.MapQScoringReadBlocker;
import com.rtg.ngs.tempstage.BinaryTempFileRecord;
import com.rtg.ngs.tempstage.TempRecordReader;
import com.rtg.ngs.tempstage.TempRecordReaderNio;
import com.rtg.reader.CgUtils;
import com.rtg.reader.FastaUtils;
import com.rtg.reader.PrereadNamesInterface;
import com.rtg.reader.SequencesReader;
import com.rtg.sam.ReadGroupUtils;
import com.rtg.sam.SamBamConstants;
import com.rtg.sam.SamUtils;
import com.rtg.util.ProgramState;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.ReferenceRegions;
import com.rtg.util.io.FileUtils;
import com.rtg.variant.sv.ReadGroupStatsCalculator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;

/**
 * Concatenates one or more intermediate SAM files and filters out just the
 * records with the best score. It also adds the <code>HI</code> and
 * <code>NH/IH</code> attributes to each output record.
 *
 * Records with the best score, but with a repeat count greater than the
 * threshold value are output into an unmapped SAM file instead of the main
 * output file.
 *
 * For speed, these methods read and parse the SAM files using large
 * arrays of bytes, rather than using the Java string methods,
 * such as String.split.
 *
 */
public abstract class AbstractSamResultsFilter {

  protected PrereadNamesInterface mNames;

  private boolean mWriteHeader = true;
  private boolean mBamOutput = false;
  protected SAMFileHeader mHeader;
  byte[] mReadBuffer;
  byte[] mQualBuffer;
  final byte[] mFlattenedQualBuffer;
  final byte[] mGQBuffer;
  int mQualStringLength;
  int mGQLength;
  final boolean mCG;
  final boolean mPaired;
  final String mReadGroupId;
  protected final boolean mUnfiltered;
  protected final long mReadIdOffset;
  protected final boolean mLegacyCigars;
  private ReadGroupStatsCalculator mStatsCalculator;
  private MapReportData mMapReportData;

  private boolean mWriteBogusBamHeader;

  /**
   *
   * @param reader1 reader for left side
   * @param reader2 reader for right side (null if single end)
   * @param readGroupId read group id
   * @param cg true if reads are from complete genomics
   * @param readIdOffset the offset to apply to <code>readId</code>s for this run
   * @param unfiltered true if we are running in <code>--allhits</code> (unfiltered) mode
   * @param legacyCigars true if legacy cigar output mode is enabled
   */
  public AbstractSamResultsFilter(SequencesReader reader1, SequencesReader reader2, String readGroupId, boolean cg, long readIdOffset, boolean unfiltered, boolean legacyCigars) {
    final int maxQualLength = getQualLength(reader1, reader2);
    if (maxQualLength > 0) {
      mQualBuffer = new byte[maxQualLength];
      mReadBuffer = new byte[maxQualLength];
      mFlattenedQualBuffer = new byte[maxQualLength];
      if (cg) {
        mGQBuffer = new byte[maxQualLength];
        mGQLength = 0;
      } else {
        mGQBuffer = null;
      }
    } else {
      mQualBuffer = null;
      mFlattenedQualBuffer = null;
      mGQBuffer = null;
      mReadBuffer = null;
    }
    mCG = cg;
    mUnfiltered = unfiltered;
    mReadGroupId = readGroupId;
    mPaired = reader2 != null;
    mReadIdOffset = readIdOffset;
    mLegacyCigars = legacyCigars;
  }

  static int getQualLength(SequencesReader reader1, SequencesReader reader2) {
    int length = 0;
    if (reader1 != null) {
      if (reader2 != null) {
        length = (int) (reader1.maxLength() > reader2.maxLength() ? reader1.maxLength() : reader2.maxLength());
      } else {
        length = (int) reader1.maxLength();
      }
    }
    return length;
  }

  /**
   * Concatenate and filter <code>inputFiles</code>, writing the output
   * into <code>destination</code>.
   * It is the callers responsibility to close <code>destination</code>.
   *
   * @param header the sam file header
   * @param destination output stream
   * @param calibrationDest output stream for calibration
   * @param referenceRegions restrict calibration to these regions
   * @param template sequences reader for the template
   * @param deleteTempFiles delete temporary files after we are done with them
   * @param inputFiles input files
   * @throws IOException if any IO errors occur
   */
  public void filterConcat(SAMFileHeader header, OutputStream destination, OutputStream calibrationDest, ReferenceRegions referenceRegions, SequencesReader template, boolean deleteTempFiles, File... inputFiles) throws IOException {
    try (SAMFileWriter samWriter = getSAMFileWriter(header, destination)) {
      if (mStatsCalculator != null) {
        mStatsCalculator.setupReadGroups(header);
      }

      final CalibratingSamRecordPopulator cp;
      if (template != null && calibrationDest != null) {
        final Calibrator cal = new Calibrator(CovariateEnum.getCovariates(CovariateEnum.DEFAULT_COVARIATES, header), referenceRegions);
        if (referenceRegions != null) {
          cal.setSequenceLengths(Calibrator.getSequenceLengthMap(template, referenceRegions));
        }
        cp = new CalibratingSamRecordPopulator(cal, template, false);
      } else {
        cp = null;
      }

      int inputRecords = 0;
      int outputRecords = 0;
      for (final File currentFile : inputFiles) {
        ProgramState.checkAbort();
        final long t0 = System.nanoTime();
        if (currentFile.length() == 0) {
          continue;
        }
        //System.out.println("Starting to filter file " + current);

        final TempRecordReader.RecordFactory fact = new TempRecordReader.RecordFactory(mPaired, mLegacyCigars, mCG, mUnfiltered);
        try (TempRecordReader recReader = new TempRecordReaderNio(FileUtils.createGzipInputStream(currentFile, false), fact)) {
          BinaryTempFileRecord rec;
          while ((rec = recReader.readRecord()) != null) {
            ++inputRecords;
            final SAMRecord filteredRecord = filterRecord(samWriter, rec, template.names());
            if (filteredRecord != null) {
              if (cp != null) {
                cp.populate(filteredRecord);
              }
              if (mStatsCalculator != null) {
                mStatsCalculator.addRecord(filteredRecord);
              }
              if (mMapReportData != null) {
                mMapReportData.processRead(filteredRecord);
              }
              ++outputRecords;
            }
          } //while
        }
        final long diff = System.nanoTime() - t0;
        final long length = currentFile.length();
        Diagnostic.developerLog("filter concat file=" + currentFile.getAbsolutePath() + " bytes=" + length + " time=" + diffDivided(diff) + "ms" + " bytes/sec=" + Utils.realFormat(bytesPerSecond(diff, length), 2));
        if (deleteTempFiles) {
          if (!currentFile.delete()) {
            Diagnostic.developerLog("can not delete temp file : " + currentFile.getPath());
          }
        }
      } //for
      // timer.stopLog();
      if (cp != null) {
        cp.calibrator().writeToStream(calibrationDest);
      }

      Diagnostic.userLog(getName() + " SAM filter outputs " + outputRecords + "/" + inputRecords + " records");
    }
  }

  private SAMFileWriter getSAMFileWriter(SAMFileHeader header, OutputStream destination) {
    final SAMFileWriter samWriter;
    if (mBamOutput && mWriteBogusBamHeader) {
      final SAMFileHeader basicHeader = new SAMFileHeader();
      basicHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
      samWriter = new SAMFileWriterFactory().makeBAMWriter(basicHeader, true, destination, true, false, true); //terminate flag irrelevant when preBlockCompressed true
    } else if (mBamOutput) {
      samWriter = new SAMFileWriterFactory().makeBAMWriter(header, true, destination, mWriteHeader, false, true); //terminate flag irrelevant when preBlockCompressed true
    } else {
      samWriter = new SAMFileWriterFactory().makeSAMWriter(header, true, destination, mWriteHeader);
    }
    return samWriter;
  }

  protected double bytesPerSecond(final long diff, final long length) {
    return length * 1.0e9 / diff;
  }

  protected long diffDivided(final long diff) {
    return diff / 1000000;
  }

  protected String processReadName(PrereadNamesInterface names, int readId, boolean paired) {
    if (names == null) {
      return Long.toString(readId + mReadIdOffset);
    }
    return SamUtils.samReadName(names.name(readId), paired);
  }

  protected String read(int readId, SequencesReader reader, boolean isReverse) throws IllegalArgumentException, IOException {
    if (reader == null) {
      return "";
    }
    final int len = reader.read(readId, mReadBuffer);
    final StringBuilder res = new StringBuilder();
    for (int i = 0; i < len; ++i) {
      res.append(DnaUtils.getBase(mReadBuffer[i]));
    }
    final String readString;
    if (isReverse) {
      readString = DnaUtils.reverseComplement(res.toString());
    } else {
      readString = res.toString();
    }
    return readString;
  }

  // result = true means the XQ field is non-empty so should be output.
  protected boolean qualField(int readId, SequencesReader reader, boolean reverse, BinaryTempFileRecord rec, boolean first, int readLen) throws IOException {
    if (reader == null || mQualBuffer == null || !reader.hasQualityData()) {
      return false;
    }
    boolean result = true;
    if (mCG && new String(rec.getSuperCigarString()).contains("B")) {
      final int length = reader.readQuality(readId, mQualBuffer);
      assert length == CgUtils.CG_RAW_READ_LENGTH || length == CgUtils.CG2_RAW_READ_LENGTH;
      if (reverse) {
        Utils.reverseInPlace(mQualBuffer);
      }
      final int discard = length - readLen;
      if (discard == 0) {
        System.arraycopy(mQualBuffer, 0, mFlattenedQualBuffer, 0, length);
        result = false; // do NOT add the XQ field, since there is no overlap
      } else {
        final boolean v2 = length == CgUtils.CG2_RAW_READ_LENGTH;
        final int overlapPos = v2 ? CgUtils.CG2_OVERLAP_POSITION : CgUtils.CG_OVERLAP_POSITION;
        final int mainChunkLength = length - (overlapPos + discard);
        final boolean overlapOnLeft = v2 && !reverse || !v2 && first == !reverse;
        if (overlapOnLeft) {
          System.arraycopy(mQualBuffer, 0, mFlattenedQualBuffer, 0, overlapPos);
          System.arraycopy(mQualBuffer, overlapPos, mGQBuffer, 0, discard);
          System.arraycopy(mQualBuffer, overlapPos + discard, mFlattenedQualBuffer, overlapPos, mainChunkLength);
        } else {
          System.arraycopy(mQualBuffer, 0, mFlattenedQualBuffer, 0, mainChunkLength);
          System.arraycopy(mQualBuffer, mainChunkLength, mGQBuffer, 0, discard);
          System.arraycopy(mQualBuffer, length - overlapPos, mFlattenedQualBuffer, mainChunkLength, overlapPos);
        }
      }
      mGQLength = discard;
      mQualStringLength = length - discard;
      for (int i = 0; i < mGQLength; ++i) {
        mGQBuffer[i] += FastaUtils.PHRED_LOWER_LIMIT_CHAR;
      }
    } else {
      final int length = reader.readQuality(readId, mQualBuffer);
      for (int i = 0; i < length; ++i) {
        mFlattenedQualBuffer[i] = mQualBuffer[reverse ? length - 1 - i : i];
      }
      mQualStringLength = length;
      result = false;
    }
    return result;
  }

  /**
   * Function to set the name sets for left and right reads
   * @param names read name set
   *
   */
  public void setReadNames(PrereadNamesInterface names) {
    mNames = names;
  }

  void setHeader(SAMFileHeader header) {
    mHeader = header;
  }

  void setWriteHeader(boolean value) {
    mWriteHeader = value;
  }

  void setWriteBogusBamHeader(boolean value) {
    mWriteBogusBamHeader = value;
  }

  void setBamOutput(boolean value) {
    mBamOutput = value;
  }

  void setStatsCalculator(ReadGroupStatsCalculator.Merger merger) {
    mStatsCalculator = merger == null ? null : merger.createReadGroupStatsCalculator();
  }

  void setMapReportData(MapReportData.Merger merger) {
    mMapReportData = merger == null ? null : merger.createMapReportData();
  }

  protected void addBaseAttributes(SAMRecord dest, BinaryTempFileRecord rec) {
    final boolean properpaired = (rec.getSamFlags() & SamBamConstants.SAM_READ_IS_MAPPED_IN_PROPER_PAIR) != 0;
    dest.setAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE, rec.getAlignmentScore());
    if (mCG) {
      if (rec.getSuperCigarString().length > 0) {
        dest.setAttribute(SamUtils.CG_SUPER_CIGAR, new String(rec.getSuperCigarString()));
      }
      if (rec.getReadDeltaString().length > 0) {
        dest.setAttribute(SamUtils.CG_READ_DELTA, new String(rec.getReadDeltaString()));
      }
    } else {
      dest.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, rec.getNumberMismatches());
    }
    if (mLegacyCigars) {
      if (rec.getMdString().length > 0) {
        dest.setAttribute(SamUtils.ATTRIBUTE_MISMATCH_POSITIONS, new String(rec.getMdString()));
      }
    }
    if (properpaired) {
      dest.setAttribute(SamUtils.ATTRIBUTE_COMBO_SCORE, rec.getComboScore());
    }
    if (mReadGroupId != null) {
      dest.setAttribute(ReadGroupUtils.RG_ATTRIBUTE, mReadGroupId);
    }
  }

  protected void addIhAttributes(SAMRecord record, int readId, boolean hits, int count, MapQScoringReadBlocker blocker) {
    if (hits) {
      record.setAttribute(SamUtils.ATTRIBUTE_IH, count);
    } else {
      record.setAttribute(SamUtils.ATTRIBUTE_IH, 1);
    }
    record.setAttribute(SamUtils.ATTRIBUTE_NH, count);
  }
  protected abstract String getName();

  protected abstract SAMRecord filterRecord(SAMFileWriter samWriter, BinaryTempFileRecord rec, PrereadNamesInterface templateNames) throws IOException;
}
