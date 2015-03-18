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

import java.io.IOException;

import com.rtg.reader.SdfReaderWrapper;
import com.rtg.reader.SdfWriterWrapper;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Abstract Output Processor which allows SDF writing
 */
public abstract class AbstractSdfOutputProcessor extends AbstractMapOutputProcessor {

  private static final String ALIGNMENTS_SDF_FILE = "alignments.sdf";
  /** Filename used for unmapped records */
  public static final String UNMAPPED_SDF_FILE = "unmapped.sdf";

  final NgsParams mParams;

  private final boolean mOutputSdf;


  /**
   * Create an abstract SDF output processor
   *
   * @param param ngs params
   * @param stats statistics container
   * @param paired whether this is a paired end read output processor
   * @param outputAlignments whether we are outputting SAM/BAM alignments
   * @throws IOException if an io exception occurs
   */
  public AbstractSdfOutputProcessor(NgsParams param, MapStatistics stats, boolean paired, boolean outputAlignments) throws IOException {
    super(param, stats, paired, outputAlignments);
    mParams = param;
    mOutputSdf = param.outputParams().sdf();
  }

  @Override
  public void finish() throws IOException {
    if (mOutputSdf) {
      Diagnostic.userLog("Starting SDF output");
      final SdfReaderWrapper reader = new SdfReaderWrapper(mSharedResources.firstReaderCopy(), mSharedResources.secondReaderCopy());
      final int maxLength = reader.maxLength();
      final byte[] dataBuffer = new byte[maxLength];
      final byte[] qualityBuffer = new byte[maxLength];
      try (SdfWriterWrapper alignments = new SdfWriterWrapper(mParams.file(ALIGNMENTS_SDF_FILE), reader, false);
           SdfWriterWrapper unmapped = new SdfWriterWrapper(mParams.file(UNMAPPED_SDF_FILE), reader, false)) {
        for (long seq = 0; seq < reader.numberSequences(); seq++) {
          if (!mUnmappedTracker.getStatus((int) seq, ReadStatusTracker.UNMAPPED_FIRST)
            || !mUnmappedTracker.getStatus((int) seq, ReadStatusTracker.UNMAPPED_SECOND)) {
            alignments.writeSequence(reader, seq, dataBuffer, qualityBuffer);
          } else {
            unmapped.writeSequence(reader, seq, dataBuffer, qualityBuffer);
          }
        }
      }
      Diagnostic.userLog("Finished SDF output");
    }
  }
}
