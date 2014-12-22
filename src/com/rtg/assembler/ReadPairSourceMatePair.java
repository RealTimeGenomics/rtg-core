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
package com.rtg.assembler;

import java.io.File;
import java.io.IOException;

import com.rtg.reader.ReaderUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.intervals.LongRange;

/**
 */
public class ReadPairSourceMatePair extends ReadPairSource454 {

  ReadPairSourceMatePair(final SequencesReader... readers) {
    super(readers);
  }

  @Override
  protected int startFlipArm() {
    // Mate pair is RF so need to flip both arms
    return 0;
  }

  static ReadPairSource makeSource(final File readDir, final LongRange region) throws IOException {
    if (ReaderUtils.isPairedEndDirectory(readDir)) {
      final SequencesReader left = SequencesReaderFactory.createDefaultSequencesReader(ReaderUtils.getLeftEnd(readDir), region);
      final SequencesReader right = SequencesReaderFactory.createDefaultSequencesReader(ReaderUtils.getRightEnd(readDir), region);
      return new ReadPairSourceMatePair(left, right);
    } else {
      final SequencesReader reader = SequencesReaderFactory.createDefaultSequencesReader(readDir, region);
      return new ReadPairSourceMatePair(reader);
    }
  }

}
