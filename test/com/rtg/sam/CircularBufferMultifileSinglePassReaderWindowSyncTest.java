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

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import com.rtg.sam.SamFilterParams.SamFilterParamsBuilder;
import com.rtg.util.Pair;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantAlignmentRecordPopulator;

/**
 */
public class CircularBufferMultifileSinglePassReaderWindowSyncTest extends CircularBufferMultifileSinglePassReaderWindowTest {

  @Override
  protected Pair<CircularBufferMultifileSinglePassReaderWindow<VariantAlignmentRecord>, RecordIterator<VariantAlignmentRecord>> getCircularBuffer(File[] samFiles, int start, int end) throws IOException {
    final List<File> list = Arrays.asList(samFiles);
    final SamRegionRestriction restriction = new SamRegionRestriction("simulatedSequence2", start, end);
    final VariantAlignmentRecordPopulator pop = new VariantAlignmentRecordPopulator("a", "b", "c");
    final RecordIterator<VariantAlignmentRecord> it = CircularBufferMultifileSinglePassReaderWindow.defaultIterator(list, new SamFilterParamsBuilder().restriction(restriction).create(), 4, pop);
    final CircularBufferMultifileSinglePassReaderWindow<VariantAlignmentRecord> ssrw =
        new CircularBufferMultifileSinglePassReaderWindowSync<>(
            it,
            pop,
            SamUtils.getUberHeader(list).getSequenceIndex("simulatedSequence2"), restriction.getStart(), Integer.MAX_VALUE
            );
    return new Pair<>(ssrw, it);
  }

}
