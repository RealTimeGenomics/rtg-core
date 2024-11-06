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

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import com.rtg.sam.SamFilterParams.SamFilterParamsBuilder;
import com.rtg.util.Pair;
import com.rtg.variant.DefaultMachineErrorChooser;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantAlignmentRecordPopulator;

/**
 */
public class CircularBufferMultifileSinglePassReaderWindowSyncTest extends CircularBufferMultifileSinglePassReaderWindowTest {

  @Override
  protected Pair<CircularBufferMultifileSinglePassReaderWindow<VariantAlignmentRecord>, RecordIterator<VariantAlignmentRecord>> getCircularBuffer(File[] samFiles, int start, int end) throws IOException {
    final List<File> list = Arrays.asList(samFiles);
    final SamRegionRestriction restriction = new SamRegionRestriction("simulatedSequence2", start, end);
    final VariantAlignmentRecordPopulator pop = new VariantAlignmentRecordPopulator(new DefaultMachineErrorChooser(), 0, "a", "b", "c");
    final RecordIterator<VariantAlignmentRecord> it = CircularBufferMultifileSinglePassReaderWindowTest.defaultIterator(list, new SamFilterParamsBuilder().restriction(restriction).create(), pop);
    final CircularBufferMultifileSinglePassReaderWindow<VariantAlignmentRecord> ssrw =
        new CircularBufferMultifileSinglePassReaderWindowSync<>(
            it,
            pop,
            SamUtils.getUberHeader(list).getSequenceIndex("simulatedSequence2"), restriction.getStart(), Integer.MAX_VALUE
            );
    return new Pair<>(ssrw, it);
  }

}
