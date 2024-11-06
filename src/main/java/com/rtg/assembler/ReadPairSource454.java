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
package com.rtg.assembler;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.rtg.assembler.graph.Graph;
import com.rtg.mode.DNA;
import com.rtg.reader.ReaderUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.intervals.LongRange;

/**
 */
public class ReadPairSource454 extends ReadPairSource {

  ReadPairSource454(SequencesReader... readers) {
    super(readers);
  }

  protected int startFlipArm() {
    // 454 is FF so flip only the second arm
    return 1;
  }

  @Override
  GraphAligner aligner(Graph graph, IntegerOrPercentage maxMismatches, GraphTraversions traversions) {
    return new GraphFlowAligner(graph, maxMismatches, traversions);
  }

  @Override
  List<byte[]> nextFragments() throws IOException {
    final List<byte[]> fragments = super.nextFragments();
    if (fragments != null) {
      for (int i = startFlipArm(); i < fragments.size(); ++i) {
        DNA.reverseComplementInPlace(fragments.get(i));
      }
    }
    return fragments;
  }

  static ReadPairSource makeSource(File readDir, LongRange region) throws IOException {
    if (ReaderUtils.isPairedEndDirectory(readDir)) {
      final SequencesReader left = SequencesReaderFactory.createDefaultSequencesReader(ReaderUtils.getLeftEnd(readDir), region);
      final SequencesReader right = SequencesReaderFactory.createDefaultSequencesReader(ReaderUtils.getRightEnd(readDir), region);
      return new ReadPairSource454(left, right);
    } else {
      final SequencesReader reader = SequencesReaderFactory.createDefaultSequencesReader(readDir, region);
      return new ReadPairSource454(reader);
    }
  }
}
