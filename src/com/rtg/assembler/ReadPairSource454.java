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
      for (int i = startFlipArm(); i < fragments.size(); i++) {
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
