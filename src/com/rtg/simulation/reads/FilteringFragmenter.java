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

package com.rtg.simulation.reads;

import java.io.IOException;

import com.rtg.reader.SequencesReader;
import com.rtg.simulation.genome.SequenceDistribution;
import com.rtg.util.intervals.ReferenceRegions;

/**
 */
public class FilteringFragmenter extends GenomeFragmenter {
  final ReferenceRegions mRegions;

  FilteringFragmenter(ReferenceRegions regions, long randomSeed, SequenceDistribution[] selectionProb, SequencesReader[] sdfs) throws IOException {
    super(randomSeed, selectionProb, sdfs);
    mRegions = regions;
  }

  @Override
  boolean emitFragment(int fragLength, int seqId, int readerId, String seqName, int fragStart) throws IOException {
    if (mRegions.overlapped(seqName, fragStart, fragStart + fragLength)) {
      return super.emitFragment(fragLength, seqId, readerId, seqName, fragStart);
    }
    return false;
  }
}
