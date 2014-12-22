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
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.List;

import com.rtg.tabix.SequenceIndex;
import com.rtg.tabix.SequenceIndexContainer;
import com.rtg.tabix.TabixIndexMerge;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.io.ByteArrayIOUtils;
import com.rtg.util.io.IOUtils;

/**
 * Methods for merging BAM index files
 */
public final class BamIndexMerge {

  private BamIndexMerge() { }

  /**
   * Merge indexes for files that will be concatenated.
   * @param output output index file
   * @param files BAM index files
   * @param dataFileSizes file size of corresponding data files
   * @throws IOException if an IO error occurs
   */
  public static void mergeBamIndexFiles(File output, List<File> files, List<Long> dataFileSizes) throws IOException {
    long pointerAdjust = 0;
    final SequenceIndex[][] indexesSquared = new SequenceIndex[files.size()][];
    final String[][] sequenceNames = new String[files.size()][];
    for (int i = 0; i < files.size(); i++) {
      final File baiFile = files.get(i);
      try (FileInputStream is = new FileInputStream(baiFile)) {
        final byte[] smallBuf = new byte[8];
        IOUtils.readFully(is, smallBuf, 0, 8);
        final int numSequences = ByteArrayIOUtils.bytesToIntLittleEndian(smallBuf, 4);
        sequenceNames[i] = new String[numSequences];
        for (int j = 0; j < numSequences; j++) {
          sequenceNames[i][j] = Integer.toString(j);
        }

        indexesSquared[i] = TabixIndexMerge.loadFileIndexes(is, numSequences, pointerAdjust);
      }
      pointerAdjust += dataFileSizes.get(i);
    }

    final List<SequenceIndex> indexes = TabixIndexMerge.collapseIndexes(indexesSquared, sequenceNames);
    final SequenceIndexContainer cont = new SequenceIndexContainer(indexes, 0);
    TabixIndexer.mergeChunks(indexes);
    try (FileOutputStream fos = new FileOutputStream(output)) {
      BamIndexer.writeIndex(cont, fos);
    }
  }
}
