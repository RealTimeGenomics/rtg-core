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
package com.rtg.tabix;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;

import com.rtg.util.io.ByteArrayIOUtils;

import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.BlockCompressedOutputStream;

/**
 * Methods for merging <code>tabix</code> index files
 */
public final class TabixIndexMerge {

  private TabixIndexMerge() { }

  /**
   * Merge indexes for files that will be concatenated.
   * @param output output index file
   * @param files <code>tabix</code> index files
   * @param dataFileSizes file size of corresponding data files
   * @throws IOException if an IO error occurs
   */
  public static void mergeTabixFiles(File output, List<File> files, List<Long> dataFileSizes) throws IOException {
    long pointerAdjust = 0;
    final SequenceIndex[][] indexesSquared = new SequenceIndex[files.size()][];
    final String[][] sequenceNames = new String[files.size()][];
    TabixHeader mergedHeader = null;
    for (int i = 0; i < files.size(); i++) {
      final File tbiFile = files.get(i);
      try (BlockCompressedInputStream bcis = new BlockCompressedInputStream(tbiFile)) {
        final TabixHeader th = TabixHeader.readHeader(bcis);
        sequenceNames[i] = th.getSequenceNamesUnpacked();
        if (mergedHeader != null) {
          mergedHeader = TabixHeader.mergeHeaders(mergedHeader, th);
        } else {
          mergedHeader = th;
        }
        indexesSquared[i] = loadFileIndexes(bcis, th.getNumSequences(), pointerAdjust);
      }
      pointerAdjust += dataFileSizes.get(i);
    }
    final List<SequenceIndex> indexes = collapseIndexes(indexesSquared, sequenceNames);
    TabixIndexer.mergeChunks(indexes);
    try (BlockCompressedOutputStream fos = new BlockCompressedOutputStream(output)) {
      TabixIndexer.writeIndex(indexes, mergedHeader.getOptions(), Arrays.asList(mergedHeader.getSequenceNamesUnpacked()), fos);
    }
  }

  /**
   * Collapses indexes from a two dimensional array into a list
   * @param indexesSquared the indexes to collapse
   * @param sequenceNames the names associated
   * @return a list
   */
  public static List<SequenceIndex> collapseIndexes(SequenceIndex[][] indexesSquared, String[][] sequenceNames) {
    final ArrayList<SequenceIndex> ret = new ArrayList<>();
    assert indexesSquared.length > 0 && sequenceNames.length == indexesSquared.length;
    SequenceIndex[] mergeMaster = indexesSquared[0];
    String[] nameMaster = sequenceNames[0];
    for (int i = 1; i < indexesSquared.length; i++) {
      final String[] mergeName = collapseNames(nameMaster, sequenceNames[i]);
      mergeMaster = collapse(mergeMaster, indexesSquared[i], mergeName, nameMaster, sequenceNames[i]);
      nameMaster = mergeName;
    }
    ret.addAll(Arrays.asList(mergeMaster));
    return ret;
  }

  private static SequenceIndex[] collapse(SequenceIndex[] master, SequenceIndex[] next, String[] jointNames, String[] masterNames, String[] nextNames) {
    assert master.length == masterNames.length && next.length == nextNames.length;
    //it is required that ordering of sequences between two instances are equivalent
    int masterArrayIndex = 0;
    int nextArrayIndex = 0;
    final ArrayList<SequenceIndex> ret = new ArrayList<>();
    for (final String currName : jointNames) {
      SequenceIndex currIndex;
      if (masterArrayIndex < masterNames.length) {
        assert masterNames[masterArrayIndex].equals(currName);
        currIndex = master[masterArrayIndex];
        if (nextArrayIndex < nextNames.length && nextNames[nextArrayIndex].equals(currName)) {
          currIndex.addChunks(next[nextArrayIndex]);
          currIndex.addLinearIndex(next[nextArrayIndex]);
          nextArrayIndex++;
        }
        masterArrayIndex++;
      } else {
        assert nextArrayIndex < nextNames.length && nextNames[nextArrayIndex].equals(currName);
        currIndex = next[nextArrayIndex];
        nextArrayIndex++;
      }
      ret.add(currIndex);
    }
    return ret.toArray(new SequenceIndex[ret.size()]);
  }

  /**
   * Collapse two arrays into one
   * @param master the first array
   * @param next the array to append
   * @return the collapsed array
   */
  public static String[] collapseNames(String[] master, String[] next) {
    final LinkedHashSet<String> ret = new LinkedHashSet<>();
    ret.addAll(Arrays.asList(master));
    ret.addAll(Arrays.asList(next));
    return ret.toArray(new String[ret.size()]);
  }

  /**
   * Load file indexes.
   * @param bcis an input stream
   * @param numberSequences the number of input sequences to read
   * @param pointerAdjust an adjustment to fix chunk positions by...
   * @return an array of sequence indexes
   * @throws IOException if an error occurs reading
   */
  public static SequenceIndex[] loadFileIndexes(InputStream bcis, int numberSequences, long pointerAdjust) throws IOException {
    final byte[] bufIn = new byte[10 * 1024];
    final ByteArrayOutputStream bufHelper = new ByteArrayOutputStream();
    int inLen;
    while ((inLen = bcis.read(bufIn, 0, bufIn.length)) != -1) {
      bufHelper.write(bufIn, 0, inLen);
    }
    final byte[] buf = bufHelper.toByteArray();
    final ArrayList<SequenceIndex> indexes = new ArrayList<>();
    long pos = 0;
    for (int seqNo = 0; seqNo < numberSequences; seqNo++) {
      final SequenceIndex bi = new SequenceIndex();
      //if (!th.getSequenceNamesUnpacked()[seqNo].equals(seqName) && seqName != null) {
      //  indexes.add(bi);
      //  bi = new BamIndex();
      //}
      //seqName = th.getSequenceNamesUnpacked()[seqNo];
      final int numBins = ByteArrayIOUtils.bytesToIntLittleEndian(buf, (int) pos);
      pos += 4;
      for (int i = 0; i < numBins; i++) {
        final int binNo = ByteArrayIOUtils.bytesToIntLittleEndian(buf, (int) pos);
        final int numChunks = ByteArrayIOUtils.bytesToIntLittleEndian(buf, (int) pos + 4);
        pos += 8;
        for (int j = 0; j < numChunks; j++) {
          final long chunkBeg;
          final long chunkEnd;
          if (binNo == TabixIndexer.META_BIN && j == 1) {
            //this refers to the second meta chunk which contains record counts instead of file pointers
            chunkBeg = ByteArrayIOUtils.bytesToLongLittleEndian(buf, (int) pos);
            chunkEnd = ByteArrayIOUtils.bytesToLongLittleEndian(buf, (int) pos + 8);
          } else {
            chunkBeg = fixChunkPosition(ByteArrayIOUtils.bytesToLongLittleEndian(buf, (int) pos), pointerAdjust);
            chunkEnd = fixChunkPosition(ByteArrayIOUtils.bytesToLongLittleEndian(buf, (int) pos + 8), pointerAdjust);
          }
          bi.addChunk(binNo, chunkBeg, chunkEnd);
          pos += 16;
        }
      }
      final int numLinear = ByteArrayIOUtils.bytesToIntLittleEndian(buf, (int) pos);
      pos += 4;
      for (int chunkNo = 0; chunkNo < numLinear; chunkNo++) {
        final long orgChunkOffset = ByteArrayIOUtils.bytesToLongLittleEndian(buf, (int) pos);
        final long chunkOffset = fixChunkPosition(orgChunkOffset, pointerAdjust);
        bi.setLinearIndex(chunkNo, chunkOffset, -1);
        pos += 8;
      }
      indexes.add(bi);
    }
    return indexes.toArray(new SequenceIndex[indexes.size()]);
  }

  private static long fixChunkPosition(long filePointer, long accumulatedFileSize) {
    return filePointer + (accumulatedFileSize << 16L);
  }
}
