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
package com.rtg.reader;

import java.io.File;
import java.io.IOException;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.LongRange;

/**
 * This stores data in memory as it appears on disk. Allowing for direct load (no decompress/re-compress)
 */
public class CompressedMemorySequencesReader2 extends AbstractSequencesReader {

  private final IndexFile mIndexFile;

  private PrereadNamesInterface mNames;
  private PrereadNamesInterface mNameSuffixes;
  private final DataInMemory mData;
  private final long mNumberSequences;
  final File mDirectory;
  private final LongRange mRegion;
  private final long mStart;
  private final long mEnd;

  /**
   * Alternative to other one with similar name
   * @param directory directory containing SDF
   * @param indexFile index file
   * @param loadNames should we load names
   * @param loadFullNames should we load full names
   * @param region region to restrict to
   * @throws IOException IO exception occurs
   */
  CompressedMemorySequencesReader2(File directory, IndexFile indexFile, boolean loadNames, boolean loadFullNames, LongRange region) throws IOException {
    mIndexFile = indexFile;
    mRegion = SequencesReaderFactory.resolveRange(indexFile, region);
    mStart = mRegion.getStart();
    mEnd = mRegion.getEnd();
    mNumberSequences = mEnd - mStart;
    mData = DataInMemory.loadDelayQuality(directory, indexFile, DataFileIndex.loadSequenceDataFileIndex(indexFile.dataIndexVersion(), directory), mStart, mEnd);
    if (mNumberSequences > Integer.MAX_VALUE) {
      throw new IllegalArgumentException("Too many sequences in region: " + region + ", maximum is: " + Integer.MAX_VALUE);
    }
    mDirectory = directory;
    if (loadNames && mIndexFile.hasNames()) {
      loadNames();
    }
    loadNameSuffixes(loadFullNames, mIndexFile.hasSequenceNameSuffixes());
  }

  //copy constructor
  CompressedMemorySequencesReader2(CompressedMemorySequencesReader2 other) {
    mIndexFile = other.mIndexFile;
    mNames = other.mNames;
    mNameSuffixes = other.mNameSuffixes;
    mData = other.mData.copy();
    mNumberSequences = other.mNumberSequences;
    mDirectory = other.mDirectory;
    mRegion = other.mRegion;
    mStart = other.mStart;
    mEnd = other.mEnd;
  }

  /**
   * Load the names if they haven't already been loaded.
   * @throws IOException if an I/O related error occurs
   */
  private void loadNames() throws IOException {
    mNames = new PrereadNames(mDirectory, mRegion, false);
    if (mIndexFile.getVersion() >= IndexFile.SEPARATE_CHECKSUM_VERSION && mRegion.getStart() == 0 && mRegion.getEnd() == mIndexFile.getNumberSequences()) {
      if (mNames.calcChecksum() != mIndexFile.getNameChecksum()) {
        throw new CorruptSdfException("Sequence names failed checksum - SDF may be corrupt: \"" + mDirectory + "\"");
      } else {
        Diagnostic.developerLog("Sequence names passed checksum");
      }
    }
  }

  private void loadNameSuffixes(boolean attemptLoad, boolean suffixExists) throws IOException {
    mNameSuffixes = attemptLoad && suffixExists ? new PrereadNames(mDirectory, mRegion, true) : new EmptyStringPrereadNames(mEnd - mStart);
    if (attemptLoad && suffixExists) {
      if (mRegion.getStart() == 0 && mRegion.getEnd() == mIndexFile.getNumberSequences()) {
        if (mNameSuffixes.calcChecksum() != mIndexFile.getNameSuffixChecksum()) {
          throw new CorruptSdfException("Sequence name suffixes failed checksum - SDF may be corrupt: \"" + mDirectory + "\"");
        } else {
          Diagnostic.developerLog("Sequence name suffixes passed checksum");
        }
      }
    }
  }

  @Override
  public IndexFile index() {
    return mIndexFile;
  }

  @Override
  public long numberSequences() {
    return mNumberSequences;
  }

  @Override
  public File path() {
    return mDirectory;
  }


  // Direct access methods

  @Override
  public int length(long sequenceIndex) {
    return mData.length((int) sequenceIndex);
  }

  @Override
  public String name(long sequenceIndex) {
    return mNames.name(sequenceIndex);
  }

  @Override
  public String nameSuffix(long sequenceIndex) {
    return mNameSuffixes.name(sequenceIndex);
  }

  @Override
  public int read(long sequenceIndex, byte[] dataOut) throws IllegalArgumentException, IOException {
    if (sequenceIndex >= mNumberSequences) {
      throw new IllegalArgumentException("Invalid sequence index: " + sequenceIndex + ", maximum is: " + mNumberSequences);
    }
    return mData.readSequence((int) sequenceIndex, dataOut, 0, Integer.MAX_VALUE);
  }

  @Override
  public int read(long sequenceIndex, byte[] dataOut, int start, int length) throws IllegalArgumentException, IOException {
    if (sequenceIndex >= mNumberSequences) {
      throw new IllegalArgumentException("Invalid sequence index: " + sequenceIndex + ", maximum is: " + mNumberSequences);
    }
    return mData.readSequence((int) sequenceIndex, dataOut, start, length);
  }

  @Override
  public int readQuality(long sequenceIndex, byte[] dest) throws IllegalArgumentException, IOException {
    if (sequenceIndex >= mNumberSequences) {
      throw new IllegalArgumentException("Invalid sequence index: " + sequenceIndex + ", maximum is: " + mNumberSequences);
    }
    return mData.readQuality((int) sequenceIndex, dest, 0, Integer.MAX_VALUE);
  }

  @Override
  public int readQuality(long sequenceIndex, byte[] dest, int start, int length) throws IllegalArgumentException, IOException {
    if (sequenceIndex >= mNumberSequences) {
      throw new IllegalArgumentException("Invalid sequence index: " + sequenceIndex + ", maximum is: " + mNumberSequences);
    }
    return mData.readQuality((int) sequenceIndex, dest, start, length);
  }

  @Override
  public void close() {
  }

  @Override
  public PrereadNamesInterface names() {
    return mNames;
  }

  @Override
  public long lengthBetween(long start, long end) {
    if (start > mNumberSequences) {
      throw new IllegalArgumentException("Invalid sequence index: " + start + ", maximum is: " + mNumberSequences);
    }
    if (end > mNumberSequences) {
      throw new IllegalArgumentException("Invalid sequence index: " + end + ", maximum is: " + mNumberSequences);
    }
    return mData.lengthBetween((int) start, (int) end);
  }

  @Override
  public int[] sequenceLengths(long start, long end) {
    if (start >= mNumberSequences) {
      throw new IllegalArgumentException("Invalid sequence index: " + start + ", maximum is: " + mNumberSequences);
    }
    if (end > mNumberSequences) {
      throw new IllegalArgumentException("Invalid sequence index: " + end + ", maximum is: " + mNumberSequences);
    }
    return mData.sequenceLengths((int) start, (int) end);
  }

  @Override
  public SequencesReader copy() {
    return new CompressedMemorySequencesReader2(this);
  }


  void initQuality() throws IOException {
    mData.initQuality();
  }
}
