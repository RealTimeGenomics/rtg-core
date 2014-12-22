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
package com.rtg.launcher;

import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Stack;

import com.rtg.reader.ReaderUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.reference.Ploidy;
import com.rtg.reference.ReferenceGenome;
import com.rtg.reference.ReferenceGenome.DefaultFallback;
import com.rtg.reference.ReferenceSequence;
import com.rtg.reference.Sex;
import com.rtg.util.Pair;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.util.StringUtils;
import com.rtg.util.array.ArrayUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.SlimException;

/**
 * Encapsulates a hashing region with clipping and padding with respect to a set of sequences.
 * Start position is exclusive.
 * End position is exclusive.
 *
 */
public class HashingRegion implements Serializable, Comparable<HashingRegion> {

  /** Indicates that a value is missing (usually meaning the range should be extended as far as possible) */
  public static final long MISSING = -1;

  /**
   * Never split workload into a chunk smaller than this.
   * This attempts to prevent multithreading from slowing down on smaller runs
   */
  public static final long DEFAULT_MIN_CHUNK_SIZE = 20000;

  /** Number of regions to create per thread */
  public static final int DEFAULT_THREAD_MULTIPLIER = 10;

  /** A singleton that treats all locations as in range. */
  public static final HashingRegion NONE = new HashingRegion() {
      @Override
      public boolean isInRange(final long sequenceId, final long position) {
        return true;
      }
      @Override
      public String toString() {
        return "[All inclusive]";
      }
    };


  private final long mStartId;
  private final long mEndId;

  private final long mStartClipPosition;
  private final long mEndClipPosition;

  private final long mPaddedStartPosition;
  private final long mPaddedEndPosition;
  private int mChunkId;

  // Used for the NONE singleton
  private HashingRegion() {
    this(MISSING, MISSING, MISSING, MISSING, MISSING, MISSING);
  }

  /**
   * @return The start reference sequence id
   */
  public long getStart() {
    return mStartId;
  }
  /**
   * @return The end reference sequence id
   */
  public long getEnd() {
    return mEndId;
  }

  public void setMapxMetaChunkId(int val) {
    mChunkId = val;
  }

  public int getChunkId() {
    return mChunkId;
  }

  /**
   * Sometimes we need the exclusive full sequence end id even when we have a clip position
   * @return the end id as if the entire last sequence was in the region
   */
  public long getExclusiveEndId() {
    return mEndClipPosition == MISSING ? mEndId : mEndId + 1;
  }

  /**
   * @return the first position in the start template
   */
  public long getStartClipPosition() {
    return mStartClipPosition;
  }

  /**
   * @return the last position in the start template
   */
  public long getEndClipPosition() {
    return mEndClipPosition;
  }

  /**
   * Constructs a region with boundaries within a sequence.
   *
   * @param startId the start reference sequence id
   * @param clipStartPosition the position on the start sequence.
   * @param endId the end reference sequence id
   * @param clipEndPosition the position on the end sequence
   * @param paddedStart position of start after padding is applied
   * @param paddedEnd position of end after padding is applied
   */
  public HashingRegion(final long startId, final long clipStartPosition, final long endId, final long clipEndPosition, final long paddedStart, final long paddedEnd) {
    mStartId = startId;
    mStartClipPosition = clipStartPosition;
    mEndId = endId;
    mEndClipPosition = clipEndPosition;
    mPaddedStartPosition = paddedStart;
    mPaddedEndPosition = paddedEnd;
    if (mEndId == MISSING && mStartId != MISSING || mStartId == MISSING && mEndId != MISSING) {
      throw new IllegalArgumentException("cannot specify wild card at only one side of start/end range");
    }
  }

  /**
   * Constructs a region using only whole sequences
   *
   * @param startId the start reference sequence id
   * @param endId the end reference sequence id
   */
  public HashingRegion(final long startId, final long endId) {
    mStartId = startId;
    mStartClipPosition = MISSING;
    mEndId = endId;
    mEndClipPosition = MISSING;
    mPaddedStartPosition = MISSING;
    mPaddedEndPosition = MISSING;
  }

  public long getStartPaddedPosition() {
    return mPaddedStartPosition;
  }

  public long getEndPaddedPosition() {
    return mPaddedEndPosition;
  }

  /**
   * Returns true if the reference sequence is within the clipping region.
   * @param sequenceId the reference sequence id to check is in range
   * @return true if the sequence is within the clip region.
   */
  public boolean isInRange(final long sequenceId) {
    if (mEndId == MISSING || mStartId == MISSING) {
      return true;
    }
    if (sequenceId < mStartId || sequenceId > mEndId) {
      return false;
    }
    if (sequenceId == mEndId) {
      if (mStartClipPosition == MISSING || mEndClipPosition == MISSING) {
        return false;
      }
    }
    return true;
  }
  /**
   * Returns true if the supplied sequence coordinate is within the clipping region.
   *
   * @param sequenceId the sequence id
   * @param position the position on the supplied sequence
   * @return true if the position is within the clip region.
   */
  public boolean isInRange(final long sequenceId, final long position) {
    if (mEndId == MISSING || mStartId == MISSING) {
      return true;
    }
    if (sequenceId < mStartId || sequenceId > mEndId) {
      return false;
    }
    if (mStartClipPosition == MISSING || mEndClipPosition == MISSING) {
        return sequenceId < mEndId;
    }

    final long tempPos = position < 0 ? 0 : position;
    if (sequenceId == mStartId && tempPos < mStartClipPosition) {
      return false;
    }
    if (sequenceId == mEndId && tempPos >= mEndClipPosition) {
      return false;
    }
    return true;
  }

  /**
   * Returns true if the supplied sequence coordinate is within the padded region.
   *
   * @param sequenceId the sequence id
   * @param position the position on the supplied sequence
   * @return -1 if position precedes padded range, 1 if position exceeds padded range, otherwise 0 for when included in padded range
   */
  public int isInPaddedRange(final long sequenceId, final long position) {
    if (mEndId == MISSING || mStartId == MISSING) {
      return 0;
    }
    if (sequenceId < mStartId) {
      return -1;
    }
    if (sequenceId > mEndId) {
      return 1;
    }
    if (mStartClipPosition == MISSING || mEndClipPosition == MISSING) {
      //endId becomes 0 based exclusive when positions are not being compared in the region
      //we know sequenceId is already <= endId
      return sequenceId != mEndId ? 0 : 1;
    }
    final long tempPos = position < 0 ? 0 : position;
    if (sequenceId == mStartId && tempPos < mPaddedStartPosition) {
      return -1;
    }
    if (sequenceId == mEndId && tempPos >= mPaddedEndPosition) {
      return 1;
    }
    return 0;
  }

  @Override
  public String toString() {
    return "[(" + Long.toString(mStartId) + ":" + Long.toString(mStartClipPosition) + "), (" + Long.toString(mEndId) + ":" + Long.toString(mEndClipPosition) + ")]";
  }

  /**
   * Splits the supplied sequences into a bunch of roughly similar sized chunks.
   * We'll use the following strategy.
   * Calculate the total length we'll process = total
   * ideal chunk size = max(total/number of chunks, minimum chunk size)
   * for each sequence split into ceiling(length/ideal chunk size) evenly sized chunks
   *
   * Chunks will not cross sequences.
   *
   * @param reader the sequences reader to split
   * @param sex the sex of the sample
   * @param startId the sequence id of the first sequence to split
   * @param endId final sequence of to split
   * @param numberChunks ideal number of chunks to split into
   * @param minChunkSize minimum size of a chunk
   * @param padding amount of padding for padded regions
   * @return an array of evenly sized regions spanning all sequences
   * @throws IOException if there is a problem with the sequences reader
   */
  public static HashingRegion[] splitWorkload(SequencesReader reader, Sex sex, final long startId, long endId
      , final int numberChunks, final long minChunkSize, final long padding) throws IOException {
    final ReferenceGenome rg;
    if (sex == null || sex == Sex.EITHER) {
      rg = new ReferenceGenome(reader, sex, DefaultFallback.DIPLOID);
    } else {
      rg = new ReferenceGenome(reader, sex);
      Diagnostic.userLog("Sex-specific Reference Genome:" + StringUtils.LS + rg.toString());
    }
    final HashingRegion[] regions = splitWorkload(reader, rg, startId, endId, numberChunks, minChunkSize, padding);

    final List<HashingRegion> excluded = excludeDuplicateRegions(rg, regions, ReaderUtils.getSequenceNameMap(reader));
    return excluded.toArray(new HashingRegion[excluded.size()]);
  }

  static List<HashingRegion> excludeDuplicateRegions(ReferenceGenome rg, HashingRegion[] regions, Map<String, Long> sequenceNameMap) {
    final Map<Long, List<Pair<Integer, Integer>>> dupMap = duplicateMap(rg, sequenceNameMap);
    final Stack<HashingRegion> stack = new Stack<>();
    for (HashingRegion r : regions) {
      stack.add(0, r);
    }
    final List<HashingRegion> excluded = new ArrayList<>();
    while (stack.size() > 0) {
      final HashingRegion current = stack.pop();
      boolean changed = false;
      for (long sequenceId = current.getStart(); sequenceId <= current.getEnd(); sequenceId++) {
        final List<Pair<Integer, Integer>> duplicates = dupMap.get(sequenceId);
        if (duplicates != null) {
          for (Pair<Integer, Integer> dup : duplicates) {
            final Integer duplicateStart = dup.getA();
            final Integer duplicateEnd = dup.getB();
            if (current.isInRange(sequenceId, duplicateEnd - 1)) {
              stack.push(new HashingRegion(sequenceId,
                  duplicateEnd, current.getEnd(), current.getEndClipPosition(),
                  duplicateEnd, current.getEndPaddedPosition()));
              changed = true;
            } else if (current.isInPaddedRange(sequenceId, duplicateEnd - 1) == 0) {
              stack.push(new HashingRegion(sequenceId, current.getStartClipPosition(), current.getEnd(),
                      current.getEndClipPosition(), duplicateEnd, current.getEndPaddedPosition()));
              changed = true;
            }
            if (current.isInRange(sequenceId, duplicateStart)) {
              stack.push(new HashingRegion(current.getStart(), current.getStartClipPosition(), sequenceId,
                  duplicateStart, current.getStartPaddedPosition(), duplicateStart));
              changed = true;
            } else if (current.isInPaddedRange(sequenceId, duplicateStart) == 0) {
              stack.push(new HashingRegion(current.getStart(), current.getStartClipPosition(), sequenceId,
                      current.getEndClipPosition(), current.getStartPaddedPosition(), duplicateStart));
            }

            if (current.getStart() == current.getEnd() && current.getStartClipPosition() >= duplicateStart && current.getEndClipPosition() <= duplicateEnd) {
              changed = true;
            }
            if (changed) {
              break;
            }
          }
          if (changed) {
            break;
          }
        }
      }
      if (!changed) {
        excluded.add(current);
      }
    }
    return excluded;
  }

  static Map<Long, List<Pair<Integer, Integer>>> duplicateMap(ReferenceGenome rg, Map<String, Long> nameMap) {
    final Map<Long, List<Pair<Integer, Integer>>> result = new HashMap<>();
    for (final ReferenceSequence sequence : rg.sequences()) {
      for (Pair<RegionRestriction, RegionRestriction> dup : sequence.duplicates()) {
        final RegionRestriction side = dup.getB();
        final Long key = nameMap.get(side.getSequenceName());
        if (!result.containsKey(key)) {
          result.put(key, new ArrayList<Pair<Integer, Integer>>());
        }
        result.get(key).add(new Pair<>(side.getStart(), side.getEnd()));
      }
    }

    return result;
  }

  /**
   * Internal split workload for testing
   *
   * @param reader the sequences reader to split
   * @param rg the reference genome info
   * @param startId the sequence id of the first sequence to split
   * @param endId final sequence to split
   * @param numberChunks ideal number of chunks to split into
   * @param minChunkSize minimum size of a chunk
   * @param padding amount of padding for padded regions
   * @return an array of evenly sized regions spanning all sequences
   * @throws IOException if there is a problem with the sequences reader
   */
  protected static HashingRegion[] splitWorkload(SequencesReader reader, ReferenceGenome rg, final long startId, long endId
      , final int numberChunks, final long minChunkSize, final long padding) throws IOException {

    long totalLengths = 0;
    final int[] seqLengths = reader.sequenceLengths(0, reader.numberSequences());
    final ArrayList<Long> sequences = new ArrayList<>();
    for (long id = startId; id < endId; id++) {
      reader.seek(id);
      if (rg == null || rg.sequence(reader.currentName()) != null && rg.sequence(reader.currentName()).ploidy() != Ploidy.NONE) {
        totalLengths += reader.currentLength();
        sequences.add(reader.currentSequenceId());
      }
    }

    final List<HashingRegion> regions = new ArrayList<>();
    long lastId = -1L;
    long splitStart = -1L;
    for (long id : sequences) {
      if (id != lastId + 1 && lastId != -1) {
        final int[] lengths = Arrays.copyOfRange(seqLengths, (int) splitStart, (int) lastId + 1);
        final long partLength = ArrayUtils.sum(lengths);
        final long chunks = Math.max(1, (long) (numberChunks / (totalLengths / (double) partLength)));
        final HashingRegion[] r = splitWorkload(lengths, splitStart, (int) chunks, DEFAULT_MIN_CHUNK_SIZE, padding);
        regions.addAll(Arrays.asList(r));
        splitStart = id;
      }
      if (splitStart == -1L) {
        splitStart = id;
      }
      lastId = id;
    }
    final int[] lengths = Arrays.copyOfRange(seqLengths, (int) splitStart, (int) lastId + 1);
    final long partLength = ArrayUtils.sum(lengths);
    final long chunks = Math.max(1, (long) (numberChunks / (totalLengths / (double) partLength)));
    final HashingRegion[] r = splitWorkload(lengths, splitStart, (int) chunks, minChunkSize, padding);
    regions.addAll(Arrays.asList(r));
    return regions.toArray(new HashingRegion[regions.size()]);

  }
  /**
   * Builds an array of clip regions with evenly sized sequence chunks
   * @param sequenceLengths lengths of the sequences to split
   * @param start the sequence id of the first sequence to split
   * @param numberChunks how many chunks to split into
   * @param minChunkSize minimum size of a chunk
   * @param padding amount of padding for padded regions
   * @return an array of evenly sized regions spanning all sequences
   */
  protected static HashingRegion[] splitWorkload(final int[] sequenceLengths, final long start, final int numberChunks, final long minChunkSize, final long padding) {
    long totalLengths = 0;
    for (final long sequenceLength : sequenceLengths) {
      totalLengths += sequenceLength;
    }
    if (numberChunks > totalLengths) {  // Mostly for unit tests, low total length
      // bizzare case where we have more threads than bases in all sequences combined
      final HashingRegion[] ranges = new HashingRegion[(int) totalLengths];
      int startId = 0;
      long startPosition = 0;
      for (int i = 0; i < totalLengths; i++, startPosition++) {
        if (startPosition > sequenceLengths[startId] - 1) {
          startId++;
          startPosition = 0;
        }
        ranges[i] = new HashingRegion(startId, startPosition, startId, startPosition + 1, Math.max(0, startPosition - padding), Math.min(startPosition + 1 + padding, sequenceLengths[startId]));
      }
      return ranges;

    } else {
      final long individualLength = Math.max(totalLengths / numberChunks, minChunkSize + (totalLengths % minChunkSize) / numberChunks);
      final List<HashingRegion> ranges = new ArrayList<>();
      int startId = 0;
      long startPosition = 0;
      long lengthSoFar = 0;
      //Split the workload for the first n - 1 threads evenly.
      while (lengthSoFar < totalLengths - individualLength * 2 + 1) {
        int endId = startId;
        long endPosition = startPosition;

        long currentLength = individualLength;
        while (currentLength > sequenceLengths[endId] - endPosition) {
          currentLength -= sequenceLengths[endId] - endPosition;
          endId++;
          endPosition = 0;
        }
        endPosition += currentLength;
        ranges.add(new HashingRegion(start + startId, startPosition, start + endId, endPosition, Math.max(0, startPosition - padding), Math.min(endPosition + padding, sequenceLengths[endId])));
        startId = endId;
        startPosition = endPosition;
        lengthSoFar += individualLength;
      }
      // last thread picks up any slack due to rounding
      final long endPosLastCase = sequenceLengths[sequenceLengths.length - 1];
      ranges.add(new HashingRegion(start + startId, startPosition, start + sequenceLengths.length - 1, endPosLastCase, Math.max(0, startPosition - padding), endPosLastCase));
      //System.out.println(Arrays.toString(ranges));
      final HashingRegion[] rangeArray = new HashingRegion[ranges.size()];

      return ranges.toArray(rangeArray);
    }
  }
  @Override
  public boolean equals(final Object obj) {
    if (this == obj) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    if (!this.getClass().equals(obj.getClass())) {
      return false;
    }
    final HashingRegion that = (HashingRegion) obj;
    if (mStartId != that.getStart()) {
      return false;
    }
    if (mStartClipPosition != that.getStartClipPosition()) {
      return false;
    }
    if (mEndId != that.getEnd()) {
      return false;
    }
    if (mEndClipPosition != that.getEndClipPosition()) {
      return false;
    }
    return true;
  }

  @Override
  public int hashCode() {
    // Auto generated by netbeans
    int hash = 3;
    hash = 59 * hash + (int) (mStartId ^ (mStartId >>> 32));
    hash = 59 * hash + (int) (mStartClipPosition ^ (mStartClipPosition >>> 32));
    hash = 59 * hash + (int) (mEndId ^ (mEndId >>> 32));
    hash = 59 * hash + (int) (mEndClipPosition ^ (mEndClipPosition >>> 32));
    return hash;
  }

  @Override
  public int compareTo(final HashingRegion other) {
    if (other == null) {
      return 1;
    }
    final long startTemplate = mStartId - other.getStart();
    if (startTemplate != 0) {
      return startTemplate > 0 ? 1 : -1;
    }
    final long startPosition = mStartClipPosition - other.getStartClipPosition();
    if (startPosition != 0) {
      return startPosition > 0 ? 1 : -1;
    }

    final long endTemplate = mEndId - other.getEnd();
    if (endTemplate != 0) {
      return endTemplate > 0 ? 1 : -1;
    }

    final long endPosition = mEndClipPosition - other.getEndClipPosition();
    if (endPosition != 0) {
      return endPosition > 0 ? 1 : -1;
    }
    return 0;

  }

  /**
   * Find the longest sub-sequence of a region
   * @param reader the location of the sequences
   * @return the length of the longest subsequence
   * @throws SlimException If a sequence is longer than allowed by the size of an integer
   * @throws IOException if an IO exception occurs
   */
  public long longestSubSequence(final SequencesReader reader) throws SlimException, IOException {
    if (this.equals(HashingRegion.NONE)) {
      final long maxLength = reader.maxLength();
      if (maxLength >= Integer.MAX_VALUE) {
        Diagnostic.error(ErrorType.SEQUENCE_TOO_LONG, maxLength + "");
        throw new SlimException();
      }
      return maxLength;
    } else {
      final long start = this.getStart();
      final long end = this.getEndClipPosition() == MISSING ? this.getEnd() : this.getEnd() + 1;
      final int[] lengths = reader.sequenceLengths(start, end);
      if (lengths.length == 1 && this.getStartClipPosition() != MISSING && this.getEndClipPosition() != MISSING) {
        lengths[0] = (int) (mPaddedEndPosition - mPaddedStartPosition);
      } else {
        if (this.getStartClipPosition() != MISSING) {
          lengths[0] -= mPaddedStartPosition;
        }
        if (this.getEndClipPosition() != MISSING) {
          lengths[lengths.length - 1] = (int) mPaddedEndPosition;
        }
      }

      long max = 0;
      for (final long l : lengths) {
        if (l > max) {
          max = l;
        }
      }
      return max;
    }
  }

  /**
   * Get the start position of a given sequence within the region
   * @param sequenceId sequence id for required sequence
   * @param padding padding on regions
   * @return the start position
   */
  public int getReferenceStart(final long sequenceId, final long padding) {
    final int startPos;
    if (this.getStartClipPosition() != MISSING && sequenceId == this.getStart()) {
      startPos = (int) Math.max(0, this.getStartClipPosition() - padding);
    } else {
      startPos = 0;
    }
    return startPos;
  }

  /**
   * Get the end position of a given sequence within the region
   * @param sequenceId sequence id for required sequence
   * @param padding padding on regions
   * @param length the length of the sequence with the supplied sequence id
   * @return the end position
   */
  public int getReferenceEnd(final long sequenceId, final long padding, final int length) {
    final int endPos;
    if (this.getEndClipPosition() != MISSING && sequenceId == this.getEnd() && this.getEndClipPosition() < length) {
      endPos = (int) Math.min(length, this.getEndClipPosition() + padding);
    } else {
      endPos = length;
    }
    return endPos;
  }
}
