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
package com.rtg.pairedend;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.Properties;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.ngs.SharedResources;
import com.rtg.reader.CgUtils;
import com.rtg.reader.PrereadType;
import com.rtg.reader.SequencesReader;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.ReferenceRegions;
import com.rtg.util.machine.MachineOrientation;


/**
 * @param <T> type of hit info subclass.
 */
@TestClass(value = "com.rtg.pairedend.SlidingWindowCollectorTest")
public abstract class AbstractSlidingWindowCollector<T extends AbstractHitInfo<T>> {

  /**
   * This should be greater than the approximate number of reads starting for each
   * reference position, to ensure the average linked list is short.
   *
   * For 30 x mapping of 35 nt long reads mapping to human reference in a single mapping run,
   * we only expect a new read starting about every 1.1 ref bases, so a factor of 30 here
   * is ample with a capital A. Memory consumption due to this is small, so there isn't a
   * downside.
   */
  private static final int HASHTABLE_FACTOR = 30;

  //this works out to roughly an extra 10 per 1 million reads on the full human reference
  private static final int READS_PP_FACTOR = 31375;
  private static final int READ_OVERLOAD_BASE_LIMIT = 10;

  //see bug #1476 for consequences of this on larger datasets
  private int mMaxHitsPerPosition;

  final int mWindowSize;
  private final int mMinFragmentLength;
  private final int mMaxFragmentLength;

  protected SequencesReader mLeftReader;
  protected SequencesReader mRightReader;

  /**
   * A hash table of opposite-end HitInfo objects.
   * Each entry is a linked list, sorted by increasing template start order.
   */
  private final AbstractHitInfo<T>[] mReadsLookup;
  private final AbstractHitInfo<T>[] mReadsLookupReverse;
  final int mReadsLookupMask;

  static final int DONT_KNOW_YET = Integer.MIN_VALUE;
  int mCurrentReferencePosition = DONT_KNOW_YET;

  final ArrayList<T>[] mReadsWindow; // the sliding window, contains unmated hits
  final int[] mReadsWindowInUse;

  final int[] mRightPairCounts = new int[5000];
  int mLeftOverloadCount = 0;
  int mRightOverloadCount = 0;

  // statistics counters
  private long mDuplicateCount = 0;
  private long mHitCount = 0;
  private long mMaxHitsExceededCount = 0;
  private final long mGenomeLength;
  private int mReferenceCount = 0;
  protected long mReferenceId = 0;
  private long mReferenceLengthTotal = 0; // sum of length of all references
  private final MachineOrientation mMachineOrientation;
  private final int mReadOverloadLimit;

  AbstractSlidingWindowCollector(int maxFragmentLength, int minFragmentLength, MachineOrientation pairOrientation, SharedResources sharedResources, ReferenceRegions bedRegions) {
    if (maxFragmentLength < minFragmentLength) {
      throw new IllegalArgumentException("Maximum window size too small: " + maxFragmentLength);
    }
    if (minFragmentLength < 0) {
      throw new IllegalArgumentException("Minimum window size too small: " + minFragmentLength);
    }

//    final int readMillions = (int) (sharedResources.firstReaderCopy().numberSequences() / 1000000);
    mGenomeLength = bedRegions != null ? genomeLength(bedRegions) : genomeLength(sharedResources.templateReaderCopy());
    if (GlobalFlags.isSet(CoreGlobalFlags.SLIDING_WINDOW_MAX_HITS_PER_POS_FLAG)) {
      setMaxHitsPerPosition(GlobalFlags.getIntegerValue(CoreGlobalFlags.SLIDING_WINDOW_MAX_HITS_PER_POS_FLAG));
    } else {
      final long numReads = sharedResources.firstReaderCopy().numberSequences();
      setMaxHitsPerPosition(1000 + calculateExtraMaxHitsPerPosition(mGenomeLength, numReads));
    }

    Diagnostic.developerLog("Setting max hits per position to: " + getMaxHitsPerPosition());
    mMinFragmentLength = minFragmentLength;
    mMaxFragmentLength = maxFragmentLength;
    mMachineOrientation = pairOrientation;

    mLeftReader = sharedResources.firstReaderCopy();
    mRightReader = sharedResources.secondReaderCopy();

    //margin needs to be added to allow for reads coming from GappedOutput in slightly wrong order, 1 RL is our best guess
    final int uncertaintyMargin = (int) Math.max(64L, Math.max(mLeftReader.maxLength(), mRightReader.maxLength()));
    mWindowSize = maxFragmentLength + (int) Math.max(64L, Math.max(mLeftReader.maxLength(), mRightReader.maxLength())) + uncertaintyMargin;
    if (GlobalFlags.isSet(CoreGlobalFlags.SLIDING_WINDOW_MAX_HITS_PER_READ_FLAG)) {
      mReadOverloadLimit = GlobalFlags.getIntegerValue(CoreGlobalFlags.SLIDING_WINDOW_MAX_HITS_PER_READ_FLAG);
    } else {
      mReadOverloadLimit = READ_OVERLOAD_BASE_LIMIT + (mWindowSize / 100); //1%
    }
    Diagnostic.developerLog("Setting max read hits to: " + mReadOverloadLimit);

    @SuppressWarnings("unchecked")
    final ArrayList<T>[] readsWindow = (ArrayList<T>[]) new ArrayList<?>[mWindowSize];
    mReadsWindow = readsWindow;
    mReadsWindowInUse = new int[mWindowSize];
    for (int i = 0; i < mWindowSize; i++) {
      mReadsWindow[i] = new ArrayList<>();
    }

    final int hashTableSize = nextPowerOfTwo(mWindowSize * HASHTABLE_FACTOR);
    mReadsLookupMask = hashTableSize - 1;
    @SuppressWarnings("unchecked")
    final AbstractHitInfo<T>[] readsLookup = (AbstractHitInfo<T>[]) new AbstractHitInfo<?>[hashTableSize];
    mReadsLookup = readsLookup;
    @SuppressWarnings("unchecked")
    final AbstractHitInfo<T>[] readsLookupReverse = (AbstractHitInfo<T>[]) new AbstractHitInfo<?>[hashTableSize];
    mReadsLookupReverse = readsLookupReverse;
  }

  static int calculateExtraMaxHitsPerPosition(double genomeLength, double numReads) {
    return (int) (numReads / genomeLength * READS_PP_FACTOR);
  }

  static long genomeLength(ReferenceRegions regions) {
    long totalLength = 0;
    for (Map.Entry<String, Integer> entry : regions.coveredLengths().entrySet()) {
      totalLength += entry.getValue();
    }
    return totalLength;
  }

  static long genomeLength(SequencesReader reader) {
    return reader.totalLength();
  }


  void setReadLookup(int i, T hit) {
    mReadsLookup[i] = hit;
  }

  protected static int nextPowerOfTwo(int x0) {
    int x = x0;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return x + 1;
  }

  // Calculate an index into <code>mReadsLookup</code> hash table.
  int readHash(long readId, boolean first) {
    return ((int) readId & mReadsLookupMask) ^ (first ? 1 : 0);
  }

  private T getHitInfo(int i) {
    final T ret;
    if (mReadsWindowInUse[i] == getMaxHitsPerPosition()) {
      mMaxHitsExceededCount++;
      if (mMaxHitsExceededCount < 5) {
        Diagnostic.userLog("Max hits per position exceeded at template: " + mReferenceId + " templateStart: " + (mReadsWindow[i].size() > 0 ? "" + mReadsWindow[i].get(0).mTemplateStart : "unknown"));
      }
      removeHits(mReadsWindow[i], mReadsWindowInUse[i]);
      mReadsWindowInUse[i] = -1; // Blacklist the position until it slides off
    }
    if (mReadsWindowInUse[i] == -1) {
      ret = null;
    } else {
      if (mReadsWindowInUse[i] < mReadsWindow[i].size()) {
        ret = mReadsWindow[i].get(mReadsWindowInUse[i]);
      } else {
        ret = createHitInfo();
        mReadsWindow[i].add(ret);
      }
      mReadsWindowInUse[i]++;
    }
    return ret;
  }

  private void returnHitInfo(int i) {
    mReadsWindowInUse[i]--;
  }


  abstract T createHitInfo();

  int windowPosition(final int i) {
    final int mod = i % mWindowSize;
    if (mod < 0) {
      return mod + mWindowSize;
    } else {
      return mod;
    }
  }

  /**
   * Record a match from one end of a pair. For a given template
   * sequence, input positions must occur in approximate monotonically
   * increasing <code>templateStart</code> coordinate.  Subsequent
   * values of <code>templateStart</code> must be greater than or
   * equal to 64 less than the maximum seen so far for the current
   * template.
   *
   * @param first true if this match report is for the read end that was first in sequencing
   * @param reverseComplement true if the match is on the reverse complement of the template
   * @param readId the 0-based read identifier
   * @param templateStart the 0-base template start position on the forward strand
   * @exception IllegalStateException if <code>nextTemplateId</code> has not been called
   * or the sequence is not monotonic.
   * @throws IOException if an IO error occurs
   */
  public void match(final boolean first, final boolean reverseComplement, int readId, final int templateStart) throws IOException {
    if (mCurrentReferencePosition == DONT_KNOW_YET) {
      mCurrentReferencePosition = templateStart;  // first template position we've seen
    }
    assert templateStart >= mCurrentReferencePosition - (mWindowSize - mMaxFragmentLength) : "Out of order template start: " + templateStart + " < " + (mCurrentReferencePosition - (mWindowSize - mMaxFragmentLength));

    //System.err.println("MATCH: " + first + " : " + readId + " : " + templateStart);
    // if new template start find prior pairs and write out
    if (templateStart > mCurrentReferencePosition) {
      // write out pairs
      flushToPosition(templateStart);
      // buffer indexes updated in flushToPosition
    }

    // add hit to current position in window
    final int windPos = windowPosition(templateStart);
    final T hit = getHitInfo(windPos);
    if (hit != null) {
      hit.setValues(first, reverseComplement, readId, templateStart);

      // add hit to lookup map
      final int hash = readHash(readId, first);

      @SuppressWarnings("unchecked")
      T list = (T) mReadsLookupReverse[hash];
      if (list == null) {
        mReadsLookup[hash] = hit;
        mReadsLookupReverse[hash] = hit;
      } else if (templateStart > list.templateStart()) {
        list.insertHit(hit);
        mReadsLookupReverse[hash] = hit;
      } else {
        boolean same;
        while (!(same = same(list, hit)) && list.prev() != null && list.prev().templateStart() >= templateStart) {
          list = list.prev();
        }
        if (same) {
          returnHitInfo(windPos); // Return it to the pool
          mDuplicateCount++;
        } else {
          if (list.prev() != null) {
            list.prev().insertHit(hit);
            if (hit.next() == null) {
              mReadsLookupReverse[hash] = hit;
            }
          } else {
            hit.setNext(list);
            list.setPrev(hit);
            mReadsLookup[hash] = hit;
          }
        }
      }
    }
    mHitCount++;
  }

  private boolean same(T first, T second) {
    return first.mTemplateStart == second.mTemplateStart
        && first.mReadId == second.mReadId
        && first.mReverseComplement == second.mReverseComplement
        && first.mFirst == second.mFirst;
  }

  /**
   * Step to the specified template identifier.  The ids should be monotonically
   * increasing. <code>Long.MAX_VALUE</code> should be used to indicate no more sequences.
   *
   * @param templateId a <code>long</code> value
   * @throws IOException if an I/O Error occurs
   */
  public void nextTemplateId(final long templateId) throws IOException {
    assert integrity();
    if (mCurrentReferencePosition != DONT_KNOW_YET) {
      mReferenceLengthTotal += mCurrentReferencePosition + 1; // add one for first position at index 0
      // flush remaining pairs from buffer
      // and reset buffer
      flushToPosition(mCurrentReferencePosition + mWindowSize);
      mCurrentReferencePosition = DONT_KNOW_YET;
    }

    mReferenceCount++;
    mReferenceId = templateId;

    // pass next template call along
    writerNextTemplateId(templateId);

    //assert mReadsLookup.size() == 0;
    Arrays.fill(mReadsLookup, null);
    Arrays.fill(mReadsLookupReverse, null);
  }

  void findNewMates(final ArrayList<T> hits, int size) throws IOException {

    nextHit:
    for (int i = 0; i < size; i++) {
      final T hit = hits.get(i);
      //      boolean mated = hit.getLeftReads() != null;

      // find right mate
      // if there is a right mate
      // send pair result to writer
      // link right to left
      final long readId = hit.readId();
      @SuppressWarnings("unchecked")
      T thisSide = (T) mReadsLookup[readHash(readId, hit.first())];
      int thisSideCount = 0;
      while (thisSide != null) {
        if (thisSide.readId() == hit.readId() && thisSide.first() == hit.first()) {
          thisSideCount++;
        }
        if (thisSideCount > mReadOverloadLimit) {
          mLeftOverloadCount++;
          continue nextHit;
        }
        thisSide = thisSide.next();
      }

      final int hash = readHash(readId, !hit.first());
      int pairCount = 0;
      @SuppressWarnings("unchecked")
      T mate = (T) mReadsLookup[hash];
      while (mate != null) {

        if (mate.readId() == readId         // the list may contain other read ids too, ignore them
            && mate.first() != hit.first()  // process only the real mates
            && orientationCorrect(hit, mate)) {

          final int mateReadLength;
          final int hitReadLength;
          //Unfortunately these read lengths are approximations of the alignment length, so the thresholding isn't based on the ultimate template length :(
          if (mLeftReader.getPrereadType() == PrereadType.CG) {
            // we adjust because the most common CG alignment size along the template depends on gap/overlap structure
            final int expectedCgReadLength = mLeftReader.maxLength() == CgUtils.CG_RAW_READ_LENGTH ? CgUtils.CG_EXPECTED_LENGTH : CgUtils.CG2_EXPECTED_LENGTH;
            mateReadLength = expectedCgReadLength;
            hitReadLength = expectedCgReadLength;
          } else {
            mateReadLength = mate.first() ? mLeftReader.length(readId) : mRightReader.length(readId);
            hitReadLength = hit.first() ? mLeftReader.length(readId) : mRightReader.length(readId);
          }
          final int fragmentLength = InsertHelper.calculateFragmentLength(mate.templateStart(), mateReadLength, hit.templateStart(), hitReadLength);
          assert mMaxFragmentLength == mWindowSize - (mWindowSize - mMaxFragmentLength);
          if (fragmentLength >= mMinFragmentLength && fragmentLength <= mMaxFragmentLength) {  // only process if length is within specified fragment bounds
            assert hit.isPair(mate);
            pairCount++;
            if (pairCount <= mReadOverloadLimit) {
              if (!checkPair(hit, mate)) {
                break;
              }
            } else {
              mRightOverloadCount++;
              break;
            }
          }
        }
        mate = mate.next();
      }
      if (pairCount < mRightPairCounts.length) {
        mRightPairCounts[pairCount]++;
      } else {
        mRightPairCounts[mRightPairCounts.length - 1]++;
      }
    }
  }

  /**
   * @param hit the hit.
   * @param mate its mate.
   * @return true iff the orientation is corrected for mated Illumina
   *
   */
  //TODO get correct for other machine types
  private boolean orientationCorrect(final T hit, final T mate) {
    if (hit.mFirst) {
      return mMachineOrientation.orientationOkay(hit.mTemplateStart, hit.mReverseComplement, mate.mTemplateStart, mate.mReverseComplement);
    } else {
      return mMachineOrientation.orientationOkay(mate.mTemplateStart, mate.mReverseComplement, hit.mTemplateStart, hit.mReverseComplement);
    }
  }

  /*
   * Size threshold to apply trimToSize. Making this higher will use more
   * memory but perhaps run slightly faster.  Reduction below 10 is probably
   * not a good idea.
   */
  private static final int ARRAY_SIZE_TRIM_LIMIT = 100;

  /**
   * Clear all the hits at a position in the window that is about to slide out.
   * Clears the slot in the window along with any other hits
   * in the hash table that have a position that will also be cleared.
   * @param hits the hit list
   * @param size the number of currently used elements
   */
  void clearHits(final ArrayList<T> hits, int size) {
    final int count = hits.size();
    for (int i = 0; i < size; i++) {
      final T hit = hits.get(i);
      final long readId = hit.readId();
      final int hash = readHash(readId, hit.first());
      @SuppressWarnings("unchecked")
      T head = (T) mReadsLookup[hash];
      final int templateStart = hit.templateStart();
      while (head != null && head.templateStart() <= templateStart) {
        head = head.next();
      }
      mReadsLookup[hash] = head;
      if (head != null) {
        head.setPrev(null);
        if (head.next() == null) {
          mReadsLookupReverse[hash] = head;
        }
      } else {
        mReadsLookupReverse[hash] = null;
      }
    }

    // Chop down memory used by the list
    if (count > ARRAY_SIZE_TRIM_LIMIT) {
      // Ideally I want cut the list to ARRAY_SIZE_TRIM_LIMIT size, but Java
      // doesn't provide a nice way to do so.  Instead we end up with a
      // capacity of 0 and another array allocation when an add is attempted
      hits.trimToSize();
    }
  }

  /**
   * Clear all the hits at an active position within the window. Clears the slot
   * in the window, along with their entries in the hash table (and no others).
   */
  private void removeHits(final ArrayList<T> hits, int size) {
    final int count = hits.size();
    for (int i = 0; i < size; i++) {
      final T hit = hits.get(i);
      final long readId = hit.readId();
      final int hash = readHash(readId, hit.first());
      @SuppressWarnings("unchecked")
      T head = (T) mReadsLookup[hash];
      T previous = null;
      final int templateStart = hit.templateStart();
      while (head != null && head.templateStart() < templateStart) {
        previous = head;
        head = head.next();
      }
      while (head != null && head.templateStart() == templateStart) {
        head = head.next();
      }
      if (previous == null) {
        mReadsLookup[hash] = head;
        if (head != null) {
          head.setPrev(null);
        } else {
          mReadsLookupReverse[hash] = null;
        }
      } else {
        previous.setNext(head);
        if (head != null) {
          head.setPrev(previous);
        } else {
          mReadsLookupReverse[hash] = previous;
        }
      }
    }
    //    System.err.println(hits.size() + " " + mReadsLookup.size());
    //hits.clear();


    // Chop down memory used by the list
    if (count > ARRAY_SIZE_TRIM_LIMIT) {
      // Ideally I want cut the list to ARRAY_SIZE_TRIM_LIMIT size, but Java
      // doesn't provide a nice way to do so.  Instead we end up with a
      // capacity of 0 and another array allocation when an add is attempted
      hits.trimToSize();
    }
  }

  /**
   * Returns internal counts about the hits that have passed through the sliding window.
   *
   * @return a <code>Properties</code> object containing counts
   */
  public Properties getStatistics() {
    final Properties stats = new Properties();

    stats.setProperty("window_size", Integer.toString(mWindowSize));

    stats.setProperty("templates", Integer.toString(mReferenceCount));
    stats.setProperty("total_hits", Long.toString(mHitCount));
    stats.setProperty("duplicate_hits", Long.toString(mDuplicateCount));
    stats.setProperty("max_pos_hits_exceeded", Long.toString(mMaxHitsExceededCount));
    stats.setProperty("max_pos_hits_threshold", Integer.toString(getMaxHitsPerPosition()));
    stats.setProperty("max_read_hits_left_exceeded", Long.toString(mLeftOverloadCount));
    stats.setProperty("max_read_hits_right_exceeded", Long.toString(mRightOverloadCount));
    stats.setProperty("max_read_hits_threshold", Integer.toString(mReadOverloadLimit));
    if (mGenomeLength > 0) {
      stats.setProperty("max_hits_exceed_percentage", String.format("%.2f", mMaxHitsExceededCount * 100.0 / mGenomeLength));
    }
    stats.setProperty("template_lengths_total", Long.toString(mReferenceLengthTotal));

    for (int i = 0; i < mRightPairCounts.length; i++) {
      if (mRightPairCounts[i] != 0) {
        Diagnostic.developerLog("Reads with " + i + " potential right side candidates: " + mRightPairCounts[i] + " effective alignments: " + (i * 2 * mRightPairCounts[i]));
      }
    }
    return stats;
  }

  /**
   * Performs an integrity check of the internal data structures.  Needs assertions enabled to function.
   *
   * @return always true
   */
  public boolean integrity() {
    assert mReadsLookupMask == mReadsLookup.length - 1;
    // check contents of arrays/hashmaps are consistent
    for (int i = 0; i < mReadsWindow.length; i++) {
      final ArrayList<T> hits = mReadsWindow[i];
      assert hits != null;
      for (int j = 0; j < mReadsWindowInUse[i]; j++) {
        final T hit = hits.get(j);
        assert hit != null;
        final int hash = readHash(hit.readId(), hit.first());
        @SuppressWarnings("unchecked")
        T others = (T) mReadsLookup[hash]; //.get(hit.readId() | (hit.first() ? (1L << 62) : 0L));

        assert others != null;
        boolean seen = false;
        while (others != null) {
          if (others == hit) {
            seen = true;
          }
          others = others.next();
        }
        assert seen;
      }
    }

    for (final AbstractHitInfo<T> hit0 : mReadsLookup) {
      @SuppressWarnings("unchecked")
      T hit = (T) hit0;
      T prevHit = null;
      while (hit != null) {
        final ArrayList<T> others = mReadsWindow[windowPosition(hit.templateStart())];
        assert others != null;
        boolean seen = false;
        for (final T other : others) {
          if (other == hit) {
            seen = true;
          }
        }
        assert seen;
        // check that they are in ascending templateStart order
        assert prevHit == null || prevHit.templateStart() <= hit.templateStart();
        prevHit = hit;
        hit = hit.next();
      }
    }

    return true;
  }

  /**
   * Compute the scores for this pair.
   * @param hit the hit
   * @param mate the mate
   * @return false if hit can never pass the thresholds (with any mate), otherwise true.
   * @throws IOException if an IOException occurs
   */
  abstract boolean checkPair(T hit, T mate) throws IOException;

  abstract void flushToPosition(final int newStart) throws IOException;

  abstract void writerNextTemplateId(long templateId) throws IOException;

  /**
   * @return Maximum number of hits at a given position in the sliding window collector
   */
  protected final int getMaxHitsPerPosition() {
    return mMaxHitsPerPosition;
  }

  void setMaxHitsPerPosition(int maxHitsPerPosition) {
    mMaxHitsPerPosition = maxHitsPerPosition;
  }
}
