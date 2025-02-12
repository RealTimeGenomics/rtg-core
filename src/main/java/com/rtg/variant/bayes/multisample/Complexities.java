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
package com.rtg.variant.bayes.multisample;

import java.util.ArrayList;
import java.util.Deque;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.util.SynchronizedLinkedList;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.bayes.multisample.ComplexRegion.RegionType;
import com.rtg.variant.bayes.multisample.population.AlleleCounts;
import com.rtg.variant.bayes.multisample.population.SiteSpecificPriors;

/**
 * Contains a single chunks complex region information.
 */
public class Complexities extends IntegralAbstract implements Iterable<ComplexRegion> {

  private static final boolean DUMP_COMPLEXITIES = GlobalFlags.getBooleanValue(CoreGlobalFlags.DUMP_COMPLEX_TRIGGER_SIGNALS);

  private static final String MER_FINDER_IMPL = GlobalFlags.getStringValue(CoreGlobalFlags.COMPLEX_REGION_SIMPLE_REPEAT_IMPL);

  private static final boolean USE_CX_SSP = true; //Boolean.parseBoolean(System.getProperty("com.rtg.complexities.ssp", "true"));

  private final RepeatMeasurer mRepeatMeasurer;
  private final int mInterestingSeparation;
  private final int mHyperComplexLength;

  private ComplexRegion mStartDangling;
  private ComplexRegion mEndDangling;

  private int mStartOfChunk;
  private int mEndOfChunk;

  private boolean mStartFixed;
  private boolean mEndFixed;

  private final String mReferenceName;

  private final SynchronizedLinkedList<ComplexRegion> mRegions;
  private final List<Variant> mOriginalCalls;


  /**
   * Assumed that any complex regions are contained within boundaries set here. Once violated the user should adjust the boundaries
   * @param chunk the variant calls for the chunk, all calls should be within bounds specified by start and end parameters
   * @param referenceName name of the reference sequence
   * @param start start of chunk (0 based inclusive)
   * @param end end of chunk (0 based exclusive)
   * @param interestingSeparation how far apart calls must be to not be determined part of same complex region
   * @param hyperComplexLength length beyond which complex regions are considered hyper complex
   * @param refNts reference sequence bytes
   * @param simpleRepeatHandling should interesting regions extend freely over simple repeats
   * @param ssp optional provider of site specific priors
   */
  public Complexities(List<Variant> chunk, String referenceName, int start, int end, int interestingSeparation, int hyperComplexLength, byte[] refNts, boolean simpleRepeatHandling, SiteSpecificPriors ssp) {
    mRepeatMeasurer = makeRepeatMeasurer(simpleRepeatHandling, refNts);
    mReferenceName = referenceName;
    mStartOfChunk = start;
    mEndOfChunk = end;
    mInterestingSeparation = interestingSeparation;
    mHyperComplexLength = hyperComplexLength;
    mOriginalCalls = chunk;

    final List<Variant> allInteresting;
    if (ssp != null && USE_CX_SSP) {
      final List<Variant> ssps = getSspVariants(ssp, referenceName, start, end, chunk);
      allInteresting = OutputUtils.merge(chunk, ssps);
    } else {
      allInteresting = chunk;
    }
    /*
    if (allInteresting.size() != chunk.size()) {
      System.err.println(chunk.toString());
      System.err.println("------------------");
      System.err.println(allInteresting.toString());
    }
     */
    // Initialisation of members should occur before these calls.
    final Deque<ComplexRegion> temp = getComplexRegions(allInteresting);
    mRegions = filter(temp);
    mStartDangling = computeStartDangle(temp);
    mEndDangling = computeEndDangle(temp);

    //System.err.println("REGIONS: " + mRegions);
  }

  private static RepeatMeasurer makeRepeatMeasurer(boolean simpleRepeatHandling, byte[] refNts) {
    final RepeatMeasurer measurer;
    if (!simpleRepeatHandling) {
      measurer = new RepeatAntiMeasurer();
    } else {
      switch (MER_FINDER_IMPL) {
        case "none":
          measurer = new RepeatAntiMeasurer();
          break;
        case "default": // Label used by GlobalFlags
        case "single-nmer":
          measurer = new SingleNMerRepeatMeasurer(refNts);
          break;
        default:
          throw new RuntimeException("Invalid repeat measurer: " + MER_FINDER_IMPL);
      }
    }
    return measurer;
  }

  // Make variants out of SSPs for complex call triggering
  private static List<Variant> getSspVariants(SiteSpecificPriors ssp, String referenceName, int start, int end, List<Variant> chunkVariants) {
    final List<Variant> ssps = new ArrayList<>();
    Variant overflow = null;
    final Iterator<Variant> it = chunkVariants.iterator();
    for (final AlleleCounts ac : ssp.getCounts(referenceName, start, end)) {
      while (it.hasNext() && (overflow == null || overflow.getLocus().getEnd() < ac.position())) {
        final Variant next = it.next();
        if (next.isOverflow()) {
          overflow = next;
          break;
        }
      }
      final int sspEnd = ac.position() + ac.refLength();
      if (overflow != null) {
        final int originalStart = overflow.getLocus().getStart();
        final int originalEnd = overflow.getLocus().getEnd();
        if ((ac.position() >= originalStart && ac.position() <= originalEnd) || (sspEnd >= originalStart && sspEnd <= originalEnd)) {
          // ignore if the SSP overlaps an extreme coverage region
          continue;
        }
      }
      // Truncate the variant trigger so that current chunk boundaries are not exceeded.
      final Variant v = new Variant(new VariantLocus(referenceName, Math.max(start, ac.position()), Math.min(sspEnd, end)));
      if (ac.isComplex()) {
        // TODO, when we evaluate indel calling along with EvidenceIndelFactory.COMPLEX_REGION_INDEL_EXTENSION
        // v.setIndel(Math.abs(ac.refLength - ac.maxLength()));
        v.setIndel(0);
      }
      v.setInteresting();
      ssps.add(v);
    }
    return ssps;
  }

  // Main routine to convert input variant triggers into a list of complex regions
  final Deque<ComplexRegion> getComplexRegions(List<? extends Variant> chunk) {
    final Deque<ComplexRegion> regions = new SynchronizedLinkedList<>(new LinkedList<ComplexRegion>());

    // These vars store state of current "in-construction" region
    int firstInteresting = -1;
    int firstIndelLength = -1;
    int lastInteresting = -1;
    int lastIndelLength = -1;
    boolean forceComplex = false;

    for (final Variant c : chunk) {
      if (c.isOverflow()) {
        // Want to write an overflow region, but need to output any current complex region first
        addRegion(regions, firstInteresting, lastInteresting, forceComplex, firstIndelLength, lastIndelLength);
        // Reset state of "in-construction" region
        firstInteresting = -1;
        firstIndelLength = -1;
        lastInteresting = -1;
        lastIndelLength = -1;
        forceComplex = false;
        addOverflowRegion(regions, c.getLocus());
      } else if (c.isInteresting()) {
        if (DUMP_COMPLEXITIES) {
          final String desc = c.isIndel() ? "indel" : c.isSoftClip() ? "soft-clip" : "generic";
          Diagnostic.developerLog(String.format("COMPLEX_CHUNK\t%s\t%d\t%d\t%s\t%d", c.getLocus().getSequenceName(), c.getLocus().getStart(), c.getLocus().getEnd(), desc, c.getIndelLength()));
        }
        assert firstInteresting > -1 || (lastInteresting == -1 && !forceComplex);
        assert firstInteresting == -1 || (firstInteresting >= mStartOfChunk && lastInteresting >= firstInteresting && lastInteresting <= mEndOfChunk);
        final int outputPosition = c.getLocus().getStart();
        final int indelLength = !c.isSoftClip() ? c.getIndelLength() : 0;

        // See whether the current variant is var enough away that we need to close out the "in-construction" region
        if (lastInteresting > 0) {
          if (!joinInteresting(lastInteresting - 1, outputPosition, lastIndelLength, indelLength)) {
            // Can't joint current variant to in-construction complex region, so finish the complex region and add it
            addRegion(regions, firstInteresting, lastInteresting, forceComplex, firstIndelLength, lastIndelLength);
            // Reset state of "in-construction" region
            firstInteresting = -1;
            firstIndelLength = -1;
            lastInteresting = -1;
            lastIndelLength = -1;
            forceComplex = false;
          }
        }

        // Add current variant info into "in-construction" state
        if (firstInteresting == -1) {
          firstInteresting = outputPosition;
          firstIndelLength = indelLength;
        }
        final int lastInterestingShadow = lastInteresting + lastIndelLength; // Where did previous indel cast to
        lastInteresting = Math.max(lastInteresting, c.getLocus().getEnd());
        lastIndelLength = Math.max(lastInterestingShadow - lastInteresting, indelLength);
        forceComplex |= c.isForceComplex();
        //System.err.println("getComplexRegions: " + outputPosition + " : " + c.getLocus().getEnd() + " : " + firstInteresting + " : " + lastInteresting);
      }
    }
    addRegion(regions, firstInteresting, lastInteresting, forceComplex, firstIndelLength, lastIndelLength);
    return regions;
  }

  // Strip out any regions corresponding to a single SNP call (since we've already called it, no complex calling is needed)
  final SynchronizedLinkedList<ComplexRegion> filter(Deque<ComplexRegion> regions) {
    final SynchronizedLinkedList<ComplexRegion> ret = new SynchronizedLinkedList<>(new LinkedList<ComplexRegion>());
    for (final ComplexRegion complex : regions) {
      if (complex.type() != ComplexRegion.RegionType.INTERESTING) {
        ret.add(complex);
      }
    }
    return ret;
  }

  final boolean joinInteresting(int positionA, int positionB) {
    return joinInteresting(positionA, positionB, 0, 0);
  }

  final boolean joinInteresting(int positionA, int positionB, int indelA, int indelB) {
    final int indelHint = Math.max(indelA, indelB);
    final int posA = positionA + indelA / 2;
    final int posB = positionB - indelB / 2;
    final int repeatTotal = mRepeatMeasurer.measureRepeats(posA, posB, indelHint);
    return posA + mInterestingSeparation + repeatTotal >= posB;
  }


  final ComplexRegion computeStartDangle(Deque<ComplexRegion> regions) {
    return regions.peekFirst();
  }

  final ComplexRegion computeEndDangle(Deque<ComplexRegion> regions) {
    return regions.peekLast();
  }

  final void addRegion(final Deque<ComplexRegion> regions, int firstInteresting, int endInteresting, boolean forceComplex, int firstIndel, int endIndel) {
    if (firstInteresting != -1) {
      assert endInteresting >= firstInteresting : endInteresting + ":" + firstInteresting;
      final ComplexRegion.RegionType regionType;
      if ((endInteresting - firstInteresting) == 1 && !forceComplex) {
        regionType = ComplexRegion.RegionType.INTERESTING;
      } else if (endInteresting - firstInteresting > mHyperComplexLength) {
        regionType = ComplexRegion.RegionType.HYPER;
      } else {
        regionType = ComplexRegion.RegionType.COMPLEX;
      }
      final ComplexRegion com = new ComplexRegion(mReferenceName, firstInteresting, endInteresting, regionType, firstIndel, endIndel);
      regions.add(com);
    }
  }

  final void addOverflowRegion(final Deque<ComplexRegion> regions, final VariantLocus locus) {
    final ComplexRegion com = new ComplexRegion(mReferenceName, locus.getStart(), locus.getEnd(), RegionType.OVERFLOW);
    regions.add(com);
  }

  static ComplexRegion mergeRegions(ComplexRegion regionA, ComplexRegion regionB, int hyperComplexLength) {
    assert regionA.getEnd() <= regionB.getStart() || (regionB.getStart() >= regionA.getStart() && regionB.getStart() <= regionA.getEnd()) : regionA.getEnd() + "-" + regionB.getStart();
    assert regionA.type() == RegionType.COMPLEX || regionA.type() == RegionType.HYPER || regionA.type() == RegionType.OVERFLOW || regionA.type() == RegionType.INTERESTING;
    assert regionB.type() == RegionType.COMPLEX || regionB.type() == RegionType.HYPER || regionB.type() == RegionType.OVERFLOW || regionB.type() == RegionType.INTERESTING;
    if ((regionA.type() == RegionType.OVERFLOW) && (regionA.type() != regionB.type())) {
      // Can't merge if only one side is an overflow region
      return null;
    }
    final int start = regionA.getStart();
    final int end = regionB.getEnd();
    final boolean isHyper = end - start > hyperComplexLength;
    final ComplexRegion.RegionType type;
    if (regionA.type() == RegionType.OVERFLOW) {
      type = RegionType.OVERFLOW;
    } else if (isHyper) {
      type = RegionType.HYPER;
    //} else if (regionA.type() == RegionType.NO_HYPOTHESES) {
    //  type = RegionType.NO_HYPOTHESES;
    } else {
      type = RegionType.COMPLEX;
    }
    return new ComplexRegion(regionA.getSequenceName(), start, end, type, regionA.getIndelStart(), regionB.getIndelEnd());
  }

  /**
   * Replace a dangling complex region at the end of the first chunk with a complete region if possible.
   * @param regionA first list of regions, we may adjust the last complex region
   * @param regionB second list of regions, we may remove the first complex region
   */
  static void fixDangling(Complexities regionA, Complexities regionB) {
    assert regionA != null || regionB != null;
    assert regionA == null || regionB == null || regionA.mInterestingSeparation == regionB.mInterestingSeparation;
    assert regionA == null || regionB == null || regionA.mHyperComplexLength == regionB.mHyperComplexLength;
    assert regionA == null || regionB == null || regionA.mEndOfChunk == regionB.mStartOfChunk;
    if (regionA == null) {
      regionB.mStartFixed = true;
      regionB.mStartDangling = null;
      return;
    }
    if (regionB == null) {
      regionA.mEndFixed = true;
      regionA.mEndDangling = null;
      return;
    }
    regionA.mEndFixed = true;
    regionB.mStartFixed = true;
    final ComplexRegion aLast = regionA.danglingEnd();
    final ComplexRegion bFirst = regionB.danglingStart();

    if (aLast == null || bFirst == null) {
      regionA.mEndDangling = null;
      regionB.mStartDangling = null;
      return;
    }

    if (!regionA.joinInteresting(aLast.getEnd() - 1, bFirst.getStart(), aLast.getIndelEnd(), bFirst.getIndelStart())) {
      regionA.mEndDangling = null;
      regionB.mStartDangling = null;
      return;
    }

    final ComplexRegion nRegion = mergeRegions(aLast, bFirst, regionA.mHyperComplexLength);
    if (nRegion != null) {

      //if they were interesting then they have already been filtered out of the deques
      if (aLast.type() != ComplexRegion.RegionType.INTERESTING) {
        regionA.mRegions.removeLast();
      }
      if (bFirst.type() != ComplexRegion.RegionType.INTERESTING) {
        regionB.mRegions.removeFirst();
      }

      // Detect situation where region spans multiple chunks
      final boolean multichunkSpanner = aLast.getStart() < regionA.startOfChunk() ||  bFirst.getEnd() >= regionB.endOfChunk();

      if (nRegion.type() == RegionType.HYPER || multichunkSpanner) {
        //split in two on boundary
        final ComplexRegion hyperLeft = new ComplexRegion(regionA.mReferenceName, nRegion.getStart(), regionA.endOfChunk(), nRegion.type());
        final ComplexRegion hyperRight = new ComplexRegion(regionA.mReferenceName, regionB.startOfChunk(), nRegion.getEnd(), nRegion.type());
        if (regionA.mStartDangling == aLast) { //evil test for very very long regions
          regionA.mStartDangling = hyperLeft;
        }
        regionA.mEndDangling = hyperLeft;
        regionA.mRegions.addLast(hyperLeft);
        if (regionB.mEndDangling == bFirst) {
          regionB.mEndDangling = hyperRight;
        }
        regionB.mStartDangling = hyperRight;
        regionB.mRegions.addFirst(hyperRight);
      } else {
        // Handle the case where a region spans an entire chunk.
        regionA.mRegions.addLast(nRegion);
        regionB.mRegions.addFirst(nRegion);
        regionA.mEndDangling = null;
        regionB.mStartDangling = null;
        if (regionA.mRegions.getLast().getEnd() >= regionA.mEndOfChunk) {
          //copy over the calls which are in the dangling part of region B into region A, required when the complex caller fails to make a call and prints originals
          for (final Variant v : regionB.mOriginalCalls) {
            if (v.getLocus().getStart() >= regionA.mEndOfChunk && v.getLocus().getStart() < regionA.mRegions.getLast().getEnd()) {
              regionA.mOriginalCalls.add(v);
            }
          }
        }
        regionA.mEndOfChunk = nRegion.getEnd();
        regionB.mStartOfChunk = nRegion.getEnd();
      }
    }
  }

  /**
   * @return start of chunk (0 based inclusive) represented by this object
   */
  public int startOfChunk() {
    return mStartOfChunk;
  }

  /**
   * @return end of chunk (0 based exclusive) represented by this object
   */
  public int endOfChunk() {
    return mEndOfChunk;
  }

  /**
   * Check that this object has had both ends correct for dangling.
   * @return true iff fixed
   */
  public boolean isFixed() {
    return mStartFixed && mEndFixed;
  }

  /**
   * @return iterator over complex regions contained in this chunk, these will be returned in coordinate order. Regions at the
   * boundary may not be represented fully or at all if this chunk hasn't had <code>fixDangling</code> invoked on that boundary.
   */
  @Override
  public Iterator<ComplexRegion> iterator() {
    if (!mStartFixed || !mEndFixed) {
      throw new IllegalStateException("Dangling boundaries haven't been fixed (left: " + mStartFixed + " right: " + mEndFixed + ")");
    }
    return iteratorInternal();
  }

  Iterator<ComplexRegion> iteratorInternal() {
    return mRegions.iterator();
  }

  /**
   * @return number of complex regions
   */
  int size() {
    return mRegions.size();
  }

  /**
   * Returns the first complex region in this chunk if it is close enough
   * to chunk boundaries that it potentially needs to be merged, or in the case
   * of a hyper complex region extends into the previous chunk.
   * @return the first region, null implies that leftmost region does not interfere with previous chunks regions
   */
  public ComplexRegion danglingStart() {
    return mStartDangling;
  }

  /**
   * Returns the last complex region in this chunk if it is close enough
   * to chunk boundaries that it potentially needs to be merged, or in the case
   * of a hyper complex region extends into the next chunk.
   * @return the last region, null implies that rightmost region does not interfere with next chunks regions
   */
  public ComplexRegion danglingEnd() {
    return mEndDangling;
  }

  /**
   * Get original calls.
   * @return Returns the original calls.
   */
  public List<Variant> getOriginalCalls() {
    return mOriginalCalls;
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append(mStartFixed ? "" : "*").append(mReferenceName).append("[").append(mStartOfChunk).append("..").append(mEndOfChunk).append(")").append(mEndFixed ? "" : "*");
    sb.append("#").append(mOriginalCalls.size());
    sb.append(mRegions);
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    final Iterator<ComplexRegion> it = mRegions.iterator();
    if (it.hasNext()) {
      final ComplexRegion cr0 = it.next();
      Exam.assertFalse(cr0.type() == ComplexRegion.RegionType.INTERESTING);
      int lastEndPosition = cr0.getEnd();
      Exam.assertTrue(lastEndPosition + ":" + mStartOfChunk, lastEndPosition == mStartOfChunk || cr0.getStart() >= mStartOfChunk);
      while (it.hasNext()) {
        final ComplexRegion cr = it.next();
        Exam.assertTrue(cr.getStart() >= lastEndPosition + mInterestingSeparation);
        lastEndPosition = cr.getEnd();
        Exam.assertFalse(cr.type() == ComplexRegion.RegionType.INTERESTING);
      }
    }
    return true;
  }

  @Override
  public boolean integrity() {
    if (mStartFixed) {
      Exam.assertTrue(mStartDangling == null || mStartDangling.type() == ComplexRegion.RegionType.HYPER);
    }
    if (mStartDangling != null && mEndDangling != null) {
      Exam.assertTrue(mStartDangling == mEndDangling || !joinInteresting(mStartDangling.getEnd() - 1, mEndDangling.getStart(), mStartDangling.getIndelEnd(), mEndDangling.getIndelStart()));
    }
    if (mEndFixed) {
      Exam.assertTrue(mEndDangling == null || mEndDangling.type() == ComplexRegion.RegionType.HYPER);
    }
    Exam.assertTrue(mStartOfChunk >= 0);
    Exam.assertTrue(mEndOfChunk >= mStartOfChunk);
    Exam.assertTrue(mHyperComplexLength > mInterestingSeparation);
    return true;
  }
}
