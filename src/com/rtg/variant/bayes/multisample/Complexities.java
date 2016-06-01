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
package com.rtg.variant.bayes.multisample;

import java.util.ArrayList;
import java.util.Deque;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.util.SynchronizedLinkedList;
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

  static final String MER_FINDER_IMPL = GlobalFlags.getStringValue(CoreGlobalFlags.COMPLEX_REGION_SIMPLE_REPEAT_IMPL);

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
        case "multi-nmer":
          measurer = new SimpleRepeatMeasurer(refNts);
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
      final Variant v = new Variant(new VariantLocus(referenceName, ac.position(), sspEnd));
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
    int lastInteresting = -1;
    boolean forceComplex = false;

    int lastInterestingShadow = -1;
    for (final Variant c : chunk) {
      if (c.isOverflow()) {
        // Want to write an overflow region, but need to output any current complex region first
        addRegion(regions, firstInteresting, lastInteresting, forceComplex);
        firstInteresting = -1;   // Reset state of "in-construction" region
        lastInteresting = -1;
        lastInterestingShadow = -1;
        forceComplex = false;
        addOverflowRegion(regions, c.getLocus());
      } else if (c.isInteresting()) {
        //        System.err.println("first: " + firstInteresting + " last: " + lastInteresting + " forced: " + forceComplex + " call: " + c + " SoC: " + mStartOfChunk + " EoC: " + mEndOfChunk);
        assert firstInteresting > -1 || (lastInteresting == -1 && !forceComplex);
        assert firstInteresting == -1 || (firstInteresting >= mStartOfChunk && lastInteresting >= firstInteresting && lastInteresting <= mEndOfChunk);
        final int outputPosition = c.getLocus().getStart();

        // See whether the current variant is var enough away that we need to close out the "in-construction" region
        if (lastInteresting > 0 && !joinInteresting(lastInterestingShadow - 1, outputPosition - c.getIndelLength())) {
          // Can't joint current variant to in-construction complex region, so finish the complex region and add it
          addRegion(regions, firstInteresting, lastInteresting, forceComplex);
          firstInteresting = -1;   // Reset state of "in-construction" region
          lastInteresting = -1;
          lastInterestingShadow = -1;
          forceComplex = false;
        }

        // Add current variant info into "in-construction" state
        if (firstInteresting == -1) {
          firstInteresting = outputPosition;
        }
        lastInteresting = Math.max(lastInteresting, c.getLocus().getEnd());
        lastInterestingShadow = Math.max(lastInterestingShadow, lastInteresting + c.getIndelLength());
        forceComplex |= c.isForceComplex();
        //System.err.println("getComplexRegions: " + outputPosition + " : " + c.getLocus().getEnd() + " : " + firstInteresting + " : " + lastInteresting);
      }
    }
    addRegion(regions, firstInteresting, lastInteresting, forceComplex);
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

  boolean joinInteresting(int positionA, int positionB) {
    final int repeatTotal = mRepeatMeasurer.measureRepeats(positionA, positionB);
    return positionA + mInterestingSeparation + repeatTotal >= positionB;
  }


  final ComplexRegion computeStartDangle(Deque<ComplexRegion> regions) {
    return regions.peekFirst();
  }

  final ComplexRegion computeEndDangle(Deque<ComplexRegion> regions) {
    return regions.peekLast();
  }

  void addRegion(final Deque<ComplexRegion> regions, int firstInteresting, int endInteresting, boolean forceComplex) {
    if (firstInteresting != -1) {
      assert endInteresting >= firstInteresting : endInteresting + ":" + firstInteresting;
      final boolean simpleInteresting = (endInteresting - firstInteresting) == 1 && !forceComplex;
      final boolean hyper = !simpleInteresting && endInteresting - firstInteresting > mHyperComplexLength;
      final ComplexRegion.RegionType regionType = hyper ? ComplexRegion.RegionType.HYPER : (simpleInteresting ? ComplexRegion.RegionType.INTERESTING : ComplexRegion.RegionType.COMPLEX);
      final ComplexRegion com = new ComplexRegion(mReferenceName, firstInteresting, endInteresting, regionType);
      regions.add(com);
    }
  }

  void addOverflowRegion(final Deque<ComplexRegion> regions, final VariantLocus locus) {
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
    return new ComplexRegion(regionA.getSequenceName(), start, end, type);
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
    final int aLastPosition = aLast.getEnd() - 1;
    final int bFirstPosition = bFirst.getStart();

    if (!regionA.joinInteresting(aLastPosition, bFirstPosition)) {
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
          for (final Variant v : regionB.mOriginalCalls) { //copy over the calls which are in the dangling region, required when the complex caller fails to make a call and prints originals
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
    sb.append(mRegions.toString());
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    final Iterator<ComplexRegion> it = mRegions.iterator();
    if (it.hasNext()) {
      final ComplexRegion cr0 = it.next();
      Exam.assertFalse(cr0.type() == ComplexRegion.RegionType.INTERESTING);
      int lastEndPosition = cr0.getEnd();
      Exam.assertTrue("" + lastEndPosition + ":" + mStartOfChunk, lastEndPosition == mStartOfChunk || cr0.getStart() >= mStartOfChunk);
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
      Exam.assertTrue(mStartDangling == mEndDangling || !joinInteresting(mStartDangling.getEnd() - 1, mEndDangling.getStart()));
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
