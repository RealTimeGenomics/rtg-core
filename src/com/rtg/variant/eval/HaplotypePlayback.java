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
package com.rtg.variant.eval;

import static com.rtg.util.StringUtils.LS;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.NoSuchElementException;
import java.util.Queue;

import com.rtg.mode.DnaUtils;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.Integrity;

/**
 */
public class HaplotypePlayback implements Integrity, Comparable<HaplotypePlayback> {

  private static final int INVALID = -1;

  /** Sorted list of variants yet to be processed. */
  private final Queue<OrientedVariant> mVariants;

  private final byte[] mTemplate;

  /** Position in template (start of current variant if one is active). 0 based.*/
  private int mTemplatePosition = -1;

  /** Position in variant. INVALID if not currently in variant. 0 based. */
  private int mPositionInVariant;

  /** Variant that currently in or next one. */
  private OrientedVariant mNextVariant;

  HaplotypePlayback(final byte[] template) {
    mVariants = new LinkedList<>();
    mTemplate = template;
    mNextVariant = null;
    mPositionInVariant = INVALID;
    assert globalIntegrity();
  }

  private HaplotypePlayback(HaplotypePlayback old) {
    mVariants = new LinkedList<>(old.mVariants);
    mTemplate = old.mTemplate;
    mNextVariant = old.mNextVariant;
    mTemplatePosition = old.mTemplatePosition;
    mPositionInVariant = old.mPositionInVariant;
    assert globalIntegrity();
  }

  /**
   * Add an <code>OrientedVariant</code> to this <code>HaplotypePlayback</code>
   * @param v the <code>OrientedVariant</code> to add
   */
  public void addVariant(OrientedVariant v) {
    //System.err.println("Adding Variant: " + v);
    if (mNextVariant == null) {
      mNextVariant = v;
    } else {
      mVariants.add(v);
    }
  }

  /**
   * Get the current nucleotide in this playback. Taken from template or variant.
   * @return the current nucleotide in this playback.
   */
  byte nt() {
    if (mPositionInVariant == INVALID) {
      return mTemplatePosition < mTemplate.length ? mTemplate[mTemplatePosition] : DnaUtils.UNKNOWN_RESIDUE;
    }
    return mNextVariant.nt(mNextVariant.isAlleleA())[mPositionInVariant];
  }

  /**
   * Test if there are more nucleotides available.
   * @return true if there are more nucleotides available.
   */
  boolean hasNext() {
    return mTemplatePosition < mTemplate.length - 1;
  }

  /**
   * Test whether a deficit of variant bases are upstream in the queue in order to perform a step.
   * A result of false indicates that no variants need to be immediately enqueued
   */
  boolean wantsFutureVariantBases() {
    if (mNextVariant == null) {
      return true;
    }
    if (mPositionInVariant != INVALID && mPositionInVariant < mNextVariant.nt(mNextVariant.isAlleleA()).length - 1) {
      return false;
    }
    for (OrientedVariant v : mVariants) {
      if (v.nt(v.isAlleleA()).length > 0) {
        return false;
      }
    }
    return true;
  }

  /**
   * Step to the next nucleotide.
   * @throws NoSuchElementException if no more nucleotides available.
   */
  void next() {
    if (!hasNext()) {
      throw new NoSuchElementException();
    }
    if (isOnTemplate()) {
      mTemplatePosition++;
      if (mNextVariant != null && mNextVariant.getStart() == mTemplatePosition) { // Position to consume the variant
        mPositionInVariant = 0;
      }
    } else {
      assert mPositionInVariant != INVALID;
      mPositionInVariant++;
    }
    assert mNextVariant != null || mPositionInVariant == INVALID;
    if (mNextVariant != null) {
      while (true) {
        //dont forget that variant may come back to back
        final byte[] norm = mNextVariant.nt(mNextVariant.isAlleleA());
        if (mPositionInVariant != norm.length) { // Haven't reached the end of the current variant.
          //this also catches the case when mPositionInVariant == INVALID
          assert integrity();
          break;
        }
        // Finished variant, so position for next baseStart consuming next variant from the queue (which may be a deletion, in
        mTemplatePosition = mNextVariant.getEnd();
        mPositionInVariant = INVALID;
        if (!mVariants.isEmpty()) {
          mNextVariant = mVariants.poll();
        } else {
          mNextVariant = null;
          assert integrity();
          break;
        }
        if (mTemplatePosition < mNextVariant.getStart()) {
          assert integrity();
          break;
        }
        mPositionInVariant = 0;
        //System.err.println("templatePosition=" + mTemplatePosition + " varPosition=" + (mNextVariant.variant().start() - 1));
        assert mTemplatePosition == mNextVariant.getStart();
      }
    }
    assert integrity();
  }

  /**
   * Test if the haplotype is currently within a variant.
   * @return true iff not in a variant.
   */
  boolean isOnTemplate() {
    return mPositionInVariant == INVALID;
  }

  /**
   * Get the current position within a variant.
   * @return the current position within the current variant.
   */
  int positionInVariant() {
    return mPositionInVariant;
  }

  /**
   * Get the current position in the template. Use <code>isOnTemplate</code> to determine
   * whether the haplotype is also within a current variant.
   * 0 based and may be equal to length of template.
   * @return the current position in the template.
   */
  int templatePosition() {
    return mTemplatePosition;
  }

  /**
   * Get the variant we are currently in
   * @return the current variant or null if we aren't in a variant
   */
  OrientedVariant currentVariant() {
    if (isOnTemplate()) {
      return null;
    }
    return mNextVariant;
  }

  /**
   * Force the template position to the first template position at or beyond "position" and the current template position which is not
   * in a variant. Force the state of any otherwise unmarked variants as <code>UNKNOWN</code>.
   * @param position to be forced to. (0 based)
   */
  void moveForward(final int position) {
    assert position >= 0 && position <= mTemplate.length;
    if (!isOnTemplate()) {
      throw new IllegalStateException();
    }
    mTemplatePosition = position - 1;
    next();
    assert templatePosition() >= position && isOnTemplate();
    assert integrity();
  }


  @Override
  public String toString() {
    return "HaplotypePlayback: position=" + templatePosition() + " inPosition=" + mPositionInVariant + LS
    + "current:" + nullToString(mNextVariant) + LS
    + "future:" + mVariants.toString() + LS;

  }

  private String nullToString(final Object obj) {
    if (obj == null) {
      return null;
    }
    return obj.toString();
  }

  @Override
  public boolean globalIntegrity() {
    Exam.assertTrue(integrity());
    /*int last = -1;
    for (final OrientedVariant dv : mVariants) {
      final int position = dv.variant().start() - 1;
      final int end = position + dv.variant().end() - 1;
      //System.err.println(dv.toVerboseString());
      //System.err.println("position=" + position + " last=" + last + " end=" + end + " length=" + mTemplate.length);
      Assert.assertTrue(dv.toString(), position > last && end <= mTemplate.length);
      last = end;
    }*/
    return true;
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mVariants);
    Exam.assertNotNull(mTemplate);
    Exam.assertTrue(-1 <= mTemplatePosition && mTemplatePosition <= mTemplate.length);
    if (mNextVariant == null) {
      Exam.assertTrue(mVariants.isEmpty());
    }
    if (mPositionInVariant == HaplotypePlayback.INVALID) {
      Exam.assertTrue(toString(), mNextVariant == null || mTemplatePosition <= mNextVariant.getStart());
    } else {
      Exam.assertTrue(mNextVariant != null && mTemplatePosition == mNextVariant.getStart());
      if (!(0 <= mPositionInVariant && mPositionInVariant < mNextVariant.nt(mNextVariant.isAlleleA()).length)) {
        System.err.println(this);
      }
      Exam.assertTrue(0 <= mPositionInVariant && mPositionInVariant < mNextVariant.nt(mNextVariant.isAlleleA()).length);
    }
    return true;
  }

  /** @return a copy of this object for use in another path */
  public HaplotypePlayback copy() {
    return new HaplotypePlayback(this);
  }

  /**
   * Performs a comparison of the variants in this path that aren't yet resolved.
   * @param that the HaplotypePlayback to compare to.
   * @return an indication that this is further advanced along its path or not.
   */
  @Override
  public int compareTo(HaplotypePlayback that) {
    final int position = this.templatePosition() - that.templatePosition();
    if (position != 0) {
      return position;
    }
    if (this.mNextVariant == null) {
      if (that.mNextVariant == null) {
        return 0;
      }
      return -1;
    }
    if (that.mNextVariant == null) {
      return 1;
    }

    final int current = this.mNextVariant.compareTo(that.mNextVariant);
    if (current != 0) {
      return current;
    }
    final int varPos = this.mPositionInVariant - that.mPositionInVariant;
    if (varPos != 0) {
      return varPos;
    }
    final Iterator<OrientedVariant> thisIt = this.mVariants.iterator();
    final Iterator<OrientedVariant> thatIt = that.mVariants.iterator();
    while (thisIt.hasNext()) {
      if (!thatIt.hasNext()) {
        return 1;
      }
      final int future = thisIt.next().compareTo(thatIt.next());
      if (future != 0) {
        return future;
      }

    }
      if (thatIt.hasNext()) {
        return -1;
      }
    return 0;
  }
  @Override
  public boolean equals(Object obj) {
    if (obj == null) {
      return false;
    }
    return compareTo((HaplotypePlayback) obj) == 0;
  }

  @Override
  public int hashCode() {
    int hash = Utils.hash(new Object[] {this.templatePosition(), this.mNextVariant, this.mPositionInVariant});
    for (final OrientedVariant v : this.mVariants) {
      hash = Utils.pairHash(hash, v.hashCode());
    }
    return hash;
  }


}
