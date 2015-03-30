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

  /** Sorted list of variations. */
  private final Queue<OrientedVariant> mVariations;

  private final byte[] mTemplate;

  /** Position in template (start of current variation if one is active). 0 based.*/
  private int mTemplatePosition = -1;

  /** Position in variation. INVALID if not currently in variation. 0 based. */
  private int mPositionInVariation;

  /** Variation that currently in or next one. */
  private OrientedVariant mNextVariant;

  HaplotypePlayback(final byte[] template) {
    mVariations = new LinkedList<>();
    mTemplate = template;
    mNextVariant = null;
    mPositionInVariation = INVALID;
    assert globalIntegrity();
  }

  private HaplotypePlayback(HaplotypePlayback old) {
    mVariations = new LinkedList<>(old.mVariations);
    mTemplate = old.mTemplate;
    mNextVariant = old.mNextVariant;
    mTemplatePosition = old.mTemplatePosition;
    mPositionInVariation = old.mPositionInVariation;
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
      mVariations.add(v);
    }
  }

  /**
   * Get the current nucleotide in this playback. Taken from template or variation.
   * @return the current nucleotide in this playback.
   */
  byte nt() {
    if (mPositionInVariation == INVALID) {
      return mTemplatePosition < mTemplate.length ? mTemplate[mTemplatePosition] : DnaUtils.UNKNOWN_RESIDUE;
    }
    return mNextVariant.variant().nt(mNextVariant.isAlleleA())[mPositionInVariation];
  }

  /**
   * Test if if there are more nucleotides available.
   * @return true if there are more nucleotides available.
   */
  boolean hasNext() {
    return mTemplatePosition < mTemplate.length - 1;
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
      if (mNextVariant != null && mNextVariant.variant().getStart() == mTemplatePosition) {
        mPositionInVariation = 0;
      }
    } else {
      assert mPositionInVariation != INVALID;
      mPositionInVariation++;
    }
    if (mNextVariant == null) {
      assert mPositionInVariation == INVALID;
    } else {
      while (true) {
        //dont forget that variations may come back to back
        final byte[] norm = mNextVariant.variant().nt(mNextVariant.isAlleleA());
        final int effectiveLength = norm.length;
        if (mPositionInVariation != effectiveLength) {
          //this also catches the case when mPositionInVariation == INVALID
          assert integrity();
          break;
        }
        mTemplatePosition = mNextVariant.variant().getEnd();
        mPositionInVariation = INVALID;
        if (!mVariations.isEmpty()) {
          mNextVariant = mVariations.poll();
        } else {
          mNextVariant = null;
          assert integrity();
          break;
        }
        if (mTemplatePosition < mNextVariant.variant().getStart()) {
          assert integrity();
          break;
        }
        mPositionInVariation = 0;
        //System.err.println("templatePosition=" + mTemplatePosition + " varPosition=" + (mNextVariant.variant().start() - 1));
        assert mTemplatePosition == mNextVariant.variant().getStart();
      }
    }
    assert integrity();
  }

  /**
   * Test if the haplotype is currently within a variant.
   * @return true iff not in a variant.
   */
  boolean isOnTemplate() {
    return mPositionInVariation == INVALID;
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
   * in a variation. Force the state of any otherwise unmarked variants as <code>UNKNOWN</code>.
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
    return "HaplotypePlayback: position=" + templatePosition() + " inPosition=" + mPositionInVariation + LS
    + "current:" + nullToString(mNextVariant) + LS
    + "future:" + mVariations.toString() + LS;

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
    for (final OrientedVariant dv : mVariations) {
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
    Exam.assertNotNull(mVariations);
    Exam.assertNotNull(mTemplate);
    Exam.assertTrue(-1 <= mTemplatePosition && mTemplatePosition <= mTemplate.length);
    if (mNextVariant == null) {
      Exam.assertTrue(mVariations.isEmpty());
    }
    if (mPositionInVariation == HaplotypePlayback.INVALID) {
      Exam.assertTrue(toString(), mNextVariant == null || mTemplatePosition <= mNextVariant.variant().getStart());
    } else {
      Exam.assertTrue(mNextVariant != null && mTemplatePosition == mNextVariant.variant().getStart());
      if (!(0 <= mPositionInVariation && mPositionInVariation < mNextVariant.variant().nt(mNextVariant.isAlleleA()).length)) {
        System.err.println(this);
      }
      Exam.assertTrue(0 <= mPositionInVariation && mPositionInVariation < mNextVariant.variant().nt(mNextVariant.isAlleleA()).length);
    }
    return true;
  }

  /** @return a copy of this object for use in another path */
  public HaplotypePlayback copy() {
    return new HaplotypePlayback(this);
  }

  /**
   * Performs a comparison of the variations in this path that aren't yet resolved.
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
    final int varPos = this.mPositionInVariation - that.mPositionInVariation;
    if (varPos != 0) {
      return varPos;
    }
    final Iterator<OrientedVariant> thisIt = this.mVariations.iterator();
    final Iterator<OrientedVariant> thatIt = that.mVariations.iterator();
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
    int hash = Utils.hash(new Object[] {this.templatePosition(), this.mNextVariant, this.mPositionInVariation});
    for (final OrientedVariant v : this.mVariations) {
      hash = Utils.pairHash(hash, v.hashCode());
    }
    return hash;
  }


}
