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

package com.rtg.util.machine;

import com.rtg.variant.sv.discord.Orientation;

import net.sf.samtools.SAMRecord;

/**
 * Enumeration of the orientations of paired reads which is technology dependent.
 */
public enum MachineOrientation {

 /** Forward Reverse orientation in paired reads. */
  FR {

    @Override
    public boolean firstOnTemplate(SAMRecord sam) {
      return !sam.getReadNegativeStrandFlag();
    }

    @Override
    public Orientation orientation(SAMRecord rec) {
      if (rec.getReadNegativeStrandFlag()) {
        if (rec.getMateNegativeStrandFlag()) {
          return Orientation.DD;
        } else {
          return Orientation.DU;
        }
      } else {
        if (rec.getMateNegativeStrandFlag()) {
          return Orientation.UD;
        } else {
          return Orientation.UU;
        }
      }
    }

    @Override
    public boolean orientationOkay(boolean firstOnReferenceReverse, boolean firstOnReferenceLeft, boolean secondOnReferenceReverse, boolean secondOnReferenceLeft) {
      assert firstOnReferenceLeft ^ secondOnReferenceLeft;
      return !firstOnReferenceReverse && secondOnReferenceReverse;
    }

    @Override
    public PairOrientation getMateOrientation(PairOrientation pairOrientation) {
      if (PairOrientation.F1.equals(pairOrientation)) {
        return PairOrientation.R2;
      } else if (PairOrientation.F2.equals(pairOrientation)) {
        return PairOrientation.R1;
      } else if (PairOrientation.R1.equals(pairOrientation)) {
        return PairOrientation.F2;
      } else if (PairOrientation.R2.equals(pairOrientation)) {
        return PairOrientation.F1;
      }
      return null;
    }

    @Override
    public boolean isMateUpstream(PairOrientation pairOrientation) {
      return PairOrientation.F1.equals(pairOrientation) || PairOrientation.F2.equals(pairOrientation);
    }
  },
  /** reverse forward orientation in paired reads */
  RF {

    @Override
    public boolean firstOnTemplate(SAMRecord sam) {
      return sam.getReadNegativeStrandFlag();
    }

    @Override
    public Orientation orientation(SAMRecord rec) {
      if (rec.getReadNegativeStrandFlag()) {
        if (rec.getMateNegativeStrandFlag()) {
          return Orientation.UU;
        } else {
          return Orientation.UD;
        }
      } else {
        if (rec.getMateNegativeStrandFlag()) {
          return Orientation.DU;
        } else {
          return Orientation.DD;
        }
      }
    }

    @Override
    public boolean orientationOkay(boolean firstOnReferenceReverse, boolean firstOnReferenceLeft, boolean secondOnReferenceReverse, boolean secondOnReferenceLeft) {
      assert firstOnReferenceLeft ^ secondOnReferenceLeft;
      return firstOnReferenceReverse && !secondOnReferenceReverse;
    }

    @Override
    public PairOrientation getMateOrientation(PairOrientation pairOrientation) {
      return FR.getMateOrientation(pairOrientation);
    }

    @Override
    public boolean isMateUpstream(PairOrientation pairOrientation) {
      return !(PairOrientation.F1.equals(pairOrientation) || PairOrientation.F2.equals(pairOrientation));
    }
  },
  /** Forward Forward orientation in paired reads. */
  TANDEM {

    @Override
    public boolean firstOnTemplate(SAMRecord sam) {
      if (sam.getFirstOfPairFlag() && !sam.getReadNegativeStrandFlag()) {
        return true;
      }
      if (sam.getSecondOfPairFlag() && sam.getReadNegativeStrandFlag()) {
        return true;
      }
      return false;
    }

    @Override
    public Orientation orientation(SAMRecord rec) {
      if (rec.getFirstOfPairFlag()) {
        if (rec.getReadNegativeStrandFlag()) {
          if (rec.getMateNegativeStrandFlag()) {
            return Orientation.DU;
          } else {
            return Orientation.DD;
          }
        } else {
          if (rec.getMateNegativeStrandFlag()) {
            return Orientation.UU;
          } else {
            return Orientation.UD; // Normal mated orientation on the forward strand
          }
        }
      } else {
        if (rec.getReadNegativeStrandFlag()) {
          if (rec.getMateNegativeStrandFlag()) {
            return Orientation.UD; // Normal mated orientation on the reverse strand
          } else {
            return Orientation.UU;
          }
        } else {
          if (rec.getMateNegativeStrandFlag()) {
            return Orientation.DD;
          } else {
            return Orientation.DU;
          }
        }
      }
    }

    @Override
    public boolean orientationOkay(boolean firstOnReferenceReverse, boolean firstOnReferenceLeft, boolean secondOnReferenceReverse, boolean secondOnReferenceLeft) {
      assert firstOnReferenceLeft ^ secondOnReferenceLeft;
      return /*F1F2*/ (firstOnReferenceLeft && !firstOnReferenceReverse && !secondOnReferenceReverse && !secondOnReferenceLeft)
            || /*R2R1*/ (firstOnReferenceReverse && !firstOnReferenceLeft && secondOnReferenceReverse && secondOnReferenceLeft);
    }

    @Override
    public PairOrientation getMateOrientation(PairOrientation pairOrientation) {
      if (PairOrientation.F1.equals(pairOrientation)) {
        return PairOrientation.F2;
      } else if (PairOrientation.F2.equals(pairOrientation)) {
        return PairOrientation.F1;
      } else if (PairOrientation.R1.equals(pairOrientation)) {
        return PairOrientation.R2;
      } else if (PairOrientation.R2.equals(pairOrientation)) {
        return PairOrientation.R1;
      }
      return null;
    }

    @Override
    public boolean isMateUpstream(PairOrientation pairOrientation) {
      return PairOrientation.F1.equals(pairOrientation) || PairOrientation.R2.equals(pairOrientation);
    }
  },

  /** Allow any orientation, many methods unsupported */
  ANY {
    @Override
    public boolean firstOnTemplate(SAMRecord sam) {
      throw new UnsupportedOperationException("Not supported.");
    }

    @Override
    public Orientation orientation(SAMRecord rec) {
      throw new UnsupportedOperationException("Not supported.");
    }

    @Override
    public boolean orientationOkay(boolean firstOnReferenceReverse, boolean firstOnReferenceLeft, boolean secondOnReferenceReverse, boolean secondOnReferenceLeft) {
      assert firstOnReferenceLeft ^ secondOnReferenceLeft;
      return true;
    }

    @Override
    public PairOrientation getMateOrientation(PairOrientation pairOrientation) {
      throw new UnsupportedOperationException();
    }

    @Override
    public boolean isMateUpstream(PairOrientation pairOrientation) {
      throw new UnsupportedOperationException();
    }
  };

  /**
   * Test if the read has a mapped mate which is in the correct orientation
   * that the two could form a normal mate on the given template.
   * This does NOT test that the insert size lies within the expected distribution.
   * @param sam the SAM record being checked.
   * @return true iff there is a mate in the correct orientation.
   */
  public final boolean hasValidMate(final SAMRecord sam) {
    assert !sam.getReadUnmappedFlag();
    if (!sam.getReadPairedFlag()) {
      return false;
    }
    if (sam.getMateUnmappedFlag()) {
      return false;
    }
    if (!sam.getReferenceName().equals(sam.getMateReferenceName())) {
      return false;
    }
    if (sam.getFirstOfPairFlag()) {
      return orientationOkay(sam.getAlignmentStart(), sam.getReadNegativeStrandFlag(), sam.getMateAlignmentStart(), sam.getMateNegativeStrandFlag());
    } else {
      return orientationOkay(sam.getMateAlignmentStart(), sam.getMateNegativeStrandFlag(), sam.getAlignmentStart(), sam.getReadNegativeStrandFlag());
    }
  }

  /**
   * Test if record would be first on template if normally mated.
   * @param sam the SAM record being checked.
   * @return true iff the record would appear first on template in normal case.
   */
  public abstract boolean firstOnTemplate(final SAMRecord sam);

  /**
   * Determine the constraint direction represented by this SAM record.
   * @param rec the sam record to determine the orientation of
   * @return the orientation of the sam record
   */
  public abstract Orientation orientation(SAMRecord rec);

  /**
   * check orientation of pairs against required orientation
   * @param aPos position of <code>a</code> read arm on reference (requires this is &quot;left&quot; arm)
   * @param aRev true if <code>a</code> read arm reverse complement  (requires this is &quot;left&quot; arm)
   * @param bPos position of <code>b</code> read arm on reference  (requires this is &quot;right&quot; arm)
   * @param bRev true if <code>b</code> read arm is reverse complement (requires this is &quot;right&quot; arm)
   * @return true if read arms obey required pair orientation
   */
  public boolean orientationOkay(int aPos, boolean aRev, int bPos, boolean bRev) {
    if (aPos <= bPos) {
      return orientationOkay(aRev, true, bRev, false);
    } else {
      return orientationOkay(bRev, false, aRev, true);
    }
  }

  /**
   * check orientation of pairs against required orientation
   * @param firstOnReferenceReverse true if smallest start position mapped read is reverse complement
   * @param firstOnReferenceLeft true if smallest start position mapped read is from &quot;left&quot arm
   * @param secondOnReferenceReverse true if greatest start position mapped read is reverse complement
   * @param secondOnReferenceLeft true if greatest start position mapped read is from &quot;left&quot arm
   * @return true if read arms obey required pair orientation
   */
  public abstract boolean orientationOkay(boolean firstOnReferenceReverse, boolean firstOnReferenceLeft, boolean secondOnReferenceReverse, boolean secondOnReferenceLeft);

  /**
   * Return the pair orientation of the mate of this pair orientation. Useful for
   * unmapped breakpoint/structural variation detection.
   * @param pairOrientation the pair orientation of the mapped read.
   * @return the pair orientation of the mate
   */
  public abstract PairOrientation getMateOrientation(PairOrientation pairOrientation);

  /**
   * Return whether the mate is up or downstream from a mapped read, assuming a concordant pair.
   * @param pairOrientation the pair orientation of the mapped read.
   * @return true if upstream (higher template start position), false if downstream (lower template start position).
   */
  public abstract boolean isMateUpstream(PairOrientation pairOrientation);
}
