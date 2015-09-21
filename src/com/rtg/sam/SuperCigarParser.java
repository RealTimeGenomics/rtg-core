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

import java.util.Arrays;

import com.rtg.mode.DNA;
import com.rtg.mode.DnaUtils;

/**
 * A super class for parsing super cigars.
 * Can also be used for standard SAM cigars, if you call <code>setStandardCigar</code>.
 *
 * To implement a new parser, just subclass this class and override some of the
 * <code>do*(...)</code> call-back methods.
 *
 */
public class SuperCigarParser {

  private static final byte N_CODE = (byte) DNA.N.ordinal();

  protected String mCigar;
  protected boolean mIsUnmapped = false;

  byte[] mReadDelta;
  private String mReadDeltaString; // just for error reporting
  int mReadDeltaPos; // current position in mReadDelta.

  /**
   * Overrides <code>mReadDelta</code>.
   * If this is non-null, then we have a standard cigar
   */
  private byte[] mRead;
  int mReadLength;
  int mReadPos;

  protected int mOperationLength;

  private byte[] mTemplate;
  private int mTemplateLength;
  private int mTemplateStart; // zero-based start position.
  private int mTemplatePos;

  /** The template position where the CG overlap region starts (-1 means no overlap). */
  private int mTemplateOverlapStart;

  /** The template position just after the end of the CG overlap region (-1 means no overlap). */
  private int mTemplateOverlapEnd;

  /**
   * The read position where the CG overlap region starts.
   * If there is an insertion at that point, then this will be after that insertion.
   */
  private int mReadOverlapStart;

  /**
   * The read position where the second half of the CG overlap region starts.
   */
  private int mReadOverlapMiddle;

  /**
   * The first read position after the end of the CG overlap.
   * If there is an insertion at that point, then this will be
   * the first read position of that insertion.
   */
  private int mReadOverlapEnd;

  private static final int TABLE_SIZE = 8; // must be > maximum CG overlap and a power of 2.
  private static final int TABLE_MASK = TABLE_SIZE - 1;

  /**
   * This remembers the mapping from the last few template positions to the largest corresponding
   * read position.  Used to calculate where on the read the CG overlap starts.
   */
  private final int[] mTemplateToReadPos = new int[TABLE_SIZE];

  /**
   * @param superCigar Extended cigar string
   * @param readDelta  The read nucleotides that are missing/modified in comparison to the template.
   */
  public void setCigar(String superCigar, String readDelta) {
    mCigar = superCigar;
    //System.err.println("cigar parser setCigar(" + superCigar + ", " + readDelta + ")");
    if (readDelta == null) {
      mReadDelta = new byte[0];
      mReadDeltaString = "";
    } else {
      mReadDelta = DnaUtils.encodeString(readDelta);
      mReadDeltaString = readDelta;
    }
    mReadOverlapStart = -1;
    mReadOverlapMiddle = -1;
    mReadOverlapEnd = -1;
    mTemplateOverlapStart = -1;
    mTemplateOverlapEnd = -1;
    Arrays.fill(mTemplateToReadPos, 0); // in case the left side of the overlap goes off the read?
    mRead = null;
    mReadLength = 0;
  }

  /**
   * @param cigar Standard SAM cigar string (must use X and =, not just M, unless template is also set).
   * @param read  The complete sequence read.  Cannot be null.  The array can be longer than necessary.
   * @param readLength length of the read
   */
  public void setStandardCigar(String cigar, byte[] read, int readLength) {
    setCigar(cigar, null);
    mRead = read;
    mReadLength = readLength;
  }

  /**
   * @param cigar Standard SAM cigar string (must use X and =, not just M, unless template is also set).
   * @param read  The complete sequence read.  Cannot be null.
   */
  public void setStandardCigar(String cigar, String read) {
    setStandardCigar(cigar, DnaUtils.encodeString(read), read.length());
  }

  /**
   * Set template.
   * @param template the template sequence.
   * @param length length of template used
   */
  public void setTemplate(byte[] template, int length) {
    mTemplate = template;
    mTemplateLength = length;
  }

  /**
   * Set template.
   * @param template the template sequence.
   */
  public void setTemplate(byte[] template) {
    if (template == null) {
      setTemplate(null, 0);
    } else {
      setTemplate(template, template.length);
    }
  }

  /**
   * Set the template start position
   * @param zeroBasedStart start position on the template. Cannot be less than zero.
   */
  public void setTemplateStart(int zeroBasedStart) {
    assert zeroBasedStart >= 0;
    mTemplateStart = zeroBasedStart;
  }

  /**
   * Get template.
   * @return Returns the template.
   */
  public byte[] getTemplate() {
    return mTemplate;
  }

//  protected byte[] getRead() {
//    return mRead;
//  }
  protected byte getReadByte(int position) {
    if (position > mReadLength) {
      throw new ArrayIndexOutOfBoundsException(position);
    }
    return mRead[position];
  }
  protected int getReadLength() {
    return mReadLength;
  }

  public int getTemplateLength() {
    return mTemplateLength;
  }

  /**
   * Get the start position of the read on the template.
   * @return Returns the zero-based start position.
   */
  public int getTemplateStartPosition() {
    return mTemplateStart;
  }

  /**
   * The current position on the template.
   * @return Returns the zero-based position.
   */
  public int getTemplatePosition() {
    return mTemplatePos;
  }

  /**
   * The current position on the unrolled read.
   * @return Returns the zero-based position.
   */
  public int getReadPosition() {
    return mReadPos;
  }

  /** @return The (largest) read position where the CG overlap starts (-1 if no overlap or cigar not parsed yet). */
  public int getReadOverlapStart() {
    return mReadOverlapStart;
  }

  /** @return The first read position in the second half of the CG overlap region (-1 if no overlap or cigar not parsed yet). */
  public int getReadOverlapMiddle() {
    return mReadOverlapMiddle;
  }

  /** @return The first read position after the CG overlap starts (-1 if no overlap or cigar not parsed yet). */
  public int getReadOverlapEnd() {
    return mReadOverlapEnd;
  }

  /** @return The position on the template where the CG overlap starts (-1 if no overlap or cigar not parsed yet). */
  public int getTemplateOverlapStart() {
    return mTemplateOverlapStart;
  }

  /** @return The first template position after the CG overlap ends (-1 if no overlap or cigar not parsed yet). */
  public int getTemplateOverlapEnd() {
    return mTemplateOverlapEnd;
  }

  /**
   * Parses and processes the cigar, calling the <code>doChunk</code> method
   * for each chunk of identical commands, which in turn calls the lower level
   * <code>do...(...)</code> methods for each individual command.
   * @throws BadSuperCigarException on bad cigar or read delta
   */
  public void parse() throws BadSuperCigarException {
    if (mIsUnmapped) {
      doChunk(SamUtils.CIGAR_UNMAPPED, mRead.length);
      return;
    }
    mTemplatePos = mTemplateStart;
    mReadPos = 0;
    mReadDeltaPos = 0;
    /* The current cigar action count */
    int cigarCount = 0;
    for (int i = 0; i < mCigar.length(); i++) {
      final char ch = mCigar.charAt(i);
      if (ch >= '0' && ch <= '9') {
        cigarCount = cigarCount * 10 + (ch - '0');
      } else {
        doChunk(ch, cigarCount);
        cigarCount = 0;
      }
    }
    if (mReadOverlapStart != -1 && mReadOverlapEnd == -1) {
      mReadOverlapEnd = mReadPos;
    }
    assert cigarCount == 0; // we should not have digits at the end.
  }

  protected byte getReadDelta(int pos) throws BadSuperCigarException {
    if (mRead != null) {
      if (mReadPos >= mReadLength) {
        throw new BadSuperCigarException("Ill-formed cigar: " + mCigar);
      }
      return mRead[mReadPos];
    } else {
      if (pos >= mReadDelta.length) {
        throw new BadSuperCigarException("Ill-formed cigar/read delta: " + mCigar + "/" + mReadDeltaString);
      }
      return mReadDelta[pos];
    }
  }

  protected byte templateNt() {
    return mTemplate != null && 0 <= mTemplatePos && mTemplatePos < mTemplateLength ? mTemplate[mTemplatePos] : DnaUtils.UNKNOWN_RESIDUE;
  }

  /**
   * Process one chunk of identical commands.
   * This method can be overridden by subclasses, but the overriding method MUST
   * call this super method to correctly increment the read and template positions.
   *
   * Note that Ns are treated as matches by this method, any subclasses need to handle Ns in the <code>doEquality</code> method.
   *
   * @param ch the command to process
   * @param count the number of times the command is repeated.
   * @throws BadSuperCigarException on bad cigar or read delta
   */
  protected void doChunk(char ch, int count) throws BadSuperCigarException {
    for (int i = 0; i < count; i++) {
      mOperationLength = count;
      final int templateNt = templateNt();
      mTemplateToReadPos[mTemplatePos & TABLE_MASK] = mReadPos;
      switch (ch) {
        case SamUtils.CIGAR_UNMAPPED:
          doUnmapped();
          advanceTemplate(true);
          mReadPos++;
          break;
        case SamUtils.CIGAR_SAME:
          final int readNt = mRead == null ? templateNt : getReadDelta(-1);
          doEquality(readNt, templateNt);
          advanceTemplate(true);
          mReadPos++;
          break;
        case SamUtils.CIGAR_MISMATCH:
          doSubstitution(getReadDelta(mReadDeltaPos), templateNt);
          mReadDeltaPos++;
          mReadPos++;
          advanceTemplate(true);
          break;
        case SamUtils.CIGAR_SAME_OR_MISMATCH:
          if (mTemplate == null) {
            throw new ReferenceSequenceRequiredException("No template has been set and legacy cigars are used.");
          }
          if (mTemplatePos < 0 || mTemplatePos >= mTemplateLength) {
            throw new BadSuperCigarException("bad cigar start=" + mTemplateStart + ", pos=" + mTemplatePos + ", len=" + mTemplateLength + ", cigar " + mCigar);
          }
          final byte readDeltaNt = getReadDelta(mReadDeltaPos);
          assert readDeltaNt >= DNA.N.ordinal() && readDeltaNt <= DNA.T.ordinal() : "ReadDelta out of 0-4 range: " + readDeltaNt; //make sure we're dealing with everything in the same number space.
          if (N_CODE != readDeltaNt && N_CODE != templateNt && readDeltaNt != templateNt) {
            doSubstitution(readDeltaNt, templateNt);
            mReadDeltaPos++;
          } else {
            doEquality(readDeltaNt, templateNt);
          }
          mReadPos++;
          advanceTemplate(true);
          break;
        case SamUtils.CIGAR_DELETION_FROM_REF:
          doTemplateOnly(templateNt);
          advanceTemplate(true);
          break;
        case SamUtils.CIGAR_INSERTION_INTO_REF:
          doReadOnly(getReadDelta(mReadDeltaPos)); // out of bounds here means readDelta is too short.
          mReadDeltaPos++;
          mReadPos++;
          break;
        case SamUtils.CIGAR_GAP_IN_READ:
          doTemplateSkip(templateNt);
          advanceTemplate(true);
          break;
        case SuperCigar.OVERLAP_IN_READ:
          if (mTemplateOverlapEnd == -1) {
            mTemplateOverlapStart = mTemplatePos - count;
            mTemplateOverlapEnd = mTemplatePos;
            mReadOverlapStart = mTemplateToReadPos[mTemplateOverlapStart & TABLE_MASK];
            mReadOverlapMiddle = mReadPos;
            // mReadOverlapEnd will be set below, as soon as we have advanced far enough along the read.
          }
          doTemplateOverlap();
          if (mTemplatePos > 0) { //this can reverse past the template start position if you have a large overlap and some inserts
            advanceTemplate(false);
          }
          break;
        case SamUtils.CIGAR_SOFT_CLIP:
          doReadSoftClip(getReadDelta(mReadDeltaPos));
          mReadDeltaPos++;
          mReadPos++;
          if (mTemplatePos > mTemplateStart) {  //this only advances the template for soft clips that aren't at the start of the read. Advancing for other cases is purely for super cigar backstep operations.
            advanceTemplate(true);
          }
          break;
        case 'H':
          doReadHardClip();
          // hard clip regions are not included in the read, so we do not increment mReadPos.
          if (mTemplatePos > mTemplateStart) {  //this only advances the template for soft clips that aren't at the start of the read. Advancing for other cases is purely for super cigar backstep operations.
            advanceTemplate(true);
          }
          break;
        case SuperCigar.UNKNOWN_READ:
          doUnknownOnRead();
          advanceTemplate(true);
          mReadPos++;
          break;
        case SuperCigar.UNKNOWN_TEMPLATE:
          doUnknownOnTemplate(getReadDelta(mReadDeltaPos), templateNt);
          mReadDeltaPos++;
          mReadPos++;
          advanceTemplate(true);
          break;
        default:
          throw new BadSuperCigarException("Illegal command in supercigar: " + mCigar);
      }
      if (mTemplatePos == mTemplateOverlapEnd
          && mReadOverlapEnd == -1
          && mReadOverlapStart != -1) {
        mReadOverlapEnd = mReadPos; // the first read position after the overlap
      }
    }
  }

  protected void advanceTemplate(boolean direction) {
    mTemplatePos = mTemplatePos + (direction ? 1 : -1);
  }

  /**
   * Process one hard-clipped nucleotide on the read.
   */
  protected void doReadHardClip() {
  }

  /**
   * Process one soft-clipped nucleotide on the read.
   * @param readNt the current read base
   */
  protected void doReadSoftClip(int readNt) {
  }

  /**
   * Skip backwards on the template by one position.
   * These commands appear only in the middle of a CG overlap region.
   */
  protected void doTemplateOverlap() {
  }

  /**
   * Skip over one nucleotide on the template.
   * These commands are used for the CG large and small gaps.
   * @param templateNt the current template nucleotide
   */
  protected void doTemplateSkip(int templateNt) {
  }

  /**
   * Process one 'I' command, which means that a nucleotide appears
   * in the read only, and is missing from the template.
   * @param readNt the current read nucleotide
   * @throws BadSuperCigarException when a bad cigar is encountered
   */
  protected void doReadOnly(int readNt) throws BadSuperCigarException {
  }

  /**
   * Process one 'D' command, which means that a nucleotide appears
   * in the template only, and is missing from the read.
   * @param templateNt the current template nucleotide
   * @throws BadSuperCigarException when a bad cigar is encountered
   */
  protected void doTemplateOnly(int templateNt) throws BadSuperCigarException {
  }

  /**
   * Process one mismatch.
   * @param readNt the current read base
   * @param templateNt the current template base
   * @throws BadSuperCigarException when a bad cigar is encountered
   */
  protected void doSubstitution(int readNt, int templateNt) throws BadSuperCigarException {
  }

  /**
   * Process one equality command.
   * @param readNt the current read base
   * @param nt nucleotide of the read and the template.
   * @throws BadSuperCigarException when a bad cigar is encountered
   */
  protected void doEquality(int readNt, int nt) throws BadSuperCigarException {
  }

  /**
   * Process one unmapped command.
   */
  protected void doUnmapped() {
  }

  /**
   * Process one Unknown nucleotide in the template
   * @param readNt the current read base
   * @param templateNt the current template base
   * @throws BadSuperCigarException when a bad cigar is encountered
   */
  protected void doUnknownOnTemplate(int readNt, int templateNt) throws BadSuperCigarException {  }

  /**
   * Process one Unknown nucleotide in the read
   * @throws BadSuperCigarException when a bad cigar is encountered
   */
  protected void doUnknownOnRead() throws BadSuperCigarException {  }
}
