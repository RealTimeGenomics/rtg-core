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
package com.rtg.protein;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.alignment.ActionsHelper;
import com.rtg.reader.PrereadNamesInterface;

/**
 * Hold a single protein alignment result.
 */
class ProteinAlignmentResult implements Comparable<ProteinAlignmentResult> {

  private static final int FRAME_COUNT = 6;
  private static final int NUCLEOTIDES_PER_CODON = 3;

  /** System dependent representation of null. */
  private static final byte[] NULL = ("" + null).getBytes();

  private static final byte[][] FRAME_BYTES = {
    "-3".getBytes(),
    "-2".getBytes(),
    "-1".getBytes(),
    NULL,
    "+1".getBytes(),
    "+2".getBytes(),
    "+3".getBytes(),
  };

  private final SharedProteinResources mResources;
  private final int mTemplateId;
  private final int mReadAndFrame;
  private final long mReadIdOffset;
  private final int[] mActions;
  private final boolean mOutputProteinSequences;

  ProteinAlignmentResult(final SharedProteinResources resources, final int templateId, final int readId, final int[] res, long offset) {
    this(resources, templateId, readId, res, offset, true);
  }

  ProteinAlignmentResult(final SharedProteinResources resources, final int templateId, final int readId, final int[] res, long offset, boolean outputProteinSequences) {
    mOutputProteinSequences = outputProteinSequences;
    mResources = resources;
    mTemplateId = templateId;
    mReadAndFrame = readId;
    mReadIdOffset = offset;
    if (res != null) {
      mActions = ActionsHelper.copy(res);
    } else {
      mActions = new int[ActionsHelper.ACTIONS_START_INDEX];
    }
  }

  int readId() {
    return mReadAndFrame / FRAME_COUNT;
  }

  private int frame() {
    return ProteinOutputProcessor.INTERNAL_ENCODED_FRAME_TO_NATURAL_FRAME[mReadAndFrame % FRAME_COUNT];
  }

  int alignmentScore() {
    return ActionsHelper.alignmentScore(mActions);
  }

  int getAlignmentStart() {
    return ActionsHelper.zeroBasedTemplateStart(mActions);
  }

  // Compares only the template position and read id/frame, not the actual alignment
  boolean equals(int templateId, int templateStart, int readAndFrame) {
    return templateId == mTemplateId && readAndFrame == mReadAndFrame && templateStart == ActionsHelper.zeroBasedTemplateStart(mActions);
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj) {
      return true;
    }
    if (!(obj instanceof ProteinAlignmentResult)) {
      return false;
    }
    final ProteinAlignmentResult other = (ProteinAlignmentResult) obj;
    return equals(other.mTemplateId, ActionsHelper.zeroBasedTemplateStart(other.mActions), other.mReadAndFrame);
  }

  @Override
  public int hashCode() {
    final long x = mTemplateId + mReadAndFrame + ActionsHelper.zeroBasedTemplateStart(mActions);
    return (int) x;
  }

  @Override
  public int compareTo(ProteinAlignmentResult other) {
    if (mTemplateId != other.mTemplateId) {
      return mTemplateId - other.mTemplateId;
    }
    final long start = ActionsHelper.zeroBasedTemplateStart(mActions);
    final long oStart = ActionsHelper.zeroBasedTemplateStart(other.mActions);
    if (start > oStart) {
      return 1;
    } else if (start < oStart) {
      return -1;
    }
    return mReadAndFrame - other.mReadAndFrame;
  }

  /**
   * Given a frame in the form +1, +2, +3, -1, -2, -3, and length of the read,
   * return the 1-based index of the leftmost nucleotide in the translated protein.
   *
   * @param frame current frame
   * @param ntLength length of read in nucleotide space
   * @return start position in nucleotide spaces
   */
  static int readNtStart(final int frame, final int ntLength) {
    return frame > 0 ? frame : 1 + (ntLength + frame - 2) % NUCLEOTIDES_PER_CODON;
  }

  /**
   * Given a 1-based start offset and length of a nucleotide sequence, compute
   * the index of the last nucleotide used in a conversion to protein.
   *
   * @param start start position
   * @param ntLength length of read in nucleotide space
   * @return end position in nucleotide space
   */
  static int readNtEnd(final int start, final int ntLength) {
    final int effectiveLength = ntLength - start + 1;
    return start + (effectiveLength / NUCLEOTIDES_PER_CODON) * NUCLEOTIDES_PER_CODON - 1;
  }

  private int templateEnd() {
    return ActionsHelper.zeroBasedTemplateStart(mActions) + ActionsHelper.actionsCount(mActions) - ActionsHelper.deletionFromReadAndOverlapCount(mActions);
  }

  static void writeBitScore(final OutputStream os, final double sc) throws IOException {
    // write nnnnn.n
    double x;
    if (sc < 0) {
      os.write('-');
      x = -sc;
    } else {
      x = sc;
    }
    x += 0.05;
    os.write(String.valueOf((int) x).getBytes());
    os.write('.');
    os.write((int) (x * 10) % 10 + '0');
  }

  static void writeEScore(final OutputStream os, final double sc) throws IOException {
    // If near 0 write 0
    // otherwise write in form d.den where d is a digit and n an integer
    if (sc < 1E-180) {
      os.write('0');
    } else {
      int exponent = (int) Math.log10(sc) - 1;
      double x = sc * Math.pow(0.1, exponent);
      x += 0.05; // rounding
      if (x > 10) {
        exponent++;
        x *= 0.1;
      }
      final int leadingDigit = (int) x;
      x -= leadingDigit;
      x *= 10; // retrieve first significant digit
      os.write(leadingDigit + '0');
      os.write('.');
      os.write((int) x + '0');
      os.write('e');
      os.write(String.valueOf(exponent).getBytes());
    }
  }

  private void writeTemplateName(final OutputStream os) throws IOException {
    final PrereadNamesInterface tnames = mResources.templateNames();
    if (tnames != null) {
      tnames.writeName(os, mTemplateId);
    } else {
      os.write(NULL);
    }
  }

  private void writeReadName(OutputStream os, int readId) throws IOException  {
    final PrereadNamesInterface rnames = mResources.readNames();
    if (rnames != null) {
      rnames.writeName(os, readId);
    } else {
      os.write(String.valueOf(readId + mReadIdOffset).getBytes());
    }
  }

  private final ByteArrayOutputStream mNullOutputStream = new ByteArrayOutputStream();

  /**
   * Write a representation of this alignment result to the given stream.
   * @param os output stream
   * @exception IOException if an I/O error occurs
   */
  public void write(final OutputStream os) throws IOException {
    final int readId = readId();
    final int frame = frame();
    final byte[] read = mResources.query(readId, frame);
    final byte[] template = mResources.template(mTemplateId);
    final int rawScore = ActionsHelper.alignmentScore(mActions);
    final int length = ActionsHelper.actionsCount(mActions);
    final int templateStart = ActionsHelper.zeroBasedTemplateStart(mActions);
    writeTemplateName(os);
    os.write('\t');
    os.write(FRAME_BYTES[frame + 3]);
    os.write('\t');
    writeReadName(os, readId);
    os.write('\t');
    os.write(String.valueOf(templateStart + 1).getBytes());
    os.write('\t');
    os.write(String.valueOf(templateEnd()).getBytes());
    os.write('\t');
    os.write(String.valueOf(mResources.templateLength(mTemplateId)).getBytes());
    os.write('\t');
    //    final int readStart = readStart(frame, mResources.queryLength(readId), length * NUCLEOTIDES_PER_CODON);
    final int ntLength = mResources.queryLength(readId);
    final int readStart = readNtStart(frame, ntLength);
    os.write(String.valueOf(readStart).getBytes());
    os.write('\t');
    os.write(String.valueOf(readNtEnd(readStart, ntLength)).getBytes());
    os.write('\t');
    os.write(String.valueOf(mResources.queryLength(readId)).getBytes());

    final int positive;
    if (mOutputProteinSequences) {
      os.write('\t');
      ActionsHelper.writeProtein(os, template, templateStart, mActions, ActionsHelper.INSERTION_INTO_REFERENCE);
      os.write('\t');
      ActionsHelper.writeProtein(os, read, 0, mActions, ActionsHelper.DELETION_FROM_REFERENCE);
      os.write('\t');
      positive = ActionsHelper.writeProteinMatches(os, read, template, mActions, mResources.proteinScoringMatrix());
    } else {
      positive = ActionsHelper.writeProteinMatches(mNullOutputStream, read, template, mActions, mResources.proteinScoringMatrix());
      mNullOutputStream.reset();
    }

    os.write('\t');
    final int matches = ActionsHelper.matchCount(mActions);
    os.write(String.valueOf(matches).getBytes());
    os.write('\t');
    os.write(String.valueOf(ScoringHelper.percentage(matches, length)).getBytes());
    os.write('\t');
    os.write(String.valueOf(positive).getBytes());
    os.write('\t');
    os.write(String.valueOf(ScoringHelper.percentage(positive, length)).getBytes());
    os.write('\t');
    os.write(String.valueOf(length - positive).getBytes());
    os.write('\t');
    os.write(String.valueOf(rawScore).getBytes());
    os.write('\t');
    writeBitScore(os, ScoringHelper.computeBitScore(rawScore, mResources.proteinScoringMatrix()));
    os.write('\t');
    writeEScore(os, ScoringHelper.computeEScore(rawScore, length, mResources.totalTemplateLength(), mResources.proteinScoringMatrix()));
  }
}
