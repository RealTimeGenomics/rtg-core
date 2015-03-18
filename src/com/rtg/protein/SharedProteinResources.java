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

import java.io.IOException;

import com.rtg.mode.Frame;
import com.rtg.mode.ProteinScoringMatrix;
import com.rtg.mode.SequenceType;
import com.rtg.reader.PrereadNamesInterface;
import com.rtg.reader.SequencesReader;

/**
 * Hold information shared between all <code>ProteinAlignmentResult</code>
 * objects.
 */
public class SharedProteinResources {

  private final ProteinScoringMatrix mProteinScoringMatrix;
  private final SequencesReader mTemplateReader;
  private final SequencesReader mQueryReader;
  private final byte[] mWorkSpace;
  private final byte[] mWorkSpaceProtein;


  private final PrereadNamesInterface mTemplateNames;
  private final PrereadNamesInterface mReadNames;

  SharedProteinResources(final ProteinScoringMatrix matrix, final SequencesReader template, final SequencesReader query, boolean readNames) throws IOException {
    mProteinScoringMatrix = matrix;
    mTemplateReader = template;

    mQueryReader = query;
    if (query != null) {
      mWorkSpace = new byte[(int) query.maxLength()];
      mWorkSpaceProtein = new byte[mWorkSpace.length / 3];
      if (query.numberSequences() > Integer.MAX_VALUE) {
        throw new IllegalArgumentException("Number of reads exceeds " + Integer.MAX_VALUE);
      }
    } else {
      mWorkSpace = null;
      mWorkSpaceProtein = null;
    }
    if (mTemplateReader != null) {
      mTemplateNames = mTemplateReader.names();
    } else {
      mTemplateNames = null;
    }
    if (readNames && mQueryReader != null) {
      mReadNames = query.names();
    } else {
      mReadNames = null;
    }
  }

  SequencesReader templateReader() {
    return mTemplateReader;
  }

  SequencesReader queryReader() {
    return mQueryReader;
  }

  ProteinScoringMatrix proteinScoringMatrix() {
    return mProteinScoringMatrix;
  }

  int templateLength(final int id) throws IOException {
    if (mTemplateReader == null) {
      return -1;
    } else {
      return mTemplateReader.length(id);
    }
  }

  PrereadNamesInterface templateNames() {
    return mTemplateNames;
  }

  PrereadNamesInterface readNames() {
    return mReadNames;
  }

  int queryLength(final int id) throws IOException {
    if (mQueryReader == null) {
      return -1;
    } else {
      return mQueryReader.length(id);
    }
  }

  long totalTemplateLength() {
    if (mTemplateReader == null) {
      return -1;
    } else {
      return mTemplateReader.totalLength();
    }
  }

  private static final byte[] EMPTY = new byte[0];

  /**
   * Return the query, resulting array may contain rubbish entries at the end.
   * Assumes the caller knows what it is doing!
   * @param id read id
   * @param frame read frame
   * @return amino acid sequence
   * @exception IOException if an I/O error occurs.
   */
  byte[] query(final int id, final int frame) throws IOException {
    if (mQueryReader == null) {
      return EMPTY;
    } else if (mQueryReader.type() == SequenceType.DNA) {
      // In this situation, use precomputed array to store intermediate nt version
      final int length = mQueryReader.read(id, mWorkSpace);
      // Convert to protein in situ
      final Frame frames = ProteinOutputProcessor.FRAMES_MAPPING[frame + 3];
      final int limit = (length - Math.abs(frame) + 1) / 3;
      for (int j = 0, i = 0; j < limit; j++, i += 3) {
        mWorkSpaceProtein[j] = frames.code(mWorkSpace, length, i);
      }
      return mWorkSpaceProtein;
    } else {
      return mQueryReader.read(id);
    }
  }

  private int mCachedTemplateId = -1;
  private byte[] mCachedTemplate = null;

  byte[] template(final int id) throws IOException {
    if (id != mCachedTemplateId) {
      if (mTemplateReader == null) {
        mCachedTemplate = EMPTY;
      } else {
        mCachedTemplate = mTemplateReader.read(id);
      }
      mCachedTemplateId = id;
    }
    return mCachedTemplate;
  }
}
