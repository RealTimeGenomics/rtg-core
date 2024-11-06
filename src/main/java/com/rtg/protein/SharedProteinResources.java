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
package com.rtg.protein;

import java.io.IOException;

import com.rtg.mode.Frame;
import com.rtg.mode.ProteinScoringMatrix;
import com.rtg.mode.SequenceType;
import com.rtg.reader.NamesInterface;
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


  private final NamesInterface mTemplateNames;
  private final NamesInterface mReadNames;

  SharedProteinResources(final ProteinScoringMatrix matrix, final SequencesReader template, final SequencesReader query, boolean readNames) throws IOException {
    if (template == null || query == null) {
      throw new NullPointerException();
    }
    if (query.numberSequences() > Integer.MAX_VALUE) {
      throw new IllegalArgumentException("Number of reads exceeds " + Integer.MAX_VALUE);
    }
    assert template.type() == SequenceType.PROTEIN;
    assert query.type() == SequenceType.PROTEIN || query.type() == SequenceType.DNA;

    mProteinScoringMatrix = matrix;
    mTemplateReader = template;
    mQueryReader = query;
    if (mQueryReader.type() == SequenceType.DNA) {
      mWorkSpace = new byte[(int) query.maxLength()];
      mWorkSpaceProtein = new byte[mWorkSpace.length / 3];
    } else {
      mWorkSpace = null;
      mWorkSpaceProtein = null;
    }
    mTemplateNames = mTemplateReader.names();
    mReadNames = readNames ? query.names() : null;
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
    return mTemplateReader.length(id);
  }

  NamesInterface templateNames() {
    return mTemplateNames;
  }

  NamesInterface readNames() {
    return mReadNames;
  }

  int queryLength(final int id) throws IOException {
    return mQueryReader.length(id);
  }

  long totalTemplateLength() {
    return mTemplateReader.totalLength();
  }

  /**
   * Return the query, resulting array may contain rubbish entries at the end.
   * Assumes the caller knows what it is doing!
   * @param id read id
   * @param frame read frame
   * @return amino acid sequence
   * @exception IOException if an I/O error occurs.
   */
  byte[] query(final int id, final int frame) throws IOException {
    if (mQueryReader.type() == SequenceType.DNA) {
      // In this situation, use precomputed array to store intermediate nt version
      final int length = mQueryReader.read(id, mWorkSpace);
      // Convert to protein in situ
      final Frame frames = ProteinOutputProcessor.FRAMES_MAPPING[frame + 3];
      final int limit = (length - Math.abs(frame) + 1) / 3;
      for (int j = 0, i = 0; j < limit; ++j, i += 3) {
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
      mCachedTemplate = mTemplateReader.read(id);
      mCachedTemplateId = id;
    }
    return mCachedTemplate;
  }
}
