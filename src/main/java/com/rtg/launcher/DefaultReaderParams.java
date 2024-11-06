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
package com.rtg.launcher;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import com.rtg.reader.IndexFile;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.Integrity;
import com.rtg.util.intervals.LongRange;

/**
 */
public class DefaultReaderParams extends ReaderParams implements Integrity {

  private final File mSequenceDir;

  private long mMaxLength = -1;

  private final boolean mUseMemSeqReader;

  private final boolean mLoadNames;

  private final boolean mLoadFullNames;

  transient SequencesReader mReader = null;

  transient int[] mLengths = null;

  private final ReaderParams mParent;

  private final LongRange mReaderRestriction;

  /**
   * @param sequenceDir directory containing sequences.
   * @param readerRestriction a region specifying the subset of the SDF to load
   * @param useMemSeqReader use {@link com.rtg.reader.CompressedMemorySequencesReader}
   * @param loadNames whether to load names from disk or not
   * @param loadFullNames whether to load full names from disk or not
   */
  public DefaultReaderParams(final File sequenceDir, LongRange readerRestriction, final boolean useMemSeqReader, final boolean loadNames, boolean loadFullNames) {
    mSequenceDir = sequenceDir;
    mParent = null;
    mUseMemSeqReader = useMemSeqReader;
    mLoadNames = loadNames;
    mLoadFullNames = loadFullNames;
    mReaderRestriction = readerRestriction;
  }

  /**
   * Create a ReaderParams corresponding to an already created ReaderParams, useful for creating a copy.
   * @param parent parent, can be null
   */
  DefaultReaderParams(final ReaderParams parent) {
    mParent = parent;
    mSequenceDir = parent.directory();
    // Rest of these make no sense in this case
    mReaderRestriction = null;
    mUseMemSeqReader = true;
    mLoadNames = true;
    mLoadFullNames = true;
  }

  /**
   * Create a ReaderParams corresponding to an already opened SequencesReader.
   * @param reader the sequence data
   */
  public DefaultReaderParams(SequencesReader reader) {
    mReader = reader;
    mSequenceDir = reader.path();
    // Rest of these make no sense in this case
    mReaderRestriction = null;
    mParent = null;
    mUseMemSeqReader = true;
    mLoadNames = true;
    mLoadFullNames = true;
  }

  /**
   * Get the directory.
   * @return the directory.
   */
  @Override
  public File directory() {
    return mSequenceDir;
  }

  /**
   * If necessary close the reader.
   * @throws IOException If an IO error occurs
   */
  @Override
  public void close() throws IOException {
    if (mReader == null) {
      return;
    }
    mReader.close();
    mReader = null;
    mLengths = null;
  }

  @Override
  public SequencesReader reader() {
    if (mReader == null) {
      if (mParent != null) {
        mReader = mParent.reader().copy();
        return mReader;
      }
      try {
        final IndexFile index = new IndexFile(mSequenceDir);
        final LongRange adjustedRegion = SequencesReaderFactory.resolveRange(index, mReaderRestriction);
        if (adjustedRegion.getLength() > Integer.MAX_VALUE) {
          throw new NoTalkbackSlimException(ErrorType.INFO_ERROR, "There are too many sequences in the specified range of \"" + mSequenceDir + "\", "
                  + "provide a range with at most " + Integer.MAX_VALUE + " sequences");
        }
        mReader = mUseMemSeqReader
            ? SequencesReaderFactory.createMemorySequencesReaderCheckEmpty(mSequenceDir, mLoadNames, mLoadFullNames, adjustedRegion)
            : SequencesReaderFactory.createDefaultSequencesReaderCheckEmpty(mSequenceDir, adjustedRegion);
      } catch (final FileNotFoundException e) {
        if (mSequenceDir.isDirectory()) {
          throw new NoTalkbackSlimException(ErrorType.SDF_INDEX_NOT_VALID, mSequenceDir.getPath());
        } else if (mSequenceDir.exists()) {
          throw new NoTalkbackSlimException(ErrorType.INFO_ERROR, "The specified file, \"" + mSequenceDir.getPath() + "\", is not an SDF.");
        } else {
          throw new NoTalkbackSlimException(ErrorType.INFO_ERROR, "The specified SDF, \"" + mSequenceDir.getPath() + "\", does not exist.");
        }
      } catch (final IOException e) {
        throw new NoTalkbackSlimException(ErrorType.SDF_INDEX_NOT_VALID, mSequenceDir.getPath());
      }
    }
    return mReader;
  }

  @Override
  public int[] lengths() throws IOException {

    //lengths are discovered lazily
    if (mLengths == null) {
      mLengths = this.reader().sequenceLengths(0, this.reader().numberSequences());
    }
    return mLengths;

  }

  /**
   * Get the length of the longest sequence.
   * @return the length of the longest sequence.
   */
  @Override
  public long maxLength() {
    if (mMaxLength == -1) {
      //being lazy makes tests much easier
      mMaxLength = reader().maxLength();
    }
    return mMaxLength;
  }

  @Override
  public String toString() {
    return "ReaderParams directory=" + mSequenceDir + " usememory=" + mUseMemSeqReader;
  }

  @Override
  public int hashCode() {
    return Utils.hash(new Object[] {mSequenceDir, mUseMemSeqReader});
  }

  @Override
  public boolean equals(final Object obj) {
    if (this == obj) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    final DefaultReaderParams that = (DefaultReaderParams) obj;
    return this.mSequenceDir.equals(that.mSequenceDir) && this.mUseMemSeqReader == that.mUseMemSeqReader;
  }

  @Override
  public boolean integrity() {
    //be careful here not to lazily open the reader
    if (mReader != null) {
      Exam.assertTrue(mMaxLength == -1 || mMaxLength == reader().maxLength());
      if (mLengths != null) {
        Exam.assertTrue(mLengths.length == mReader.numberSequences());
      }
    } else {
      Exam.assertTrue(mLengths == null);
      Exam.assertTrue(mMaxLength >= -1);
    }
    return true;
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    return true;
  }

}
