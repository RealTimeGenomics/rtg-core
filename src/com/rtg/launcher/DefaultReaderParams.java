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
package com.rtg.launcher;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import com.rtg.mode.SequenceMode;
import com.rtg.reader.IndexFile;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.Integrity;

/**
 */
public class DefaultReaderParams extends ReaderParams implements Integrity {

  private final File mSequenceDir;

  private final SequenceMode mMode;

  private long mMaxLength = -1;

  private boolean mUseMemSeqReader;

  private boolean mLoadNames;

  private boolean mLoadFullNames;

  transient SequencesReader mReader = null;

  transient int[] mLengths = null;

  private final ReaderParams mParent;

  private final LongRange mReaderRestriction;

  private LongRange mAdjustedRegion;

  /**
   * @param sequenceDir directory containing preread sequences.
   * @param mode the sequence mode.
   * @param readerRestriction a region specifying the subset of the SDF to load
   * @param parent parent, can be null
   * @param useMemSeqReader use {@link com.rtg.reader.CompressedMemorySequencesReader}
   * @param loadNames whether to load names from disk or not
   * @param loadFullNames whether to load full names from disk or not
   */
  public DefaultReaderParams(final File sequenceDir, LongRange readerRestriction, final SequenceMode mode, final ReaderParams parent, final boolean useMemSeqReader, final boolean loadNames, boolean loadFullNames) {
    mSequenceDir = sequenceDir;
    mMode = mode;
    mParent = parent;
    mUseMemSeqReader = useMemSeqReader;
    mLoadNames = loadNames;
    mLoadFullNames = loadFullNames;
    mReaderRestriction = readerRestriction;
  }

  /**
   * Create a reader not backed by a directory.
   * @param reader the sequence data
   * @param readerRestriction a region specifying the subset of the SDF to load
   * @param mode the sequence mode.
   */
  public DefaultReaderParams(SequencesReader reader, LongRange readerRestriction, SequenceMode mode) {
    mReader = reader;
    mMode = mode;
    mReaderRestriction = readerRestriction;
    mParent = null;
    mUseMemSeqReader = true;
    mLoadNames = true;
    mLoadFullNames = true;
    mSequenceDir = reader.path();
  }

  /**
   * Get the mode.
   * @return the mode.
   */
  @Override
  public SequenceMode mode() {
    return mMode;
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
        final LongRange initialRange = mReaderRestriction;
        mAdjustedRegion = SequencesReaderFactory.resolveRange(index, initialRange);
        if (mAdjustedRegion.getLength() > Integer.MAX_VALUE) {
          throw new NoTalkbackSlimException(ErrorType.INFO_ERROR, "There are too many sequences in the specified range of \"" + mSequenceDir + "\", "
                  + "provide a range with at most " + Integer.MAX_VALUE + " sequences");
        }
        mReader = mUseMemSeqReader
            ? SequencesReaderFactory.createMemorySequencesReaderCheckEmpty(mSequenceDir, mLoadNames, mLoadFullNames, mAdjustedRegion)
            : SequencesReaderFactory.createDefaultSequencesReaderCheckEmpty(mSequenceDir, mAdjustedRegion);
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
    return "SequenceParams mode=" + mMode + " directory=" + mSequenceDir + " usememory=" + mUseMemSeqReader;
  }

  @Override
  public int hashCode() {
    return Utils.hash(new Object[] {mMode, mSequenceDir, mUseMemSeqReader});
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
    return this.mMode.equals(that.mMode) && this.mSequenceDir.equals(that.mSequenceDir) && this.mUseMemSeqReader == that.mUseMemSeqReader;
  }

  @Override
  public boolean integrity() {
    //be careful here not to lazily open the reader
    if (mReader != null) {
      if (reader().type() != mMode.type()) {
        throw new NoTalkbackSlimException(ErrorType.SDF_FILETYPE_ERROR, reader().type().toString(), mMode.type().toString());
      }
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

  @Override
  public LongRange adjustedRegion() {
    reader(); // To force the reader creation
    return mAdjustedRegion;
  }
}
