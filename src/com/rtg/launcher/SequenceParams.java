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
import java.io.IOException;

import com.rtg.mode.SequenceMode;
import com.rtg.reader.SequencesReader;
import com.rtg.reference.Sex;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.Integrity;

/**
 */
public final class SequenceParams implements ISequenceParams, Integrity {

  /**
   * Creates a SequenceParams builder.
   * @return the builder.
   */
  public static SequenceParamsBuilder builder() {
    return new SequenceParamsBuilder();
  }

  /**
   * A builder class for <code>SequenceParams</code>.
   */
  public static class SequenceParamsBuilder {

    protected File mSequenceDir;

    protected SequenceMode mMode = SequenceMode.BIDIRECTIONAL;

    protected HashingRegion mRegion = HashingRegion.NONE;

    protected LongRange mReaderRestriction = LongRange.NONE;

    protected ReaderParams mReaderParams = null;

    protected boolean mUseMemReader = false;

    protected boolean mLoadNames = false;

    protected boolean mLoadFullNames = false;

    protected Sex mSex = null;

    /**
     * Sets the directory containing the sequence source.
     *
     * @param sequenceDir a <code>File</code> value.
     * @return this builder, so calls can be chained.
     */
    public SequenceParamsBuilder directory(File sequenceDir) {
      mSequenceDir = sequenceDir;
      return this;
    }

    /**
     * Sets the sequence mode used when reading the sequence
     * source. The default value is
     * <code>SequenceMode.BIDIRECTIONAL</code>.
     *
     * @param mode the <code>SequenceMode</code> to use.
     * @return this builder, so calls can be chained.
     */
    public SequenceParamsBuilder mode(SequenceMode mode) {
      mMode = mode;
      return this;
    }

    /**
     * Sets the sex of the sequences to read.
     *
     * @param sex the sex for this params
     * @return this builder, so calls can be chained.
     */
    public SequenceParamsBuilder sex(Sex sex) {
      mSex = sex;
      return this;
    }

    /**
     * Sets the region of the sequences to process during search.
     *
     * @param region the region this params
     * @return this builder, so calls can be chained.
     */
    public SequenceParamsBuilder region(HashingRegion region) {
      mRegion = region;
      return this;
    }

    /**
     * Sets the region of the sequences to read.
     *
     * @param readerRestriction the region to restrict readers to
     * @return this builder, so calls can be chained.
     */
    public SequenceParamsBuilder readerRestriction(LongRange readerRestriction) {
      mReaderRestriction = readerRestriction;
      return this;
    }

    /**
     * Sets the parent {@link ReaderParams} object, null for new first instance
     *
     * @param param null or existing {@link ReaderParams}
     * @return this builder, so calls can be chained.
     */
    public SequenceParamsBuilder readerParam(ReaderParams param) {
      mReaderParams = param;
      return this;
    }

    /**
     * Sets if {@link ReaderParams} object should use {@link com.rtg.reader.CompressedMemorySequencesReader}
     *
     * @param useMemReader true if {@link com.rtg.reader.CompressedMemorySequencesReader} to be used
     * @return this builder, so calls can be chained.
     */
    public SequenceParamsBuilder useMemReader(boolean useMemReader) {
      mUseMemReader = useMemReader;
      return this;
    }

    /**
     * Sets if {@link ReaderParams} object should load sequence names
     *
     * @param loadNames true if names should be loaded
     * @return this builder, so calls can be chained.
     */
    public SequenceParamsBuilder loadNames(boolean loadNames) {
      mLoadNames = loadNames;
      return this;
    }

    /**
     * Sets if {@link ReaderParams} object should load full sequence names
     *
     * @param loadFullNames true if full names should be loaded
     * @return this builder, so calls can be chained.
     */
    public SequenceParamsBuilder loadFullNames(boolean loadFullNames) {
      mLoadFullNames = loadFullNames;
      return this;
    }

    /**
     * Creates a SequenceParams using the current builder
     * configuration.
     * @return the new SequenceParams
     */
    public SequenceParams create() {
      return new SequenceParams(this);
    }
  }

  private final File mSequenceDir;

  private final SequenceMode mMode;

  private final ReaderParams mReaderParams;

  private final Sex mSex;
  private final HashingRegion mRegion;
  private final LongRange mReaderRestriction;

  private SequenceParams(SequenceParamsBuilder builder) {
    mSequenceDir = builder.mSequenceDir;
    mMode = builder.mMode;
    final boolean useMemoryReader = builder.mUseMemReader;
    final boolean loadNames = builder.mLoadNames;
    final boolean loadFullNames = builder.mLoadFullNames;
    mSex = builder.mSex;
    mReaderRestriction = builder.mReaderRestriction;
    mReaderParams = new DefaultReaderParams(mSequenceDir, mReaderRestriction, mMode, builder.mReaderParams, useMemoryReader, loadNames, loadFullNames);
    if (builder.mRegion == HashingRegion.NONE) {
      mRegion = new HashingRegion(0, mReaderParams.reader().numberSequences());
    } else {
      mRegion = builder.mRegion;
    }
    integrity();
  }

  @Override
  public ISequenceParams subSequence(HashingRegion region) {
     return new SequenceParams.SequenceParamsBuilder().directory(mSequenceDir)
                                               .sex(mSex)
                                               .mode(mMode)
                                               .region(region)
                                               .readerParam(mReaderParams)
                                               .create();
  }

  @Override
  public SequenceMode mode() {
    return mReaderParams.mode();
  }

  @Override
  public File directory() {
    return mReaderParams.directory();
  }

  @Override
  public Sex sex() {
    return mSex;
  }

  @Override
  public HashingRegion region() {
    return mRegion;
  }

  @Override
  public LongRange readerRestriction() {
    return mReaderRestriction;
  }

  @Override
  public long numberSequences() {
    return mRegion.getEnd() - mRegion.getStart();
  }

  @Override
  public SequencesReader reader() {
    return mReaderParams.reader();
  }

  @Override
  public ReaderParams readerParams() {
    return mReaderParams;
  }

  @Override
  public long maxLength() {
    return mReaderParams.maxLength();
  }

  @Override
  public void close() throws IOException {
    mReaderParams.close();
  }

  @Override
  public int hashCode() {
    return Utils.pairHash(Utils.pairHash(mReaderParams.hashCode(), (int) mRegion.getStart()), (int) mRegion.getEnd());
  }

  @Override
  public boolean equals(final Object obj) {
    if (this == obj) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    final SequenceParams that = (SequenceParams) obj;
    if (!this.mReaderParams.equals(that.mReaderParams)) {
      return false;
    }
    if (!this.region().equals(that.region())) {
      return false;
    }
    return true;
  }

  @Override
  public String toString() {
    return "SequenceParams mode=" + mReaderParams.mode() + " region=" + mRegion.toString() + " directory=" + mReaderParams.directory()
    + (mSex != null ? " sex=" + mSex : "");
  }

  @Override
  public boolean integrity() {
    if (mReaderParams == null) {
      throw new NullPointerException();
    }
    Exam.assertTrue(0 <= mRegion.getStart() && mRegion.getStart() <= mRegion.getEnd());
    return true;
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    Exam.integrity(mReaderParams);
    return true;
  }

}

