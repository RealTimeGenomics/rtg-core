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
import java.io.IOException;

import com.rtg.mode.SequenceMode;
import com.rtg.reader.SequencesReader;
import com.rtg.reference.Sex;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.Integrity;
import com.rtg.util.intervals.LongRange;

/**
 * Encapsulates parameters relating to how a SequencesReader should be obtained and how it will later be accessed.
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

  private final SequenceMode mMode;

  private final ReaderParams mReaderParams;

  private final Sex mSex;
  private final HashingRegion mRegion;
  private final LongRange mReaderRestriction;

  private SequenceParams(SequenceParamsBuilder builder) {
    mMode = builder.mMode;
    mSex = builder.mSex;
    mReaderRestriction = builder.mReaderRestriction;
    if (builder.mReaderParams != null) {
      mReaderParams = new DefaultReaderParams(builder.mReaderParams);
    } else {
      final boolean useMemoryReader = builder.mUseMemReader;
      final boolean loadNames = builder.mLoadNames;
      final boolean loadFullNames = builder.mLoadFullNames;
      mReaderParams = new DefaultReaderParams(builder.mSequenceDir, mReaderRestriction, useMemoryReader, loadNames, loadFullNames);
    }
    if (builder.mRegion == HashingRegion.NONE) {
      mRegion = new HashingRegion(0, mReaderParams.reader().numberSequences());
    } else {
      mRegion = builder.mRegion;
    }
    integrity();
  }

  @Override
  public ISequenceParams subSequence(HashingRegion region) {
     return new SequenceParams.SequenceParamsBuilder().directory(directory())
                                               .sex(mSex)
                                               .mode(mode())
                                               .region(region)
                                               .readerParam(mReaderParams)
                                               .create();
  }

  @Override
  public SequenceMode mode() {
    return mMode;
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
    return Utils.pairHash(Utils.pairHash(Utils.pairHash(mReaderParams.hashCode(), (int) mRegion.getStart()), (int) mRegion.getEnd()), mMode.hashCode());
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
    if (!this.mode().equals(that.mode())) {
      return false;
    }
    return true;
  }

  @Override
  public String toString() {
    return "SequenceParams mode=" + mode() + " region=" + mRegion + " directory=" + mReaderParams.directory()
    + (mSex != null ? " sex=" + mSex : "");
  }

  @Override
  public boolean integrity() {
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

