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
package com.rtg.similarity;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.List;

import com.rtg.index.params.CountParams;
import com.rtg.index.params.CreateParams;
import com.rtg.launcher.BuildParams;
import com.rtg.launcher.ModuleParams;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.ProgramMode;
import com.rtg.util.Pair;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.Integrity;

/**
 * Holds all the parameters needed for doing a build and a search.
 */
public class BuildSearchParams extends ModuleParams implements Integrity {


  /**
   * Creates a BuildSearchParams builder.
   * @return the builder.
   */
  public static BuildSearchParamsBuilder builder() {
    return new BuildSearchParamsBuilder();
  }

  /**
   * A builder class for <code>BuildSearchParams</code>.
   */
  public static class BuildSearchParamsBuilder {

    protected String mName = "BuildSearchParams";

    protected ProgramMode mProgramMode;

    protected BuildParams mBuildParams;

    protected CreateParams mIndexParams;

    protected List<Pair<String, List<SequenceParams>>> mSequenceParams;

    protected CountParams mCountParams;

    protected boolean mUniqueWords;

    /**
     * Sets the program mode.
     * @param mode program mode.
     * @return this builder, so calls can be chained.
     */
    public BuildSearchParamsBuilder mode(ProgramMode mode) {
      mProgramMode = mode;
      return this;
    }

    /**
     * Count unique words only
     * @param uniqueWords true to count unique words, false to count all words
     * @return this builder, so calls can be chained.
     */
    public BuildSearchParamsBuilder uniqueWords(boolean uniqueWords) {
      mUniqueWords = uniqueWords;
      return this;
    }

    /**
     * Sets the parameters for doing build.
     * @param build the build parameters.
     * @return this builder, so calls can be chained.
     */
    public BuildSearchParamsBuilder build(BuildParams build) {
      mBuildParams = build;
      return this;
    }

    /**
     * Sets the parameters for doing create.
     * @param index the index parameters.
     * @return this builder, so calls can be chained.
     */
    public BuildSearchParamsBuilder index(CreateParams index) {
      mIndexParams = index;
      return this;
    }

    /**
     * Sets the parameters for doing builds.
     * @param sequences the sequences parameters.
     * @return this builder, so calls can be chained.
     */
    public BuildSearchParamsBuilder sequences(List<Pair<String, List<SequenceParams>>> sequences) {
      mSequenceParams = sequences;
      return this;
    }

    /**
     * Sets the parameters for creating hit count collector.
     * @param count the count parameters.
     * @return this builder, so calls can be chained.
     */
    public BuildSearchParamsBuilder count(CountParams count) {
      mCountParams = count;
      return this;
    }

    /**
     * Sets the application name.
     * @param name the application name.
     * @return this builder, so calls can be chained.
     */
    public BuildSearchParamsBuilder name(final String name) {
      mName = name;
      return this;
    }

    /**
     * Creates a BuildSearchParams using the current builder
     * configuration.
     * @return the new BuildSearchParams.
     */
    public BuildSearchParams create() {
      return new BuildSearchParams(this);
    }
  }

  private final ProgramMode mProgramMode;

  private final BuildParams mBuildParams;

  private final CreateParams mIndexParams;

  private final List<Pair<String, List<SequenceParams>>> mSequenceParams;

  private final CountParams mCountParams;

  private final boolean mUniqueWords;

  /**
   * Create a set of parameters to use from the builder.
   * @param builder the builder object.
   */
  public BuildSearchParams(final BuildSearchParamsBuilder builder) {
    super(builder.mName);
    mProgramMode = builder.mProgramMode;
    mBuildParams = builder.mBuildParams;
    mSequenceParams = builder.mSequenceParams;
    mCountParams = builder.mCountParams;
    mUniqueWords = builder.mUniqueWords;
    mIndexParams = builder.mIndexParams;
  }

  /**
   * @return the program mode.
   */
  public ProgramMode mode() {
    return mProgramMode;
  }

  /**
   * @return the build parameters.
   */
  public BuildParams build() {
    return mBuildParams;
  }

  /**
   * @return the create parameters
   */
  public CreateParams index() {
    return mIndexParams;
  }

  /**
   * @return the sequences parameters.
   */
  public List<Pair<String, List<SequenceParams>>> sequences() {
    return mSequenceParams;
  }

  /**
   * @return the hit count parameters.
   */
  public CountParams countParams() {
    return mCountParams;
  }

  /**
   *  Get the length required for the shared subject/query buffer.
   * @return the length required for the shared subject/query buffer.
   */
  public long bufferLength() {
    if (mBuildParams.directory() != null) {
      return mBuildParams.sequences().maxLength();
    } else {
      long ret = 0;
      for (Pair<String, List<SequenceParams>> mSequenceParam : mSequenceParams) {
        for (int j = 0; j < mSequenceParam.getB().size(); ++j) {
          ret = Math.max(ret, mSequenceParam.getB().get(j).maxLength());
        }
      }
      return ret;
    }
  }

  @Override
  public void close() throws IOException {
    mBuildParams.close();
  }

  @Override
  public int hashCode() {
    return Utils.pairHash(
        mCountParams.hashCode(),
        Utils.pairHash(
            mProgramMode.hashCode(),
            mBuildParams.hashCode()
            )
        );
  }

  @Override
  public boolean equals(final Object obj) {
    if (this == obj) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    final BuildSearchParams that = (BuildSearchParams) obj;
    return this.mProgramMode.equals(that.mProgramMode)
        && this.mBuildParams.equals(that.mBuildParams)
        && this.mCountParams.equals(that.mCountParams);
  }

  @Override
  public String toString() {
    final String linePrefix = com.rtg.util.StringUtils.LS + "..";
    return name() + " mode=" + mProgramMode + linePrefix
        + " hits={" + countParams() + "} " + linePrefix
        + " index={" + index() + "}" + linePrefix
        + " build={" + build() + "}"
        + com.rtg.util.StringUtils.LS;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mProgramMode != null);
    Exam.assertNotNull(mBuildParams);
    Exam.assertNotNull(mCountParams);
    if (mProgramMode != null) {
      Exam.assertEquals(mProgramMode.subjectMode(), mBuildParams.sequences().mode());
    }
    return true;
  }

  @Override
  public File file(String name) {
    return mCountParams.file(name);
  }

  /**
   * Get a stream to the output file.
   * @param name file name
   * @return the stream.
   * @throws IOException if an I/O error occurs.
   */
  public OutputStream outStream(String name) throws IOException {
    return mCountParams.outStream(name);
  }

  @Override
  public File directory() {
    return mCountParams.directory();
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    Exam.globalIntegrity(mBuildParams);
    Exam.globalIntegrity(mCountParams);
    return true;
  }

  /**
   * Only record unique words.  This is used by the probe finding code.
   * @return true if only unique words should be retained
   */
  public boolean uniqueWords() {
    return mUniqueWords;
  }

}

