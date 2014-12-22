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
package com.rtg.position.output;



import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.index.params.CreateParams;
import com.rtg.launcher.BuildParams;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.OutputDirParams;
import com.rtg.mode.ProgramMode;
import com.rtg.ngs.NgsParams;
import com.rtg.util.ObjectParams;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.Integrity;
import com.rtg.util.io.FileUtils;

/**
 * Holds all the parameters needed for doing a build and a search.
 */
public class PositionParams extends ObjectParams implements OutputDirParams, Integrity {


  /**
   * Name of the file where unmapped sequences are written.
   */
  public static final String UNMAPPED_SUFFIX = "unmapped";

  static final String OUT_SUFFIX = "out";


  /**
   * Creates a PositionParams builder.
   * @return the builder.
   */
  public static PositionParamsBuilder builder() {
    return new PositionParamsBuilder();
  }

  /**
   * A builder class for <code>PositionParams</code>.
   */
  public static class PositionParamsBuilder {

    protected ProgramMode mProgramMode = ProgramMode.SLIMN;

    protected Integer mHashCountThreshold = 1000;

    protected boolean mProgress = false;

    protected int mNumberThreads = 1;

    protected BuildParams mBuildParams = null;

    protected BuildParams mBuildSecondParams = null;

    protected BuildParams mSearchParams = null;

    protected PositionOutputParams mOutputParams = null;

    protected NgsParams mNgsParams = null;

    /**
     * Sets the program mode.
     *
     * @param mode the program mode. The default is <code>ProgramMode.SLIMN</code>.
     * @return this builder, so calls can be chained.
     */
    public PositionParamsBuilder mode(final ProgramMode mode) {
      mProgramMode = mode;
      return this;
    }

    /**
     * Sets the maximum hash count threshold.
     *
     * @param threshold the threshold. The default is 1000.
     * @return this builder, so calls can be chained.
     */
    public PositionParamsBuilder hashCountThreshold(final int threshold) {
      mHashCountThreshold = threshold;
      return this;
    }

    /**
     * Sets the number of threads to use.
     *
     * @param numThreads the number of threads to employ. The default is 1.
     * @return this builder, so calls can be chained.
     */
    public PositionParamsBuilder numberThreads(final int numThreads) {
      mNumberThreads = numThreads;
      return this;
    }

    /**
     * Sets whether progress should be output.
     *
     * @param progress true if progress should be output.
     * @return this builder, so calls can be chained.
     */
    public PositionParamsBuilder progress(final boolean progress) {
      mProgress = progress;
      return this;
    }

    /**
     * Sets the build parameters.
     *
     * @param params the parameters used during build.
     * @return this builder, so calls can be chained.
     */
    public PositionParamsBuilder buildParams(final BuildParams params) {
      mBuildParams = params;
      return this;
    }

    /**
     * Sets the build parameters.
     *
     * @param params the parameters used during build.
     * @return this builder, so calls can be chained.
     */
    public PositionParamsBuilder buildSecondParams(final BuildParams params) {
      mBuildSecondParams = params;
      return this;
    }

    /**
     * Sets the search parameters.
     *
     * @param params the parameters used during search.
     * @return this builder, so calls can be chained.
     */
    public PositionParamsBuilder searchParams(final BuildParams params) {
      mSearchParams = params;
      return this;
    }

    /**
     * Sets the output parameters.
     *
     * @param params the parameters used for output.
     * @return this builder, so calls can be chained.
     */
    public PositionParamsBuilder outputParams(final PositionOutputParams params) {
      mOutputParams = params;
      return this;
    }

    /**
     * Sets the output parameters.
     *
     * @param params the <code>NgsParams</code> used to create these position params
     * @return this builder, so calls can be chained.
     */
    public PositionParamsBuilder ngsParams(final NgsParams params) {
      mNgsParams = params;
      return this;
    }

    /**
     * Creates a PositionParams using the current builder
     * configuration.
     * @return the new PositionParams
     */
    public PositionParams create() {
      return new PositionParams(this);
    }
  }



  private final ProgramMode mProgramMode;

  private final BuildParams mBuildParams;
  private final BuildParams mBuildSecondParams;

  private final BuildParams mSearchParams;

  private final Integer mHashCountThreshold;

  private final boolean mProgress;

  private final int mNumberThreads;

  private final PositionOutputParams mOutputParams;

  private final NgsParams mNgsParams;


  /**
   * @param builder the builder object.
   */
  public PositionParams(final PositionParamsBuilder builder) {
    mProgramMode = builder.mProgramMode;
    mNumberThreads = builder.mNumberThreads;
    mHashCountThreshold = builder.mHashCountThreshold;
    mBuildParams = builder.mBuildParams;
    mBuildSecondParams = builder.mBuildSecondParams;
    mSearchParams = builder.mSearchParams;
    mOutputParams = builder.mOutputParams;
    mProgress = builder.mProgress;
    mNgsParams = builder.mNgsParams;
    mObjects = new Object[] {mOutputParams, mProgramMode, mBuildParams, mBuildSecondParams, mSearchParams, mHashCountThreshold, mProgress, mNumberThreads, mNgsParams};
  }

  /**
   * @param programMode program mode.
   * @param threshold maximum repeat frequency threshold.
   * @param buildParams parameters for doing create and build.
   * @param searchParams parameters for doing search.
   * @param outParams used for logging and other auxiliary output.
   * @param progress flag. (true iff progress is to be output).
   * @param numberThreads number of threads available to be used in parallel
   */
  private PositionParams(
      final ProgramMode programMode,
      final Integer threshold,
      final BuildParams buildParams,
      final BuildParams buildSecondParams,
      final BuildParams searchParams,
      final PositionOutputParams outParams,
      final boolean progress,
      final int numberThreads,
      final NgsParams ngsParams
  ) {
    mProgramMode = programMode;
    mBuildParams = buildParams;
    mBuildSecondParams = buildSecondParams;
    mSearchParams = searchParams;
    mOutputParams = outParams;
    mHashCountThreshold = threshold;
    mProgress = progress;
    mNumberThreads = numberThreads;
    mNgsParams = ngsParams;
    mObjects = new Object[] {mOutputParams, mProgramMode, mBuildParams, mBuildSecondParams, mSearchParams, mHashCountThreshold, mProgress, mNumberThreads, mNgsParams};
  }

  /**
   * Create a new <code>PositionParams</code> object that is identical to this except
   * that the search range is set ot a subset of the current range.
   * @param region of subsequence to construct sub params for
   * @return the new <code>PositionParams</code> object.
   * @throws IOException If an I/O error occurs
   */
  public PositionParams subSearch(final HashingRegion region) throws IOException {
    final BuildParams search = search();
    //System.err.println(search.getClass().getName() + "search.class " + search);
    final BuildParams build = build();
    final ISequenceParams sequences = build.sequences();
    //TODO improve efficiency of duplicating build
    return new PositionParams(
        mProgramMode, hashCountThreshold(), build.subSequence(sequences.region()), buildSecond() != null ? buildSecond().subSequence(sequences.region()) : null, search.subSequence(region), output(), progress(), numberThreads(), ngsParams()
    );
  }

  /**
   * Create a new <code>PositionParams</code> object that is identical to this except
   * that the search range is set ot a subset of the current range.
   * @param region of subsequence to construct sub params for
   * @return the new <code>PositionParams</code> object.
   * @throws IOException If an I/O error occurs
   */
  public PositionParams subBuild(final HashingRegion region) throws IOException {
    final BuildParams search = search();
    //System.err.println(search.getClass().getName() + "search.class " + search);
    final BuildParams build = build();
    //TODO improve efficiency of duplicating build
    return new PositionParams(
        mProgramMode, hashCountThreshold(), build.subSequence(region), buildSecond() != null ? buildSecond().subSequence(region) : null, search.subSequence(search.sequences().region()), output(), progress(), numberThreads(), ngsParams()
    );
  }

  /**
   * Get the maximum repeat frequency threshold.
   * @return the maximum repeat frequency threshold.
   */
  public Integer hashCountThreshold() {
    return mHashCountThreshold;
  }

  /**
   * Params used when in <code>ngsMode</code>
   * @return the params
   */
  public NgsParams ngsParams() {
    return mNgsParams;
  }

  /**
   * Get the progress flag.
   * @return the progress flag. (true iff progress is to be output).
   */
  public boolean progress() {
    return mProgress;
  }

  /**
   * @return the program mode.
   */
  public ProgramMode mode() {
    return mProgramMode;
  }

  /**
   * @return the build and create parameters.
   */
  public BuildParams build() {
    return mBuildParams;
  }

  /**
   * @return the build and create parameters.
   */
  public BuildParams buildSecond() {
    return mBuildSecondParams;
  }

  /**
   * @return the search parameters.
   */
  public BuildParams search() {
    return mSearchParams;
  }

  /**
   * @return the output parameters.
   */
  public PositionOutputParams output() {
    return mOutputParams;
  }

  @Override
  public File directory() {
    return output().directory();
  }

  /**
   * @return the output directory.
   */
  public File outputDir() {
    return mOutputParams.directory();
  }

  /**
   * Return create parameters.
   * @return the parameters
   */
  public CreateParams indexParams() {
    if (mBuildSecondParams != null) {
      return new CreateParams(mBuildParams.size() + mBuildSecondParams.size(), mBuildParams.hashBits(), mBuildParams.hashBits(), mNgsParams.compressHashes(), false, false);
    } else {
      return mBuildParams;
    }
  }

  /**
   *  Get the length required for the shared subject/query buffer.
   * @return the length required for the shared subject/query buffer.
   */
  public long bufferLength() {
    if (mBuildSecondParams != null) {
      return Math.max(mBuildSecondParams.sequences().maxLength(), Math.max(mBuildParams.sequences().maxLength(), mSearchParams.sequences().maxLength()));

    }
    return Math.max(mBuildParams.sequences().maxLength(), mSearchParams.sequences().maxLength());
  }

  /**
   * Get the number of threads available for parallel execution.
   * @return the number of threads available for parallel execution.
   */
  public int numberThreads() {
    return mNumberThreads;
  }

  @Override
  public void close() throws IOException {
    mBuildParams.close();
    if (mBuildSecondParams != null) {
      mBuildSecondParams.close();
    }
    mSearchParams.close();
  }

  /**
   * Check if queries are closed.
   * @return true iff the queries are currently closed.
   */
  @Override
  public boolean closed() {
    return mBuildParams.closed() && mSearchParams.closed();
  }

  /**
   * Get the name of the logfile.
   * @return the name of the logfile.
   */
  public String logFile() {
    return new File(outputDir(), FileUtils.LOG_SUFFIX).getPath();
  }

  /**
   * Compute smallest l such that n &lt;= 10^l.
   * @param n value being checked.
   * @return the log.
   */
  public static int log10(final int n) {
    if (n <= 0) {
      throw new IllegalArgumentException();
    }
    int i = 10;
    int l = 1;
    while (n > i) {
      i = i * 10;
      if (i < 0) {
        return l - 1;
      }
      l++;
    }
    return l;
  }

  /**
   * Get a number formatted so that all positive numbers less than n will be
   * have the same length (using leading 0's if necessary to pad them out).
   * @param i number to be formatted.
   * @param n maximum value.
   * @return formatted version of i.
   */
  public static String zeroFormat(final int i, final int n) {
    final int l = log10(n);
    if (n == 1) {
      return "";
    }
    final String formatString;
    formatString = "%0" + l + "d";
    //
    return String.format(formatString, i);
  }

  /**
   * Get a stream to the output file.
   * @return the stream.
   * @throws IOException if an I/O error occurs.
   */
  public OutputStream outStream() throws IOException {
    return stream(OUT_SUFFIX);
  }

  /**
   * Get a stream to the output file.
   * @param thread the current thread (starts at 0).
   * @param numberThreads total number of threads.
   * @return the stream.
   * @throws IOException if an I/O error occurs.
   */
  public OutputStream outStream(final int thread, final int numberThreads) throws IOException {
    return stream(OUT_SUFFIX + zeroFormat(thread, numberThreads));
  }

  /**
   * Get a stream to the output file with a suffix attached.
   * @param suffix to apply to the name.
   * @return the stream.
   * @throws IOException if an I/O error occurs.
   */
  OutputStream stream(final String suffix) throws IOException {
    final String name = mOutputParams.zip() ? suffix + FileUtils.GZ_SUFFIX : suffix;
    //System.err.println("name=" + name);
    if (!outputDir().exists() && !outputDir().mkdirs()) {
      throw new FileNotFoundException();
    }
    return FileUtils.createOutputStream(new File(outputDir(), name), mOutputParams.zip(), false);
  }

  @Override
  public String toString() {
    final String linePrefix = LS + "..";
    return "PositionParams mode=" + mProgramMode + " threshold=" + mHashCountThreshold + " progress=" + mProgress + " number threads=" + mNumberThreads + linePrefix
        + " output={" + mOutputParams + "}" + linePrefix
        + " build={" + build() + "}" + linePrefix
        + " search={" + search() + "}" + LS;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mProgramMode != null);
    Exam.assertNotNull(mSearchParams);
    Exam.assertNotNull(mBuildParams);
    Exam.assertNotNull(mOutputParams);
    if (mProgramMode != null) {
      Exam.assertEquals(mProgramMode.subjectMode(), mBuildParams.sequences().mode());
      Exam.assertEquals(mProgramMode.queryMode(), mSearchParams.sequences().mode());
    }
    Exam.assertTrue(mHashCountThreshold != null && mHashCountThreshold >= 1);
    Exam.assertTrue(mNumberThreads >= 1);
    return true;
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    mBuildParams.globalIntegrity();
    mSearchParams.globalIntegrity();
    mOutputParams.globalIntegrity();
    return true;
  }
}

