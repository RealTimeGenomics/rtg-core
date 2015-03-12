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
package com.rtg.ngs;

import static com.rtg.util.StringUtils.LS;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.Collections;

import com.rtg.launcher.SequenceParams;
import com.rtg.mode.SequenceMode;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.usage.UsageMetric;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.ListenerType;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.Assert;

/**
 */
public final class NgsTestUtils {

  /** Prevent instantiation. */
  private NgsTestUtils() { }

  static final String HEADER = "#template-name\tframe\tread-id\ttemplate-start\tscore\tscore-indel" + LS;

  /**
   * class to reduce parameter list
   */
  public static class ParamsParams {
    String mSubjects;
    String mQueries;
    int mErrorLimit;
    boolean mProgress;
    boolean mZip;
    //int mStepSize;

    /**
     * Constructor to set the parameters
     * @param subjects subjects - probably the reads
     * @param queries queries - probably the reference
     * @param errorLimit error limit
     * @param progress flag if progressed
     * @param zip flag if output should be zipped
     */
    public ParamsParams(final String subjects, final String queries, final int errorLimit, final boolean progress, final boolean zip) {
      mSubjects = subjects;
      mQueries = queries;
      mErrorLimit = errorLimit;
      mProgress = progress;
      mZip = zip;
    }

  }

  /**
   * class for the test data, reference and query sequence and expected
   */
  public static class TestPairedEndParams extends TestParams {
    final String mSubjectsSecond;

    final Integer mMinInsertSize;

    final Integer mMaxInsertSize;

    /**
     * @param subjects left read
     * @param subjectsSecond right reads
     * @param queries queries template
     * @param expected expected results
     * @param minInsertSize minimum insert size
     * @param maxInsertSize maxiumum insert size
     * @param usage expected usage count.
     */
    public TestPairedEndParams(final String subjects, final String subjectsSecond, final String queries, final String expected, final Integer minInsertSize, final Integer maxInsertSize, final Long usage) {
      super(subjects, queries, expected, usage);
      //System.err.println("Expected\n" + expected);

      mSubjectsSecond = subjectsSecond;
      mMinInsertSize = minInsertSize;
      mMaxInsertSize = maxInsertSize;
    }
  }

  /**
   * class for the test data, reference and query sequence and expected
   */
  public static class TestParams {
    static final boolean D_USELONGREADMAPPING = false;
    static final int D_STEPSIZE = -1;
    static final ListenerType D_LISTENER = ListenerType.NULL;
    static final boolean D_PROGRESS = false;
    static final String[] D_ARGSTRING = null;
    static final long D_USAGE = 0;

    final boolean mUseLongReadMapping;
    final int mStepSize;
    final String mSubjects;
    final String mQueries;
    final String mExpected;
    final ListenerType mListener;
    final boolean mProgress;
    final String[] mArgsString;
    final long mUsage;

    /**
     * set simple three test parameters
     * @param subjects template sequences
     * @param queries query sequences
     * @param expected flag expected
     * @param usage expected usage count.
     */
    public TestParams(final String subjects, final String queries, final String expected, long usage) {
      mUseLongReadMapping = D_USELONGREADMAPPING;
      mSubjects = subjects;
      mQueries = queries;
      mExpected = expected;
      mListener = D_LISTENER;
      mProgress = D_PROGRESS;
      mStepSize = D_STEPSIZE;
      mArgsString = D_ARGSTRING;
      mUsage = usage;
    }

    /**
     * set simple three test parameters
     * @param subjects template sequences
     * @param queries query sequences
     * @param expected flag expected
     * @param argsString argument string
     * @param usage expected usage count.
     */
    public TestParams(final String subjects, final String queries, final String expected, final String[] argsString, final long usage) {
      mUseLongReadMapping = D_USELONGREADMAPPING;
      mStepSize = D_STEPSIZE;
      mSubjects = subjects;
      mQueries = queries;
      mExpected = expected;
      mListener = D_LISTENER;
      mProgress = D_PROGRESS;
      mArgsString = new String[argsString.length];
      System.arraycopy(argsString, 0, mArgsString, 0, mArgsString.length);
      mUsage = usage;
    }

    /**
     * set test parameter, and enforce-longreadmapping is possible to be set
     * @param useLongReadMapping flag to enforce long or short readmapping
     * @param stepSize step size for long read mapping
     * @param subjects the template
     * @param queries the reads
     * @param expected the expected result
     * @param usage expected usage count.
     */
    public TestParams(final boolean useLongReadMapping, final int stepSize, final String subjects, final String queries, final String expected, long usage) {
      mUseLongReadMapping = useLongReadMapping;
      mStepSize = stepSize;
      mSubjects = subjects;
      mQueries = queries;
      mExpected = expected;
      mListener = D_LISTENER;
      mProgress = D_PROGRESS;
      mArgsString = D_ARGSTRING;
      mUsage = usage;
    }
    /**
     * set test parameter, and enforce-longreadmapping is possible to be set
     * @param useLongReadMapping flag to enforce long or short readmapping
     * @param stepSize step size for long read mapping
     * @param subjects the template
     * @param queries the reads
     * @param expected the expected result
     * @param argsString argument string
     */
    public TestParams(final boolean useLongReadMapping, final int stepSize, final String subjects, final String queries, final String expected, final String[] argsString) {
      mUseLongReadMapping = useLongReadMapping;
      mStepSize = stepSize;
      mSubjects = subjects;
      mQueries = queries;
      mExpected = expected;
      mListener = D_LISTENER;
      mProgress = D_PROGRESS;
      mArgsString = new String[argsString.length];
      System.arraycopy(argsString, 0, mArgsString, 0, mArgsString.length);
      mUsage = D_USAGE;
    }
  }

  /**
   * short version of NgsFilterParams,
   * exclude and usid are fixed
   */
  public static class NgsFilterPartlyParams {
    static final boolean D_ZIP = false;
    static final boolean D_EXCLUDEREPEATS = false;
    static final boolean D_USESEQUENCEIDS = false;

    private final OutputFilter mOutputFilter;
    private final boolean mZip;
    private final int mTopN;
    private final boolean mExcludeRepeats;
    private final boolean mUseSequenceIds;
    private final int mErrorLimit;
    //private int mDistance;

    /**
     * sets part of the filter parameters, rest default
     * @param outputFilter output filter
     * @param topn num for topn filter
     */
    public NgsFilterPartlyParams(final OutputFilter outputFilter, final int topn) {
      mOutputFilter = outputFilter; //param
      mZip = D_ZIP;
      mTopN = topn; //param
      mExcludeRepeats = D_EXCLUDEREPEATS;
      mUseSequenceIds = D_USESEQUENCEIDS;
      mErrorLimit = 10; //D_ERRORLIMIT;
      //mDistance = 1;
    }

    NgsFilterParams getFilterParams() {
      return NgsFilterParams.builder().outputFilter(mOutputFilter).zip(mZip).topN(mTopN).exclude(mExcludeRepeats).useids(mUseSequenceIds).errorLimit(mErrorLimit).create();
    }
  }

  /** the better check */
  static void checkShort(final NgsMaskParams mask, final TestParams test, final NgsFilterPartlyParams filterParams) throws Exception {
    final ByteArrayOutputStream out = new ByteArrayOutputStream();
    try {
      final File tmpDir = FileHelper.createTempDirectory();
      try {
        final NgsParams params = getParamsShort(out, mask, test, filterParams, tmpDir);
        execNgs(params);
      } finally {
        Assert.assertTrue(FileHelper.deleteAll(tmpDir));
      }
    } finally {
      out.close();
    }
    if (!test.mExpected.equals("")) {
      Assert.assertEquals(NgsTestUtils.HEADER + test.mExpected, out.toString());
    } else {
      Assert.assertEquals(test.mExpected, out.toString());
    }
  }

  /** the better getParams */
  static NgsParams getParamsShort(final OutputStream out, final NgsMaskParams mask, final TestParams test, final NgsFilterPartlyParams filter, File tmpDir) throws IOException {
    final File subjectsDir = FileHelper.createTempDirectory(tmpDir);
    final File queriesDir = FileHelper.createTempDirectory(tmpDir);
    final File hitsDir = FileHelper.createTempDirectory(tmpDir);
    //System.err.println("hitsDir=" + hitsDir);
    ReaderTestUtils.getReaderDNA(test.mSubjects, subjectsDir, null).close();
    ReaderTestUtils.getReaderDNA(test.mQueries, queriesDir, null).close();
    final SequenceParams subjectParams = SequenceParams.builder().directory(subjectsDir).mode(SequenceMode.UNIDIRECTIONAL).create();
    final SequenceParams queryParams = SequenceParams.builder().directory(queriesDir).create();
    final NgsFilterParams filterParams = filter.getFilterParams();
    final NgsOutputParams outputParams = new OverriddenNgsOutputParams(OverriddenNgsOutputParams.builder().outStream(out).progress(test.mProgress).outputDir(hitsDir).filterParams(filterParams));

    return NgsParams.builder().useLongReadMapping(test.mUseLongReadMapping).stepSize(test.mStepSize).buildFirstParams(subjectParams).searchParams(queryParams).outputParams(outputParams).maskParams(mask).listeners(Collections.singleton(test.mListener)).create();
  }

  /**
   * Override the output stream.
   */
  public static class OverriddenNgsOutputParams extends NgsOutputParams {

    /** Override the builder */
    public static class OverriddenNgsOutputParamsBuilder extends NgsOutputParamsBuilder {
      protected OutputStream mOut = new ByteArrayOutputStream();
      protected OutputStream mRepeats = new ByteArrayOutputStream();
      protected OutputStream mUnmapped = new ByteArrayOutputStream();

      /**
       * @param stream override unmapped output stream
       * @return this builder, so calls can be chained.
       */
      public OverriddenNgsOutputParamsBuilder outStream(final OutputStream stream) {
        mOut = stream;
        return this;
      }
    }

    /**
     * Creates a NgsOutputParams builder.
     * @return the builder.
     */
    public static OverriddenNgsOutputParamsBuilder builder() {
      return new OverriddenNgsOutputParamsBuilder();
    }
    private OutputStream mOut;
    private OutputStream mRepeats;
    private OutputStream mUnmapped;


    @Override
    public boolean equals(Object obj) {
      if (obj instanceof OverriddenNgsOutputParams) {
        final OverriddenNgsOutputParams op = (OverriddenNgsOutputParams) obj;
        return super.equals(obj)
            && (mOut == null || mOut.equals(op.mOut))
            && (mRepeats == null || mRepeats.equals(mRepeats))
            && (mUnmapped == null || mUnmapped.equals(mUnmapped));

      }
      return false;
    }

    @Override
    public int hashCode() {
      return super.hashCode()
          + (mOut == null ? 0 : mOut.hashCode())
          + (mRepeats == null ? 0 : mRepeats.hashCode())
          + (mUnmapped == null ? 0 : mUnmapped.hashCode());
    }


    /**
     * Constructor
     * @param builder other parameters
     */
    public OverriddenNgsOutputParams(final NgsOutputParamsBuilder builder) {
      super(builder);
      if (builder instanceof OverriddenNgsOutputParamsBuilder) {
        mOut = ((OverriddenNgsOutputParamsBuilder) builder).mOut;
        mRepeats = ((OverriddenNgsOutputParamsBuilder) builder).mRepeats;
        mUnmapped = ((OverriddenNgsOutputParamsBuilder) builder).mUnmapped;
      }
    }

    @Override
    public OutputStream outStream() throws IOException {
      return mOut != null ? mOut : super.outStream();
    }
    @Override
    public OutputStream repeatsStream() throws IOException {
      return mRepeats != null ? mRepeats : super.repeatsStream();
    }
    @Override
    public OutputStream unmappedStream() throws IOException {
      return mUnmapped != null ? mUnmapped : super.unmappedStream();
    }
  }

  /**
   * Test long reads
   * @param mask mask used
   * @param test test parameters
   * @param filterParams filter parameters
   * @param containsOnly flag if contains only
   * @throws Exception exception
   */
  public static void checkLong(final NgsMaskParams mask, final TestParams test, final NgsFilterPartlyParams filterParams, final boolean containsOnly) throws Exception {
    final File hitsDir = FileHelper.createTempDirectory();
    try {
      final NgsParams params = getParamsLong(hitsDir, mask, test, filterParams);
      final NgsTask task = execNgs(params);
      final String output = FileUtils.fileToString(new File(hitsDir, "out"));
      if (test.mExpected.equals("")) {
        Assert.assertEquals(test.mExpected, output);
      } else {
        if (!containsOnly) {
          Assert.assertEquals(HEADER + test.mExpected, output);
        } else {
          final String exp = HEADER + test.mExpected;
          final String[] outputArray = output.split(LS);
          containsAllWithoutReplacement(exp, outputArray);
        }
      }
      Assert.assertEquals(test.mUsage, task.usage());
    } finally {
      Assert.assertTrue(FileHelper.deleteAll(hitsDir));
    }
  }

  static void checkLongThreads(final NgsMaskParams mask, final TestParams test, final NgsFilterPartlyParams filterParams, final boolean containsOnly) throws Exception {
    final File hitsDir = FileHelper.createTempDirectory();
    try {
      final NgsParams params = getParamsLongThreads(hitsDir, mask, test, filterParams, 4);
      execNgs(params);
      final String output = FileUtils.fileToString(new File(hitsDir, "out"));
      if (test.mExpected.equals("")) {
        Assert.assertEquals(test.mExpected, output);
      } else {
        if (!containsOnly) {
          if (!test.mExpected.equals("")) {
            final String exps = TestUtils.sortLines(NgsTestUtils.HEADER + test.mExpected);
            final String acts = TestUtils.sortLines(output);
            Assert.assertEquals("expected:" + LS + exps + LS + "actual:" + LS + acts, exps, acts);
          } else {
            Assert.assertEquals(test.mExpected, output);
          }
        } else {
          final String exp = HEADER + test.mExpected;
          final String[] outputArray = output.split(LS);
          containsAllWithoutReplacement(exp, outputArray);
        }
      }
    } finally {
      Assert.assertTrue(FileHelper.deleteAll(hitsDir));
    }
  }

  static NgsParams getParamsLong(final File hitsDir, final NgsMaskParams mask, final TestParams test, final NgsFilterPartlyParams filter) throws IOException {
    return getParamsLongThreads(hitsDir, mask, test, filter, 1);
  }

  static NgsParams getParamsLongThreads(final File hitsDir, final NgsMaskParams mask, final TestParams test, final NgsFilterPartlyParams filter, final int threads) throws IOException {
    final File subjectsDir = FileHelper.createTempDirectory(hitsDir);
    final File queriesDir = FileHelper.createTempDirectory(hitsDir);
    //System.err.println("hitsDir=" + hitsDir);
    ReaderTestUtils.getReaderDNA(test.mSubjects, subjectsDir, null).close();
    ReaderTestUtils.getReaderDNA(test.mQueries, queriesDir, null).close();
    final SequenceParams subjectParams = SequenceParams.builder().directory(subjectsDir).mode(SequenceMode.UNIDIRECTIONAL).useMemReader(true).create();
    final SequenceParams queryParams = SequenceParams.builder().directory(queriesDir).create();
    final NgsFilterParams filterParams = filter.getFilterParams();
    //  public NgsOutputParams(final NgsMaskParams windowSelector, final boolean progress, final File outputDir, final NgsFilterParams filterParams) {

    final NgsOutputParams outputParams = NgsOutputParams.builder().progress(test.mProgress).outputDir(hitsDir).filterParams(filterParams).create();

    return NgsParams.builder().useLongReadMapping(test.mUseLongReadMapping).stepSize(test.mStepSize).numberThreads(threads).buildFirstParams(subjectParams).searchParams(queryParams).outputParams(outputParams).maskParams(mask).listeners(Collections.singleton(test.mListener)).create();
  }


  static NgsTask getNgs(final NgsParams ngs)  {
    return new NgsTask(ngs, null, new UsageMetric());
  }

  static NgsTask execNgs(final NgsParams params) throws Exception {
    final NgsTask ngs = getNgs(params);
    ngs.exec();
    params.close();
    return ngs;
  }

  /**
   * Check that a number of strings are contained in str is used for comparison of result
   * output;
   * string contained only once and no strings are left
   * Substitutes lines with "xxx" string
   * @param str string to be checked.
   * @param subs all the items that should be in the string.
   */
  public static void containsAllWithoutReplacement(final String str, final String... subs) {
    boolean ok = true;
    final StringBuilder sb = new StringBuilder();
    String newString = str;
    //int len = subs.length;
    //int count = 0;
    for (final String sub : subs) {
      if (!str.contains(sub)) {
        sb.append("'").append(sub).append("' was not contained in:").append(str).append(StringUtils.LS);
        ok = false;
      } else {
        newString = newString.replaceFirst(sub, "xxx");
      }
      //count++;
      Assert.assertTrue(sb.toString(), ok);
    }

    final String[] newStrings = newString.split(StringUtils.LS);
    for (final String xxx : newStrings) {
      if (!xxx.equals("xxx")) {
        sb.append("'").append(xxx).append("' was not delivered in results:").append(Arrays.toString(subs)).append(StringUtils.LS);
        ok = false;
      }
      Assert.assertTrue(sb.toString(), ok);
    }
  }
}

