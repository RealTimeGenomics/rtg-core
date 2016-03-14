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
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

import com.rtg.index.hash.ngs.HashFunctionFactory;
import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.NgsHashLoop;
import com.rtg.index.hash.ngs.NgsHashLoopImpl;
import com.rtg.index.hash.ngs.ReadCall;
import com.rtg.index.hash.ngs.ReadEncoder;
import com.rtg.index.hash.ngs.ReadHashFunction;
import com.rtg.index.hash.ngs.TemplateCall;
import com.rtg.index.hash.ngs.TemplateHashFunction;
import com.rtg.index.params.CreateParams;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.HashingRegion;
import com.rtg.ngs.NgsTestUtils.TestPairedEndParams;
import com.rtg.usage.UsageMetric;
import com.rtg.util.MathUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ListenerType;
import com.rtg.util.io.MemoryPrintStream;

/**
 */
public class NgsTaskTest extends NgsPairedEndTest {

  protected static final String SEQ_DNA_Y1 = ">x" + LS + "actg" + LS;
  protected static final String SEQ_DNA_Y2 = ">u" + LS + "actg" + LS + ">v" + LS + "actg" + LS;

  private CreateParams makeIndexParams(NgsParams params) {
    final HashFunctionFactory factory = params.maskParams().maskFactory((int) params.getMaxReadLength());
    final long numSeqs = params.buildFirstParams().numberSequences() + (params.paired() ? params.buildSecondParams().numberSequences() : 0);
    return new CreateParams(numSeqs, factory.hashBits(), factory.windowBits(), 31, params.compressHashes(), true, false, false);
  }

  public void testThreads1Log() throws Exception {
    final ByteArrayOutputStream ba = new ByteArrayOutputStream();
    final PrintStream pr = new PrintStream(ba);
    try {
      Diagnostic.setLogStream(pr);
      try {
        final int numberThreads = 1;
        final NgsTaskFunctionalTest ntt = new NgsTaskFunctionalTest();
        ntt.mDir = mDir;
        try (NgsParams params = ntt.getParams(ba, new NgsMaskParamsGeneral(4, 0, 0, 1), new NgsTestUtils.ParamsParams(SEQ_DNA_Y2, SEQ_DNA_Y1, 10, true, false), ListenerType.NULL, OutputFilter.NONE, 1, numberThreads)) {
          final int threadBits = MathUtils.ceilPowerOf2Bits(numberThreads - 1);
          assertEquals(8, NgsTask.indexThenSearchShortReads(params, new NgsHashLoopImpl(params.buildFirstParams().numberSequences(), params.outputParams().progress(), 0x3FFFFL, ((0x1FFFFL + 1L) << threadBits) - 1L), null, makeIndexParams(params)));
        }
        pr.flush();
        final String logs = ba.toString();
        //deprecated timers to be removed
        TestUtils.containsAll(logs, " Timer Index_initialization ", " Timer Index_sort ", " Timer Index_pointer ", " Timer Index_position ", " Timer Index_bitVector ");
        //index statistics
        TestUtils.containsAll(logs, "Index[", "] statistics", "] search performance");
        assertTrue(logs.contains("Start create job 0"));
        assertTrue(logs.contains("Finish create job 0"));
        assertFalse(logs.contains("Start create job 1"));
        assertFalse(logs.contains("Finish create job 1"));
      } finally {
        Diagnostic.setLogStream();
        pr.close();
      }
    } finally {
      ba.close();
    }
  }

  private class NgsHashLoopDummyImpl implements NgsHashLoop {
    boolean mSecondSeen = false;
    int mReverseCounts = 0;

    @Override
    public long readLoop(ISequenceParams params, ReadHashFunction hashFunction, ReadEncoder encoder, boolean reverse) {
      mSecondSeen = encoder == ReadEncoder.PAIRED_SECOND;
      if (reverse) {
        mReverseCounts++;
      }
      return 43;
    }

    @Override
    public void templateLoop(ISequenceParams params, TemplateHashFunction hashFunction) {
    }

    @Override
    public void templateLoopMultiCore(ISequenceParams params, NgsHashFunction hf, int numberThreads, int threadMultiplier) throws IOException {
      hf.threadClone(HashingRegion.NONE);
    }

  }

  private static class NgsBlahHashFunction implements NgsHashFunction {
    long mReads = -1;
    private TemplateCall mTemplateCall;
    private void setTemplateCall(TemplateCall tc) {
      mTemplateCall = tc;
    }
    @Override
    public void templateSet(long name, int length) { }
    @Override
    public void templateReverse(int endPosition) { }
    @Override
    public void templateForward(int endPosition) { }
    @Override
    public void templateBidirectional(int endPosition) { }
    @Override
    public void logStatistics() { }
    @Override
    public int indelScore(int readId) {
      return 0;
    }
    @Override
    public int fastScore(int readId) {
      return 0;
    }
    @Override
    public void endSequence() { }
    @Override
    public int windowSize() {
      return 0;
    }
    @Override
    public void reset() { }
    @Override
    public int readLength() {
      return 0;
    }
    @Override
    public int numberWindows() {
      return 0;
    }
    @Override
    public void hashStep() { }
    @Override
    public void hashStep(byte code) { }
    @Override
    public void setValues(int id2, boolean reverse) { }
    @Override
    public void setReadSequences(long numberReads) {
      mReads = numberReads;
    }
    @Override
    public void readAll(int readId, boolean reverse) { }
    @Override
    public void threadFinish() { }
    @Override
    public NgsHashFunction threadClone(HashingRegion region) throws IOException {
      mTemplateCall.threadClone(region).threadFinish();
      return null;
    }
  }

  private static final NgsBlahHashFunction HASH_FUNCT = new NgsBlahHashFunction();

  private static final HashFunctionFactory FACTORY = new HashFunctionFactory() {
    @Override
    public NgsHashFunction create(final ReadCall readCall, final TemplateCall templateCall) {
      HASH_FUNCT.setTemplateCall(templateCall);
      return HASH_FUNCT;
    }
    @Override
    public int hashBits() {
      return 1;
    }
    @Override
    public int windowBits() {
      return 1;
    }
    @Override
    public int numberWindows() {
      return 0;
    }
    @Override
    public int windowSize() {
      return 0;
    }
  };

  private class DummyMaskParams extends NgsMaskParamsGeneral {
    DummyMaskParams(int wordSize, int substitutions, int indels, int indelLength) {
      super(wordSize, substitutions, indels, indelLength);
    }
    @Override
    public HashFunctionFactory maskFactory(int readLength) {
      return FACTORY;
    }
  }

  public void testThreads2Log() throws Exception {

    final ByteArrayOutputStream ba = new ByteArrayOutputStream();
    try (PrintStream pr = new PrintStream(ba)) {
      Diagnostic.setLogStream(pr);
      try {
        final int numberThreads = 2;
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        final NgsTaskFunctionalTest ntt = new NgsTaskFunctionalTest();
        ntt.mDir = mDir;
        try (NgsParams params = ntt.getParams(out, new NgsMaskParamsGeneral(4, 0, 0, 1), new NgsTestUtils.ParamsParams(SEQ_DNA_Y2, SEQ_DNA_Y1, 10, true, false), ListenerType.NULL, OutputFilter.NONE, 1, numberThreads)) {
          final NgsHashLoopDummyImpl hashdummy = new NgsHashLoopDummyImpl();
          try {
            assertEquals(43, NgsTask.indexThenSearchShortReads(params, hashdummy, null, makeIndexParams(params)));
          } finally {
            out.close();
          }
          assertFalse(hashdummy.mSecondSeen);
          assertEquals(0, hashdummy.mReverseCounts);
        }

        pr.flush();
        final String logs = ba.toString();
        //deprecated timers to be removed
        TestUtils.containsAll(logs, " Timer Index_initialization ", " Timer Index_sort ", " Timer Index_pointer ", " Timer Index_position ", " Timer Index_bitVector ");
        //index statistics
        TestUtils.containsAll(logs, "Index[0] statistics", "] search performance");
        //thread stats
        TestUtils.containsAll(logs, "Start create job 0", "Finish create job 0", "Start freeze job 0", "Finish freeze job 0");
      } finally {
        Diagnostic.setLogStream();
      }
    } finally {
      ba.close();
    }
  }
  public void testBuildQueryPairedLog() throws Exception {
    final MemoryPrintStream pr = new MemoryPrintStream();
    Diagnostic.setLogStream(pr.printStream());
    try {
      final String cgdata = "@test\nacgtacgtacgtacgtacgtacgtacgtacgtacg\n+\n###################################";
      try (DummyMaskParams dmp = new DummyMaskParams(3, 1, 1, 1)) {
        NgsParams params = getParamsPairedEnd(pr.outputStream(), dmp, new TestPairedEndParams(cgdata, cgdata, ">t\na", "", 1, 1, 0L), true);
        try {
          final NgsHashLoopDummyImpl hashdummy = new NgsHashLoopDummyImpl();
          assertEquals(2 * 43, NgsTask.indexThenSearchShortReads(params, hashdummy, null, makeIndexParams(params)));

          assertTrue(hashdummy.mSecondSeen);
          assertEquals(2, hashdummy.mReverseCounts);
          assertEquals(2, HASH_FUNCT.mReads);

          final String fastadata = ">test\nacgtacgtacgtacgtacgtacgtacgtacgtacg\n>test2\nactttttttt\n";
          final String fastadata2 = ">test\nacgtacgtacgtacgtacgtacgtacgtacgtacg\n";
          params = getParamsPairedEnd(pr.outputStream(), dmp, new TestPairedEndParams(fastadata, fastadata2, ">t\na", "", 1, 1, 0L), false);
          assertEquals(2 * 43, NgsTask.indexThenSearchShortReads(params, hashdummy, null, makeIndexParams(params)));
          assertTrue(hashdummy.mSecondSeen);
          assertEquals(2, hashdummy.mReverseCounts);
          assertEquals(3, HASH_FUNCT.mReads);
        } finally {
          params.close();
        }
      }
    } finally {
      Diagnostic.setLogStream();
    }
  }

  public void testConstructor() throws Exception {
    final OutputStream os = new ByteArrayOutputStream();
    TestPairedEndParams params = new NgsTestUtils.TestPairedEndParams(">test\nACGT\n>test2\nGCTA\n", ">test\nAGCT\n", ">t\na", "", 1, 1, 0L);
    try {
      new NgsTask(getParamsPairedEnd(os, new NgsMaskParamsGeneral(4, 0, 0, 1), params, false), null, new UsageMetric());
      fail();
    } catch (final RuntimeException re) {
      assertEquals("Number of reads in first and second read sets must be equal.", re.getMessage());
    }
    params = new NgsTestUtils.TestPairedEndParams(">test\nACGTT\n", ">test\nAGCT\n", ">t\na", "", 1, 1, 0L);
    try {
      new NgsTask(getParamsPairedEnd(os, new NgsMaskParamsGeneral(4, 0, 0, 1), params, false), null, new UsageMetric());
      fail();
    } catch (final RuntimeException re) {
      assertEquals("Length of reads in first and second read sets must be equal.", re.getMessage());
    }
    final NgsParams params2 = getParamsPairedEnd(os, new NgsMaskParamsGeneral(4, 0, 0, 1), new NgsTestUtils.TestPairedEndParams(">test\nACGT\n", ">test\nAGCT\n", ">t\na", "", 1, 1, 0L), false);
    final NgsTask ngst = new NgsTask(params2, os, new UsageMetric());
    assertEquals(params2, ngst.parameters());
    assertNotNull(ngst.hashCode());
    assertTrue(ngst.getStatistics() instanceof PairedEndMapStatistics);
  }

}
