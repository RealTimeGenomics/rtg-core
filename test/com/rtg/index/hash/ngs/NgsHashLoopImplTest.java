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
package com.rtg.index.hash.ngs;

import static com.rtg.util.StringUtils.LS;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.io.StringWriter;
import java.util.ArrayList;

import com.rtg.launcher.DefaultReaderParamsTest;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.MockSequenceParams;
import com.rtg.launcher.ReaderParams;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.DNAFastaSymbolTable;
import com.rtg.mode.SequenceMode;
import com.rtg.reader.FastaSequenceDataSource;
import com.rtg.reader.PrereadType;
import com.rtg.reader.ReaderLongMock;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesWriter;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.CliDiagnosticListener;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.SlimException;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public final class NgsHashLoopImplTest extends TestCase {

  private File mDir;

  @Override
  public void setUp() throws Exception {
    Diagnostic.setLogStream();
    mDir = FileHelper.createTempDirectory();
  }

  @Override
  public void tearDown() {
    assertTrue(FileHelper.deleteAll(mDir));
    mDir = null;
  }

  private ReaderParams getReads(final String inputDnaSequence) throws IOException {
    final File dir = FileHelper.createTempDirectory(mDir);
    final ArrayList<InputStream> inputStreams = new ArrayList<>();
    inputStreams.add(new ByteArrayInputStream(inputDnaSequence.getBytes()));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(inputStreams,
            new DNAFastaSymbolTable());
    final SequencesWriter sequenceWriter = new SequencesWriter(ds, dir, 20, PrereadType.UNKNOWN, false);
    sequenceWriter.processSequences();
    return DefaultReaderParamsTest.createDefaultReaderParams(dir, LongRange.NONE, SequenceMode.UNIDIRECTIONAL);
  }

  private ReaderParams getTemplate(final String inputDnaSequence) throws IOException {
    final File dir = FileHelper.createTempDirectory(mDir);
    final ArrayList<InputStream> inputStreams = new ArrayList<>();
    inputStreams.add(new ByteArrayInputStream(inputDnaSequence.getBytes()));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(inputStreams,
            new DNAFastaSymbolTable());
    final SequencesWriter sequenceWriter = new SequencesWriter(ds, dir, 20, PrereadType.UNKNOWN, false);
    sequenceWriter.processSequences();
    return DefaultReaderParamsTest.createDefaultReaderParams(dir, LongRange.NONE, SequenceMode.BIDIRECTIONAL);
  }

  private void check(final String readstr, final String template, final String expected) throws IOException {
    final ReaderParams tem = getTemplate(template);
    final ReaderParams reads = getReads(readstr);
    final NgsHashLoopImpl hl = new NgsHashLoopImpl(reads.reader().numberSequences(), false);
    hl.integrity();
    assertEquals("HashLoop", hl.toString());

    final StringWriter sb = new StringWriter();
    final ByteArrayOutputStream byteOut = new ByteArrayOutputStream();
    final PrintStream progressOut = new PrintStream(byteOut);

    final CliDiagnosticListener clir = new CliDiagnosticListener(progressOut);
    Diagnostic.addListener(clir);

    final TemplateHashFunction thf = new HashFunctionMock(sb);
    try {
      hl.templateLoop(new MockSequenceParams(tem, 0, tem.reader().numberSequences()), thf);
      fail();
    } catch (final RuntimeException e) {
      assertEquals("Read sequences not defined", e.getMessage());
    }

    final ReadHashFunction rhf = new HashFunctionMock(sb);

    sb.append("bulding").append(LS);
    try (ISequenceParams sp1 = new MockSequenceParams(reads, 0, reads.reader().numberSequences())) {
      rhf.setReadSequences(reads.reader().numberSequences());
      hl.readLoop(sp1, rhf, ReadEncoder.SINGLE_END, false);
    }
    sb.append("searching").append(LS);
    try (ISequenceParams sp2 = new MockSequenceParams(tem, 0, tem.reader().numberSequences())) {
      hl.templateLoop(sp2, thf);
    }
    Diagnostic.removeListener(clir);
    assertEquals(expected, sb.toString());

  }

  private static void append(final Appendable ap, final String str) {
    try {
      ap.append(str);
    } catch (final IOException e) {
      fail();
    }
  }

  static class HashFunctionMock implements NgsHashFunction, Cloneable {

    final Appendable mOut;
    private long[] mReadSequences;
    private int mSoFar = 0;
    int mClones = 0;

    HashFunctionMock(final Appendable out) {
      mOut = out;
    }

    @Override
    public void templateForward(final int endPosition) {
      if (mSoFar >= readLength()) {
        append(mOut, "pf=" + endPosition + LS);
      }
    }

    @Override
    public void templateBidirectional(final int endPosition) {
      templateForward(endPosition);
      templateReverse(endPosition);
    }

    @Override
    public void templateReverse(final int endPosition) {
      if (mSoFar >= 64) {
        append(mOut, "pr=" + (endPosition - (64 - readLength())) + LS);
      }
    }

    @Override
    public void endSequence() {
      // do nothing
    }

    @Override
    public void templateSet(final long name, final int length) {
      append(mOut, name + " " + length + LS);
    }

    @Override
    public void hashStep(final byte code) {
      append(mOut, "code " + code + LS);
      hashStep();
    }

    @Override
    public void hashStep() {
      mSoFar++;
    }

    @Override
    public int numberWindows() {
      return 0;
    }

    @Override
    public void reset() {
      append(mOut, "reset" + LS);
      mSoFar = 0;
    }

    @Override
    public int fastScore(final int readId) {
      return 0;
    }

    @Override
    public int indelScore(final int readId) {
      return 0;
    }

    @Override
    public void setReadSequences(final long numberReads) {
      mReadSequences = new long[(int) (2 * numberReads)];
    }

    @Override
    public int windowSize() {
      return 0;
    }

    @Override
    public void readAll(final int readId, final boolean reverse) {
      append(mOut, "readId " + readId + LS);
    }

    @Override
    public int readLength() {
      return 4;
    }

    @Override
    public void setValues(final int id2, final boolean reverse) {
      if (mReadSequences == null) {
        throw new NullPointerException();
      }
      if (!(id2 >= 0 && id2 < mReadSequences.length)) {
        throw new RuntimeException("id2=" + id2 + " length=" + mReadSequences.length);
      }
      append(mOut, "set " + id2 + LS);
    }

    /**
     */
    @Override
    public NgsHashFunction threadClone(final HashingRegion region) {
      mClones++;
      try {
        return (HashFunctionMock) super.clone();
      } catch (final CloneNotSupportedException e) {
        throw new RuntimeException(e);
      }
    }

    @Override
    public void threadFinish() {
    }

    /**
     * Standard semantics of clone.
     */
    @Override
    public HashFunctionMock clone() throws CloneNotSupportedException {
      return (HashFunctionMock) super.clone();
    }

    @Override
    public void logStatistics() {
      // do nothing
    }
  }
  private static final String READS_1 = "" + ">r1" + LS + "actg" + LS;
  private static final String TEM_1 = "" + ">t1" + LS + "actg" + LS;
  private static final String EXPECTED_1 = "" + "bulding" + LS + "reset" + LS + "code 0" + LS + "code 1" + LS + "code 3" + LS + "code 2" + LS + "readId 0" + LS + "set 0" + LS + "reset" + LS + "searching" + LS + "reset" + LS + "0 4" + LS + "code 0" + LS + "code 1" + LS + "code 3" + LS + "code 2" + LS + "pf=3" + LS + StringUtils.repeat("code 0" + LS, 60) + "pr=3" + LS;

  public void test1() throws IOException {
    check(READS_1, TEM_1, EXPECTED_1);
  }
  private static final String READS_2 = "" + ">r1" + LS + "actg" + LS + ">r2" + LS + "actg" + LS + ">r3" + LS + "actg" + LS + ">r4" + LS + "actg" + LS + ">r5" + LS + "actg" + LS + ">r6" + LS + "actg" + LS + ">r7" + LS + "actg" + LS + ">r8" + LS + "actg" + LS + ">r9" + LS + "actg" + LS + ">r10" + LS + "actg" + LS;
  private static final String TEM_2 = "" + ">t1" + LS + "actg" + LS + ">t2" + LS + "accg" + LS + ">t3" + LS + "acgt" + LS + ">t4" + LS + "actg" + LS;
  private static final String EXPECTED_2 = "" + "bulding" + LS + "set 0" + LS + "set 1" + LS + "set 2" + LS + "set 3" + LS + "set 4" + LS + "set 5" + LS + "set 6" + LS + "set 7" + LS + "set 8" + LS + "set 9" + LS + "reset" + LS + "code 0" + LS + "code 1" + LS + "code 3" + LS + "code 2" + LS + "readId 0" + LS + "set 0" + LS + "reset" + LS + "reset" + LS + "code 0" + LS + "code 1" + LS + "code 3" + LS + "code 2" + LS + "readId 1" + LS + "set 1" + LS + "reset" + LS + "reset" + LS + "code 0" + LS + "code 1" + LS + "code 3" + LS + "code 2" + LS + "readId 2" + LS + "set 2" + LS + "reset" + LS + "reset" + LS + "code 0" + LS + "code 1" + LS + "code 3" + LS + "code 2" + LS + "readId 3" + LS + "set 3" + LS + "reset" + LS + "reset" + LS + "code 0" + LS + "code 1" + LS + "code 3" + LS + "code 2" + LS + "readId 4" + LS + "set 4" + LS + "reset" + LS + "reset" + LS + "code 0" + LS + "code 1" + LS + "code 3" + LS + "code 2" + LS + "readId 5" + LS + "set 5" + LS + "reset" + LS + "reset" + LS + "code 0" + LS + "code 1" + LS + "code 3" + LS + "code 2" + LS + "readId 6" + LS + "set 6" + LS + "reset" + LS + "reset" + LS + "code 0" + LS + "code 1" + LS + "code 3" + LS + "code 2" + LS + "readId 7" + LS + "set 7" + LS + "reset" + LS + "reset" + LS + "code 0" + LS + "code 1" + LS + "code 3" + LS + "code 2" + LS + "readId 8" + LS + "set 8" + LS + "reset" + LS + "reset" + LS + "code 0" + LS + "code 1" + LS + "code 3" + LS + "code 2" + LS + "readId 9" + LS + "set 9" + LS + "reset" + LS + "searching" + LS + "reset" + LS + "1 4" + LS + "code 0" + LS + "code 1" + LS + "code 1" + LS + "code 2" + LS + "pf=3" + LS + StringUtils.repeat("code 0" + LS, 60) + "pr=3" + LS + "reset" + LS + "2 4" + LS + "code 0" + LS + "code 1" + LS + "code 2" + LS + "code 3" + LS + "pf=3" + LS + StringUtils.repeat("code 0" + LS, 60) + "pr=3" + LS;

  /**
   * Tests progress and start end values.
   */
  public void test2() throws IOException {
    Diagnostic.setLogStream();
    try (ReaderParams tem = getTemplate(TEM_2)) {
      final ReaderParams reads = getReads(READS_2);
      try {
        final NgsHashLoopImpl hl = new NgsHashLoopImpl(reads.reader().numberSequences(), true, 0L, 0x3L);
        hl.integrity();

        final StringWriter sb = new StringWriter();
        final ByteArrayOutputStream byteOut = new ByteArrayOutputStream();
        final PrintStream progressOut = new PrintStream(byteOut);

        final CliDiagnosticListener clir = new CliDiagnosticListener(progressOut);
        Diagnostic.addListener(clir);

        final TemplateHashFunction thf = new HashFunctionMock(sb);
        try {
          hl.templateLoop(new MockSequenceParams(tem, 0, tem.reader().numberSequences()), thf);
          fail();
        } catch (final RuntimeException e) {
          assertEquals("Read sequences not defined", e.getMessage());
        }

        final ReadHashFunction rhf = new HashFunctionMock(sb);
        final ReadHashFunction rhf36 = new HashFunctionMock(sb) {

          @Override
          public int readLength() {
            return 36;
          }
        };

        sb.append("bulding").append(LS);
        rhf36.setReadSequences(reads.reader().numberSequences());
        hl.readLoop(new MockSequenceParams(reads, 0, reads.reader().numberSequences()), rhf36, ReadEncoder.SINGLE_END, false);
        rhf.setReadSequences(reads.reader().numberSequences());
        hl.readLoop(new MockSequenceParams(reads, 0, reads.reader().numberSequences()), rhf, ReadEncoder.SINGLE_END, false);
        sb.append("searching").append(LS);
        hl.templateLoop(new MockSequenceParams(tem, 1, tem.reader().numberSequences() - 1), thf);
        Diagnostic.removeListener(clir);
        assertEquals(sb.toString(), EXPECTED_2, sb.toString());
        progressOut.close();
        final String pr = byteOut.toString();
        //System.err.println(pr.replaceAll("\r", "\\r").replaceAll("\n", "\\n"));
        assertTrue(pr.contains("For read \"r1\" the length (4) does not agree with the length expected for all reads (36). It will be ignored."));
        assertTrue(pr.contains("For read \"r2\" the length (4) does not agree with the length expected for all reads (36). It will be ignored."));
      } finally {
        reads.close();
      }
    }
  }
  private static final String EXPECTED_3 = "" + "reset" + LS + "code 0" + LS + "code 0" + LS + "code 0" + LS + "code 0" + LS + "readId 0" + LS + "set 0" + LS + "reset" + LS + "reset" + LS + "code 0" + LS + "code 0" + LS + "code 0" + LS + "code 3" + LS + "readId 1" + LS + "set 1" + LS + "reset" + LS + "reset" + LS + "code 0" + LS + "code 0" + LS + "code 1" + LS + "code 1" + LS + "readId 2" + LS + "set 2" + LS + "reset" + LS + "reset" + LS + "code 0" + LS + "code 0" + LS + "code 1" + LS + "code 3" + LS + "readId 3" + LS + "set 3" + LS + "reset" + LS;

  /**
   * Tests progress and start end values.
   */
  public void testPcr() throws IOException {
    Diagnostic.setLogStream();
    final long[] hashes = {0, 3, 5, 7};
    final NgsHashLoopImpl hl = new NgsHashLoopImpl(hashes.length, true, 0L, 0x3L);

    hl.integrity();

    final StringWriter sb = new StringWriter();
    final ReadHashFunction rhf = new HashFunctionMock(sb);

    hl.pcrLoop(hashes, rhf);
    assertEquals(EXPECTED_3, sb.toString());
  }

  public void testMakeBuffer() throws IOException {
    Diagnostic.setLogStream();
    final ByteArrayOutputStream byteOut = new ByteArrayOutputStream();
    final PrintStream progressOut = new PrintStream(byteOut);

    final CliDiagnosticListener clir = new CliDiagnosticListener(progressOut);
    Diagnostic.addListener(clir);
    try {
      NgsHashLoopImpl.makeBuffer(new ReaderLongMock(Integer.MAX_VALUE));
      fail();
    } catch (final SlimException e) {
      //expected
    }
    Diagnostic.removeListener(clir);
    progressOut.close();
    //    assertEquals(""
    //                 + "There is a sequence which is too long to process. Its length is \"2147483647\" bytes. See the preread output for the name of the sequence." + LS
    //                 + "RTG has encountered a difficulty, please contact support@realtimegenomics.com" + LS,
    //                 byteOut.toString()
    //                 );
    assertTrue(byteOut.toString().contains("There is a sequence which is too long to process. Its length is \"2147483647\" bytes"));
  }

  public void testContigSplitting() throws Exception {
    Diagnostic.setLogStream();
    final File dir = FileUtils.createTempDir("proteinhashloop", "test");
    try {
      final File readDir = new File(dir, "read");
      final File templateDir = new File(dir, "template");
      final String read = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
      final String read2 = "cccccccccccccccccccccccccccccccccc";
      final String template = "tttaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaccccccccccccccccccccccccccccccccccggg";
      ReaderTestUtils.getReaderDNA(">a\n" + read + "\n>b\n" + read2, readDir, null).close();
      ReaderTestUtils.getReaderProtein(">b\n" + template + "\n>c\n" + template + "\n>d\n" + template + "\n>e\n" + template, templateDir).close();
      final NgsHashLoopImpl loop = new NgsHashLoopImpl(1, false);
      loop.mMinChunkSize = 2;
      //final FakeProteinMask mask = new FakeProteinMask(new Skeleton(adjLength, adjLength, 0, 0, 1), new FakeReadCall(), new ImplementHashFunctionTest.TemplateCallMock());

      final StringWriter sb = new StringWriter();
      final TemplateHashFunction thf = new HashFunctionMock(sb);
      final ReadHashFunction rhf = new HashFunctionMock(sb);
      rhf.setReadSequences(2);

      try (ISequenceParams readParams = SequenceParams.builder().directory(readDir).mode(SequenceMode.BIDIRECTIONAL).create()) {
        loop.readLoop(readParams, rhf, ReadEncoder.SINGLE_END, false);
      }
      try (ISequenceParams templateParams = SequenceParams.builder().directory(templateDir).mode(SequenceMode.BIDIRECTIONAL).create()) {
        loop.templateLoopMultiCore(templateParams, (NgsHashFunction) thf, 2, 10);
      }
      assertTrue(((HashFunctionMock) thf).mClones >= 4);
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testTemplateLoopMultiCore() throws Exception {
    Diagnostic.setLogStream();
    final File dir = FileUtils.createTempDir("proteinhashloop", "test");
    try {
      final File readDir = new File(dir, "read");
      final File templateDir = new File(dir, "template");
      final String read = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
      final String template = "KWRKNRKSKKNQRNYNHDAADA";
      ReaderTestUtils.getReaderDNA(">a\n" + read, readDir, null).close();
      ReaderTestUtils.getReaderProtein(">b\n" + template + "\n>c\n" + template + "\n>d\n" + template + "\n>e\n" + template, templateDir).close();
      final NgsHashLoopImpl loop = new NgsHashLoopImpl(1, false);
      loop.mMinChunkSize = 2;
      //final FakeProteinMask mask = new FakeProteinMask(new Skeleton(adjLength, adjLength, 0, 0, 1), new FakeReadCall(), new ImplementHashFunctionTest.TemplateCallMock());

      final StringWriter sb = new StringWriter();
      final TemplateHashFunction thf = new HashFunctionMock(sb);
      final ReadHashFunction rhf = new HashFunctionMock(sb);
      rhf.setReadSequences(1);

      try (ISequenceParams readParams = SequenceParams.builder().directory(readDir).mode(SequenceMode.TRANSLATED).create()) {
        loop.readLoop(readParams, rhf, ReadEncoder.SINGLE_END, false);
      }
      try (ISequenceParams templateParams = SequenceParams.builder().directory(templateDir).mode(SequenceMode.PROTEIN).create()) {
        loop.templateLoopMultiCore(templateParams, (NgsHashFunction) thf, 4, HashingRegion.DEFAULT_THREAD_MULTIPLIER);
      }
      assertTrue("" + ((HashFunctionMock) thf).mClones, ((HashFunctionMock) thf).mClones > 4);
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  //We don't seem to need hash priming...
  /*
  public void testHashPriming() {
    // Check that priming the hash prevents non mid region starts from generating invalid hits.
    Diagnostic.setLogStream();
    final File dir = FileUtils.createTempDir("proteinhashloop", "test");
    try {
      final File templateDir = new File(dir, "template");
      final String template = "ccccccccctccccccccctgctctcccccccgcccccccccccccccccccccccccccccccaccccccccccacccaccccgggggggggtttttttttttccccctttttttggggg";
      ReaderTestUtils.getReaderDNA(">b\n" + template, templateDir).close();
      final NgsHashLoopImpl loop = new NgsHashLoopImpl(1, false, null);
      loop.mMinChunkSize = 2;
      //final FakeProteinMask mask = new FakeProteinMask(new Skeleton(adjLength, adjLength, 0, 0, 1), new FakeReadCall(), new ImplementHashFunctionTest.TemplateCallMock());

      final ReadCall readCall = new ReadCallMock();
      final TemplateCallMock singleTemplateCall = new TemplateCallMock();
      final TemplateHashFunction singleThf = new CGMaska1b1(readCall, singleTemplateCall);

      final ISequenceParams singleTemplateParams = SequenceParams.builder().directory(templateDir).mode(SequenceMode.BIDIRECTIONAL).create();
      try {
        final NgsHashLoopImpl.SequenceLoop sequenceLoop = new NgsHashLoopImpl.SequenceLoop(loop, singleTemplateParams, (NgsHashFunction) singleThf, new HashingRegion(0, 0, 0, template.length() + 1), "asdf", 0);
        sequenceLoop.run();
      } finally {
        singleTemplateParams.close();
      }

      final NgsHashLoopImpl loop2 = new NgsHashLoopImpl(1, false, null);
      loop2.mMinChunkSize = 2;
      final TemplateCallMock templateCall = new TemplateCallMock();
      final TemplateHashFunction thf = new CGMaska1b1(readCall, templateCall);

      final ISequenceParams templateParams = SequenceParams.builder().directory(templateDir).mode(SequenceMode.BIDIRECTIONAL).create();
      try {
        final NgsHashLoopImpl.SequenceLoop sequenceLoop = new NgsHashLoopImpl.SequenceLoop(loop2, templateParams, (NgsHashFunction) thf, new HashingRegion(0, 18, 0, template.length() + 1), "asdf", 0);
        sequenceLoop.run();
      } finally {
        templateParams.close();
      }
      for (final Called c : templateCall.mCalls) {
        assertTrue(singleTemplateCall.mCalls.contains(c));
      }
      for (final Called c : singleTemplateCall.mCalls) {
        if (c.mEndPosition > 18) {
          assertTrue(templateCall.mCalls.contains(c));
        }
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  private class ReadCallMock implements ReadCall {

    @Override
    public void readCall(final int id, final long hash, final int index) {
    }

  }


  private static class Called {
    int mEndPosition;
    long mHash;
    int mIndex;
    public Called(final int pos, final long hash, final int index) {
      mEndPosition = pos;
      this.mHash = hash;
      this.mIndex = index;
    }

    @Override
    public boolean equals(final Object other) {
      if (other == null) {
        return false;
      }
      if (!getClass().equals(other.getClass())) {
        return false;
      }
      final Called call = (Called) other;
      return mEndPosition == call.mEndPosition && mHash == call.mHash && mIndex == call.mIndex;
    }

    @Override
    public int hashCode() {
      int hash = 3;
      hash = 89 * hash + this.mEndPosition;
      hash = 89 * hash + (int) (this.mHash ^ (this.mHash >>> 32));
      hash = 89 * hash + this.mIndex;
      return hash;
    }

  }*/
  /*
  private class TemplateCallMock implements TemplateCall, Cloneable {

    List<Called> mCalls = new ArrayList<>();

    @Override
    public boolean globalIntegrity() {
      throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean integrity() {
      throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void dumpString(final StringBuilder sb) {
      throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void toString(final StringBuilder sb) {
      throw new UnsupportedOperationException("Not supported yet.");
    }

    /*
    @Override
    public TemplateCall clone() {
      return (TemplateCall) super.clone();
    }

    @Override
    public void done() {
    }

    public void endSequence() {
    }

    @Override
    public boolean isReverse() {
      return true;
    }

    @Override
    public void logStatistics() {
    }

    @Override
    public void set(final long name, final int length) {
    }

    @Override
    public void setHashFunction(final NgsHashFunction hashFunction) {
    }

    @Override
    public void setReverse(final boolean reverse) {
    }

    @Override
    public void templateCall(final int endPosition, final long hash, final int index) {
      mCalls.add(new Called(endPosition, hash, index));
    }

    @Override
    public TemplateCall threadClone(final HashingRegion region) {
      return this;
    }

    @Override
    public void threadFinish() {
    }

  }
  */
}

