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
package com.rtg.sam;

import static com.rtg.sam.SharedSamConstants.SAMHEADER1;
import static com.rtg.util.StringUtils.LS;
import static com.rtg.util.StringUtils.TAB;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;

import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.MockSequenceParams;
import com.rtg.launcher.NoStatistics;
import com.rtg.mode.SequenceMode;
import com.rtg.mode.SequenceType;
import com.rtg.reader.MockSequencesReader;
import com.rtg.reader.ReaderUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.CliDiagnosticListener;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.DiagnosticListener;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LogFile;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;
import junit.framework.TestCase;

/**
 */
public class SamIteratorTaskTest extends TestCase {

  static final class DummySingleMappedParams extends SingleMappedParams {
    protected DummySingleMappedParams(Collection<File> mapped, ISequenceParams genome, int threads) {
      super(new SingleMappedParamsTest.MockSingleMappedParamsBuilder().name("DummyBuilder").mapped(mapped).genome(genome).ioThreads(threads).execThreads(1));
    }
    @Override
    public File file(String name) {
      return null;
    }
    @Override
    public File directory() {
      return null;
    }
  }

  static class DummySamIteratorTask extends SamIteratorTask<DummySingleMappedParams, NoStatistics> {
    protected DummySamIteratorTask(DummySingleMappedParams params) {
      super(params, null, null, SamFilterParams.builder().create());
      assertEquals(-1, mPreviousStart);
      assertEquals(-2, mTemplateLength);
      assertNull(mTemplateName);
    }
    @Override
    protected void finalPostFlush() {
    }
    @Override
    public int flush(int start, int last) {
      return 0;
    }
    @Override
    protected boolean processRecord(SAMRecord rec) {
      mTemplateName = rec.getReferenceName();
      mPreviousStart = rec.getAlignmentStart();
      mTemplateLength = 100000;
      return !rec.getReferenceName().equals("processwarning");
    }
  }

  File mTempDir;

  @Override
  public void setUp() throws Exception {
    Diagnostic.setLogStream();
    mTempDir = FileUtils.createTempDir("samitertask", "test");
  }

  @Override
  public void tearDown() {
    Diagnostic.setLogStream();
    assertTrue(FileHelper.deleteAll(mTempDir));
    mTempDir = null;
  }

  static final String SAMINVALID_IH = ""
    + SAMHEADER1
    + "0" + TAB + "0" + TAB + "g1" + TAB +  "3" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACTG" + TAB + "````````" + TAB + "AS:i:0" + TAB + "IH:i:0" + LS
    + "0" + TAB + "0" + TAB + "g1" + TAB +  "3" + TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "ATCGACTG" + TAB + "````````" + TAB + "AS:i:0" + TAB + "IH:i:1" + LS;
  public void testStaticMethods() throws IOException {
    //validateRecord
    final SamReader samReader = SamUtils.makeSamReader(new ByteArrayInputStream(SAMINVALID_IH.getBytes()));
    final CloseableIterator<SAMRecord> records = samReader.iterator();
    DummySamIteratorTask task = new DummySamIteratorTask(new DummySingleMappedParams(new ArrayList<File>(), null, 0));
    assertFalse(task.validateRecord(records.next()));
    assertTrue(task.validateRecord(records.next()));

    final SequencesReader reader = new MockSequencesReader(SequenceType.DNA, 10, 10);
    final Map<String, Long> names = ReaderUtils.getSequenceNameMap(reader);

    //getTemplate
    reader.seek(0);
    assertEquals(1, DummySamIteratorTask.getTemplate(reader).length);

    //getSequenceBytes
    reader.seek(-1);
    assertNull(DummySamIteratorTask.getSequenceBytes(reader, names, "blah"));
    assertEquals(1, DummySamIteratorTask.getSequenceBytes(reader, names, "seq2").length);
    assertEquals(2, reader.currentSequenceId());

    //getSequenceLength
    reader.seek(-1);
    assertEquals(-1, DummySamIteratorTask.getSequenceLength(reader, names, "blah"));
    assertEquals(1, DummySamIteratorTask.getSequenceLength(reader, names, "seq3"));
    assertEquals(3, reader.currentSequenceId());
  }

  static final String SAMHD = "@HD" + TAB + "VN:1.0" + TAB + "SO:coordinate" + LS;
  static final String SAMSQS = "@SQ" + TAB + "SN:";
  static final String SAMSQE = TAB + "LN:";
  static final String SAMFLAGS = TAB + "0" + TAB;
  static final String SAMREST = TAB + "255" + TAB + "8M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "TTCAGCTA" + TAB + "````````" + TAB + "AS:i:1" + TAB + "IH:i:1" + LS;

  static final String SAMSEQUENCES1 = ""
    + SAMHD
    + SAMSQS + "g1" + SAMSQE + (DummySamIteratorTask.PROGRESS_NT * 5) + LS
    + SAMSQS + "gempty" + SAMSQE + "0" + LS
    + SAMSQS + "processwarning" + SAMSQE + "10" + LS
    + "0" + SAMFLAGS + "g1" + TAB + "1" + SAMREST
    + "1" + SAMFLAGS + "g1" + TAB + "2" + SAMREST
    + "2" + SAMFLAGS + "g1" + TAB + "3" + SAMREST
    + "3" + SAMFLAGS + "g1" + TAB + "4" + SAMREST
    + "4" + SAMFLAGS + "g1" + TAB + "5" + SAMREST
    + "5" + SAMFLAGS + "g1" + TAB + (DummySamIteratorTask.PROGRESS_NT >> 1) + SAMREST
    + "6" + SAMFLAGS + "g1" + TAB + (DummySamIteratorTask.PROGRESS_NT + 2) + SAMREST
    + "7" + SAMFLAGS + "g1" + TAB + (DummySamIteratorTask.PROGRESS_NT * 3 + 5) + SAMREST
    + "8" + SAMFLAGS + "g1" + TAB + (DummySamIteratorTask.PROGRESS_NT * 3 + 6) + SAMREST
    + "filter" + SAMFLAGS + "g1" + TAB + (DummySamIteratorTask.PROGRESS_NT * 3 + 7) + SAMREST
    + "0" + SAMFLAGS + "processwarning" + TAB + "1" + SAMREST
    ;
  static final String SAMSEQUENCES2 = ""
    + SAMSEQUENCES1
    + "1" + SAMFLAGS + "processwarning" + TAB + "2" + SAMREST
    + "2" + SAMFLAGS + "processwarning" + TAB + "3" + SAMREST
    + "3" + SAMFLAGS + "processwarning" + TAB + "4" + SAMREST
    + "4" + SAMFLAGS + "processwarning" + TAB + "5" + SAMREST
    + "5" + SAMFLAGS + "processwarning" + TAB + "6" + SAMREST
    + "6" + SAMFLAGS + "processwarning" + TAB + "7" + TAB + "255" + TAB + "8MM" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "TTCAGCTA" + TAB + "````````" + TAB + "AS:i:1" + TAB + "IH:i:1" + LS
    ;

  public void testExec() throws IOException {
    final Collection<File> mapped = new ArrayList<>();
    final File samFile = FileUtils.stringToFile(SAMHEADER1 + "blah blah" + LS, new File(mTempDir, "sam.sam"));
    mapped.add(samFile);
    DummySamIteratorTask task = new DummySamIteratorTask(new DummySingleMappedParams(mapped, null, 0));

    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try {
      try {
        task.exec();
        fail();
      } catch (final NoTalkbackSlimException e) {
        e.logException();
        assertEquals(ErrorType.SAM_BAD_FORMAT, e.getErrorType());
        assertTrue(mps.toString(), mps.toString().contains("Error: SAM record has an irrecoverable problem in file " + samFile.getPath() + ". Error parsing text SAM file. Not enough fields; Line 6"));
      }
    } finally {
      Diagnostic.setLogStream();
    }
    assertTrue(samFile.delete());
    final File outputFile = new File(mTempDir, "output");
    final File progressFile = new File(mTempDir, "progress");
    Diagnostic.setLogStream(new LogFile(outputFile));
    FileUtils.stringToFile(SAMSEQUENCES2, samFile);
    task = new DummySamIteratorTask(new DummySingleMappedParams(mapped, null, 0));
    task.exec();
    //System.err.println(FileUtils.fileToString(progressFile));
    //System.err.println(FileUtils.fileToString(outputFile));
    TestUtils.containsAll(FileUtils.fileToString(outputFile),
        "threads = 0 : " + com.rtg.sam.MultifileIterator.class.getName(),
        " 7 records skipped because of SAM format problems.",
        "SAM record is invalid.",
        "Invalid record:",
        "processwarning\t1",
        "processwarning\t2",
        "processwarning\t3",
        "processwarning\t4",
        "processwarning\t5",
        " 10 records processed");
    assertFalse(FileUtils.fileToString(outputFile).contains("processwarning\t6"));
    TestUtils.containsAll(FileUtils.fileToString(progressFile),
        "Processed " + DummySamIteratorTask.PROGRESS_NT + " NT of g1",
        "Processed " + (DummySamIteratorTask.PROGRESS_NT * 3) + " NT of g1",
        "Reference g1 completed (1)",
        "Reference processwarning completed (2)");
    assertFalse(FileUtils.fileToString(progressFile).contains("g1: " + (DummySamIteratorTask.PROGRESS_NT * 2) + " Reference NT Processed"));
    assertFalse(FileUtils.fileToString(progressFile).contains("/"));
    assertEquals(4, FileUtils.fileToString(progressFile).split(LS).length);
    Diagnostic.setLogStream();
    assertTrue(outputFile.delete());
    assertTrue(progressFile.delete());
    assertTrue(samFile.delete());
    FileUtils.stringToFile(SAMSEQUENCES1, samFile);
    final MemoryPrintStream listenerBuffer = new MemoryPrintStream();
    final DiagnosticListener listener = new CliDiagnosticListener(listenerBuffer.printStream());
    Diagnostic.addListener(listener);
    task = new DummySamIteratorTask(new DummySingleMappedParams(mapped, null, 0));
    task.exec();
    Diagnostic.removeListener(listener);
    assertTrue(listenerBuffer.toString(), listenerBuffer.toString().contains("1 records skipped because of SAM format problems."));
  }

  static final String SAMSEQUENCES3 = ""
    + SAMHD
    + SAMSQS + "seq0" + SAMSQE + "10000" + LS
    + SAMSQS + "seq1" + SAMSQE + "2300" + LS
    + "0" + SAMFLAGS + "seq0" + TAB + "124" + SAMREST
    + "1" + SAMFLAGS + "seq0" + TAB + "369" + SAMREST
    + "0" + SAMFLAGS + "seq1" + TAB + "123" + SAMREST
    + "1" + SAMFLAGS + "seq1" + TAB + "137" + SAMREST
    + "2" + SAMFLAGS + "seq1" + TAB + "160" + SAMREST
    ;

  static final String SAMSEQUENCES4 = ""
    + SAMHD
    + SAMSQS + "seq0" + SAMSQE + "10000" + LS
    + SAMSQS + "seq1" + SAMSQE + "2300" + LS
    + SAMSQS + "seqfake" + SAMSQE + "2300" + LS
    ;

  public void testExecWithGenome() throws IOException {
    final Collection<File> mapped = new ArrayList<>();
    final File samFile = FileUtils.stringToFile(SAMSEQUENCES3, new File(mTempDir, "sam.sam"));
    mapped.add(samFile);
    final MockSequencesReader reader = new MockSequencesReader(SequenceType.DNA, 2, 12300);
    reader.setLengths(new int[] {10000, 2300});
    final ISequenceParams genome = new MockSequenceParams(reader, SequenceMode.BIDIRECTIONAL);
    final File outputFile = new File(mTempDir, "output");
    final File progressFile = new File(mTempDir, "progress");
    Diagnostic.setLogStream(new LogFile(outputFile));
    final DummySamIteratorTask task = new DummySamIteratorTask(new DummySingleMappedParams(mapped, genome, 1));
    task.exec();
    //System.err.println(FileUtils.fileToString(progressFile));
    //System.err.println(FileUtils.fileToString(outputFile));
    TestUtils.containsAll(FileUtils.fileToString(outputFile),
        "threads = 1 : " +  com.rtg.sam.ThreadedMultifileIterator.class.getName(),
        " 0 records skipped because of SAM format problems.",
        " 5 records processed");
    TestUtils.containsAll(FileUtils.fileToString(progressFile),
        "Processed 1% of seq0",
        "Processed 1% of total",
        "Processed 3% of seq0",
        "Processed 3% of total",
        "Reference seq0 completed (1/2)",
        "Processed 5% of seq1",
        "Processed 82% of total",
        "Processed 6% of seq1",
        "Reference seq1 completed (2/2)"
        );


    // Should fail, as we have an @SQ that is not in the genome
    final Collection<File> mapped2 = new ArrayList<>();
    final File samFile2 = FileUtils.stringToFile(SAMSEQUENCES4, new File(mTempDir, "sam.sam"));
    mapped2.add(samFile2);
    final DummySamIteratorTask task2 = new DummySamIteratorTask(new DummySingleMappedParams(mapped2, genome, 1));
    try {
      task2.exec();
      fail();
    } catch (NoTalkbackSlimException e) {
      // Expected
    }

  }
}
