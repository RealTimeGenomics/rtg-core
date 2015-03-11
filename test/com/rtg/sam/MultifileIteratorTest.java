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

import static com.rtg.util.StringUtils.LS;
import static com.rtg.util.StringUtils.TAB;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.CliDiagnosticListener;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import junit.framework.TestCase;

/**
 */
public class MultifileIteratorTest extends TestCase {

  @Override
  public void setUp() throws IOException {
    Diagnostic.setLogStream();
  }
  private static final String OUT_SAM = "alignments";

  static File samFile(File outDir, String label) {
    return new File(outDir, label + ".sam");
  }

  private void check(final String errExpected, String[] orderExpected, final String... samContents)
  throws IOException {
    final File alignmentsDir = FileUtils.createTempDir("sammultifileiterator", "check");
    Diagnostic.setLogStream();
    Diagnostic.clearListeners();
    final MemoryPrintStream out = new MemoryPrintStream();
    final CliDiagnosticListener listener = new CliDiagnosticListener(out.printStream());
    try {
      Diagnostic.addListener(listener);
      final ArrayList<File> files = new ArrayList<>();
      try {
        int i = 0;
        for (final String content : samContents) {
          final File ff = samFile(alignmentsDir, OUT_SAM + i);
          FileUtils.stringToFile(content, ff);
          files.add(ff);
          i++;
        }
        final RecordIterator<SAMRecord> smfi = getIterator(files);
        try {
          int ii = 0;
          while (smfi.hasNext()) {
            final SAMRecord sr = smfi.next();
            assertTrue("readName: " + sr.getReadName() + " orderExpected: " + orderExpected[ii], sr.getReadName().equals(orderExpected[ii]));
            ii++;
          }
          assertTrue(errExpected == null);
        } finally {
          smfi.close();
          assertEquals(Math.abs(smfi.getInvalidRecordsCount()), smfi.getInvalidRecordsCount());
        }
      } catch (final NoTalkbackSlimException e) {
        final String str = out.toString();
        //TestUtils.containsAll(errmsg, errExpected);

        //System.err.println("error \n" + str);
        assertTrue(str.contains(errExpected));
      } finally {
        Diagnostic.clearListeners();
      }
    } finally {
      FileHelper.deleteAll(alignmentsDir);
    }
    Diagnostic.clearListeners();
  }

  RecordIterator<SAMRecord> getIterator(Collection<File> files) throws IOException {
    return new MultifileIterator(new SamReadingContext(files, 1, new SamFilterParams.SamFilterParamsBuilder().create(), SamUtils.getUberHeader(files), null));
  }

  public void testSamMultifileIteratorTest() throws IOException {
    final File alignmentsDir = FileUtils.createTempDir("sammultifileiterator", "check");
    try {
      Diagnostic.setLogStream();

      final ArrayList<File> files = new ArrayList<>();

      final File ffa = samFile(alignmentsDir, OUT_SAM + "a");
      FileUtils.stringToFile(SAM_HEAD1 + SAM_REC_OK5 + SAM_REC_OK6, ffa);
      files.add(ffa);
      final File ffb = samFile(alignmentsDir, OUT_SAM + "b");
      FileUtils.stringToFile(SAM_HEAD1 + SAM_REC_OK3 + SAM_REC_OK4, ffb);
      files.add(ffb);
      final RecordIterator<?> smfi = getIterator(files);
      final SAMSequenceDictionary dict = smfi.header().getSequenceDictionary();
      assertEquals(0, dict.getSequenceIndex("gi"));
      try {
        assertTrue(smfi.hasNext());
        smfi.remove();
        fail();
      } catch (final UnsupportedOperationException e) {
        // expected

      } finally {
        smfi.close();
      }
      assertEquals(0, smfi.getInvalidRecordsCount());
    } finally {
      FileHelper.deleteAll(alignmentsDir);
    }
  }

  static final String SAM_HEAD1 = ""
    + "@HD" + TAB + "VN:1.0" + TAB + "SO:coordinate\n"
    + "@SQ" + TAB + "SN:gi" + TAB + "LN:30\n";

  static final String SAM_HEADWRONG = ""
    + "@HD" + TAB + "VN:1.0" + TAB + "SO:coordinate\n"
    + "@SQ" + TAB + "SN:wrong" + TAB + "LN:30\n";

  static final String SAM_HEAD_LONGER = ""
    + "@HD" + TAB + "VN:1.0" + TAB + "SO:coordinate\n"
    + "@SQ" + TAB + "SN:gi" + TAB + "LN:30\n"
    + "@SQ" + TAB + "SN:gj" + TAB + "LN:60\n";

  static final String SAM_REC_OK1 = ""
    + "read1" + TAB + "0" + TAB + "gi" + TAB + "2" + TAB + "255" + TAB + "10M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AAAAAAAAAA" + TAB +  "IB7?*III<I" + TAB + "AS:i:0" + TAB + "IH:i:1" + LS;
  static final String SAM_REC_OK2 = ""
    + "read2" + TAB + "0" + TAB + "gi" + TAB + "5" + TAB + "255" + TAB + "10M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AAAAAAAAAA" + TAB +  "IB7?*III<I" + TAB + "AS:i:0" + TAB + "IH:i:1" + LS;
  static final String SAM_REC_OK3 = ""
    + "read3" + TAB + "0" + TAB + "gi" + TAB + "7" + TAB + "255" + TAB + "10M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AAAAAAAAAA" + TAB +  "IB7?*III<I" + TAB + "AS:i:0" + TAB + "IH:i:1" + LS;
  static final String SAM_REC_OK4 = ""
    + "read4" + TAB + "0" + TAB + "gi" + TAB + "9" + TAB + "255" + TAB + "10M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AAAAAAAAAA" + TAB +  "IB7?*III<I" + TAB + "AS:i:0" + TAB + "IH:i:1" + LS;
  static final String SAM_REC_OK5 = ""
    + "read5" + TAB + "0" + TAB + "gi" + TAB + "11" + TAB + "255" + TAB + "10M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AAAAAAAAAA" + TAB +  "IB7?*III<I" + TAB + "AS:i:0" + TAB + "IH:i:1" + LS;
  static final String SAM_REC_OK6 = ""
    + "read6" + TAB + "0" + TAB + "gi" + TAB + "13" + TAB + "255" + TAB + "10M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AAAAAAAAAA" + TAB +  "IB7?*III<I" + TAB + "AS:i:0" + TAB + "IH:i:1" + LS;
  static final String SAM_REC_OK2WRONG = ""
    + "read2" + TAB + "0" + TAB + "wrong" + TAB + "5" + TAB + "255" + TAB + "10M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AAAAAAAAAA" + TAB +  "IB7?*III<I" + TAB + "AS:i:0" + TAB + "IH:i:1" + LS;

  //flags say unpaired but we supplied mate reference name
  static final String SAM_REC_BAD = ""
    + "badread" + TAB + "16" + TAB + "gi" + TAB + "3" + TAB + "255" + TAB + "10M" + TAB + "40" + TAB + "0" + TAB + "0" + TAB + "AAAAAAAAAA" + TAB + "GC@=I3IIII" + TAB + "AS:i:0" + TAB + "IH:i:0\n";

  public void testOk1() throws Exception {
    // no warning
    final String sam1 = SAM_HEAD1 + SAM_REC_OK1;
    final String sam2 = SAM_HEAD1 + SAM_REC_OK2;
    final String sam3 = SAM_HEAD1 + SAM_REC_OK3;
    final String[] order = {"read1", "read2", "read3"};
    check(null, order, sam1, sam2, sam3);
    final String sam4 = SAM_HEAD1 + SAM_REC_OK5 + SAM_REC_OK6;
    final String sam5 = SAM_HEAD1 + SAM_REC_OK3 + SAM_REC_OK4;
    final String sam6 = SAM_HEAD1 + SAM_REC_OK1 + SAM_REC_OK2;
    final String[] order2 = {"read1", "read2", "read3", "read4", "read5", "read6"};
    check(null, order2, sam4, sam5, sam6);
    final String sam7 = SAM_HEAD1 + SAM_REC_OK1 + SAM_REC_OK6;
    final String sam8 = SAM_HEAD1 + SAM_REC_OK3 + SAM_REC_OK4;
    final String sam9 = SAM_HEAD1 + SAM_REC_OK2 + SAM_REC_OK5;
    check(null, order2, sam7, sam8, sam9);
    final String sam10 = SAM_HEAD1 + SAM_REC_OK1 + SAM_REC_OK4 + SAM_REC_OK6;
    final String sam11 = SAM_HEAD1;
    final String sam12 = SAM_HEAD1 + SAM_REC_OK2 + SAM_REC_OK3 + SAM_REC_OK5;
    check(null, order2, sam10, sam11, sam12);
  }

  public void testSomeBad() throws Exception {
    // no warning
    final String[] order = {"read1", "read2", "read3", "read4", "read5", "read6"};
    final String sam1 = SAM_HEAD1 + SAM_REC_BAD + SAM_REC_OK1 + SAM_REC_OK4 + SAM_REC_OK6;
    final String sam2 = SAM_HEAD1;
    final String sam3 = SAM_HEAD1 + SAM_REC_OK2 + SAM_REC_OK3 + SAM_REC_BAD + SAM_REC_OK5;
    check(null, order, sam1, sam2, sam3);
  }

  public void testCheckHeader() throws Exception {
    // no warning
    final String sam1 = SAM_HEAD1 + SAM_REC_OK1;
    final String sam2 = SAM_HEADWRONG + SAM_REC_OK2WRONG;
    final File alignmentsDir = FileUtils.createTempDir("sammultifileiterator", "check");
    Diagnostic.setLogStream();
    Diagnostic.clearListeners();
    final MemoryPrintStream out = new MemoryPrintStream();
    final CliDiagnosticListener listener = new CliDiagnosticListener(out.printStream());
    try {
      Diagnostic.addListener(listener);
      final ArrayList<File> files = new ArrayList<>();
      final File file1 = samFile(alignmentsDir, OUT_SAM + 1);
      FileUtils.stringToFile(sam1, file1);
      files.add(file1);
      final File file2 = samFile(alignmentsDir, OUT_SAM + 2);
      FileUtils.stringToFile(sam2, file2);
      files.add(file2);
      try {
        getIterator(files);
        fail();
      } catch (final NoTalkbackSlimException e) {
        e.printErrorNoLog();
        final String str = out.toString();
        //System.err.println("error \n" + str);
        TestUtils.containsAll(str, file1.getPath(), file2.getPath(), "The SAM file", "cannot be merged with the SAM file", "because their headers are incompatible.", "Execution stopping because there were \"1\" SAM files with incompatible headers.");
      } finally {
        Diagnostic.clearListeners();
      }
    } finally {
      FileHelper.deleteAll(alignmentsDir);
    }
    Diagnostic.clearListeners();
  }

  public void testCheckHeader2() throws Exception {
    // no warning
    final String sam1 = SAM_HEAD_LONGER + SAM_REC_OK1;
    final String sam2 = SAM_HEAD1 + SAM_REC_OK1;
    final File alignmentsDir = FileUtils.createTempDir("sammultifileiterator", "check");
    Diagnostic.setLogStream();
    Diagnostic.clearListeners();
    final MemoryPrintStream out = new MemoryPrintStream();
    final CliDiagnosticListener listener = new CliDiagnosticListener(out.printStream());
    try {
      Diagnostic.addListener(listener);
      final ArrayList<File> files = new ArrayList<>();
      final File file1 = samFile(alignmentsDir, OUT_SAM + 1);
      FileUtils.stringToFile(sam1, file1);
      files.add(file1);
      final File file2 = samFile(alignmentsDir, OUT_SAM + 2);
      FileUtils.stringToFile(sam2, file2);
      files.add(file2);
      try {
        getIterator(files);
        fail();
      } catch (final NoTalkbackSlimException e) {
        e.printErrorNoLog();
        final String str = out.toString();
        //System.err.println("error \n" + str);
        final String[] exp = {""
            + "The SAM file \"" + file1.getPath()
            + "\" cannot be merged with the SAM file \"" + file2.getPath()
            + "\" because their headers are incompatible." + LS,
            "Execution stopping because there were \"1\" SAM files with incompatible headers."
            };
        TestUtils.containsAll(str, exp);
      } finally {
        Diagnostic.clearListeners();
      }
    } finally {
      FileHelper.deleteAll(alignmentsDir);
    }
    Diagnostic.clearListeners();
  }

  public void testCheckHeaderTwoWrong() throws Exception {
    // no warning
    final String sam1 = SAM_HEAD_LONGER + SAM_REC_OK1;
    final String sam2 = SAM_HEAD1 + SAM_REC_OK1;
    final String sam3 = SAM_HEAD_LONGER + SAM_REC_OK1;
    final String sam4 = SAM_HEADWRONG + SAM_REC_OK1;
    final File alignmentsDir = FileUtils.createTempDir("sammultifileiterator", "check");
    Diagnostic.setLogStream();
    Diagnostic.clearListeners();
    final MemoryPrintStream out = new MemoryPrintStream();
    final CliDiagnosticListener listener = new CliDiagnosticListener(out.printStream());
    try {
      Diagnostic.addListener(listener);
      final ArrayList<File> files = new ArrayList<>();
      final File file1 = samFile(alignmentsDir, OUT_SAM + 1);
      FileUtils.stringToFile(sam1, file1);
      files.add(file1);
      final File file2 = samFile(alignmentsDir, OUT_SAM + 2);
      FileUtils.stringToFile(sam2, file2);
      files.add(file2);
      final File file3 = samFile(alignmentsDir, OUT_SAM + 3);
      FileUtils.stringToFile(sam3, file3);
      files.add(file3);
      final File file4 = samFile(alignmentsDir, OUT_SAM + 4);
      FileUtils.stringToFile(sam4, file4);
      files.add(file4);
      try {
        getIterator(files);
        fail();
      } catch (final NoTalkbackSlimException e) {
        e.printErrorNoLog();
        final String str = out.toString();
        //System.err.println("error \n" + str);
        final String[] exp = {
          "The SAM file \"",
          "\" cannot be merged with the SAM file \"",
          "\" because their headers are incompatible." + LS,
          "Execution stopping because there were \"1\" SAM files with incompatible headers."
        };
        TestUtils.containsAll(str, exp);
      } finally {
        Diagnostic.clearListeners();
      }
    } finally {
      FileHelper.deleteAll(alignmentsDir);
    }
    Diagnostic.clearListeners();
  }

  public void testCheckHeaderFiveWrong() throws Exception {
    // no warning
    final String sam1 = SAM_HEAD_LONGER + SAM_REC_OK1;
    final String sam2 = SAM_HEAD1 + SAM_REC_OK1;
    final String sam3 = SAM_HEAD1 + SAM_REC_OK1;
    final String sam4 = SAM_HEADWRONG + SAM_REC_OK1;
    final String sam5 = SAM_HEAD1 + SAM_REC_OK1;
    final String sam6 = SAM_HEADWRONG + SAM_REC_OK1;
    final String sam7 = SAM_HEAD1 + SAM_REC_OK1;
     final File alignmentsDir = FileUtils.createTempDir("sammultifileiterator", "check");
    Diagnostic.setLogStream();
    Diagnostic.clearListeners();
    final MemoryPrintStream out = new MemoryPrintStream();
    final CliDiagnosticListener listener = new CliDiagnosticListener(out.printStream());
    try {
      Diagnostic.addListener(listener);
      final ArrayList<File> files = new ArrayList<>();
      final File file1 = samFile(alignmentsDir, OUT_SAM + 1);
      FileUtils.stringToFile(sam1, file1);
      files.add(file1);
      final File file2 = samFile(alignmentsDir, OUT_SAM + 2);
      FileUtils.stringToFile(sam2, file2);
      files.add(file2);
      final File file3 = samFile(alignmentsDir, OUT_SAM + 3);
      FileUtils.stringToFile(sam3, file3);
      files.add(file3);
      final File file4 = samFile(alignmentsDir, OUT_SAM + 4);
      FileUtils.stringToFile(sam4, file4);
      files.add(file4);
      final File file5 = samFile(alignmentsDir, OUT_SAM + 5);
      FileUtils.stringToFile(sam5, file5);
      files.add(file5);
      final File file6 = samFile(alignmentsDir, OUT_SAM + 6);
      FileUtils.stringToFile(sam6, file6);
      files.add(file6);
      final File file7 = samFile(alignmentsDir, OUT_SAM + 7);
      FileUtils.stringToFile(sam7, file7);
      files.add(file7);
      try {
        getIterator(files);
        fail();
      } catch (final NoTalkbackSlimException e) {
        e.printErrorNoLog();
        final String str = out.toString();
        //System.err.println("error \n" + str);
        final String[] exp = {
          "The SAM file \"",
          "\" cannot be merged with the SAM file \"",
          "\" because their headers are incompatible." + LS,
          "Execution stopping because there were \"1\" SAM files with incompatible headers."
        };
        TestUtils.containsAll(str, exp);
      } finally {
        Diagnostic.clearListeners();
      }
    } finally {
      FileHelper.deleteAll(alignmentsDir);
    }
    Diagnostic.clearListeners();
  }


  static final String SAM_REC_RG1 = ""
    + "read1" + TAB + "0" + TAB + "gi" + TAB + "2" + TAB + "255" + TAB + "10M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AAAAAAAAAA" + TAB +  "IB7?*III<I" + TAB + "AS:i:0" + TAB + "IH:i:1" + TAB + "RG:Z:1" + LS;
  static final String SAM_REC_RG2 = ""
    + "read1" + TAB + "0" + TAB + "gi" + TAB + "2" + TAB + "255" + TAB + "10M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + "AAAAAAAAAA" + TAB +  "IB7?*III<I" + TAB + "AS:i:0" + TAB + "IH:i:1" + TAB + "RG:Z:2" + LS;

  public void testMerging2ReadGroups() throws IOException {
    final File top = FileUtils.createTempDir("sam-multi-itr", "merge");
    try {
      final File first = new File(top, "first");
      FileUtils.stringToFile(SAM_HEAD1 + "@RG\tID:1\tSM:sample1\tPL:ILLUMINA" + LS + SAM_REC_RG1, first);
      final File second = new File(top, "second");
      FileUtils.stringToFile(SAM_HEAD1 + "@RG\tID:2\tSM:sample1\tPL:ILLUMINA" + LS + SAM_REC_RG2, second);
      final List<File> fileList = new ArrayList<>();
      fileList.add(first);
      fileList.add(second);
      try (RecordIterator<?> it = getIterator(fileList)) {
        final List<SAMReadGroupRecord> readGroups = it.header().getReadGroups();
        assertEquals(2, readGroups.size());
        boolean foundFirst = false;
        boolean foundSecond = false;
        for (final SAMReadGroupRecord r : readGroups) {
          if (r.getReadGroupId().equals("1")) {
            foundFirst = true;
          }
          if (r.getReadGroupId().equals("2")) {
            foundSecond = true;
          }
        }
        assertTrue(foundFirst && foundSecond);
      }
    } finally {
      assertTrue(FileHelper.deleteAll(top));
    }
  }

  public void testUnmappedRecords() throws IOException {
    try (final TestDirectory tmpDir = new TestDirectory()) {
      final File f = new File(tmpDir, "f.sam");
      FileUtils.stringToFile(SAM_HEAD1 + "r1\t77\t*\t0\t0\t*\t*\t0\t0\tTCAATNGTAG\tCCCFF#2=CF", f);
      final List<File> fileList = new ArrayList<>();
      fileList.add(f);
      final RecordIterator<?> it = new MultifileIterator(new SamReadingContext(fileList, 1, SamFilterParams.builder().excludeUnmapped(true).create(), SamUtils.getUberHeader(fileList), null));

      try {
        assertFalse(it.hasNext());
      } finally {
        it.close();
      }
      assertEquals(0, it.getInvalidRecordsCount());
      assertEquals(1, it.getFilteredRecordsCount());
      assertEquals(0, it.getOutputRecordsCount());
      assertEquals(1, it.getTotalRecordsCount());
    }
  }
}
