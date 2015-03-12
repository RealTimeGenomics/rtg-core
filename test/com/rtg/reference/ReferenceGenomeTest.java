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

package com.rtg.reference;

import static com.rtg.util.StringUtils.LS;
import static com.rtg.util.StringUtils.TAB;

import java.io.File;
import java.io.IOException;
import java.io.Reader;
import java.io.StringReader;
import java.util.Iterator;
import java.util.Map;

import com.rtg.mode.SequenceType;
import com.rtg.reader.MockSequencesReader;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.ReaderUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.reference.ReferenceGenome.DefaultFallback;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class ReferenceGenomeTest extends TestCase {

  public void testNames() throws IOException {
    Diagnostic.setLogStream();
    final File out = FileUtils.createTempDir("reference", "test1");
    try {
      final SequencesReader sr = ReaderTestUtils.getReaderDNA(">s1" + LS + "acgt" + LS + ">s2" + LS + "ac" + LS, out, null);
      final Map<String, Integer> names = ReaderUtils.getSequenceLengthMap(sr);
      assertEquals(2, names.size());
      assertEquals(4, (int) names.get("s1"));
      assertEquals(2, (int) names.get("s2"));
      assertEquals(null, names.get("xx"));
      sr.close();
    } finally {
      assertTrue(FileHelper.deleteAll(out));
    }
  }

  public void testNamesEmpty() throws IOException {
    Diagnostic.setLogStream();
    final File out = FileUtils.createTempDir("reference", "test1");
    try {
      final SequencesReader sr = ReaderTestUtils.getReaderDNA("", out, null);
      final Map<String, Integer> names = ReaderUtils.getSequenceLengthMap(sr);
      assertEquals(0, names.size());
      sr.close();
    } finally {
      assertTrue(FileHelper.deleteAll(out));
    }
  }

  public void test() throws IOException {
    final MemoryPrintStream ps = new MemoryPrintStream();
    Diagnostic.setLogStream(ps.printStream());
    try (TestDirectory out = new TestDirectory("referencetest1")) {
      try (SequencesReader sr = ReaderTestUtils.getReaderDNA(">s1" + LS + "acgt" + LS
          + ">s3" + LS + "ac" + LS
          + ">s2" + LS + "ac" + LS, out, null)) {
        final String refStr = ""
            + "#comment" + LS
            + "version" + TAB + "0" + LS
            + LS
            + "either" + TAB + "def" + TAB + "diploid" + TAB + "linear" + LS
            + "female" + TAB + "seq" + TAB + "s1" + TAB + "diploid" + TAB + "circular" + LS
            + "male" + TAB + "seq" + TAB + "s1" + TAB + "haploid" + TAB + "circular" + LS
            + "male" + TAB + "seq" + TAB + "s2" + TAB + "haploid" + TAB + "circular" + LS
            + LS
            + "male" + TAB + "dup" + TAB + "s1:3-4" + TAB + "s2:1-2";
        final Reader ref = new StringReader(refStr);
        final ReferenceGenome rg = new ReferenceGenome(sr, ref, Sex.MALE);
        final String s1 = ""
            + "s1 HAPLOID circular 4" + LS
            + "    s1:3-4  s2:1-2" + LS;
        final String s2 = ""
            + "s2 HAPLOID circular 2" + LS
            + "    s1:3-4  s2:1-2" + LS;
        final String s3 = ""
            + "s3 DIPLOID linear 2" + LS;
        final String exp = s1 + s3 + s2;

        assertEquals(exp, rg.toString());
        //System.err.println(ps.toString());
        assertFalse(ps.toString().contains("reference file"));
        assertEquals(s1, rg.sequence("s1").toString());
        assertTrue(rg.sequence("s1").isSpecified());
        assertFalse(rg.sequence("s3").isSpecified());
        assertEquals(null, rg.sequence("xx"));

        assertEquals(3, rg.sequences().size());
      } finally {
        ps.close();
      }
    }
  }

  public void testConstructor() throws IOException {
    final MemoryPrintStream ps = new MemoryPrintStream();
    Diagnostic.setLogStream(ps.printStream());
    final File dir = ReaderTestUtils.getDNADir(">s1" + LS + "acgt" + LS + ">s2" + LS + "acgt" + LS + ">s3" + LS + "acgt" + LS);
    try {
      try {
        new ReferenceGenome(new MockSequencesReader(SequenceType.DNA), Sex.EITHER);
      } catch (IOException e) {
        assertEquals("No reference file found", e.getMessage());
      }
      try (SequencesReader sr = SequencesReaderFactory.createDefaultSequencesReader(dir)) {
        try {
          new ReferenceGenome(sr, Sex.EITHER);
        } catch (IOException e) {
          assertEquals("No reference file contained in: " + dir.getPath(), e.getMessage());
        }
        ReferenceGenome rg = new ReferenceGenome(sr, Sex.EITHER, DefaultFallback.DIPLOID);
        final String[] rgStrings = new String[3];
        final StringBuilder rgString = new StringBuilder();
        Iterator<ReferenceSequence> it = rg.sequences().iterator();
        for (int i = 0; i < 3; i++) {
          rgStrings[i] = "s" + (i + 1) + " DIPLOID linear 4" + LS;
          rgString.append(rgStrings[i]);
          assertTrue(it.hasNext());
          assertEquals(rgStrings[i], it.next().toString());
        }
        assertEquals(rgString.toString(), rg.toString());
        assertFalse(it.hasNext());

        rg = new ReferenceGenome(sr, Sex.EITHER, DefaultFallback.HAPLOID);
        rgString.setLength(0);
        it = rg.sequences().iterator();
        for (int i = 0; i < 3; i++) {
          rgStrings[i] = "s" + (i + 1) + " HAPLOID linear 4" + LS;
          rgString.append(rgStrings[i]);
          assertTrue(it.hasNext());
          assertEquals(rgStrings[i], it.next().toString());
        }
        assertEquals(rgString.toString(), rg.toString());
        assertFalse(it.hasNext());
      }
    } finally {
      ps.close();
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testFromFile() throws IOException {
    final MemoryPrintStream ps = new MemoryPrintStream();
    Diagnostic.setLogStream(ps.printStream());
    final File dir = ReaderTestUtils.getDNADir(">s1" + LS + "acgt" + LS + ">s2" + LS + "ac" + LS + ">s3" + LS + "ac" + LS);
    try {
      try (SequencesReader sr = SequencesReaderFactory.createDefaultSequencesReader(dir)) {
        final String refStr = ""
            + "#comment" + LS
            + "version" + TAB + "0" + LS
            + LS
            + "either" + TAB + "def" + TAB + "diploid" + TAB + "linear" + LS
            + "female" + TAB + "seq" + TAB + "s1" + TAB + "diploid" + TAB + "circular" + LS
            + "male" + TAB + "seq" + TAB + "s1" + TAB + "haploid" + TAB + "circular" + LS
            + "male" + TAB + "seq" + TAB + "s2" + TAB + "haploid" + TAB + "circular" + LS
            + LS
            + "male" + TAB + "dup" + TAB + "s1:3-4" + TAB + "s2:1-2";
        FileUtils.stringToFile(refStr, new File(dir, "reference.txt"));
        final ReferenceGenome rg = new ReferenceGenome(sr, Sex.MALE);
        final String s1 = ""
            + "s1 HAPLOID circular 4" + LS
            + "    s1:3-4  s2:1-2" + LS;
        final String s2 = ""
            + "s2 HAPLOID circular 2" + LS
            + "    s1:3-4  s2:1-2" + LS;
        final String s3 = ""
            + "s3 DIPLOID linear 2" + LS;
        final String exp = s1 + s2 + s3;

        assertEquals(exp, rg.toString());
        //System.err.println(ps.toString());
        assertFalse(ps.toString().contains("reference file"));
        assertEquals(s1, rg.sequence("s1").toString());
        assertEquals(null, rg.sequence("xx"));

        assertEquals(3, rg.sequences().size());
      }
    } finally {
      ps.close();
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  //leaves default out
  public void testBad() throws IOException {
    final MemoryPrintStream ps = new MemoryPrintStream();
    Diagnostic.setLogStream(ps.printStream());
    try (TestDirectory out = new TestDirectory("referencetest1")) {
      try (SequencesReader sr = ReaderTestUtils.getReaderDNA(">s1" + LS + "acgt" + LS + ">s2" + LS + "ac" + LS + ">s3" + LS + "ac" + LS, out, null)) {
        final String refStr = ""
            + "#comment" + LS
            + "version" + TAB + "0" + LS
            + LS
            + "female" + TAB + "seq" + TAB + "s1" + TAB + "diploid" + TAB + "circular" + LS
            + "male" + TAB + "seq" + TAB + "s1" + TAB + "haploid" + TAB + "circular" + LS
            + "male" + TAB + "seq" + TAB + "s2" + TAB + "haploid" + TAB + "circular" + LS
            + LS
            + "male" + TAB + "dup" + TAB + "s1:3-4" + TAB + "s2:1-2";
        final Reader ref = new StringReader(refStr);
        try {
          new ReferenceGenome(sr, ref, Sex.MALE);
          fail();
        } catch (final NoTalkbackSlimException e) {
          assertEquals("Invalid reference file (see earlier warning messages).", e.getMessage());
        }
        final String pstr = ps.toString();
        assertTrue(pstr.contains("No default specified but required for sequence:"));
      } finally {
        ps.close();
      }
    }
  }


}
