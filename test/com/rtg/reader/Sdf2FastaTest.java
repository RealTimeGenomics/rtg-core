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
package com.rtg.reader;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Locale;
import java.util.ResourceBundle;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.mode.DNAFastaSymbolTable;
import com.rtg.mode.ProteinFastaSymbolTable;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.TestCFlags;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.MockEventListener;

import junit.framework.Assert;

/**
 * Test class for corresponding class.
 */
public class Sdf2FastaTest extends AbstractCliTest {


  @Override
  public void setUp() throws IOException {
    super.setUp();
    Diagnostic.setLogStream();
  }

  @Override
  public void tearDown() throws IOException {
    super.tearDown();
    Diagnostic.setLogStream();
  }

  @Override
  protected AbstractCli getCli() {
    return new Sdf2Fasta();
  }

  private static final String JUNITOUT = ".junitout";

  public void testCliFlags() throws Exception {
    final CFlags flags = new CFlags();
    new Sdf2Fasta().initFlags(flags);
    assertNotNull(flags);
    final ResourceBundle resource = ResourceBundle.getBundle(FormatCli.PREREAD_RESOURCE_BUNDLE, Locale.getDefault());
    assertNotNull(flags.getFlag(resource.getString("FOUTPUT_FLAG")));
    //assertNotNull(flags.getFlag(resource.getString("CHUNKS_FLAG")));
    assertNotNull(flags.getFlag(resource.getString("LINE_FLAG")));
    final Flag inFlag = flags.getFlag(resource.getString("FINPUT_FLAG"));
    assertNotNull(inFlag);
    assertEquals(1, inFlag.getMinCount());
    assertEquals(1, inFlag.getMaxCount());
    final File testDir = FileUtils.createTempDir("PrereadToFastaTest", null);
    final File xx = new File(testDir, "xxnnxxnnxx");
    xx.deleteOnExit();
    try {
      assertTrue(xx.mkdir());
      flags.setFlags("-o", "xxdog", "-i", xx.getPath());
      assertEquals("xxdog", (String) flags.getValue(resource.getString("FOUTPUT_FLAG")));
      assertEquals(0, ((Integer) flags.getValue(resource.getString("LINE_FLAG"))).intValue());
      assertEquals(1, ((Integer) flags.getValue(resource.getString("CHUNKS_FLAG"))).intValue());
      assertEquals(xx.getPath(), ((File) inFlag.getValue()).getPath());
      assertTrue(xx.delete());
      assertTrue(xx.mkdir());
      flags.setFlags("-c", "2", "-l", "5", "-o", "xxdog", "-i", xx.getPath());
      assertEquals("xxdog", (String) flags.getValue(resource.getString("FOUTPUT_FLAG")));
      assertEquals(5, ((Integer) flags.getValue(resource.getString("LINE_FLAG"))).intValue());
      assertEquals(2, ((Integer) flags.getValue(resource.getString("CHUNKS_FLAG"))).intValue());
      assertEquals(xx.getPath(), ((File) inFlag.getValue()).getPath());
    } finally {
      assertTrue(xx.delete());
      FileHelper.deleteAll(testDir);
    }
  }

  private void compareToFile(final String str, final File f) throws IOException {
    final String main = FileUtils.fileToString(f);
    assertEquals(str, main);
  }

  private void checkContent1(File x) throws Exception {
    assertTrue(x.exists());
    //final BufferedReader r = new BufferedReader(new FileReader(x));
    try {
      compareToFile(">test" + StringUtils.LS + "ACGT" + StringUtils.LS + ">bob" + StringUtils.LS + "TAGTACCC" + StringUtils.LS + ">cat" + StringUtils.LS + "CAT" + StringUtils.LS + ">dog" + StringUtils.LS + "CCC" + StringUtils.LS, x);
    } finally {
      assertTrue(x.delete());
    }
  }

  private void checkContent2(File x) throws Exception {
    assertTrue(x.exists());
    try {
      compareToFile(">test" + StringUtils.LS + "ACG" + StringUtils.LS + "T" + StringUtils.LS + ">bob" + StringUtils.LS + "TAG" + StringUtils.LS + "TAC" + StringUtils.LS + "CC" + StringUtils.LS + ">cat" + StringUtils.LS + "CAT" + StringUtils.LS + ">dog" + StringUtils.LS + "CCC" + StringUtils.LS, x);
    } finally {
      assertTrue(x.delete());
    }
  }

  public void testWorks() throws Exception {
    final File xd = FileHelper.createTempDirectory();
    final File fastaDirectory = FileHelper.createTempDirectory();
    final File jUnitMy = new File(fastaDirectory, "junitmy");
    try {
      final ArrayList<InputStream> al = new ArrayList<>();
      al.add(new ByteArrayInputStream(">test\nacgt\n>bob\ntagt\naccc\n>cat\ncat\n>dog\nccc".getBytes()));
      final FastaSequenceDataSource ds = new FastaSequenceDataSource(al, new DNAFastaSymbolTable());
      final SequencesWriter sw = new SequencesWriter(ds, xd, 300000, PrereadType.UNKNOWN, false);
      sw.processSequences();
      File x;
      SequencesReader sr = SequencesReaderFactory.createDefaultSequencesReader(xd);
      try {
        x = new File(fastaDirectory, "junitmy.fasta");
        checkMainInitOk("-i", xd.toString(), "-o", x.toString(), "-Z");
        //Sdf2Fasta.doPrereadToFasta(sr, jUnitMy.getPath() ".fasta", 1, 0, false, false);
        checkContent1(x);
      } finally {
        sr.close();
      }
      sr = SequencesReaderFactory.createDefaultSequencesReader(xd);
      try {
        checkMainInitOk("-i", xd.toString(), "-o", x.toString(), "-l", "3", "-Z");
        //Sdf2Fasta.doPrereadToFasta(sr, jUnitMy.getPath(), ".fasta", 1, 3, false, false);
        checkContent2(x);
      } finally {
        sr.close();
      }
      sr = SequencesReaderFactory.createDefaultSequencesReader(xd);
      try {
        checkMainInitOk("-i", xd.toString(), "-o", jUnitMy.toString() + ".fasta", "-c", "2", "-Z");
        //Sdf2Fasta.doPrereadToFasta(sr, jUnitMy.getPath(), ".fasta", 2, 0, false, false);
        x = new File(fastaDirectory, "junitmy0.fasta");
        assertTrue(x.exists());
        x.deleteOnExit();
        try {
          compareToFile(">test" + StringUtils.LS + "ACGT" + StringUtils.LS + ">cat" + StringUtils.LS + "CAT" + StringUtils.LS, x);
        } finally {
          assertTrue(x.delete());
        }
        x = new File(fastaDirectory, "junitmy1.fasta");
        assertTrue(x.exists());
        x.deleteOnExit();
        try {
          compareToFile(">bob" + StringUtils.LS + "TAGTACCC" + StringUtils.LS + ">dog" + StringUtils.LS + "CCC" + StringUtils.LS, x);
        } finally {
          assertTrue(x.delete());
        }
      } finally {
        sr.close();
      }
      sr = SequencesReaderFactory.createDefaultSequencesReader(xd);
      try {
        x = new File(fastaDirectory, "junitmx0.fasta");
        checkMainInitOk("-i", xd.toString(), "-o", new File(fastaDirectory, "junitmx").toString(), "-c", "3", "-Z");
        //Sdf2Fasta.doPrereadToFasta(sr, new File(fastaDirectory, "junitmx").getPath(), ".fasta", 3, 0, false, false);
        assertTrue(x.exists());
        x.deleteOnExit();
        try {
          compareToFile(">test" + StringUtils.LS + "ACGT" + StringUtils.LS + ">dog" + StringUtils.LS + "CCC" + StringUtils.LS, x);
        } finally {
          assertTrue(x.delete());
        }
        x = new File(fastaDirectory, "junitmx1.fasta");
        assertTrue(x.exists());
        x.deleteOnExit();
        try {
          compareToFile(">bob" + StringUtils.LS + "TAGTACCC" + StringUtils.LS, x);
        } finally {
          assertTrue(x.delete());
        }
        x = new File(fastaDirectory, "junitmx2.fasta");
        assertTrue(x.exists());
        x.deleteOnExit();
        try {
          compareToFile(">cat" + StringUtils.LS + "CAT" + StringUtils.LS, x);
        } finally {
          assertTrue(x.delete());
        }
      } finally {
        sr.close();
      }
    } finally {
      assertTrue(FileHelper.deleteAll(xd));
      assertTrue(FileHelper.deleteAll(fastaDirectory));
    }
  }

  public void testBadArgs() throws Exception {
    final File xd = FileHelper.createTempDirectory();
    try {
      final ArrayList<InputStream> al = new ArrayList<>();
      final FastaSequenceDataSource ds = new FastaSequenceDataSource(al, new DNAFastaSymbolTable());
      final SequencesWriter sw = new SequencesWriter(ds, xd, 300000, PrereadType.UNKNOWN, false);
      sw.processSequences();
      checkHandleFlagsErr("-i", xd.toString(), "-o", "junithi", "-c", "0");
      checkHandleFlagsErr("-i", xd.toString(), "-o", "junithi", "-l", "-1");
    } finally {
      //      assertTrue(new File("junithi.fasta").delete());
      assertTrue(FileHelper.deleteAll(xd));
    }
  }

  public void testWorksColorspace() throws Exception {
    final File xd = FileHelper.createTempDirectory();
    final File fastaDirectory = FileHelper.createTempDirectory();
    try {
      final ArrayList<InputStream> al = new ArrayList<>();
      al.add(new ByteArrayInputStream(">test\nacgt\n>bob\ntagt\naccc\n>cat\ncat\n>dog\nccc".getBytes()));
      final FastaSequenceDataSource ds = new FastaSequenceDataSource(al, new DNAFastaSymbolTable());
      final SequencesWriter sw = new SequencesWriter(ds, xd, 300000, PrereadType.UNKNOWN, false);
      sw.processSequences();
      final File jUnitMy = new File(fastaDirectory, "junitmy");
      File x;
      SequencesReader sr = SequencesReaderFactory.createDefaultSequencesReader(xd);
      try {
        x = new File(fastaDirectory, "junitmy.fasta");
        checkMainInitOk("-i", xd.toString(), "-o", jUnitMy.getPath(), "--Xcolorspace", "-Z");
        assertTrue(x.exists());
        try {
          compareToFile(">test" + StringUtils.LS + "A131" + StringUtils.LS + ">bob" + StringUtils.LS + "T3213100" + StringUtils.LS + ">cat" + StringUtils.LS + "C13" + StringUtils.LS + ">dog" + StringUtils.LS + "C00" + StringUtils.LS, x);
        } finally {
          assertTrue(x.delete());
        }
      } finally {
        sr.close();
      }
      sr = SequencesReaderFactory.createDefaultSequencesReader(xd);
      try {
        x = new File(fastaDirectory, "junitmy.fasta");
        checkMainInitOk("-i", xd.toString(), "-o", jUnitMy.getPath(), "-l", "3", "--Xcolorspace", "-Z");
        //Sdf2Fasta.doPrereadToFasta(sr, jUnitMy.getPath(), ".fasta", 1, 3, true, false);
        assertTrue(x.exists());
        try {
          compareToFile(">test" + StringUtils.LS + "A13" + StringUtils.LS + "1" + StringUtils.LS + ">bob" + StringUtils.LS + "T32" + StringUtils.LS + "131" + StringUtils.LS + "00" + StringUtils.LS + ">cat" + StringUtils.LS + "C13" + StringUtils.LS + ">dog" + StringUtils.LS + "C00" + StringUtils.LS, x);
        } finally {
          assertTrue(x.delete());
        }
      } finally {
        sr.close();
      }
    } finally {
      assertTrue(FileHelper.deleteAll(xd));
      assertTrue(FileHelper.deleteAll(fastaDirectory));
    }
  }

  public void testFlags() {
    final Appendable err = new StringWriter();
    final Appendable out = new StringWriter();
    final CFlags flags = new CFlags("", out, err);
    new Sdf2Fasta().initFlags(flags);
    TestCFlags.check(flags,
            "SDF containing sequences",
            "output filename (extension added if not present)",
            "maximum number of nucleotides or amino acids");
  }

  private void checkContent(final String name, final String content) throws Exception {
    final File f = new File(name);
    assertTrue(f.exists());
    try (BufferedReader r = new BufferedReader(new FileReader(f))) {
      assertEquals(">x", r.readLine());
      assertEquals(content, r.readLine());
      assertNull(r.readLine());
    }
  }

  private void checkString(final String s) {
    if (s.length() != 0) {
      assertEquals("The current environment (operating system, JVM or machine) has not been tested. There is a risk of performance degradation or failure." + StringUtils.LS, s);
    }
  }

  private void runCommandWithNamedOutput(final String name, final String pathpr, final String content) throws Exception {
    final File tempDir = FileUtils.createTempDir("PrereadToFasta", null);
    try {
      final ByteArrayOutputStream bos = new ByteArrayOutputStream();
      final File output = new File(tempDir, name);
      try {
        try (PrintStream err = new PrintStream(bos)) {
          final ByteArrayOutputStream out = new ByteArrayOutputStream();
          assertEquals(0, new Sdf2Fasta().mainInit(new String[]{"-o", output.getPath(), "-i", pathpr, "-Z"}, out, err));
          assertEquals(0, out.toString().length());
          if (name.toLowerCase(Locale.getDefault()).endsWith(".fa") || name.toLowerCase(Locale.getDefault()).endsWith(".fasta")) {
            checkContent(new File(tempDir, name).getPath(), content);
          } else {
            checkContent(new File(tempDir, JUNITOUT + ".fasta").getPath(), content);
          }
        }
      } finally {
        bos.close();
      }
      checkString(bos.toString());
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  private void runCommandWithNamedOutput(final String name, final String pathpr, final String contentLeft, String contentRight) throws Exception {
    final File tempDir = FileUtils.createTempDir("PrereadToFasta", null);
    try {
      final ByteArrayOutputStream bos = new ByteArrayOutputStream();
      final File output = new File(tempDir, name);
      try {
        try (PrintStream err = new PrintStream(bos)) {
          final ByteArrayOutputStream out = new ByteArrayOutputStream();
          assertEquals(0, new Sdf2Fasta().mainInit(new String[]{"-o", output.getPath(), "-i", pathpr, "-Z"}, out, err));
          assertEquals(0, out.toString().length());
          if (name.toLowerCase(Locale.getDefault()).endsWith(".fa") || name.toLowerCase(Locale.getDefault()).endsWith(".fasta")) {
            final String ext = name.substring(name.lastIndexOf('.'));
            checkContent(new File(tempDir, JUNITOUT + "_1" + ext).getPath(), contentLeft);
            checkContent(new File(tempDir, JUNITOUT + "_2" + ext).getPath(), contentRight);
          } else {
            checkContent(new File(tempDir, JUNITOUT + "_1.fasta").getPath(), contentLeft);
            checkContent(new File(tempDir, JUNITOUT + "_2.fasta").getPath(), contentRight);
          }
        }
      } finally {
        bos.close();
      }
      checkString(bos.toString());
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  private static final String COL_FLAG = "--" + Sdf2Fasta.COLORSPACE_FLAG;

  private void runCommandColorspace2(final String pathpr) throws Exception {
    final File tempDir = FileUtils.createTempDir("PrereadToFasta", null);
    try {
      final ByteArrayOutputStream bos = new ByteArrayOutputStream();
      try {
        try (PrintStream err = new PrintStream(bos)) {
          final ByteArrayOutputStream out = new ByteArrayOutputStream();
          assertEquals(0, new Sdf2Fasta().mainInit(new String[]{"-o", new File(tempDir, JUNITOUT).getPath(), "-i", pathpr, COL_FLAG, "-Z"}, out, err));
          assertEquals(0, out.toString().length());
          final File f = new File(tempDir, JUNITOUT + ".fasta");
          assertTrue(f.exists());
          compareToFile(">x" + StringUtils.LS + "A121" + StringUtils.LS, f);
        }
      } finally {
        bos.close();
      }
      checkString(bos.toString());
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  private void runCommandColorspace(final String pathpr) throws Exception {
    final File tempDir = FileUtils.createTempDir("PrereadToFasta", null);
    try {
      final ByteArrayOutputStream bos = new ByteArrayOutputStream();
      try {
        try (PrintStream err = new PrintStream(bos)) {
          final ByteArrayOutputStream out = new ByteArrayOutputStream();
          assertEquals(0, new Sdf2Fasta().mainInit(new String[]{"-o", new File(tempDir, JUNITOUT).getPath(), "-i", pathpr, COL_FLAG, "-Z"}, out, err));
          assertEquals(0, out.toString().length());
          final File f = new File(tempDir, JUNITOUT + ".fasta");
          assertTrue(f.exists());
          final BufferedReader r = new BufferedReader(new FileReader(f));
          try {
            assertNull(r.readLine());
          } finally {
            r.close();
          }
        }
      } finally {
        bos.close();
      }
      final String s = bos.toString();
      assertTrue(s.contains("Colorspace conversion failed for \"1\" sequences due to presence of unknown nucleotides"));
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  private void runCommandLineLength2(final String pathpr) throws Exception {
    final File tempDir = FileUtils.createTempDir("PrereadToFasta", null);
    try {
      final ByteArrayOutputStream bos = new ByteArrayOutputStream();
      try {
        try (PrintStream err = new PrintStream(bos)) {
          final ByteArrayOutputStream out = new ByteArrayOutputStream();
          assertEquals(0, new Sdf2Fasta().mainInit(new String[]{"-o", new File(tempDir, JUNITOUT).getPath(), "-i", pathpr, "-l", "2", "-Z"}, out, err));
          assertEquals(0, out.toString().length());
          final File f = new File(tempDir, JUNITOUT + ".fasta");
          assertTrue(f.exists());
          compareToFile(">x" + StringUtils.LS + "AC" + StringUtils.LS + "TG" + StringUtils.LS + "N" + StringUtils.LS, f);
        }
      } finally {
        bos.close();
      }
      checkString(bos.toString());
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  private void runCommandChunks2(final String pathpr) throws Exception {
    final File testDir = FileUtils.createTempDir("PrereadToFastaTest", null);
    try {
      final ByteArrayOutputStream bos = new ByteArrayOutputStream();
      try {
        try (PrintStream err = new PrintStream(bos)) {
          final ByteArrayOutputStream out = new ByteArrayOutputStream();
          assertEquals(0, new Sdf2Fasta().mainInit(new String[]{"-o", new File(testDir, JUNITOUT).getPath(), "-i", pathpr, "-c", "2", "-Z"}, out, err));
          assertEquals(0, out.toString().length());
          final File f = new File(testDir, JUNITOUT + "0.fasta");
          assertTrue(f.exists());
          assertFalse(new File(testDir, JUNITOUT + "1.fasta").exists());
          assertFalse(new File(testDir, JUNITOUT + ".fasta").exists());
          compareToFile(">x" + StringUtils.LS + "ACTGN" + StringUtils.LS, f);
        } finally {
          assertTrue(new File(testDir, JUNITOUT + "0.fasta").delete());
        }
      } finally {
        bos.close();
      }
      checkString(bos.toString());
    } finally {
      assertTrue(FileHelper.deleteAll(testDir));
    }
  }

  private void createPreread(final String s, final File dir) throws IOException {
    final ArrayList<InputStream> al = new ArrayList<>();
    al.add(new ByteArrayInputStream(s.getBytes()));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(al, new DNAFastaSymbolTable());
    new SequencesWriter(ds, dir, 100000, PrereadType.UNKNOWN, false).processSequences();
  }

  private void createPrereadProtein(final File dir) throws IOException {
    final ArrayList<InputStream> al = new ArrayList<>();
    al.add(new ByteArrayInputStream((">x" + StringUtils.LS + "X*ARNDCQEGHILKMFPSTWYV" + StringUtils.LS).getBytes()));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(al, new ProteinFastaSymbolTable());
    new SequencesWriter(ds, dir, 100000, PrereadType.UNKNOWN, false).processSequences();
  }

  public void testValidUse() throws Exception {
    final File dir = FileHelper.createTempDirectory();
    try {
      final File outDir = FileHelper.createTempDirectory();
      try {
        createPreread(">x" + StringUtils.LS + "actgn" + StringUtils.LS, dir);
        final String pathpr = dir.getPath();
        try {
          runCommandWithNamedOutput(JUNITOUT, pathpr, "ACTGN");
          runCommandWithNamedOutput(JUNITOUT + ".FA", pathpr, "ACTGN");
          runCommandWithNamedOutput(JUNITOUT + ".fasta", pathpr, "ACTGN");
          runCommandColorspace(pathpr);
          runCommandLineLength2(pathpr);
          runCommandChunks2(pathpr);
        } finally {
          FileHelper.deleteAll(dir);
        }
        createPrereadProtein(dir);
        try {
          runCommandWithNamedOutput(JUNITOUT, pathpr, "X*ARNDCQEGHILKMFPSTWYV");
        } finally {
          FileHelper.deleteAll(dir);
        }
        createPreread(">x" + StringUtils.LS + "actg" + StringUtils.LS, dir);
        try {
          runCommandColorspace2(pathpr);
        } finally {
          FileHelper.deleteAll(dir);
          FileHelper.deleteAll(outDir);
        }
        assertFalse(Sdf2Fasta.toColorspace(new byte[] {0}, 1));
      } finally {
      assertTrue(!outDir.exists() || FileHelper.deleteAll(outDir));
      }
    } finally {
      assertTrue(!dir.exists() || FileHelper.deleteAll(dir));
    }
  }

  public void testValidUse2() throws Exception {
    final File dir = FileHelper.createTempDirectory();
    createPreread(">x" + StringUtils.LS + "actgn" + StringUtils.LS, new File(dir, "left"));
    createPreread(">x" + StringUtils.LS + "actgn" + StringUtils.LS, new File(dir, "right"));
    final String pathpr = dir.getPath();
    try {
      runCommandWithNamedOutput(JUNITOUT, pathpr, "ACTGN", "ACTGN");
      runCommandWithNamedOutput(JUNITOUT + ".FA", pathpr, "ACTGN", "ACTGN");
      runCommandWithNamedOutput(JUNITOUT + ".fasta", pathpr, "ACTGN", "ACTGN");
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  /**
   * @throws Exception
   */
  public void testGetStream() throws Exception {
    final File tempDir = FileUtils.createTempDir("PrereadToFasta", null);
    try {
      try (BufferedWriter ps = Sdf2Fasta.getStream(new File(tempDir, JUNITOUT).getPath(), false, false)) {
        final char[] buf = (char[]) TestUtils.getField("cb", ps);
        Assert.assertEquals(FileUtils.BUFFERED_STREAM_SIZE, buf.length);
      }
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  public void testChunkFlag() throws IOException {
    //Diagnostic.setLogStream();
    final MockEventListener ev = new MockEventListener();
    Diagnostic.addListener(ev);
    final Appendable sb = new StringWriter();
    final Appendable out = new StringWriter();
    final CFlags flags = new CFlags("", sb, out);
    new Sdf2Fasta().initFlags(flags);
    final File que = FileUtils.createTempDir("chunk", "p2f");
    final String[] args = {
      "-o", "temp",
      "-i", que.getPath(),
      "-c", "-40"
    };
    assertFalse(flags.setFlags(args));
    assertTrue(ev.compareErrorMessage("Error: Expected a positive integer for parameter \"Xchunks\"."));
    Diagnostic.removeListener(ev);
    assertTrue(FileHelper.deleteAll(que));
  }

  public void testLineFlag() throws IOException {
    Diagnostic.setLogStream();
    final MockEventListener ev = new MockEventListener();
    Diagnostic.addListener(ev);
    final StringWriter sb = new StringWriter();
    final StringWriter out = new StringWriter();
    final CFlags flags = new CFlags("", sb, out);
    new Sdf2Fasta().initFlags(flags);
    final File que = FileUtils.createTempDir("chunk", "p2f");
    final String[] args = {
      "-o", "testFile",
      "-i", que.getPath(),
      "-l", "-5"
    };
    assertFalse(flags.setFlags(args));
    assertTrue(ev.compareErrorMessage("Error: Expected a nonnegative integer for parameter \"line-length\"."));
    Diagnostic.removeListener(ev);
    assertTrue(FileHelper.deleteAll(que));
  }

  public void testInputAsFile() throws IOException {
    Diagnostic.setLogStream();
    final MockEventListener ev = new MockEventListener();
    Diagnostic.addListener(ev);
    final StringWriter sb = new StringWriter();
    final StringWriter out = new StringWriter();
    final CFlags flags = new CFlags("", sb, out);
    new Sdf2Fasta().initFlags(flags);
    final File que = File.createTempFile("p2f", "flag");
    //assertTrue(que.createNewFile());
    final String[] args = {
      "-o", "testFile",
      "-i", que.getPath()
    };
    assertFalse(flags.setFlags(args));
    assertTrue(ev.getEvent().getMessage(), ev.compareErrorMessage("Error: The specified file, \"" + que.getPath() + "\", is not an SDF."));
    Diagnostic.removeListener(ev);
    assertTrue(FileHelper.deleteAll(que));
  }


  private static final String FULL_NAME_DATA = ""
          + ">name suffix" + StringUtils.LS
          + "ACGTCG" + StringUtils.LS
          + ">second suffix" + StringUtils.LS
          + "ACGGGT" + StringUtils.LS;

  public void testFullName() throws IOException {
    final File dir = FileUtils.createTempDir("testsdf2fasta", "fullname");
    try {
      final File sdf = ReaderTestUtils.getDNADir(FULL_NAME_DATA, new File(dir, "sdf"));
      final File fasta = new File(dir, "fs.fasta.gz");
      final MemoryPrintStream out = new MemoryPrintStream();
      final MemoryPrintStream err = new MemoryPrintStream();
      final int code = new Sdf2Fasta().mainInit(new String[] {"-i", sdf.getPath(), "-o", fasta.getPath()}, out.outputStream(), err.printStream());
      assertEquals(err.toString(), 0, code);
      assertEquals(FULL_NAME_DATA, FileHelper.gzFileToString(fasta));
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }
}

