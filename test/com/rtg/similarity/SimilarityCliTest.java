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

import static com.rtg.util.StringUtils.LS;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;

import com.rtg.index.Finder;
import com.rtg.index.FinderHashValue;
import com.rtg.index.Index;
import com.rtg.index.hash.HashLoop;
import com.rtg.index.params.CountParams;
import com.rtg.index.params.CreateParams;
import com.rtg.index.similarity.IndexSimilarity;
import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.launcher.BuildParams;
import com.rtg.launcher.GlobalFlags;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.MockReaderParams;
import com.rtg.launcher.MockSequenceParams;
import com.rtg.launcher.ReaderParams;
import com.rtg.launcher.SequenceLengthMessageBuildSearchParams;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.DNAFastaSymbolTable;
import com.rtg.mode.ProgramMode;
import com.rtg.mode.SequenceMode;
import com.rtg.mode.SequenceType;
import com.rtg.reader.CompressedMemorySequencesReader;
import com.rtg.reader.FastaSequenceDataSource;
import com.rtg.reader.MockArraySequencesReader;
import com.rtg.reader.PrereadType;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SdfWriter;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesWriter;
import com.rtg.similarity.SimilarityCli.SimilarityTask;
import com.rtg.usage.UsageMetric;
import com.rtg.util.MockAppendable;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.DiagnosticEvent;
import com.rtg.util.diagnostic.DiagnosticListener;
import com.rtg.util.diagnostic.ErrorEvent;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.SlimException;
import com.rtg.util.integrity.Exam;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LogRecord;
import com.rtg.util.io.LogStream;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;
public class SimilarityCliTest extends AbstractCliTest {

  private static final String APP_NAME = "rtg similarity";

  private static final String TREEEXP2 = ""
      + "(" + LS
      + " (" + LS
      + "  baa:0.2500," + LS
      + "  bar:0.2500" + LS
      + " ):0.1250," + LS
      + " foo:0.1250" + LS
      + ")" + LS;

  private static final String SIMEXP2 = ""
      + "#rtg " + null + LS + "\tfoo\tbar\tbaa" + LS
      + "foo\t4\t4\t4" + LS
      + "bar\t4\t4\t4" + LS
      + "baa\t4\t4\t4" + LS;

  private static final String EXPTREE_1 = ""
      + "(" + LS
      + " baaa:0.5000," + LS
      + " bar:0.5000" + LS
      + ")" + LS;

  private static final String EXPSIMI_1 = ""
      + "#rtg " + null + LS + "\tbar\tbaaa" + LS
      + "bar\t1\t0" + LS
      + "baaa\t0\t1" + LS;

  /**
   * @param inputSequence the sequence
   * @param dir the dir
   * @throws IOException if badness
   */
  public static void writeDNA(final String inputSequence, final File dir) throws IOException {
    final ArrayList<InputStream> inputStreams = new ArrayList<>();
    inputStreams.add(new ByteArrayInputStream(inputSequence.getBytes()));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(inputStreams,
        new DNAFastaSymbolTable());
    final SequencesWriter sequenceWriter = new SequencesWriter(ds, dir, 20, PrereadType.UNKNOWN, true);
    sequenceWriter.processSequences();
  }

  /**
   * Make sure there are no logs in the working directory.
   */
  public static void checkNoLogs() {
    final File wd = new File(System.getProperty("user.dir"));
    final File[] logs = wd.listFiles(new FilenameFilter() {
      @Override
      public boolean accept(final File dir, final String name) {
        return name.startsWith("log") && name.endsWith(".log");
      }
    });
    if (logs != null && logs.length > 0) {
      System.err.println("Found " + logs.length + " logs:");
      for (final File log : logs) {
        System.err.println(log);
      }
      fail();
    }
  }

  /**
   * Take an array of strings which may contain nulls and return one with
   * the nulls elided.
   * Useful when generating command line arguments.
   * @param args0 to be processed.
   * @return args0 without the nulls.
   */
  public static String[] elideNulls(final String[] args0) {
    int count = 0;
    for (final String str : args0) {
      if (str != null) {
        count++;
      }
    }
    final String[] args = new String[count];
    int i = 0;
    for (final String str : args0) {
      if (str != null) {
        args[i] = str;
        i++;
      }
    }
    return args;
  }

  @Override
  public final void testApplicationName() {
    assertEquals(APP_NAME, new SimilarityCli().applicationName() + " " + new SimilarityCli().moduleName());
  }

  public final void testInitFlags() {
    checkHelp("Produces a similarity matrix and nearest neighbor tree from the input sequences or reads.",
        "-o,",
        "--output=DIR",
        "directory for output",
        "-i,",
        "--input=SDF",
        "SDF containing subject dataset",
        "-I,",
        "--input-list-file",
        "file containing a labeled list of SDF files (1 label and file per line format:[label][space][file])",
        "-s,",
        "--step=INT",
        "step size (Default is 1)",
        "-w,",
        "--word=INT",
        "word size (Default is 20)",
        "--unique-words",
        "count only unique words",
        "--max-reads=INT",
        "maximum number of reads to use from each input SDF");
  }


  private SequencesReader getReaderDNA(final String inputDnaSequence) throws IOException {
    final ArrayList<InputStream> inputStreams = new ArrayList<>();
    inputStreams.add(new ByteArrayInputStream(inputDnaSequence.getBytes()));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(inputStreams, new DNAFastaSymbolTable());

    return CompressedMemorySequencesReader.createSequencesReader(ds);
    /*
    final SequencesWriter sequenceWriter = new SequencesWriter(ds, dir, 20);
    sequenceWriter.processSequences();
    final DefaultSequencesReader dsr = DefaultSequencesReader.createSequencesReader(dir);
    return dsr;*/
  }

  public final void check(final String subjects, final int windowSize, final int stepSize, final String expSimi, final String expTree, long expUsage) throws Exception {
    checkTask(subjects, windowSize, stepSize, expSimi, expTree);
    checkExecErr(subjects, windowSize, stepSize, expSimi, expTree);
    checkMain(subjects, windowSize, stepSize, expSimi, expTree, expUsage);
  }

  private void runPhylogeny(final BuildSearchParams bsp, UsageMetric usageMetric) throws IOException {
    final SimilarityTask phy = new SimilarityTask(bsp, TestUtils.getNullOutputStream(), usageMetric);
    phy.run();
  }

  public final void checkTask(final String subjects, final int windowSize, final int stepSize, final String expSimi, final String expTree) throws Exception {
    final ProgramMode mode = ProgramMode.PHYLOGENY;
    final File outputDir = FileHelper.createTempDirectory();
    try {
      try (final SequencesReader sr = getReaderDNA(subjects)) {
        final ISequenceParams subjectParams = new MockSequenceParams(sr, mode.subjectMode());
        try (final BuildParams buildParams = BuildParams.builder().windowSize(windowSize).stepSize(stepSize).compressHashes(true).sequences(subjectParams).create()) {
          //try (final ISequenceParams queryParams = new MockSequenceParams(sr, mode.queryMode())) {
          final CountParams countParams = new CountParams(outputDir, 1, 1, false);
          final BuildSearchParams bsp = BuildSearchParams.builder()
              .mode(mode).build(buildParams)
              .count(countParams).create();
          final UsageMetric usageMetric = new UsageMetric();
          runPhylogeny(bsp, usageMetric);
          final String actualSimi = FileUtils.fileToString(new File(outputDir, SimilarityCli.SIMILARITY_SUFFIX));
          final String actualTree = FileUtils.fileToString(new File(outputDir, SimilarityCli.TREE_SUFFIX));
          assertEquals(expSimi, actualSimi);
          if (expTree != null) {
            assertEquals(expTree, actualTree);
          }
        }
      }
    } finally {
      assertTrue(FileHelper.deleteAll(outputDir));
    }
  }

  public final String checkExecErr(final String subjects, final int windowSize, final int stepSize, final String expSimi, final String expTree) throws Exception {
    final String err;
    final File outputDir = FileHelper.createTempDirectory();
    try {
      FileHelper.deleteAll(outputDir);
      err = checkExec(subjects, windowSize, stepSize, expSimi, expTree, outputDir);
    } finally {
      assertTrue(FileHelper.deleteAll(outputDir));
    }
    return err;
  }

  public final String checkExecLog(final String subjects, final int windowSize, final int stepSize, final String expSimi, final String expTree) throws Exception {
    final String log;
    final File outputDir = FileHelper.createTempDirectory();
    try {
      FileHelper.deleteAll(outputDir);
      checkExec(subjects, windowSize, stepSize, expSimi, expTree, outputDir);
      log = FileUtils.fileToString(new File(outputDir, SimilarityCli.MODULE_NAME + ".log"));
    } finally {
      assertTrue(FileHelper.deleteAll(outputDir));
    }
    return log;
  }

  private String checkExec(String subjects, int windowSize, int stepSize, String expSimi, String expTree, File outputDir) throws IOException {
    final File subjectsDir = FileHelper.createTempDirectory();
    final String errStr;
    try {
      writeDNA(subjects, subjectsDir);
      try (final ByteArrayOutputStream err = new ByteArrayOutputStream()) {
        try (final PrintStream errs = new PrintStream(err)) {
          try (final ByteArrayOutputStream out = new ByteArrayOutputStream()) {
            final String[] args0 = {
                "-i", subjectsDir.getPath(),
                "-o", outputDir.getPath(),
                "-w", windowSize + "",
                "-s", stepSize + ""
            };
            final String[] args = elideNulls(args0);
            GlobalFlags.resetAccessedStatus();
            final int ret = new SimilarityCli().mainInit(args, out, errs);
            assertEquals(0, ret);

            final String logStr = FileUtils.fileToString(new File(outputDir, SimilarityCli.MODULE_NAME + ".log"));
            //  System.err.println("log stream:");
            //System.err.println(logStr);
            checkLog(logStr);
          }
        }
        errStr = err.toString();
      }
      final String actualSimi = FileUtils.fileToString(new File(outputDir, SimilarityCli.SIMILARITY_SUFFIX));
      final String actualTree = FileUtils.fileToString(new File(outputDir, SimilarityCli.TREE_SUFFIX));
      assertTrue(new File(outputDir, SimilarityCli.XML_SUFFIX).isFile());
      //System.err.println("similarity" + LS + actualSimi);
      //System.err.println("tree" + LS + actualTree);
      assertEquals(expSimi, actualSimi);
      if (expTree != null) {
        assertEquals(expTree, actualTree);
      }
    } finally {
      assertTrue(FileHelper.deleteAll(subjectsDir));
    }
    return errStr;
  }

  private void checkLog(final String logStr) {
    final String versionString = "java.version = ";
    TestUtils.containsAll(logStr, versionString, "user.timezone = ", "Parameters:" + LS + "BuildSearchParams",
      "Usage of memory", " Timer ", "Shared_buffer",
      "mode=PHYLOGENY", "Ph_similarity_matrix ", "Ph_similarity_neighbor ", "Ph_similarity_out ",
      " Memory performance");
  }

  public final String checkMain(final String subjects, final int windowSize, final int stepSize, final String expSimi, final String expTree, final long expUsage) throws Exception {
    final String errStr;
    final File subjectsDir = FileHelper.createTempDirectory();
    try {
      writeDNA(subjects, subjectsDir);
      final File outputDir = FileHelper.createTempDirectory();
      try {
        FileHelper.deleteAll(outputDir);
        try (final ByteArrayOutputStream err = new ByteArrayOutputStream()) {
          try (final PrintStream errs = new PrintStream(err)) {
            try (final ByteArrayOutputStream out = new ByteArrayOutputStream()) {
              checkNoLogs();
              final String[] args0 = {
                  "-i", subjectsDir.getPath(),
                  "-o", outputDir.getPath(),
                  "-w", windowSize + "",
                  "-s", stepSize + ""
              };
              final String[] args = elideNulls(args0);
              final SimilarityCli cli = new SimilarityCli();
              cli.mainInit(args, out, errs);
              checkNoLogs();

              final String logStr = FileUtils.fileToString(new File(outputDir, SimilarityCli.MODULE_NAME + ".log"));
              checkLog(logStr);
              final String usageLog = cli.usageLog();
              //System.err.println(usageLog);
              TestUtils.containsAll(usageLog,
                  "[Usage beginning module=similarity runId=", "Usage end module=similarity runId=", " metric=" + expUsage + " success=true]");
            }
          }
          errStr = err.toString();
        }

        final String actualSimi = FileUtils.fileToString(new File(outputDir, SimilarityCli.SIMILARITY_SUFFIX));
        final String actualTree = FileUtils.fileToString(new File(outputDir, SimilarityCli.TREE_SUFFIX));

        //System.err.println("similarity" + LS + actualSimi);
        //System.err.println("tree" + LS + actualTree);
        //System.err.println(FileUtils.fileToString(new File(outputDir, "similarity.log")));
        assertEquals(expSimi, actualSimi);
        if (expTree != null) {
          assertEquals(expTree, actualTree);
        }
      } finally {
        assertTrue(!outputDir.exists() || FileHelper.deleteAll(outputDir));
      }
    } finally {
      assertTrue(FileHelper.deleteAll(subjectsDir));
    }
    return errStr;
  }

  /* Single sequence. */
  private static final String SEQ_0 = ""
      + ">baa" + LS + "aaaa";

  public void test0() throws Exception {
    final String simExp = ""
        + "#rtg " + null + LS + "\tbaa" + LS
        + "baa\t1" + LS
        ;
    final String treeExp = ""
        + "baa" + LS
        ;
    check(SEQ_0, 4, 4, simExp, treeExp, 4);
  }

  private static final String SEQ_1 = ""
      + ">bar" + LS  + "acgt" + LS
      + ">baaa" + LS + "aaaa";

  public void test1() throws Exception {
    final String simExp = EXPSIMI_1
        ;
    final String treeExp = EXPTREE_1
        ;
    check(SEQ_1, 4, 4, simExp, treeExp, 8);
  }

  private static final String SEQ_2 = ""
      + ">foo" + LS + "acgt" + LS
      + ">bar" + LS + "acgt" + LS
      + ">baa" + LS + "acgt";

  public void test2a() throws Exception {
    final String simExp = ""
        + "#rtg " + null + LS + "\tfoo\tbar\tbaa" + LS
        + "foo\t1\t1\t1" + LS
        + "bar\t1\t1\t1" + LS
        + "baa\t1\t1\t1" + LS
        ;
    final String treeExp = TREEEXP2
        ;
    check(SEQ_2, 4, 4, simExp, treeExp, 12);
  }

  /* Step size 2. */
  public void test2b() throws Exception {
    check(SEQ_2, 2, 2, SIMEXP2, TREEEXP2, 12);
  }

  /* Various progress options - explicitly do no progress. */
  public void test2bProgress1() throws Exception {
    checkMain(SEQ_2, 2, 2, SIMEXP2, TREEEXP2, 12);
  }

  /* Various progress options - only to error stream. */
  public void test2bProgress2() throws Exception {
    checkMain(SEQ_2, 2, 2, SIMEXP2, TREEEXP2, 12);
  }

  /* Various progress options - only to file. */
  public void test2bProgress3() throws Exception {
    final String errStr = checkMain(SEQ_2, 2, 2, SIMEXP2, TREEEXP2, 12);
    assertEquals("", errStr);
  }

  /* Various progress options - to file and error stream. */
  public void test2bProgress4() throws Exception {
    checkMain(SEQ_2, 2, 2, SIMEXP2, TREEEXP2, 12);
  }

  public void test1Log() throws Exception {
    try (final ByteArrayOutputStream ba = new ByteArrayOutputStream()) {
      try (final PrintStream pr = new PrintStream(ba)) {
        Diagnostic.setLogStream(pr);
        final String s = checkExecLog(SEQ_1, 4, 4, EXPSIMI_1, EXPTREE_1);
        //      System.err.println(s);
        checkLog(s);
        assertTrue(s.contains("Hash counts\t0\t1\t2"));
        assertTrue(s.contains("Bucket counts\t0\t1\t2"));
        assertTrue(s.contains("Usage of memory"));
        assertTrue(s.contains("Memory performance"));
      } finally {
        Diagnostic.setLogStream();
      }
    }
  }

  //testing longer word sizes
  private static final String LONG_SEQ = //128 nt
      "actgactgactgactgactgactgactgactgactgactgactgactgactgactgactgactgactgactgactgactgactgactgactgactgactgactgactgactgactgactgactgactg";

  private static final String SEQ_DNA_LONG = ""
      + ">x" + LS
      + LONG_SEQ + LS;

  public void testL16() throws Exception {
    assertEquals(128, LONG_SEQ.length());
    final String simExp = ""
        + "#rtg " + null + LS + "\tx" + LS
        + "x\t64" + LS
        ;
    check(SEQ_DNA_LONG, 16, 16, simExp, null, 128);
  }

  public void testL17() throws Exception {
    assertEquals(128, LONG_SEQ.length());
    final String simExp = ""
        + "#rtg " + null + LS + "\tx" + LS
        + "x\t13" + LS
        ;
    check(SEQ_DNA_LONG, 17, 17, simExp, null, 128);
  }

  public void testL32() throws Exception {
    assertEquals(128, LONG_SEQ.length());
    final String simExp = ""
        + "#rtg " + null + LS + "\tx" + LS
        + "x\t16" + LS
        ;
    check(SEQ_DNA_LONG, 32, 32, simExp, null, 128);
  }

  private static class MyListener implements DiagnosticListener {
    @Override
    public void handleDiagnosticEvent(final DiagnosticEvent<?> event) {
      Exam.assertTrue(event instanceof ErrorEvent);
      final ErrorEvent e = (ErrorEvent) event;
      Exam.assertEquals(ErrorType.SEQUENCE_TOO_LONG, e.getType());
      Exam.assertEquals("2147483648", e.getParams()[0]);
    }
    @Override
    public void close() {
    }
  }

  public void testSequenceLengthMessage() throws Exception {
    final ProgramMode mode = ProgramMode.SLIMN;
    final File hitsDir = FileHelper.createTempDirectory();
    try {
      try (final SequencesReader sr = getReaderDNA(SEQ_1)) {
        final ISequenceParams subjectParams = new MockSequenceParams(sr, mode.subjectMode());
        try (final BuildParams buildParams = BuildParams.builder().windowSize(4).stepSize(1).sequences(subjectParams).create()) {
          final CountParams countParams = new CountParams(hitsDir, 20, 1, false);
          final BuildSearchParams bsp =
              new SequenceLengthMessageBuildSearchParams(BuildSearchParams.builder().mode(mode)
                  .build(buildParams).count(countParams));
          final MyListener listener = new MyListener();
          try {
            Diagnostic.addListener(listener);
            SimilarityCli.buildAndSearch(bsp);
            fail();
          } catch (final SlimException e) {
            // ok
          } finally {
            Diagnostic.removeListener(listener);
          }
        }
      }
    } finally {
      assertTrue(FileHelper.deleteAll(hitsDir));
    }
  }

  public void testValidUse() throws Exception {
    final File tempDir = FileHelper.createTempDirectory();
    try {
      final File dir = new File(tempDir, "dir");
      final File outDir = new File(tempDir, "outDir");
      final ArrayList<InputStream> al = new ArrayList<>();
      al.add(new ByteArrayInputStream((">x" + LS + "actg" + LS).getBytes()));
      final FastaSequenceDataSource ds = new FastaSequenceDataSource(al, new DNAFastaSymbolTable());
      final SequencesWriter sw = new SequencesWriter(ds, dir, 100000, PrereadType.UNKNOWN, false);
      sw.processSequences();
      final String pathpr = dir.getPath();
      try (final ByteArrayOutputStream bos = new ByteArrayOutputStream()) {
        try (final ByteArrayOutputStream bout = new ByteArrayOutputStream()) {
          try (final PrintStream err = new PrintStream(bos)) {
            assertEquals(0, new SimilarityCli().mainInit(new String[] {"-o", outDir.getPath(), "-i", pathpr}, bout, err));
            assertTrue(outDir.exists());
            assertTrue(outDir.toString(), new File(outDir, SimilarityCli.TREE_SUFFIX).exists());
            assertTrue(outDir.toString(), new File(outDir, SimilarityCli.XML_SUFFIX).exists());
            assertTrue(outDir.toString(), new File(outDir, SimilarityCli.SIMILARITY_SUFFIX).exists());
            assertTrue(outDir.toString(), FileUtils.fileToString(new File(outDir, SimilarityCli.TREE_SUFFIX)).length() > 0);
            assertTrue(outDir.toString(), FileUtils.fileToString(new File(outDir, SimilarityCli.XML_SUFFIX)).length() > 0);
            assertTrue(outDir.toString(), FileUtils.fileToString(new File(outDir, SimilarityCli.SIMILARITY_SUFFIX)).length() > 0);
          }
        }
        final String s = bos.toString();
        if (s.length() != 0) {
          assertEquals("The current environment (operating system, JVM or machine) has not been tested. There is a risk of performance degradation or failure." + LS, bos.toString());
        }
      }
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  private static final String EXPECTED_ERROR = ""
      + "Error: Unknown flag --no-such-option" + LS
      + LS
      + "Usage: rtg similarity [OPTION]... -o DIR -i SDF" + LS
      + "                      [OPTION]... -o DIR -I FILE" + LS + LS
      + "Try '--help' for more information" + LS
      ;

  public void testBadUse() throws Exception {
    try (final ByteArrayOutputStream bos = new ByteArrayOutputStream()) {
      try (final ByteArrayOutputStream bout = new ByteArrayOutputStream()) {
        try (final PrintStream err = new PrintStream(bos)) {
          assertEquals(1, new SimilarityCli().mainInit(new String[] {"--no-such-option"}, bout, err));
        }
      }
      final String e = bos.toString();
      assertEquals(e, EXPECTED_ERROR, e);
    }
  }

  /* Single sequence. */
  private static final String NT_DEFU = "actgactgactgactgactgactga";
  private static final String SEQ_DEFU = ""
      + ">baa" + LS + NT_DEFU;

  public void testDefaultUsage() throws Exception {
    final String simExp = ""
        + "#rtg " + null + LS + "\tbaa" + LS
        + "baa\t10" + LS
        ;
    final String treeExp = ""
        + "baa" + LS
        ;
    check(SEQ_DEFU, SimilarityCli.DEFAULT_WORD_SIZE, SimilarityCli.DEFAULT_STEP_SIZE, simExp, treeExp, NT_DEFU.length());

    checkDefTask(simExp, treeExp); //TODO
    checkDefMain(simExp, treeExp);
  }

  private void checkDefTask(final String simExp, final String treeExp) throws IOException {
    Diagnostic.setLogStream();
    final ProgramMode mode = ProgramMode.PHYLOGENY;
    final File outputDir = FileHelper.createTempDirectory();
    try {
      try (final SequencesReader sr = getReaderDNA(SEQ_DEFU)) {
        final ISequenceParams subjectParams = new MockSequenceParams(sr, mode.subjectMode());
        try (final BuildParams buildParams = BuildParams.builder()
            .windowSize(SimilarityCli.DEFAULT_WORD_SIZE)
            .stepSize(SimilarityCli.DEFAULT_STEP_SIZE)
            .compressHashes(true)
            .sequences(subjectParams)
            .create()) {
          final CountParams countParams = new CountParams(outputDir, 1, 1, false);
          final BuildSearchParams bsp = BuildSearchParams.builder()
              .mode(mode).build(buildParams)
              .count(countParams).create();
          final UsageMetric usageMetric = new UsageMetric();
          runPhylogeny(bsp, usageMetric);
          final String actualSimi = FileUtils.fileToString(new File(outputDir, SimilarityCli.SIMILARITY_SUFFIX));
          final String actualTree = FileUtils.fileToString(new File(outputDir, SimilarityCli.TREE_SUFFIX));
          assertEquals(simExp, actualSimi);
          if (treeExp != null) {
            assertEquals(treeExp, actualTree);
          }
        }
      }
    } finally {
      Diagnostic.setLogStream();
      assertTrue(FileHelper.deleteAll(outputDir));
    }
  }

  private void checkDefMain(final String simExp, final String treeExp)
      throws IOException {
    final File subjectsDir = FileHelper.createTempDirectory();
    try {
      writeDNA(SEQ_DEFU, subjectsDir);
      final File outputDir = FileHelper.createTempDirectory();
      try {
        FileHelper.deleteAll(outputDir);
        final ByteArrayOutputStream ba = new ByteArrayOutputStream();
        final PrintStream err = new PrintStream(ba);
        final ByteArrayOutputStream bout = new ByteArrayOutputStream();
        checkNoLogs();
        GlobalFlags.resetAccessedStatus();
        new SimilarityCli().mainInit(new String[]
            {"-i", subjectsDir.getPath(), "-o", outputDir.getPath()},
            bout, err
            );
        checkNoLogs();

        final String logStr = FileUtils.fileToString(new File(outputDir, SimilarityCli.MODULE_NAME + ".log"));
        checkLog(logStr);

        final String actualSimi = FileUtils.fileToString(new File(outputDir, SimilarityCli.SIMILARITY_SUFFIX));
        final String actualTree = FileUtils.fileToString(new File(outputDir, SimilarityCli.TREE_SUFFIX));
        //System.err.println("similarity" + LS + actualSimi);
        //System.err.println("tree" + LS + actualTree);
        assertEquals(simExp, actualSimi);
        if (treeExp != null) {
          assertEquals(treeExp, actualTree);
        }
      } finally {
        assertTrue(!outputDir.exists() || FileHelper.deleteAll(outputDir));
      }
    } finally {
      assertTrue(FileHelper.deleteAll(subjectsDir));
    }

  }

  public void testIOException() throws IOException {
    final CreateParams cp = new CreateParams(0, 1, 1, false, true, false);
    final IndexSimilarity index = new IndexSimilarity(cp, 1000, false, 1000, 0, false, 1);
    index.freeze();

    final SequencesReader sr = new MockArraySequencesReader(SequenceType.DNA, new int[] {1}, new String[] {"name"});
    final ReaderParams rp = new MockReaderParams(sr, SequenceMode.BIDIRECTIONAL);
    final ISequenceParams sp = new MockSequenceParams(rp);
    final BuildSearchParams bsp = BuildSearchParams.builder().build(BuildParams.builder().stepSize(1).windowSize(1).sequences(sp).create()).create();
    final Appendable ok = new StringWriter();
    final Appendable ok2 = new StringWriter();
    final Appendable bad = new MockAppendable();

    final LogStream log1 = new LogRecord();
    Diagnostic.setLogStream(log1);
    try {
      SimilarityCli.similarity(index, 1, null, bsp, "outDir", ok, ok2, bad, null);
      fail();
    } catch (final SlimException e) {
      assertTrue(e.getMessage().contains(java.io.IOException.class.getName()));
      e.logException();
    }
    final String logStr1 = log1.toString();
    //System.err.println(logStr);
    assertTrue(logStr1.contains(java.io.IOException.class.getName()));
    assertTrue(logStr1.contains("Error: A problem has occurred while attempting to write information into the directory \"outDir\". Check that this directory exists and that you have permission to write it."));


    final LogStream log2 = new LogRecord();
    Diagnostic.setLogStream(log2);
    try {
      SimilarityCli.similarity(index, 1, null, bsp, "dirOut", bad, ok2, ok, null);
      fail();
    } catch (final SlimException e) {
      assertTrue(e.getMessage().contains(java.io.IOException.class.getName()));
    }
    final String logStr2 = log1.toString();
    //System.err.println(logStr);
    assertTrue(logStr2.contains(java.io.IOException.class.getName()));
    assertTrue(logStr2.contains("Error: A problem has occurred while attempting to write information into the directory \"outDir\". Check that this directory exists and that you have permission to write it."));
    Diagnostic.setLogStream();
  }

  public void testMakeNames() throws IOException {
    final SequencesReader reader = new MockArraySequencesReader(SequenceType.DNA, 4);
    final ArrayList<String> na = SimilarityCli.makeNames(reader);
    assertEquals("[seq0, seq1, seq2, seq3]", Arrays.toString(na.toArray()));
  }

  public void testMakeBuild() throws IOException {
    final SequenceParams dummySubjectParams = SequenceParams.builder().region(new HashingRegion(0, 0)).mode(SequenceMode.BIDIRECTIONAL).create();
    BuildParams buildParams = BuildParams.builder().windowSize(33).stepSize(1).size(0).sequences(dummySubjectParams).create();
    final DummyIndex index = new DummyIndex() {
      @Override
      public void add(long hash, long id) {
        if (hash == -1) {
          assertEquals(67, id);
        } else if (hash == -2) {
          assertEquals(39, id);
        } else if (hash == -3) {
          assertEquals(15, id);
        } else {
          fail("hash: " + hash + " id: " + id);
        }
      }
    };
    try {
      SimilarityCli.makeBuild(index, buildParams, -1, null);
      fail();
    } catch (final SlimException e) {
      assertEquals("Word size > 32", e.getMessage());
    }
    buildParams = BuildParams.builder().windowSize(31).stepSize(1).size(0).sequences(dummySubjectParams).create();
    HashLoop hl = SimilarityCli.makeBuild(index, buildParams, -1, null);
    assertNotNull(hl);
    hl.hashCallBidirectional(1, -1, 1, 67);
    hl.hashCallBidirectional(-2, 2, 1, 39);

    hl = SimilarityCli.makeBuild(index, buildParams, 15, null);
    hl.hashCallBidirectional(3, -3, 1, 76);
  }

  public void testMakeParams() throws IOException {
    final File tempDir = FileUtils.createTempDir("makeParams", "test");
    try {
      final File input = new File(tempDir, "inputReads");
      ReaderTestUtils.getReaderDNA(REPEAT_TEMPLATE, input, null).close();
      final File output = new File(tempDir, "output");
      String error = checkMainInitBadFlags(
          "-i", input.getPath(),
          "-o", output.getPath(),
          "-s", "26"
          );
      TestUtils.containsAll(error, "The specified flag \"--step\" has invalid value \"26\". It should be less than or equal to \"--word\".");
      final File inputList = new File(tempDir, "inputList.txt");
      assertTrue(inputList.createNewFile());
      error = checkMainInitBadFlags(
          "-I", inputList.getPath(),
          "-o", output.getPath()
          );
      TestUtils.containsAll(error, "The input list file contained no target files.");
      assertTrue(inputList.delete());
      StringBuilder inputListBuilder = new StringBuilder();
      for (int i = 0; i < 11; i++) {
        inputListBuilder.append(i).append(" ").append(inputList.getPath()).append(LS);
      }
      FileUtils.stringToFile(inputListBuilder.toString(), inputList);
      error = checkMainInitBadFlags(
          "-I", inputList.getPath(),
          "-o", output.getPath()
          );
      TestUtils.containsAll(error, "There were 11 incorrect input list lines.", "The specified file, \"" + inputList.getPath() + "\", is not an SDF.");
      assertEquals(11, error.split("The specified file, \".+?\", is not an SDF.").length);

      assertTrue(inputList.delete());
      final File nonExistant = new File(tempDir, "nonExistant");
      inputListBuilder = new StringBuilder();
      for (int i = 0; i < 12; i++) {
        inputListBuilder.append(i).append(" ").append(nonExistant.getPath()).append(LS);
      }
      FileUtils.stringToFile(inputListBuilder.toString(), inputList);
      error = checkMainInitBadFlags(
          "-I", inputList.getPath(),
          "-o", output.getPath()
          );
      TestUtils.containsAll(error, "There were 12 incorrect input list lines.", "The specified SDF, \"" + nonExistant.getPath() + "\", does not exist.");
      assertEquals(11, error.split("The specified SDF, \".+?\", does not exist.").length);
    } finally {
      FileHelper.deleteAll(tempDir);
    }
  }

  private static final String REPEAT_TEMPLATE = ">A"  + LS
      + "acgtacgtacgtacgtacgtacgtacgt" + LS
      + ">B" + LS
      + "acgtacgtacgtacgtacgtacgtacgt" + LS;

  private static final String REPEAT_EXPECTED_1 = ""
      + "#rtg " + null + LS
      + "\tA\tB" + LS
      + "A\t49\t49" + LS
      + "B\t49\t49" + LS;

  private static final String REPEAT_EXPECTED_2 = ""
      + "#rtg " + null + LS
      + "\tA\tB" + LS
      + "A\t1\t1" + LS
      + "B\t1\t1" + LS;

  public void testEnd2End() throws Exception {
    final File f = FileUtils.createTempDir("phylogeny", "end2end");
    try {
      final File input = new File(f, "Reads");
      ReaderTestUtils.getReaderDNA(REPEAT_TEMPLATE, input, null).close();

      final File output1 = new File(f, "output1");
      final File output2 = new File(f, "output2");

      final String[] args = {
          "-i", input.getPath(),
          "-w", "4",
          "-s", "4",
          "-o", output1.getPath()
      };

      final SimilarityCli p = new SimilarityCli();
      assertEquals(0, p.mainInit(args, TestUtils.getNullOutputStream(), System.err));
      final String sim = FileUtils.fileToString(new File(output1, "similarity.tsv"));
      assertEquals(REPEAT_EXPECTED_1, sim);
      final String progress = FileUtils.fileToString(new File(output1, "progress"));
      TestUtils.containsAll(progress, "Input for index starting"
          , "Input for post-freeze index starting"
          , "Input for index finished. Starting indexing"
          , "Indexing finished. Starting similarity."
          , "Similarity finished."
          , "Finished successfully in "
          );

      final String[] args2 = {
          "-i", input.getPath(),
          "-w", "4",
          "-s", "4",
          "-o", output2.getPath(),
          "--unique-words"
      };

      assertEquals(0, p.mainInit(args2, TestUtils.getNullOutputStream(), System.err));
      final String sim2 = FileUtils.fileToString(new File(output2, "similarity.tsv"));
      assertEquals(REPEAT_EXPECTED_2, sim2);

    } finally {
      assertTrue(FileHelper.deleteAll(f));
    }
  }

  public void testEnd2EndMultifile() throws Exception {
    final File tempDir = FileUtils.createTempDir("multiFile", "test");
    try {
      final SimilarityCli p = new SimilarityCli();

      final File emptyInput = new File(tempDir, "empty");
      assertTrue(emptyInput.mkdir());
      final File emptyInputLeft = new File(emptyInput, "left");
      final File emptyInputRight = new File(emptyInput, "right");
      new SdfWriter(emptyInputLeft, 2147483648L, PrereadType.UNKNOWN, false, true, false, SequenceType.DNA).close();
      new SdfWriter(emptyInputRight, 2147483648L, PrereadType.UNKNOWN, false, true, false, SequenceType.DNA).close();

      final File pairedInput = new File(tempDir, "paired");
      assertTrue(pairedInput.mkdir());
      final File inputLeft = new File(pairedInput, "left");
      final File inputRight = new File(pairedInput, "right");
      ReaderTestUtils.getReaderDNA(REPEAT_TEMPLATE, inputLeft, null).close();
      ReaderTestUtils.getReaderDNA(REPEAT_TEMPLATE, inputRight, null).close();
      final File inputList = new File(tempDir, "input.txt");
      FileUtils.stringToFile("A " + inputLeft.getPath() + LS
          + "B " + inputRight.getPath() + LS
          + "E " + emptyInputLeft.getAbsolutePath(), inputList);
      final File output = new File(tempDir, "output");
      final String[] args = {
          "-I", inputList.getPath(),
          "-w", "4",
          "-s", "4",
          "-o", output.getPath(),
      };

      MemoryPrintStream ps = new MemoryPrintStream();
      assertEquals(0, p.mainInit(args, ps.outputStream(), ps.printStream()));
      assertEquals("The SDF \"" + emptyInputLeft.getPath() + "\" contains no sequences" + LS, ps.toString());
      String result = FileUtils.fileToString(new File(output, "similarity.tsv"));
      assertEquals("#rtg " + null + LS
          + "\tA\tB\tE" + LS
          + "A\t196\t196\t0" + LS
          + "B\t196\t196\t0" + LS
          + "E\t0\t0\t0" + LS, result);
      final String progress = FileUtils.fileToString(new File(output, "progress"));
      TestUtils.containsAll(progress, "Input for index starting"
          , "Input for index finished. Starting indexing"
          , "Indexing finished. Starting similarity."
          , "Similarity finished."
          , "Input for index starting label \"A\" (1/3)"
          , "Input for index label \"A\" starting directory (1/1)"
          , "Input for index starting label \"B\" (2/3)"
          , "Input for index label \"B\" starting directory (1/1)"
          , "Input for index starting label \"E\" (3/3)"
          , "Input for post-freeze index starting label \"A\" (1/3)"
          , "Input for post-freeze index label \"A\" starting directory (1/1)"
          , "Input for post-freeze index starting label \"B\" (2/3)"
          , "Input for post-freeze index label \"B\" starting directory (1/1)"
          , "Input for post-freeze index starting label \"E\" (3/3)"
          , "Finished successfully in "
          );

      assertTrue(inputList.delete());
      assertTrue(FileHelper.deleteAll(output));
      ps = new MemoryPrintStream();
      FileUtils.stringToFile("A " + pairedInput.getPath() + LS
          + " BLARGH" + LS
          + "C " + LS
          + LS
          + "E " + emptyInput.getAbsolutePath(), inputList);
      assertEquals(0, p.mainInit(args, ps.outputStream(), ps.printStream()));
      assertEquals("The SDF \"" + emptyInput.getPath() + "\" contains no sequences" + LS, ps.toString());
      result = FileUtils.fileToString(new File(output, "similarity.tsv"));
      assertEquals("#rtg " + null + LS
          + "\tA\tE" + LS
          + "A\t784\t0" + LS
          + "E\t0\t0" + LS
          , result);
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
      Diagnostic.setLogStream();
    }
  }

  public void testGetAdjustedRegion() throws IOException {
    final File tempDir = FileUtils.createTempDir("adjustedRegion", "test");
    try {
      final MemoryPrintStream ps = new MemoryPrintStream();
      Diagnostic.setLogStream(ps.printStream());
      try {
        SimilarityCli.resolveRange(tempDir, tempDir, null);
        fail();
      } catch (final IOException e) {
        // Expected
      }
      final File emptyInput = new File(tempDir, "empty");
      new SdfWriter(emptyInput, 2147483648L, PrereadType.UNKNOWN, false, true, false, SequenceType.DNA).close();
      ps.reset();
      assertNull(SimilarityCli.resolveRange(emptyInput, emptyInput, null));
      assertTrue(ps.toString(), ps.toString().contains("The SDF \"" + emptyInput.getPath() + "\" contains no sequences"));
      final File validInput = new File(tempDir, "valid");
      ReaderTestUtils.getReaderDNA(REPEAT_TEMPLATE, validInput, null).close();
      ps.reset();
      final LongRange inRegion = new LongRange(-1, 10);
      final LongRange ret = SimilarityCli.resolveRange(validInput, validInput, inRegion);
      assertNotNull(ret);
      assertEquals(2, ret.getEnd());
      assertTrue(ps.toString(), ps.toString().contains("The end sequence id"));
      ps.reset();
      assertEquals(ret.toString(), SimilarityCli.resolveRange(validInput, validInput, ret).toString());
      assertTrue(ps.toString().isEmpty());
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
      Diagnostic.setLogStream();
    }
  }

  private static class DummyIndex implements Index {

    @Override
    public void add(long hash, long id) { }
    @Override
    public long bytes() {
      return 0;
    }
    @Override
    public void dumpValues(PrintStream out) {
    }
    @Override
    public void freeze() {
    }
    @Override
    public long getHash(long found) {
      return 0;
    }
    @Override
    public long getValue(long found) {
      return 0;
    }
    @Override
    public String infoString() {
      return null;
    }
    @Override
    public int maxHashCount() {
      return 0;
    }
    @Override
    public long numberEntries() {
      return 0;
    }
    @Override
    public long numberHashes() {
      return 0;
    }
    @Override
    public String perfString() {
      return null;
    }
    @Override
    public void search(long hash, Finder finder) throws IOException, IllegalStateException {
    }
    @Override
    public void scan(FinderHashValue finder) throws IOException, IllegalStateException {
    }
    @Override
    public long search(long hash) {
      return 0;
    }
    @Override
    public int searchCount(long hash) {
      return 0;
    }
  }

  @Override
  protected AbstractCli getCli() {
    return new SimilarityCli();
  }
}

