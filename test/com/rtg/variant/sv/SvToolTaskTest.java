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
package com.rtg.variant.sv;

import static com.rtg.sam.SharedSamConstants.OUT_SAM;
import static com.rtg.sam.SharedSamConstants.REF_SEQS;
import static com.rtg.sam.SharedSamConstants.SAM_UNSORTED;
import static com.rtg.util.StringUtils.FS;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;

import com.rtg.Slim;
import com.rtg.launcher.GlobalFlags;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.Resources;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.SimpleArchive;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class SvToolTaskTest extends TestCase {

  public void testUnsortedReadStats() throws IOException {
    check(REF_SEQS, SAM_UNSORTED, "RG1\t200000\t56893\t1991255\t69693925\t55590\t22261070\t9022416526\t55590\t22261070\t9022416526\t5\t55590\t236\t1067", new String[] {}, null, new String[] {"is not sorted in coordinate order."}, "", 1, true, null);
  }

  protected static void check(String refSeq, String sam, String readStats, String[] args0, String expFile,
      String[] errorMsgs, String logMsg, int errCode, boolean gzip, String[] erlines) throws IOException {
    final File output = FileUtils.createTempDir("testcheck", "sv_out");
    FileHelper.deleteAll(output);
    final File input = FileUtils.createTempDir("testcheck", "sv_in");
    try {
      FileUtils.stringToFile(sam, new File(input, OUT_SAM));
      final File readStatsFile = new File(input, "readstats");
      FileUtils.stringToFile(readStats, readStatsFile);
      final String outn = output.getPath();
      final String inn = input.getPath();
      final File templ = ReaderTestUtils.getDNADir(refSeq);
      try {
        final String[] svArgs = {
            "sv",
            "-t", templ.getPath(),
            "-o", outn,
            "-r", readStatsFile.toString(),
            inn + FS + OUT_SAM,
            "-s", "1", "-f", "1", "--simple-signals"
        };
        final String[] args = Utils.append(svArgs, args0);
        final ByteArrayOutputStream berr = new ByteArrayOutputStream();
        final String errStr;
        final int intMain;
        try (PrintStream err = new PrintStream(berr); ByteArrayOutputStream out = new ByteArrayOutputStream()) {
          intMain = new Slim().intMain(args, out, err);
        }
        errStr = berr.toString();
        //System.err.println("err\n" + errStr);
        assertEquals("Error" + errStr, errCode, intMain);
        //System.err.println(errStr);
        TestUtils.containsAll(errStr, errorMsgs);
        final String lout = FileUtils.fileToString(new File(outn, "sv.log"));
        if (logMsg != null && !logMsg.equals("")) {
          assertTrue(lout.contains(logMsg));
        }
        if ((errorMsgs.length == 0 || (errorMsgs.length == 1 && errorMsgs[0].equals(""))) && expFile != null) {
          checkFile(output, SvToolParams.NAME_SIMPLE, gzip, expFile);
          checkFile(output, SvToolParams.NAME_BAYESIAN, gzip, expFile + "_bayes");
        }
        if (erlines != null) {
          // check sequencer.errors output file
          final File sefile = new File(output, "sequencer.errors");
          final String result = FileUtils.fileToString(sefile);
          TestUtils.containsAll(result, erlines);
        }
      } finally {
        FileHelper.deleteAll(templ);
      }
    } finally {
      FileHelper.deleteAll(output);
      FileHelper.deleteAll(input);
    }
  }

  static void checkFile(File output, String fileName, boolean gzip, String expFile) throws IOException {
    final File file = new File(output, fileName + (gzip ? ".gz" : ""));
    final String result = gzip ? FileHelper.gzFileToString(file) : FileUtils.fileToString(file);
    //System.out.println(SvToolParams.NAME + ".txt" + " contains:\n" + result);
    final String res = "com/rtg/variant/sv/resources/svtasktest" + expFile + ".txt";

    final InputStream stream = Resources.getResourceAsStream(res);
    assertNotNull("Cant find:" + res, stream);
    try {
      final String expected = FileUtils.streamToString(stream);
      assertTrue(TestUtils.startLines(expected.replaceAll("\t", " "), result.replaceAll("\t", " "), false));
    } finally {
      stream.close();
    }
  }

  public void testNormalize() {
    final double[] n1 = SvToolTask.normalize(new double[] {1.0, 2.0, 3.0});
    assertEquals("[0.298, -0.489, -1.005]", Utils.realFormat(n1, 3));
  }

  public void testNormalize1() {
    final double inc = 1000000.0;
    final double[] n1 = SvToolTask.normalize(new double[] {1.0 + inc, 2.0 + inc, 3.0 + inc});
    assertEquals("[0.298, -0.489, -1.005]", Utils.realFormat(n1, 3));
  }

  public void testNormalize2() {
    final double[] n1 = SvToolTask.normalize(new double[] {1.0, 12.0, 22.0});
    assertEquals("[4.777, -4.777, -9.120]", Utils.realFormat(n1, 3));
  }

  public void testNormalize3() {
    final double[] n1 = SvToolTask.normalize(new double[] {1.0, 50.0, 100.0});
    assertEquals("[21.280, -21.280, -42.995]", Utils.realFormat(n1, 3));
  }

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
    GlobalFlags.resetAccessedStatus();
  }

  public void testInteresting() throws Exception {
    final File tempDir = FileUtils.createTempDir("testInteresting", "test");
    try {
      final File sam = FileHelper.resourceToFile("com/rtg/variant/sv/resources/smallsvMappings.sam.gz", new File(tempDir, "mappings.sam.gz"));
      FileHelper.resourceToFile("com/rtg/variant/sv/resources/smallsvMappings.sam.gz.tbi", new File(tempDir, "mappings.sam.gz.tbi"));
      final File rgstats = FileHelper.resourceToFile("com/rtg/variant/sv/resources/rgstats.txt", new File(tempDir, "rgstats.tsv"));
      final File output = new File(tempDir, "output");
      final File template = new File(tempDir, "template");
      assertTrue(template.mkdir());
      SimpleArchive.unpackArchive(FileHelper.resourceToFile("com/rtg/variant/sv/resources/smallsvTemplate.dwa", new File(tempDir, "template.dwa")), template);

      final MemoryPrintStream mps = new MemoryPrintStream();

      final String[] args = {
        "sv",
        "--Xheterozygous",
        "--simple-signals",
        "-r", rgstats.getPath(),
        "-t", template.getPath(),
        "-o", output.getPath(),
        "-s", "10",
        "--region", "simulatedSequence1:47400-58600",
        "--fine-step", "1",
        sam.getPath()
      };
      final int result = new Slim().intMain(args, mps.outputStream(), mps.printStream());
      assertEquals(mps.toString(), 0, result);

      final String bayes = StringUtils.grep(FileHelper.gzFileToString(new File(output, "sv_bayesian.tsv.gz")), "^[^#]");

      final File expBayesFile = FileHelper.resourceToFile("com/rtg/variant/sv/resources/smallsvBayesExp.txt.gz", new File(tempDir, "expBayes.txt.gz"));
      final String expBayes = StringUtils.grep(FileHelper.gzFileToString(expBayesFile), "^[^#]");
      assertEquals(expBayes, bayes);

      final String interesting = FileHelper.gzFileToString(new File(output, "sv_interesting.bed.gz"));

      TestUtils.containsAll(interesting,
                            "#chr\tstart\tend\tareas\tmaxscore\taverage",
                            "simulatedSequence1\t47400\t47694\t3\t1.0889\t0.3722",
                            "simulatedSequence1\t49712\t52717\t5\t2.5766\t1.1907",
                            "simulatedSequence1\t52747\t52758\t1\t0.0000\t-0.0721",
                            "simulatedSequence1\t52766\t52767\t1\t0.0000\t-0.0630",
                            "simulatedSequence1\t52769\t52789\t1\t0.1322\t0.1099",
                            "simulatedSequence1\t52793\t52807\t1\t0.0000\t-0.0258",
                            "simulatedSequence1\t52810\t53134\t1\t1.8786\t1.1216",
                            "simulatedSequence1\t55032\t55119\t1\t0.4204\t0.0312",
                            "simulatedSequence1\t55190\t55246\t1\t0.1123\t-0.0294",
                            "simulatedSequence1\t55779\t56265\t4\t7.5133\t1.7407",
                            "simulatedSequence1\t58291\t99999\t1\t1.7161\t0.5013"
                            );

      final String simple = StringUtils.grep(FileHelper.gzFileToString(new File(output, "sv_simple.tsv.gz")), "^[^#]");

      final File expSimpleFile = FileHelper.resourceToFile("com/rtg/variant/sv/resources/smallsvSimpleExp.txt.gz", new File(tempDir, "expSimple.txt.gz"));
      final String expSimple = StringUtils.grep(FileHelper.gzFileToString(expSimpleFile), "^[^#]");
      assertEquals(expSimple, simple);

    } finally {
      FileHelper.deleteAll(tempDir);
    }

  }
}
