/*
 * Copyright (c) 2017. Real Time Genomics Limited.
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

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.TestUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class PairedEndTrimCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new PairedEndTrimCli();
  }

  public void testHelp() {
    checkHelp("rtg petrim"
      , "-l", "--left=FILE", "left input FASTQ file (AKA R1)"
      , "-r", "--right=FILE", "right input FASTQ file (AKA R2)"
      , "-o", "--output=FILE", "output filename prefix. Use '-' to write to standard output"
      , "-q", "--quality-format=FORMAT", "quality data encoding", "Allowed values are [sanger, solexa, illumina] (Default is sanger)"
      , "-Z", "--no-gzip", "do not gzip the output"
      , "-T", "--threads=INT", "number of threads (Default is the number of available cores)"
        , "--aligner-band-width=FLOAT", "aligner indel band width scaling factor, fraction of read length allowed as an indel (Default is 0.5)"
    , "--gap-extend-penalty=INT", "penalty for a gap extension during alignment (Default is 1)"
      , "--gap-open-penalty=INT", "penalty for a gap open during alignment (Default is 19)"
      , "--mismatch-penalty=INT", "penalty for a mismatch during alignment (Default is 9)"
    , "--soft-clip-distance=INT", "soft clip alignments if indels occur INT bp from either end (Default is 5)"
      , "--unknowns-penalty=INT", "penalty for unknown nucleotides during alignment (Default is 5)"
    , "--min-identity=INT", "minimum percent identity in overlap"
    , "--min-overlap-length=INT", "minimum number of bases in overlap"
    , "--left-probe-length=INT", "trim R2 bases that overlap"
      , "--interleave", "interleave paired data into a single output file. Default is to split to separate output files"
    );
    checkExtendedHelp("rtg petrim",
      "--Xbatch-size=INT",  "number of pairs to process per batch"
    );
  }

  public void testFlags() throws IOException {
    assertParseMessage("You must provide a value for -o FILE", "-l", "foo" , "-r", "bar");
    assertParseMessage("You must provide a value for -r FILE", "-l", "foo" , "-o", "bar");
    assertParseMessage("You must provide a value for -l FILE", "-r", "foo" , "-o", "bar");
    try (TestDirectory tmp = new TestDirectory()) {
      final File[] files = {
        new File(tmp, "foo"),
        new File(tmp, "bar"),
        new File(tmp, "bazout.fastq.gz"),
      };
      for (File f : files) {
        assertTrue(f.createNewFile());
      }

      assertParseMessage("must be in the range [1, 100]", "--min-identity", "101", "-l", files[0].getPath(), "-r", files[1].getPath(), "-o", files[2].getPath());
      assertParseMessage("must be at least 1", "--min-overlap-length", "0", "-l", files[0].getPath(), "-r", files[1].getPath(), "-o", files[2].getPath());
      assertParseMessage("must be at least 1", "--Xbatch-size", "0", "-l", files[0].getPath(), "-r", files[1].getPath(), "-o", files[2].getPath());
      assertParseMessage("must be at least 0", "--left-probe-length", "-1", "-l", files[0].getPath(), "-r", files[1].getPath(), "-o", files[2].getPath());
      assertParseMessage("non-interleaved paired-end data to stdout is not supported", "-l", files[0].getPath(), "-r", files[1].getPath(), "-o", "-");
      checkHandleFlags("-l", files[0].getPath(), "-r", files[1].getPath(), "--interleave", "-o", files[2].getPath());
      assertParseMessage("already exists", "-l", files[0].getPath(), "-r", files[1].getPath(), "--interleave", "-o", files[2].getPath());
      checkHandleFlags("-l", files[0].getPath(), "-r", files[1].getPath(), "--interleave", "-o", files[2].getPath(), "--Xforce");
    }
  }

  private void assertParseMessage(String expected, String... strings) {
    TestUtils.containsAllUnwrapped(checkHandleFlagsErr(strings), expected);
  }

  public void testEndToEnd() throws IOException {
    try (TestDirectory dir = new TestDirectory("petrim")) {
      final File left = new File(dir, "reads50_R1.gz");
      final File right = new File(dir, "reads50_R2.gz");
      FileHelper.resourceToGzFile("com/rtg/reader/resources/reads50_R1.fastq", left);
      FileHelper.resourceToGzFile("com/rtg/reader/resources/reads50_R2.fastq", right);

      final File out1 = new File(dir, "readsout1.fastq");
      checkMainInitOk("-l", left.toString(), "-r", right.toString(), "-o", out1.toString(), "--interleave", "-Z", "--left-probe-length=30", "--right-probe-length=30");
      mNano.check("petrim-e2e-probes.fastq", FileHelper.fileToString(out1));
      mNano.check("petrim-e2e-probes.summary.txt", FileHelper.fileToString(new File(dir, "readsout1.summary.txt")));
      mNano.check("petrim-e2e-probes.overlap-lengths.tsv", FileHelper.fileToString(new File(dir, "readsout1.overlap-lengths.tsv")));
      mNano.check("petrim-e2e-probes.left-read-lengths.tsv", FileHelper.fileToString(new File(dir, "readsout1.left-read-lengths.tsv")));
      mNano.check("petrim-e2e-probes.right-read-lengths.tsv", FileHelper.fileToString(new File(dir, "readsout1.right-read-lengths.tsv")));

      final File out2 = new File(dir, "readsout2.fastq");
      checkMainInitOk("-l", left.toString(), "-r", right.toString(), "-o", out2.toString(), "--interleave", "-Z", "--left-probe-length=30", "--right-probe-length=30", "--min-read-length=30");
      mNano.check("petrim-e2e-min-length.fastq", FileHelper.fileToString(out2));

      final File out3 = new File(dir, "readsout3.fastq");
      checkMainInitOk("-l", left.toString(), "-r", right.toString(), "-o", out3.toString(), "--interleave", "-Z", "--left-probe-length=30", "--right-probe-length=30", "--min-read-length=30", "--discard-empty-pairs");
      mNano.check("petrim-e2e-drop-pairs.fastq", FileHelper.fileToString(out3));

      final File out4 = new File(dir, "readsout4.fastq");
      checkMainInitOk("-l", left.toString(), "-r", right.toString(), "-o", out4.toString(), "--interleave", "-Z", "--left-probe-length=30", "--right-probe-length=30", "--min-read-length=30", "--discard-empty-reads");
      mNano.check("petrim-e2e-drop-reads.fastq", FileHelper.fileToString(out4));
    }
  }
}
