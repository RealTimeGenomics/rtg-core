/*
 * Copyright (c) 2018. Real Time Genomics Limited.
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

import static com.rtg.sam.SharedSamConstants.SAMHEADER1;
import static com.rtg.util.StringUtils.LS;
import static com.rtg.util.StringUtils.TAB;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;

/**
 * Tests for corresponding class.
 */
public class SoftClip2FastqTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new SoftClip2Fastq();
  }

  public void testHelp() {
    checkHelp("containing long soft clips",
      "output only the portion",
      "minimum length of soft clip"
    );
  }

  public void testValidator() throws IOException {
    try (TestDirectory tmp = new TestDirectory()) {
      final File bam = new File(tmp, "blah.bam");
      TestUtils.containsAllUnwrapped(checkMainInitBadFlags(bam.getPath(), "-o", "-"),
        "Error: File not found");
      assertTrue(bam.createNewFile());
      TestUtils.containsAllUnwrapped(checkMainInitBadFlags(bam.getPath(), "-o", "-"),
        "Error: Sending non-interleaved paired-end");
      TestUtils.containsAllUnwrapped(checkMainInitBadFlags(bam.getPath(), "-o", "-", "--interleave", "--min-soft-clip-length", "-3"),
        "min-soft-clip-length must be at least 1");
    }
  }

  static final String READ =  "ATCGACTGATCGACTGATCGACTGATCGACTGATCGACTGAAAATTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTT";
  static final String READQ = "AABBAABBAABBAABBAABBAABB````````````````AAAAAAAABBBBBBBBCCCCCCCC````````````````";
  static final String SAM = SAMHEADER1
    + "0" + TAB + "0" + TAB + "g1" + TAB +  "3" + TAB + "255" + TAB + "80M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + READ + TAB + READQ + TAB + "RG:Z:RG1" + LS
    + "1" + TAB + "0" + TAB + "g1" + TAB +  "3" + TAB + "255" + TAB + "50S30M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + READ + TAB + READQ + TAB + "RG:Z:RG1" + LS
    + "2" + TAB + "16" + TAB + "g1" + TAB +  "3" + TAB + "255" + TAB + "40M40S" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + READ + TAB + READQ + TAB + "RG:Z:RG1" + LS
    + "3" + TAB + "65" + TAB + "g1" + TAB +  "3" + TAB + "255" + TAB + "80M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + READ + TAB + READQ + TAB + "RG:Z:RG1" + LS
    + "4" + TAB + "65" + TAB + "g1" + TAB +  "3" + TAB + "255" + TAB + "50S30M" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + READ + TAB + READQ + TAB + "RG:Z:RG1" + LS
    + "5" + TAB + "129" + TAB + "g1" + TAB +  "3" + TAB + "255" + TAB + "40M40S" + TAB + "*" + TAB + "0" + TAB + "0" + TAB + READ + TAB + READQ + TAB + "RG:Z:RG1" + LS
    ;

  public void testSimple() throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      final File sam = FileUtils.stringToFile(SAM, new File(dir, "small.sam"));
      final File fastq = new File(dir, "small.fastq");
      checkMainInitOk(sam.getPath(), "-Z", "--interleave", "-o", fastq.getPath());
      mNano.check("softclip2fastq.fq", FileUtils.fileToString(fastq));

      checkMainInitOk(sam.getPath(), "-Z", "-o", fastq.getPath());
      mNano.check("softclip2fastq_1.fq", FileUtils.fileToString(new File(dir, "small_1.fastq")));
      mNano.check("softclip2fastq_2.fq", FileUtils.fileToString(new File(dir, "small_2.fastq")));

      checkMainInitOk(sam.getPath(), "-Z", "--soft-clipped-only", "-o", fastq.getPath());
      mNano.check("softclip2fastq_1.fq", FileUtils.fileToString(fastq));

      final File fastqrf = new File(dir, "smallrf.fastq");
      checkMainInitOk(sam.getPath(), "-Z", "--interleave", "--orientation", "rf", "-o", fastqrf.getPath());
      mNano.check("softclip2fastq_rf.fq", FileUtils.fileToString(fastqrf));
    }
  }
}
