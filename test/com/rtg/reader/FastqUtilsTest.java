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

import java.io.File;
import java.io.IOException;
import java.io.StringWriter;

import com.rtg.mode.DnaUtils;

import junit.framework.TestCase;

/**
 *
 */
public class FastqUtilsTest extends TestCase {

  public void testBaseFile() {
    checkBaseFile("input", ".fastq", true, FastqUtils.baseFile(new File("input"), true));
    checkBaseFile("input", ".fastq", true, FastqUtils.baseFile(new File("input.fastq"), true));
    checkBaseFile("input", ".fastq", true, FastqUtils.baseFile(new File("input.fastq.gz"), true));
    checkBaseFile("input", ".fq", true, FastqUtils.baseFile(new File("input.fq"), true));
    checkBaseFile("input", ".fq", true, FastqUtils.baseFile(new File("input.fq.gz"), true));
    checkBaseFile("input.fast", ".fastq", true, FastqUtils.baseFile(new File("input.fast.gz"), true));

    final FastqUtils.BaseFile bf = FastqUtils.baseFile(new File("input.fastq.gz"), true);
    assertEquals("input.fastq.gz", bf.suffixedFile("").getName());
    assertEquals("input_moo.fastq.gz", bf.suffixedFile("_moo").getName());
    final FastqUtils.BaseFile bf2 = FastqUtils.baseFile(new File("input.fastq.gz"), false);
    assertEquals("input.fastq", bf2.suffixedFile("").getName());
    assertEquals("input_moo.fastq", bf2.suffixedFile("_moo").getName());
  }

  public void testWriteSequence() throws IOException {
    try (final StringWriter sw = new StringWriter()) {
      final String seqString = "GATCAGGTAGTT";
      final byte[] seqData = DnaUtils.encodeArray(seqString.getBytes());
      final String qualString = "%^@%#^@#*HDA";
      final byte[] qualData = FastaUtils.asciiToRawQuality(qualString);
      FastqUtils.writeFastqSequence(sw, "thefirstsequence", seqData, qualData);
      assertEquals("@thefirstsequence\n" + seqString + "\n+\n" + qualString + "\n", sw.toString());
    }
  }

  private void checkBaseFile(String expectedBase, String expExtension, boolean gz, FastqUtils.BaseFile res) {
    checkBaseFile(new File(expectedBase), expExtension, gz, res);
  }
  private void checkBaseFile(File expectedBase, String expExtension, boolean gz, FastqUtils.BaseFile res) {
    assertEquals(expectedBase, res.getBaseFile());
    assertEquals(expExtension, res.getExtension());
    assertEquals(gz, res.isGzip());
  }
}