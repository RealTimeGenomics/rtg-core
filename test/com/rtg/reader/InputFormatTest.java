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

import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class InputFormatTest extends TestCase {

  /**
   */
  public InputFormatTest(String name) {
    super(name);
  }

  public void testEnum() {
    TestUtils.testEnum(InputFormat.class, "[FASTA, FASTQ, SOLEXA, SOLEXA1_3, CG, SAM_SE, SAM_PE]");
  }
}
