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

package com.rtg.assembler;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;

/**
 */
public class CorrectReadsCliTest extends AbstractCliTest {
  @Override
  protected AbstractCli getCli() {
    return new CorrectReadsCli();
  }
  public void testInitParams() {
    checkHelp("correctreads [OPTION]... -i SDF -k INT -o DIR",
        "-i, --input", "read SDF to correct",
        "-k, --kmer-size=INT", "size of kmer to use in correction",
        "-o, --output", "output directory",
        "-c, --threshold", "override the calculated frequency threshold"
    );
  }
}
