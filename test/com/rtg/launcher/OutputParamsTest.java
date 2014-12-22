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
package com.rtg.launcher;

import java.io.File;
import java.io.StringWriter;

import com.rtg.util.ObjectParamsTest;
import com.rtg.util.cli.CFlags;

/**
 */
public class OutputParamsTest extends ObjectParamsTest {

  @Override
  protected OutputParams getParams(final String[] args) {
    final Appendable out = new StringWriter();
    final CFlags flags = new CFlags("testOutputParams", out, null);
    OutputParams.initFlags(flags);
    flags.setFlags(args);
    final OutputParams params = new OutputParams((File) flags.getValue(CommonFlags.OUTPUT_FLAG), flags.isSet(BuildCommon.PROGRESS_FLAG), !flags.isSet(CommonFlags.NO_GZIP));
    assertTrue(params.closed());
    return params;
  }

}
