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
package com.rtg.position.output;

import java.io.File;
import java.io.IOException;
import java.io.StringWriter;

import com.rtg.launcher.BuildParams;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.MockReaderParams;
import com.rtg.launcher.MockSequenceParams;
import com.rtg.launcher.ReaderParams;
import com.rtg.mode.SequenceMode;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 * Test the enumeration OutputType.
 */
public class OutputFormatTypeTest extends TestCase {

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }

  public void test() {
    assertEquals("[SEGMENT, NGS]", OutputFormatType.values().toString());
  }


  private PositionParams makeParams(final OutputFormatType type) throws IOException {
    final ReaderParams srp = new MockReaderParams(0, 0, SequenceMode.UNIDIRECTIONAL);
    final ISequenceParams subjectParams = new MockSequenceParams(srp, 0, 0);
    final BuildParams buildParams = BuildParams.builder().windowSize(5).stepSize(2).sequences(subjectParams).create();

    final ReaderParams qrp = new MockReaderParams(0, 0, SequenceMode.BIDIRECTIONAL);
    final ISequenceParams queryParams = new MockSequenceParams(qrp, 0, 0);
    final BuildParams searchParams = BuildParams.builder().windowSize(5).stepSize(1).sequences(queryParams).create();
    final PositionDistributionParams distr = new PositionDistributionParams(0.001, 0.009, 20, 0);
    final PositionOutputParams outParams = new PositionOutputParams(new File(""), type, distr, null, false, 10);
    return PositionParams.builder().hashCountThreshold(1000).buildParams(buildParams).searchParams(searchParams).outputParams(outParams).create();
  }

  public void testOutput() throws IOException {
    final Appendable out = new StringWriter();
    final Appendable unmappedOut = new StringWriter();
    for (final OutputFormatType type : OutputFormatType.values()) {
      //System.err.println(type);
      final PositionParams params = makeParams(type);
      final GapBucketsInfo bucketInfo = new GapBucketsInfo(params, 1);
      final PositionOutput output = type.output(params, bucketInfo, null, out, unmappedOut, null, params.hashCountThreshold());
      assertNotNull(output);
    }
  }

  /** Check works ok with CFlags. */
  /**
   We don't seem to be doing this anywhere anymore
  public void testCFlags() {
    for (final OutputFormatType<?> type : OutputFormatType.values()) {
      checkFlag(type.toString(), type);
      checkFlag(type.toString().toLowerCase(Locale.getDefault()), type);
    }
    checkFlag("wOrD", OutputFormatType.WORD);
  }

  private void checkFlag(final String output, final OutputFormatType<?> expected) {
    final CFlags flags = new CFlags();
    flags.registerOptional('o', "output", OutputFormatType.class, "", "");
    assertTrue(flags.setFlags(new String[] {"-o", output}));
    final OutputFormatType<?> res = (OutputFormatType<?>) flags.getValue("output");
    assertEquals(expected, res);
  }
   */
}
