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
package com.rtg.position;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;

import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.MockReaderParams;
import com.rtg.launcher.MockSequenceParams;
import com.rtg.launcher.ReaderParams;
import com.rtg.launcher.HashingRegion;
import com.rtg.mode.SequenceMode;
import com.rtg.mode.SequenceType;
import com.rtg.ngs.DefaultOutputProcessorSynch;
import com.rtg.ngs.NgsOutputParams;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.NgsTestUtils;
import com.rtg.position.PositionWriterFactory.NgsPositionWriterFactory;
import com.rtg.position.output.OutputProcessorWrapper;
import com.rtg.position.output.PositionWriter;
import com.rtg.reader.MockArraySequencesReader;
import com.rtg.reader.SequencesReader;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;


/**
 */
public class PositionWriterFactoryTest extends TestCase {

  public void test2() throws IOException {
    final File dir = FileHelper.createTempDirectory();
    try {
      try (ByteArrayOutputStream b = new ByteArrayOutputStream()) {
        //final PositionOutputParams params = new PositionOutputParams(new File("foobar"), OutputFormatType.MAP, null, Double.valueOf(0.0), false, 5);
        final NgsOutputParams outparams = new NgsTestUtils.OverriddenNgsOutputParams(NgsTestUtils.OverriddenNgsOutputParams.builder().outStream(b).outputDir(dir));
        try (DefaultOutputProcessorSynch dops = new DefaultOutputProcessorSynch(NgsParams.builder().outputParams(outparams).create())) {
          final NgsPositionWriterFactory fact = new NgsPositionWriterFactory(dops);

          final SequencesReader reader = new MockArraySequencesReader(SequenceType.DNA, new int[]{30, 31, 48});
          final ReaderParams srp = new MockReaderParams(reader, SequenceMode.BIDIRECTIONAL);
          final ISequenceParams subjectParams = new MockSequenceParams(srp, 0, reader.numberSequences());

          final PositionWriter pw = fact.makeNgs(subjectParams, null, subjectParams, HashingRegion.NONE);
          assertNotNull(pw);
          assertTrue(pw instanceof OutputProcessorWrapper);
        }
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }
}
