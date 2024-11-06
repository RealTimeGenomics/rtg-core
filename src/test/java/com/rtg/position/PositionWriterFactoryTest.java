/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
          final ReaderParams srp = new MockReaderParams(reader);
          final ISequenceParams subjectParams = new MockSequenceParams(srp, SequenceMode.BIDIRECTIONAL, 0, reader.numberSequences());

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
