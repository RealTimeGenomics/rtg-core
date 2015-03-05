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
package com.rtg.util;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.util.io.FileUtils;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;

/**
 */
public final class TestSamThroughput {

  private TestSamThroughput() { }

  /**
   * @param args command line arguments
   * @throws IOException if an IO error occurs
   */
  public static void main(String[] args) throws IOException {
    try (SAMFileReader sr = new SAMFileReader(new File(args[0]))) {
      final File output = File.createTempFile("samor", "bam", new File("."));
      SAMFileWriter writer;
      final OutputStream out;
      if (args[1].equals("bam")) {
        out = FileUtils.createOutputStream(output, true);
        writer = new SAMFileWriterFactory().makeBAMWriter(sr.getFileHeader(), true, out, true);
      } else {
        out = FileUtils.createOutputStream(output, true);
        writer = new SAMFileWriterFactory().makeSAMWriter(sr.getFileHeader(), true, out);
      }
      try {
        for (SAMRecord aSr : sr) {
          writer.addAlignment(aSr);
        }
      } finally {
        writer.close();
      }
    }
  }
}
