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
package com.rtg.ngs;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.List;

import com.rtg.sam.SamUtils;
import com.rtg.util.PortableRandom;
import com.rtg.util.diagnostic.Timer;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMValidationError;

/**
 * Test memory and limits of SAM format.
 *
 */
public final class TestSamOutput {

  private TestSamOutput() { }

  // Realistic human names and lengths
  private static final String[] HUMAN_NAMES = {
    "gi|89161185|ref|NC_000001.9|NC_000001",
    "gi|89161187|ref|NC_000010.9|NC_000010",
    "gi|42406306|ref|NC_000019.8|NC_000019",
    "gi|89161199|ref|NC_000002.10|NC_000002",
    "gi|51511747|ref|NC_000020.9|NC_000020",
    "gi|51511750|ref|NC_000021.7|NC_000021",
    "gi|89161203|ref|NC_000022.9|NC_000022",
    "gi|89161205|ref|NC_000003.10|NC_000003",
    "gi|89161207|ref|NC_000004.10|NC_000004",
    "gi|51511721|ref|NC_000005.8|NC_000005",
    "gi|89161210|ref|NC_000006.10|NC_000006",
    "gi|89161213|ref|NC_000007.12|NC_000007",
    "gi|51511727|ref|NC_000011.8|NC_000011",
    "gi|51511724|ref|NC_000008.9|NC_000008",
    "gi|89161216|ref|NC_000009.10|NC_000009",
    "gi|89161218|ref|NC_000023.9|NC_000023",
    "gi|89161220|ref|NC_000024.8|NC_000024",
    "gi|89161190|ref|NC_000012.10|NC_000012",
    "gi|51511729|ref|NC_000013.9|NC_000013",
    "gi|51511730|ref|NC_000014.7|NC_000014",
    "gi|51511731|ref|NC_000015.8|NC_000015",
    "gi|51511732|ref|NC_000016.8|NC_000016",
    "gi|51511734|ref|NC_000017.9|NC_000017",
    "gi|51511735|ref|NC_000018.8|NC_000018",
  };

  private static final int[] HUMAN_LENGTHS = {
    247249759,
    135374777,
    63811691,
    242951190,
    62436004,
    46944363,
    49691472,
    199501868,
    191273104,
    180857906,
    170900033,
    158821465,
    134452424,
    146274866,
    140273293,
    154913794,
    57772994,
    132349575,
    114143020,
    106368625,
    100338955,
    88827294,
    78774782,
    76117193,
  };

  private static final String[] CIGARS = {
    "36M",
    "3I2M1D22M",
    "30M6I",
    "30M6D",
    "32M",
    "1I2D3M",
  };

  private static SAMFileHeader createHeader(final boolean sorted) {
    final SAMFileHeader header = new SAMFileHeader();
    header.setSortOrder(sorted ? SAMFileHeader.SortOrder.coordinate : SAMFileHeader.SortOrder.unsorted);
    for (int k = 0; k < HUMAN_NAMES.length; ++k) {
      final SAMSequenceRecord record = new SAMSequenceRecord(HUMAN_NAMES[k], HUMAN_LENGTHS[k]);
      header.addSequence(record);
    }
    return header;
  }

  private static final String BASES = "ACGTN";

  private static String randomRead(final PortableRandom rnd, final int len) {
    final StringBuilder sb = new StringBuilder();
    for (int k = 0; k < len; ++k) {
      sb.append(BASES.charAt(rnd.nextInt(5)));
    }
    return sb.toString();
  }

  private static byte[] randomQuality(final PortableRandom rnd, final int len) {
    final byte[] r = new byte[len];
    for (int k = 0; k < r.length; ++k) {
      r[k] =  (byte) rnd.nextInt(63);
    }
    return r;
  }

  private static boolean valid(final SAMRecord rec) {
    final List<SAMValidationError> errors = rec.isValid();
    if (errors == null) {
      return true;
    }
    System.err.println(errors.size() + " errors were found");
    for (final SAMValidationError e : errors) {
      System.err.println(e.toString());
    }
    return false;
  }

  private static SAMRecord randomRecord(final SAMFileHeader header, final PortableRandom r, final int reads) {
    final SAMRecord rec = new SAMRecord(header);
    final boolean rc = r.nextBoolean();
    rec.setReadName(String.valueOf(r.nextInt(reads)));
    final int t = r.nextInt(HUMAN_NAMES.length);
    rec.setReferenceName(HUMAN_NAMES[t]);
    rec.setInferredInsertSize(0);
    rec.setMateAlignmentStart(0);
    rec.setMateReferenceName("*");
    final boolean pairedEnd = r.nextBoolean();
    rec.setReadPairedFlag(pairedEnd);
    rec.setReadNegativeStrandFlag(rc);
    if (pairedEnd) {
      rec.setProperPairFlag(false);
      final boolean first = r.nextBoolean();
      rec.setFirstOfPairFlag(first);
      rec.setSecondOfPairFlag(!first);
      rec.setMateUnmappedFlag(false);
      rec.setInferredInsertSize(r.nextInt(2000));
      rec.setProperPairFlag(true);
      rec.setMateNegativeStrandFlag(r.nextBoolean());
      rec.setMateAlignmentStart(1 + r.nextInt(HUMAN_LENGTHS[t] - 500));
      rec.setMateReferenceName(HUMAN_NAMES[t]);
      rec.setAttribute("MQ", 255);
    }
    rec.setAlignmentStart(1 + r.nextInt(HUMAN_LENGTHS[t] - 500));
    rec.setMappingQuality(255);
    rec.setCigarString(CIGARS[r.nextInt(CIGARS.length)]);
    rec.setReadString(randomRead(r, 36));
    rec.setBaseQualities(randomQuality(r, 36));
    rec.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, r.nextInt(15));
    assert valid(rec);
    return rec;
  }

  private static long du(final String dir) throws IOException {
    final Process p = new ProcessBuilder("du", "-s", dir).start();
    try {
      try (BufferedReader r = new BufferedReader(new InputStreamReader(p.getInputStream()))) {
        final String line = r.readLine();
        if (line == null) {
          throw new IOException();
        }
        return Long.parseLong(line.split("\\s+")[0]);
      }
    } finally {
      p.getErrorStream().close();
      p.getInputStream().close();
      p.getOutputStream().close();
    }
  }

  private static void test(final File file, final int records, final boolean sorted) throws IOException {
    final SAMFileHeader header = createHeader(sorted);
    final SAMFileWriter w = new SAMFileWriterFactory().makeSAMWriter(header, false, file);
    final PortableRandom r = new PortableRandom();
    final long initial = du("/tmp");
    final Timer t0 = new Timer("Adding");
    for (int k = 0; k < records; ++k) {
      t0.start();
      w.addAlignment(randomRecord(header, r, records)); // randomness will give some same id
      t0.stop();
    }
    final long after = du("/tmp");
    final Timer t1 = new Timer("Closing");
    t1.start();
    w.close();
    t1.stop();
    System.out.println(t0);
    System.out.println(t1);
    System.out.println("Disk usage: " + (after - initial) / 1024);
  }

  /**
   * Simulate SAM output.
   *
   * @param args number of records
   * @exception Exception if an error occurs.
   */
  public static void main(final String[] args) throws Exception {
    if (args.length == 0) {
      System.err.println("Usage: TestSamOutput records [unsorted]");
      return;
    }
    final int records = Integer.parseInt(args[0]);
    final File f = new File("sam.out");
    test(f, records, args.length == 1);
  }

}
