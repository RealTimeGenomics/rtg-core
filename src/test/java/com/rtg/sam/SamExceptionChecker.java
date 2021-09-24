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
package com.rtg.sam;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.text.DecimalFormat;
import java.util.HashSet;
import java.util.Random;

import htsjdk.samtools.SAMFormatException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;

/**
 * Test the quality of Picard exception handling when parsing damaged SAM files.
 */
public final class SamExceptionChecker {

  private SamExceptionChecker() { }

  private static final ValidationStringency[] MODES = {
    ValidationStringency.SILENT,
    ValidationStringency.LENIENT,
    ValidationStringency.STRICT
  };

  private static final byte[] SAM_INPUT = (""
    + "@HD\tVN:1.0\tSO:coordinate\n"
    + "@PG\tID:abc\tVN:test\tCL:map\n"
    + "@CO\tREAD-SDF-ID:4f083ce3f12bbb4b\n"
    + "@SQ\tSN:ss1\tLN:10000000\n"
    + "@SQ\tSN:ss2\tLN:10000000\n"
    + "4\t99\tss1\t74719\t255\t40=\t=\t74903\t225\tCTATTTACCCGACCACCGCGAGAAATTAAGGCTCCGGCCC\t*\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:2\tIH:i:1\n"
    + "4\t147\tss1\t74903\t255\t16=1D24=\t=\t74719\t-225\tGGAAGGAGACTCACCAGAGGGACTTGTTGACTGGTCCTGT\t*\tAS:i:2\tNM:i:1\tMQ:i:255\tXA:i:2\tIH:i:1\n"
    + "9\t99\tss1\t3565772\t255\t1=1X24=1I13=\t=\t3565967\t235\tTTGAGGACATTTCGCTGCGGATCACGTTTTGTGAATAATT\t*\tAS:i:3\tNM:i:2\tMQ:i:255\tXA:i:4\tIH:i:1\n"
    + "9\t147\tss1\t3565967\t255\t16=1X23=\t=\t3565772\t-235\tGATAAAATATGTCTCTTTTGCTGAACGAAAAAAGGCCGAA\t*\tAS:i:1\tNM:i:1\tMQ:i:255\tXA:i:4\tIH:i:1\n"
    + "5\t163\tss1\t5385195\t255\t15=1I24=\t=\t5385369\t214\tCGGTACCTCAGCCCGTCGGTGATAGAGTGGTTTTCGCTTC\t*\tAS:i:2\tNM:i:1\tMQ:i:255\tXA:i:2\tIH:i:1\n"
    + "5\t83\tss1\t5385369\t255\t40=\t=\t5385195\t-214\tGAGTCTAACCTAGCGAGAAGTAAGACCACCCCCGGTCGCA\t*\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:2\tIH:i:1\n"
    + "6\t99\tss1\t5421632\t255\t40=\t=\t5421815\t223\tATCCAGTTGTCACCATCATCTCCATACGTGATGACTAAGT\t*\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tIH:i:1\n"
    + "6\t147\tss1\t5421815\t255\t40=\t=\t5421632\t-223\tGCGAACTCGTTACCTTCATGGCGCCCTGCGGGACCGTAGT\t*\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tIH:i:1\n"
    + "8\t99\tss1\t5500417\t255\t35=1X4=\t=\t5500607\t230\tTGGGCTCAAATCCAACAGACCGTACACACGGCAACTGCTT\t*\tAS:i:1\tNM:i:1\tMQ:i:25\tXA:i:1\tIH:i:1\n"
    + "8\t147\tss1\t5500607\t255\t40=\t=\t5500417\t-230\tGGCACGTCTGACGGTACCACGTGCGAAACTCACAGGGCTT\t*\tAS:i:0\tNM:i:0\tMQ:i:55\tXA:i:1\tIH:i:1\n"
    + "3\t163\tss1\t6721963\t255\t40M\t=\t6722136\t213\tCGTCTATTACCTCGCTCGAAGGTCTACATGTCAAAGCCCC\t*\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tIH:i:1\n"
    + "3\t83\tss1\t6722136\t255\t40=\t=\t6721963\t-213\tCACCTGCAAGGCAGTCATTATAAGGCACCGTGCCTATATG\t*\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:0\tIH:i:1\n"
    + "7\t99\tss1\t8734083\t255\t40=\t=\t8734276\t233\tAATCAGCATTCTGCCACACCTTTAGTGCCTGTGATAATGG\t*\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:3\tIH:i:1\n"
    + "7\t147\tss1\t8734276\t255\t1X14=2X23=\t=\t8734083\t-233\tGGACGGCCACCCATGCGTCTACTTATGAGCTTACGCGATG\t550#55##0###0F5#5F#F0#FF0F505#FF0#5#50F5\tAS:i:3\tNM:i:3\tMQ:i:255\tXA:i:3\tIH:i:1\n"
    + "0\t163\tss1\t9047338\t255\t21=1X18=\t=\t9047516\t219\tACGCTCTGCATATGCTATCAAAGTGCGAAGAATAACTTGT\t*\tAS:i:1\tNM:i:1\tMQ:i:255\tXA:i:3\tIH:i:1\n"
    + "0\t83\tss1\t9047516\t255\t37=1D3=\t=\t9047338\t-219\tTGCATGTTCCACTTAGTGTAGAGCTGTTCGGAGTATAGAG\t*\tAS:i:2\tNM:i:1\tMQ:i:255\tXA:i:3\tIH:i:1\n"
    + "2\t99\tss2\t1927565\t255\t40=\t=\t1927738\t213\tTAACTTGTTCCGTGTGTCGTTAAGCCCCGGAGGGAGTGAC\t*\tAS:i:0\tNM:i:0\tMQ:i:255\tXA:i:1\tIH:i:1\n"
    + "2\t147\tss2\t1927738\t255\t30=1X9=\t=\t1927565\t-213\tCCTAAGGCCGCGTCCTTCCCGTCACTCATACCGGTGTATG\t*\tAS:i:1\tNM:i:1\tMQ:i:255\tXA:i:1\tIH:i:1\n").getBytes();


  static final int DEFAULT_LIMIT = 10;

  private static String fix(final String n) {
    if (n.charAt(1) == '.') {
      return " " + n;
    }
    return n;
  }

  private static void recordTrace(final HashSet<String> traces, final Exception e) throws IOException {
    final StringWriter sw = new StringWriter();
    try {
      try (PrintWriter pw = new PrintWriter(sw)) {
        e.printStackTrace(pw);
      }
    } finally {
      sw.close();
    }
    // Shorten and remove variable pieces from stack traces, to get a
    // better representative list of problems.
    traces.add(sw.toString()
               .replaceFirst("(:[^:]*:).*", "$1")
               .replaceFirst("(SAMFileHeader\\$SortOrder).*", "$1")
               .replaceAll(".*\\.SamExceptionChecker\\..*", "")
               .trim());
  }

  /**
   * Check what happens when we randomly modify a valid SAM input.  This is
   * to collect information about unexpected exceptions coming out of the
   * Picard implementation of SAM reading.
   */
  static void check(final int limit) throws IOException {

    // Record distinct stack traces, so they can be output later
    final HashSet<String> traces = new HashSet<>();

    // Replace default error stream because LENIENT writes messages to screen,
    // for this checking we just ignore such messages.
    final PrintStream oldErr = System.err;
    try {
      try (PrintStream newErr = new PrintStream(new ByteArrayOutputStream())) {
        System.setErr(newErr);

        // Select a seed, so that same sequence of changes can be applied
        // to each of the stringency modes.
        final long seed = System.nanoTime();
        for (final ValidationStringency stringency : MODES) {
          if (limit != DEFAULT_LIMIT) {
            System.out.println();
            System.out.println("Mode: " + stringency.toString());
          }
          final SamReaderFactory factory = SamReaderFactory.make().validationStringency(stringency);
          final Random rnd = new Random(seed);
          final long startTime = System.currentTimeMillis();

          // Counters for various types of parsing events
          int format = 0;
          int illegal = 0;
          int array = 0;
          int runtime = 0;
          int validation = 0;
          int headerOk = 0;

          for (int k = 0; k < limit; ++k) {

            // Introduce a mutation into the SAM input
            final int size = SAM_INPUT.length;
            final byte[] mutant = new byte[size];
            System.arraycopy(SAM_INPUT, 0, mutant, 0, size);
            final int placeToMutate = rnd.nextInt(size);
            final int current = mutant[placeToMutate];
            int replace;
            do {
              replace = rnd.nextInt(256);
            } while ((byte) replace == (byte) current);
            mutant[placeToMutate] = (byte) replace;

            // Now try and iterate over its records using the SAM reader
            try {
              ++headerOk;
              try (SamReader r = factory.open(SamInputResource.of(new ByteArrayInputStream(mutant)))) {
                try (CloseableIterator<SAMRecord> it = r.iterator()) {
                  while (it.hasNext()) {
                    it.next();
                  }
                }
              }
            } catch (final SAMFormatException e) {
              ++format;
            } catch (final IllegalArgumentException e) {
              ++illegal;
              recordTrace(traces, e);
            } catch (final ArrayIndexOutOfBoundsException e) {
              ++array;
              recordTrace(traces, e);
            } catch (final RuntimeException e) {
              if (e.getMessage().startsWith("SAM validation error:")) {
                ++validation;
              } else {
                ++runtime;
                recordTrace(traces, e);
              }
            }
          }
          if (limit != DEFAULT_LIMIT) {
            final DecimalFormat nf = new DecimalFormat("0.00");
            System.out.println("Total runs: " + limit);
            System.out.println(fix(nf.format(100.0 * format / limit)) + "% SAMFormatException: " + format + " (this is an acceptable class)");
            System.out.println(fix(nf.format(100.0 * illegal / limit)) + "% IllegalArgumentException: " + illegal);
            System.out.println(fix(nf.format(100.0 * array / limit)) + "% ArrayIndexOutOfBoundsException: " + array);
            System.out.println(fix(nf.format(100.0 * validation / limit)) + "% RuntimeException: SAM validation error: " + validation);
            System.out.println(fix(nf.format(100.0 * runtime / limit)) + "% Other RuntimeException: " + runtime);
            final int undetected = limit - format - illegal - array - runtime - validation;
            System.out.println(fix(nf.format(100.0 * undetected / limit)) + "% Undetected change: " + undetected + " (this is an acceptable class)");
            System.out.println();
            System.out.println(fix(nf.format(100.0 * headerOk / limit)) + "% Number of times header successfully parsed: " + headerOk);
            final double duration = (System.currentTimeMillis() - startTime) / 1000.0;
            System.out.println("Duration: " + duration + "s");

          }
        }

        if (traces.size() > 0 && limit != DEFAULT_LIMIT) {
          System.out.println("Observed traces not extending SAMFormatException");
          for (final String t : traces) {
            System.out.println();
            System.out.println(t);
          }
        }
      }
    } finally {
      System.setErr(oldErr);
    }
  }

  /**
   * Run a mutation check against a valid SAM input.
   *
   * @param args number of iterations
   * @exception IOException if an I/O error occurs.
   */
  public static void main(final String[] args) throws IOException {
    check(args == null || args.length == 0 ? DEFAULT_LIMIT : Integer.parseInt(args[0]));
  }
}
