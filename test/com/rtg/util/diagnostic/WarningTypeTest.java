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
package com.rtg.util.diagnostic;

import java.io.ObjectStreamException;

import com.rtg.util.TestUtils;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Tests the corresponding class.
 *
 */
public class WarningTypeTest extends TestCase {

  /**
   */
  public WarningTypeTest(final String name) {
    super(name);
  }
  public static Test suite() {
    return new TestSuite(WarningTypeTest.class);
  }
  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }

  public void testEnum() {
    TestUtils.testPseudoEnum(WarningType.class, ""
        + "[BAD_TIDE, NUMBER_OF_BAD_TIDE, NO_NAME, SEQUENCE_TOO_LONG, DIRECTORY_DELETION_FAILED, "
        + "BAD_PATH, MEMORY_ALERT, EMPTY_SEQUENCE_INPUT, NO_SEQUENCE, SEQUENCE_LABEL_TOO_LONG, SEQUENCE_LABEL_MISMATCH, COLORSPACE_WARNING, "
        + "ENVIRONMENT_MAY_NOT_BE_USEABLE, INCORRECT_LENGTH, SCORE_NOT_SNP, TOPN_NOT_TOPEQUAL, EXCLUDED_SEQUENCES, "
        + "MASK_FLAGS_SET, FILE_DOES_NOT_EXIST, POSSIBLY_NOT_PROTEIN, NUMBER_OF_INCORRECT_LENGTH, NOT_FASTA_FILE, "
        + "POSSIBLY_SOLEXA, FILE_CORRUPTION, BAD_CHAR_WARNINGS, OVERFLOW_RECORDS, PARAMS_WARNING, "
        + "SAM_IGNORED_RECORDS_CIGAR, SAM_IGNORED_RECORDS, SAM_INCOMPATIBLE_HEADERS, SAM_BAD_FORMAT_WARNING1, SAM_BAD_FORMAT_WARNING, "
        + "CIGAR_BAD_FORMAT_WARNING, INFO_WARNING, NOT_FASTQ_FILE]"
    );
    assertEquals(3, WarningType.BAD_TIDE.getNumberOfParameters());
    assertEquals(1, WarningType.NUMBER_OF_BAD_TIDE.getNumberOfParameters());
    assertEquals(1, WarningType.NO_NAME.getNumberOfParameters());
    assertEquals(1, WarningType.SEQUENCE_TOO_LONG.getNumberOfParameters());
    assertEquals(1, WarningType.DIRECTORY_DELETION_FAILED.getNumberOfParameters());
    assertEquals(1, WarningType.BAD_PATH.getNumberOfParameters());
    assertEquals(0, WarningType.MEMORY_ALERT.getNumberOfParameters());
    assertEquals(1, WarningType.NO_SEQUENCE.getNumberOfParameters());
    assertEquals(0, WarningType.EMPTY_SEQUENCE_INPUT.getNumberOfParameters());
    assertEquals(1, WarningType.SEQUENCE_LABEL_TOO_LONG.getNumberOfParameters());
    assertEquals(2, WarningType.SEQUENCE_LABEL_MISMATCH.getNumberOfParameters());
    assertEquals(1, WarningType.COLORSPACE_WARNING.getNumberOfParameters());
    assertEquals(0, WarningType.ENVIRONMENT_MAY_NOT_BE_USEABLE.getNumberOfParameters());
    assertEquals(3, WarningType.INCORRECT_LENGTH.getNumberOfParameters());
    assertEquals(1, WarningType.SCORE_NOT_SNP.getNumberOfParameters());
    assertEquals(1, WarningType.TOPN_NOT_TOPEQUAL.getNumberOfParameters());
    assertEquals(1, WarningType.EXCLUDED_SEQUENCES.getNumberOfParameters());
    assertEquals(0, WarningType.MASK_FLAGS_SET.getNumberOfParameters());
    assertEquals(1, WarningType.FILE_DOES_NOT_EXIST.getNumberOfParameters());
    assertEquals(2, WarningType.POSSIBLY_NOT_PROTEIN.getNumberOfParameters());
    assertEquals(1, WarningType.NUMBER_OF_INCORRECT_LENGTH.getNumberOfParameters());
    assertEquals(1, WarningType.NOT_FASTA_FILE.getNumberOfParameters());
    assertEquals(0, WarningType.POSSIBLY_SOLEXA.getNumberOfParameters());
    assertEquals(1, WarningType.FILE_CORRUPTION.getNumberOfParameters());
    assertEquals(1, WarningType.BAD_CHAR_WARNINGS.getNumberOfParameters());
    assertEquals(1, WarningType.OVERFLOW_RECORDS.getNumberOfParameters());
    assertEquals(1, WarningType.PARAMS_WARNING.getNumberOfParameters());
    assertEquals(1, WarningType.SAM_IGNORED_RECORDS_CIGAR.getNumberOfParameters());
    assertEquals(2, WarningType.SAM_IGNORED_RECORDS.getNumberOfParameters());
    assertEquals(2, WarningType.SAM_INCOMPATIBLE_HEADERS.getNumberOfParameters());
    assertEquals(1, WarningType.SAM_BAD_FORMAT_WARNING1.getNumberOfParameters());
    assertEquals(2, WarningType.SAM_BAD_FORMAT_WARNING.getNumberOfParameters());
    assertEquals(1, WarningType.CIGAR_BAD_FORMAT_WARNING.getNumberOfParameters());
    assertEquals(1, WarningType.INFO_WARNING.getNumberOfParameters());
    assertEquals(1, WarningType.NOT_FASTQ_FILE.getNumberOfParameters());
  }

  public void testReadResolve() throws ObjectStreamException {
    for (WarningType t : WarningType.values()) {
      assertEquals(t, t.readResolve());
    }
  }

  public void testPrefix() {
    assertEquals("", WarningType.INFO_WARNING.getMessagePrefix());
  }
}

