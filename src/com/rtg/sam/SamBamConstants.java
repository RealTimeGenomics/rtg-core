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

import com.reeltwo.jumble.annotations.TestClass;

/**
 * Class containing useful constants
 */
@TestClass(value = {"com.rtg.sam.BamReaderTest"})
public final class SamBamConstants {

  private SamBamConstants() {
  }

  /** These flags are defined in the SAM Tools specification. */
  public static final int
      SAM_READ_IS_PAIRED = 0x0001,
      SAM_READ_IS_MAPPED_IN_PROPER_PAIR = 0x0002,
      SAM_READ_IS_UNMAPPED = 0x0004,
      SAM_MATE_IS_UNMAPPED = 0x0008,
      SAM_READ_IS_REVERSE = 0x0010,
      SAM_MATE_IS_REVERSE = 0x0020,
      SAM_READ_IS_FIRST_IN_PAIR = 0x0040,
      SAM_READ_IS_SECOND_IN_PAIR = 0x0080,
      SAM_SECONDARY_ALIGNMENT = 0x0100,
      SAM_PCR_OR_OPTICAL_DUPLICATE = 0x0400;

  /**
   * The fixed field positions in SAM records.
   * We use these constants to access fields in BAM records too.
   */
  public static final int
      QNAME_FIELD = 0, // Query pair NAME if paired; or Query NAME if unpaired
      FLAG_FIELD = 1,  // bitwise FLAG
      RNAME_FIELD = 2, // Reference sequence NAME
      POS_FIELD = 3,   // 1-based leftmost POSition/coordinate of the clipped sequence
      MAPQ_FIELD = 4,  // MAPping Quality - phred-scaled posterior probability that position is incorrect
      CIGAR_FIELD = 5, // extended CIGAR string
      MRNM_FIELD = 6,  // Mate Reference sequence NaMe
      MPOS_FIELD = 7,  // 1-based leftmost Mate POSition of the clipped sequence
      ISIZE_FIELD = 8, // inferred Insert SIZE 5
      SEQ_FIELD = 9,   // query SEQuence; "=" for a match to the reference; n/N/. for ambiguity
      QUAL_FIELD = 10, // query QUALity; ASCII-33 gives the Phred base quality
      ATTRIBUTES_FIELD = 11; // optional TAG:TYPE:VALUE attribute fields start here.


}
