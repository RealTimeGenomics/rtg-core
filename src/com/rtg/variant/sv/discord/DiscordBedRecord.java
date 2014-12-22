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

package com.rtg.variant.sv.discord;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.bed.BedRecord;

/**
 * Encapsulates a BED format record with filtered status.
 */
@TestClass("com.rtg.variant.sv.discord.SmartBedWriterTest")
public class DiscordBedRecord extends BedRecord {

  private boolean mIsFiltered = false;

  /**
   * Create a new BED record.
   * @param sequence the name of the reference sequence
   * @param start the zero-based start position
   * @param end the zero-based end position, exclusive
   * @param annotations extra string valued columns
   */
  public DiscordBedRecord(String sequence, int start, int end, String... annotations) {
    super(sequence, start, end, annotations);
  }

  /**
   * Sets this bed record as being filtered
   * @return this, for call chaining
   */
  public DiscordBedRecord setFiltered() {
    mIsFiltered = true;
    return this;
  }

  /**
   * @return true if this bed record is filtered
   */
  public boolean isFiltered() {
    return mIsFiltered;
  }

}
