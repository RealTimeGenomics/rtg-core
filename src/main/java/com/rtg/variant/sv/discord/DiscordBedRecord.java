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
