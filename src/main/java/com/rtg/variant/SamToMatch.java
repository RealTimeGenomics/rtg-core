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
package com.rtg.variant;

import com.rtg.sam.ReaderRecord;

/**
 */
public interface SamToMatch {

  /**
   * Process a single alignment record by passing along to appropriate matchers (which
   * are likely to be supplied in constructors).
   * @param templateBytes the nucleotides in the template.
   * @param var alignment record to be processed.
   * @return true iff the record is processed without error.
   */
  boolean process(byte[] templateBytes, VariantAlignmentRecord var);

  /**
   * Get the effective start position of the alignment record. That is, the earliest
   * reference nucleotide position that might be affected by the record. This
   * may differ from the start position specified in the alignment record in places
   * like all paths processing.
   * @param var record whose start position is to be extracted.
   * @return start position (0 based).
   */
  int start(ReaderRecord<?> var);

  /**
   * Get the effective end position of the alignment record. That is, the latest
   * reference nucleotide position that might be affected by the record. This
   * may differ from the end position specified in the alignment record in places
   * like all paths processing.
   * @param var record whose end position is to be extracted.
   * @return end position (0 based).
   */
  int end(VariantAlignmentRecord var);

}
