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
package com.rtg.index.hash.ngs;

import java.io.IOException;

/**
 * Common functions to both <code>ReadHashFunction</code> and <code>TemplateHashFunction</code>.
 */
public interface ReadHashFunction extends CommonHashFunction {


  /**
   * Force creation of an array to hold the read bit sequences from the <code>HashLoop</code> supplied to the <code>HashFunction</code>
   * @param numberReads number of read sequences.
   */
  void setReadSequences(long numberReads);

  /**
   * Process all windows for the specified read.
   * @param readId number of the read (&gt;=0).
   * @param reverse if true store reverse complement.
   * @throws IOException If an I/O error occurs
   */
  void readAll(final int readId, final boolean reverse) throws IOException;

  /**
   * Store the current values in <code>readSequences</code>.
   * @param id2 read id.
   * @param reverse if true use reverse complement.
   */
  void setValues(int id2, boolean reverse);
}

