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

package com.rtg.visualization;

import com.rtg.sam.BadSuperCigarException;

import htsjdk.samtools.SAMRecord;

/**
 * Uses a SAM record and unfolds it into strings that can be displayed.
 */
public interface SamAssistance {

  /**
   * Use information from SAM record (including read, cigar and GS/GC fields) and <code>templateInserts</code> to construct one or more strings that can be
   * inserted into a display.
   * @param sam SAM format record.
   * @param templateStr the template with inserts already applied (this is in screen co-ordinates and is a subset of the original with modifications added).
   * @param templateBytes the template as a byte array (the original with no modifications applied)
   * @param readStart where the read starts on the template (in screen co-ordinates).
   * @param displayDots if lines match display dot
   * @param displaySoftClip if soft clipped bases should be displayed
   * @return array of strings each string is padded with spaces at start (but not end). There will usually be one string or 2 in case of CG reads.
   * @throws BadSuperCigarException on bad cigar
   */
  String[] samToReads(SAMRecord sam, String templateStr, byte[] templateBytes, int readStart, boolean displayDots, boolean displaySoftClip) throws BadSuperCigarException;
}
