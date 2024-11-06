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
package com.rtg.sam.probe;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.rtg.sam.SamUtils;
import com.rtg.util.intervals.Interval;
import com.rtg.util.intervals.RangeList;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

/**
 *
 */
class PosChecker extends PositionAndStrandChecker {
  PosChecker(int tolerance) {
    super(tolerance);
  }

  @Override
  public boolean checkStrand(SAMRecord record) {
    return !record.getReadNegativeStrandFlag();
  }

  @Override
  public boolean checkPosition(SAMRecord record, Interval data) {
    assert checkStrand(record);
    final int alignmentStart = record.getAlignmentStart() - 1;
    if (data.getStart() >= alignmentStart - mTolerance && data.getStart() <= alignmentStart + mTolerance) {
//                    System.err.println(record.getSAMString() + " strip forward to: " + data.getEnd() + " (" + data.getStart() + " : " + data.getEnd() + ")");
      return true;
    }
    return false;
  }

  @Override
  public int getStartDataIndex(SAMRecord record, RangeList<?> list) {
    final int alignmentStart = record.getAlignmentStart() - 1;
    return list.findFullRangeIndex(alignmentStart - mTolerance);
  }

  @Override
  void stripRecord(SAMRecord record, SAMRecord mate, Interval data) {
    final int diff = record.getAlignmentStart() - 1 - data.getStart();
    mPosDiffStats[mTolerance + diff]++;
    setAlignmentStart(record, mate, data.getEnd());
  }

  void setAlignmentStart(SAMRecord record, SAMRecord mate, int alignmentStart) {
    int readStart = 0;
    int refStart = record.getAlignmentStart() - 1;
    final List<CigarElement> cigarElements = new ArrayList<>();
    for (CigarElement e : record.getCigar().getCigarElements()) {
      final CigarOperator operator = e.getOperator();
      if (alignmentStart > refStart) {
        final int consume = operator.consumesReferenceBases() ? Math.min(alignmentStart - refStart, e.getLength()) : e.getLength();
        if (operator.consumesReferenceBases()) {
          refStart += consume;
        }
        if (operator.consumesReadBases()) {
          readStart += consume;
        }
        updateStrippedStats(operator, consume);
        if (e.getLength() - consume > 0) {
          cigarElements.add(new CigarElement(e.getLength() - consume, operator));
        }
      } else {
        cigarElements.add(e);
      }
    }
    final byte[] readBases = record.getReadBases();
    record.setReadBases(Arrays.copyOfRange(readBases, readStart, readBases.length));
    final byte[] baseQualities = record.getBaseQualities();
    record.setBaseQualities(Arrays.copyOfRange(baseQualities, readStart, readBases.length));
    record.setCigar(new Cigar(cigarElements));
    record.setAlignmentStart(alignmentStart + 1);
    if (mate != null) {
      updateTlenAndMateStart(record, mate);
    }
    record.setAttribute(DeProbeCli.ATTRIBUTE_STRIPPED_PROBE, new String(readBases, 0, readStart));
    record.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, null);
    mBasesTrimmed += readStart;
  }

}
