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
import java.util.Collections;
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
class NegChecker extends PositionAndStrandChecker {

  NegChecker(int tolerance) {
    super(tolerance);
  }

  @Override
  public boolean checkStrand(SAMRecord record) {
    return record.getReadNegativeStrandFlag();
  }

  @Override
  public boolean checkPosition(SAMRecord record, Interval data) {
    assert checkStrand(record);
    final int alignmentEnd = record.getAlignmentEnd();
    if (data.getEnd() >= alignmentEnd - mTolerance && data.getEnd() <= alignmentEnd + mTolerance) {
//                    System.err.println(record.getSAMString() + "strip back to: " + data.getStart() + " (" + data.getStart() + " : " + data.getEnd() + ")");
      return true;
    }
    return false;
  }

  @Override
  public int getStartDataIndex(SAMRecord record, RangeList<?> list) {
    return list.findFullRangeIndex(record.getAlignmentEnd() - 1 - mTolerance);
  }

  @Override
  void stripRecord(SAMRecord record, SAMRecord mate, Interval data) {
    final int diff = record.getAlignmentEnd() - data.getEnd();
    mPosDiffStats[mTolerance + diff]++;
    setAlignmentEnd(record, mate, data.getStart());
  }

  void setAlignmentEnd(SAMRecord record, SAMRecord mate, int alignmentEnd) {
    //end positions are 1 based exclusive
    int readEnd = record.getReadLength();
    int refEnd = record.getAlignmentEnd();
    final List<CigarElement> newCigarElements = new ArrayList<>();
    final List<CigarElement> cigarElements = record.getCigar().getCigarElements();
    for (int i = cigarElements.size() - 1; i >= 0; --i) {
      final CigarElement e = cigarElements.get(i);
      final CigarOperator operator = e.getOperator();
      if (alignmentEnd < refEnd) {
        final int consume = operator.consumesReferenceBases() ? Math.min(refEnd - alignmentEnd, e.getLength()) : e.getLength();
        if (operator.consumesReferenceBases()) {
          refEnd -= consume;
        }
        if (operator.consumesReadBases()) {
          readEnd -= consume;
        }
        updateStrippedStats(operator, consume);
        if (e.getLength() - consume > 0) {
          newCigarElements.add(new CigarElement(e.getLength() - consume, operator));
        }
      } else {
        newCigarElements.add(e);
      }
    }
    Collections.reverse(newCigarElements);
    final byte[] readBases = record.getReadBases();
    record.setReadBases(Arrays.copyOfRange(readBases, 0, readEnd));
    final int trimmed = readBases.length - readEnd;
    final byte[] baseQualities = record.getBaseQualities();
    record.setBaseQualities(Arrays.copyOfRange(baseQualities, 0, readEnd));
    record.setCigar(new Cigar(newCigarElements));
    if (mate != null) {
      updateTlenAndMateStart(record, mate);
    }
    record.setAttribute(DeProbeCli.ATTRIBUTE_STRIPPED_PROBE, new String(readBases, readEnd, trimmed));
    record.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, null);
    mBasesTrimmed += trimmed;
  }
}
