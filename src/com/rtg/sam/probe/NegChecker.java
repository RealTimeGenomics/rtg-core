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
package com.rtg.sam.probe;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import com.rtg.sam.SamUtils;
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
  public boolean check(SAMRecord record, RangeList.RangeData<?> data) {
    if (record.getReadNegativeStrandFlag()) {
      final int alignmentEnd = record.getAlignmentEnd();
      if (data.getEnd() > alignmentEnd - mTolerance && data.getEnd() < alignmentEnd + mTolerance) {
//                    System.err.println(record.getSAMString() + "strip back to: " + data.getStart() + " (" + data.getStart() + " : " + data.getEnd() + ")");
        return true;
      }
    }
    return false;
  }

  @Override
  public int getStartDataIndex(SAMRecord record, RangeList<?> list) {
    return list.findFullRangeIndex(record.getAlignmentEnd() - mTolerance);
  }

  @Override
  void stripRecord(SAMRecord record, SAMRecord mate, RangeList.RangeData<?> data) {
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
    record.setAttribute("XP", new String(readBases, readEnd, trimmed));
    record.setAttribute(SamUtils.ATTRIBUTE_NUM_MISMATCHES, null);
    mBasesTrimmed += trimmed;
  }
}
