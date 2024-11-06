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
package com.rtg.variant.bayes.multisample;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import com.rtg.sam.CigarFormatter;
import com.rtg.variant.CalibratedMachineErrorChooser;
import com.rtg.variant.DefaultMachineErrorChooser;
import com.rtg.variant.MachineErrorChooserInterface;
import com.rtg.variant.ReadGroupMachineErrorChooser;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantParams;
import com.rtg.variant.match.AlignmentMatch;

/**
 */
public final class MultisampleUtils {

  private MultisampleUtils() { }

  /**
   * @param params command line parameters for variant.
   * @return global machine error chooser.
   * @throws IOException when reading calibration files.
   */
  public static MachineErrorChooserInterface chooser(final VariantParams params) throws IOException {
    final String errorName = params.machineErrorName();
    if (errorName == null) {
      if (params.calibrator() != null) {
        return new CalibratedMachineErrorChooser(params.calibrator());
      } else {
        return new ReadGroupMachineErrorChooser(params.uberHeader());
      }
    } else {
      return new DefaultMachineErrorChooser(errorName);
    }
  }

  /**
   * Obtain matches for use in complex region hypothesis creation
   * @param start 0 based inclusive start position of complex region
   * @param end 0 based exclusive end position of complex region
   * @param records records covering complex region
   * @param chooser machine error chooser
   * @param params parameters
   * @return match objects for region, or null if an overflow record was encountered
   */
  public static List<AlignmentMatch> intersectSet(final int start, final int end, final Iterator<VariantAlignmentRecord> records, final MachineErrorChooserInterface chooser, final VariantParams params) {
    final ArrayList<AlignmentMatch> ret = new ArrayList<>();
    while (records.hasNext()) {
      final VariantAlignmentRecord rec = records.next();
      //TODO rec.endPosition correct? may need a toMatch
      if (rec.getStart() + rec.getLength() >= start && rec.getStart() <= end) {
        if (rec.isOverflow()) {
          //System.err.println("Intersect encountered overflow record, aborting");
          return null;
        }

        //System.err.println("start=" + start + " end=" + end + " rec.start=");
        final AlignmentMatch match = CigarFormatter.cigarSubsequence(rec, chooser, start, end, params);
        if (match != null) {
          ret.add(match);
        //} else {
        //  System.err.println("In Range, No Match [" + rec.getStart() + "-" + rec.getEnd() + ") " + rec.getCigar());
        }
      //} else {
      //  System.err.println("Out Of Range [" + rec.getStart() + "-" + rec.getEnd() + ") " + rec.getCigar());
      }
    }
    //System.err.println(ret);
    return ret;
  }


}
