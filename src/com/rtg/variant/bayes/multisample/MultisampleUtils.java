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
        return new ReadGroupMachineErrorChooser(params.mapped());
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
