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
package com.rtg.simulation.reads;

import java.io.File;
import java.util.Locale;

import com.rtg.launcher.CommonFlags;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.Validator;
import com.rtg.util.machine.MachineType;

/**
 * Flag validator for the read simulator.
 */
class ReadSimCliValidator implements Validator {

  private static boolean okProbability(final Object o) {
    if (!(o instanceof Double)) {
      return false;
    }
    final double d = (Double) o;
    return d >= 0 && d <= 1;
  }

  @Override
  public boolean isValid(final CFlags cflags) {

    if (cflags.isSet(ReadSimCli.MNP_EVENT_RATE) && !okProbability(cflags.getValue(ReadSimCli.MNP_EVENT_RATE))) {
      cflags.setParseMessage("--" + ReadSimCli.MNP_EVENT_RATE + " must be a value in [0,1]");
      return false;
    }
    if (cflags.isSet(ReadSimCli.INS_EVENT_RATE) && !okProbability(cflags.getValue(ReadSimCli.INS_EVENT_RATE))) {
      cflags.setParseMessage("--" + ReadSimCli.INS_EVENT_RATE + " must be a value in [0,1]");
      return false;
    }
    if (cflags.isSet(ReadSimCli.DEL_EVENT_RATE) && !okProbability(cflags.getValue(ReadSimCli.DEL_EVENT_RATE))) {
      cflags.setParseMessage("--" + ReadSimCli.DEL_EVENT_RATE + " must be a value in [0,1]");
      return false;
    }

    if (!(cflags.isSet(ReadSimCli.COVERAGE) ^ cflags.isSet(ReadSimCli.READS))) {
      cflags.setParseMessage("Exactly one of --" + ReadSimCli.COVERAGE + " or --" + ReadSimCli.READS + " must be specified");
      return false;
    }
    if (cflags.isSet(ReadSimCli.READS)) {
      if ((Long) cflags.getValue(ReadSimCli.READS) <= 0) {
        cflags.setParseMessage("Number of reads should be greater than 0");
        return false;
      }
    } else if (cflags.isSet(ReadSimCli.COVERAGE)) {
      final double coverage = (Double) cflags.getValue(ReadSimCli.COVERAGE);
      if (coverage <= 0.0) {
        cflags.setParseMessage("Coverage should be positive");
        return false;
      } else if (coverage > 1000000.0) {
        cflags.setParseMessage("Coverage cannot be greater than 1000000.0");
        return false;
      }
    }
    if (!CommonFlags.validateOutputDirectory((File) cflags.getValue(ReadSimCli.OUTPUT))) {
      return false;
    }
    if (cflags.isSet(ReadSimCli.DISTRIBUTION) && cflags.isSet(ReadSimCli.TAXONOMY_DISTRIBUTION)) {
      cflags.setParseMessage("At most one of --" + ReadSimCli.DISTRIBUTION + " or --" + ReadSimCli.TAXONOMY_DISTRIBUTION + " should be set");
      return false;
    }
    if (cflags.isSet(ReadSimCli.DISTRIBUTION)) {
      final File distFile = (File) cflags.getValue(ReadSimCli.DISTRIBUTION);
      if (!distFile.exists() || distFile.isDirectory()) {
        cflags.setParseMessage("File: " + distFile + " does not exist");
        return false;
      }
    }
    if (!cflags.isSet(ReadSimCli.TAXONOMY_DISTRIBUTION) && (cflags.isSet(ReadSimCli.ABUNDANCE) || cflags.isSet(ReadSimCli.DNA_FRACTION))) {
      cflags.setParseMessage("--" + ReadSimCli.ABUNDANCE + " and --" + ReadSimCli.DNA_FRACTION + " are only applicable if using --" + ReadSimCli.TAXONOMY_DISTRIBUTION);
      return false;
    }
    if (cflags.isSet(ReadSimCli.TAXONOMY_DISTRIBUTION) && !(cflags.isSet(ReadSimCli.ABUNDANCE) || cflags.isSet(ReadSimCli.DNA_FRACTION))) {
      cflags.setParseMessage("either --" + ReadSimCli.ABUNDANCE + " or --" + ReadSimCli.DNA_FRACTION + " must be set if using --" + ReadSimCli.TAXONOMY_DISTRIBUTION);
      return false;

    }
    if (cflags.isSet(ReadSimCli.ABUNDANCE) && cflags.isSet(ReadSimCli.DNA_FRACTION)) {
      cflags.setParseMessage("At most one of --" + ReadSimCli.ABUNDANCE + " or --" + ReadSimCli.DNA_FRACTION + " should be set");
      return false;
    }
    if (cflags.isSet(ReadSimCli.TAXONOMY_DISTRIBUTION)) {
      final File distFile = (File) cflags.getValue(ReadSimCli.TAXONOMY_DISTRIBUTION);
      if (!distFile.exists() || distFile.isDirectory()) {
        cflags.setParseMessage("File: " + distFile + " does not exist");
        return false;
      }
    }

    if (cflags.isSet(ReadSimCli.QUAL_RANGE)) {
      final String range = (String) cflags.getValue(ReadSimCli.QUAL_RANGE);
      final String[] vals = range.split("-");
      if (vals.length != 2) {
        cflags.setParseMessage("Quality range is not of form qualmin-qualmax");
        return false;
      } else {
        try {
          final int l = Integer.parseInt(vals[0]);
          if ((l < 0) || (l > 63)) {
            cflags.setParseMessage("Minimum quality value must be between 0 and 63");
            return false;
          }
          final int u = Integer.parseInt(vals[1]);
          if ((u < 0) || (u > 63)) {
            cflags.setParseMessage("Maximum quality value must be between 0 and 63");
            return false;
          }
          if (l > u) {
            cflags.setParseMessage("Minimum quality value cannot be greater than maximum quality value");
            return false;
          }
        } catch (final NumberFormatException e) {
          cflags.setParseMessage("Quality range is not of form qualmin-qualmax");
          return false;
        }
      }
    }
    final File f = (File) cflags.getValue(ReadSimCli.INPUT);
    if (!f.exists()) {
      cflags.setParseMessage("The specified SDF, \"" + f.getPath() + "\", does not exist.");
      return false;
    }
    if (!f.isDirectory()) {
      cflags.setParseMessage("The specified file, \"" + f.getPath() + "\", is not an SDF.");
      return false;
    }
    if (cflags.isSet(ReadSimCli.TWIN_INPUT)) {
      final File tf = (File) cflags.getValue(ReadSimCli.TWIN_INPUT);
      if (!tf.exists()) {
        cflags.setParseMessage("The specified SDF, \"" + tf.getPath() + "\", does not exist.");
        return false;
      }
      if (!tf.isDirectory()) {
        cflags.setParseMessage("The specified file, \"" + tf.getPath() + "\", is not an SDF.");
        return false;
      }
      if (tf.equals(f)) {
        cflags.setParseMessage("The --" + ReadSimCli.TWIN_INPUT + " SDF cannot be the same as that given with --" + ReadSimCli.INPUT);
        return false;
      }
    }
    if ((Integer) cflags.getValue(ReadSimCli.MIN_FRAGMENT) > (Integer) cflags.getValue(ReadSimCli.MAX_FRAGMENT)) {
      cflags.setParseMessage("--" + ReadSimCli.MAX_FRAGMENT + " should not be smaller than --" + ReadSimCli.MIN_FRAGMENT);
      return false;
    }

    final File bed = (File) cflags.getValue(ReadSimCli.BED_FILE);
    if (cflags.isSet(ReadSimCli.BED_FILE)) {
      if (!bed.exists()) {
        cflags.setParseMessage("The --" + ReadSimCli.BED_FILE + " specified file doesn't exist: " + bed.getPath());
        return false;
      }
      if (cflags.isSet(ReadSimCli.COVERAGE)) {
        cflags.setParseMessage("--" + ReadSimCli.BED_FILE + " is incompatible with --" + ReadSimCli.COVERAGE);
        return false;
      }
    }

    return checkMachines(cflags);
  }

  protected boolean checkMachines(CFlags cflags) {
    final MachineType mt = MachineType.valueOf(cflags.getValue(ReadSimCli.MACHINE_TYPE).toString().toLowerCase(Locale.getDefault()));
    if (mt == MachineType.ILLUMINA_SE) {
      if (!cflags.checkRequired(ReadSimCli.READLENGTH)) {
        return false;
      }
      if ((Integer) cflags.getValue(ReadSimCli.READLENGTH) <= 1) {
        cflags.setParseMessage("Read length is too small");
        return false;
      }
      if ((Integer) cflags.getValue(ReadSimCli.READLENGTH) > (Integer) cflags.getValue(ReadSimCli.MIN_FRAGMENT)) {
        cflags.error("Read length is too large for selected fragment size");
        return false;
      }
      if (!cflags.checkBanned(ReadSimCli.LEFT_READLENGTH, ReadSimCli.RIGHT_READLENGTH, ReadSimCli.MIN_TOTAL_454_LENGTH, ReadSimCli.MAX_TOTAL_454_LENGTH, ReadSimCli.MIN_TOTAL_IONTORRENT_LENGTH, ReadSimCli.MAX_TOTAL_IONTORRENT_LENGTH)) {
        return false;
      }
    } else if (mt == MachineType.ILLUMINA_PE) {
      if (!cflags.checkRequired(ReadSimCli.LEFT_READLENGTH, ReadSimCli.RIGHT_READLENGTH)) {
        return false;
      }
      if (((Integer) cflags.getValue(ReadSimCli.LEFT_READLENGTH) <= 1)
        || ((Integer) cflags.getValue(ReadSimCli.RIGHT_READLENGTH) <= 1)) {
        cflags.error("Read length is too small");
        return false;
      }
      if (((Integer) cflags.getValue(ReadSimCli.LEFT_READLENGTH) > (Integer) cflags.getValue(ReadSimCli.MIN_FRAGMENT))
        || ((Integer) cflags.getValue(ReadSimCli.RIGHT_READLENGTH) > (Integer) cflags.getValue(ReadSimCli.MIN_FRAGMENT))) {
        cflags.error("Read length is too large for selected fragment size");
        return false;
      }
      if (!cflags.checkBanned(ReadSimCli.READLENGTH, ReadSimCli.MIN_TOTAL_454_LENGTH, ReadSimCli.MAX_TOTAL_454_LENGTH)) {
        return false;
      }

    } else if (mt == MachineType.COMPLETE_GENOMICS || mt == MachineType.COMPLETE_GENOMICS_2) {
      if (!cflags.checkBanned(ReadSimCli.READLENGTH, ReadSimCli.LEFT_READLENGTH, ReadSimCli.RIGHT_READLENGTH, ReadSimCli.MIN_TOTAL_454_LENGTH, ReadSimCli.MAX_TOTAL_454_LENGTH, ReadSimCli.MIN_TOTAL_IONTORRENT_LENGTH, ReadSimCli.MAX_TOTAL_IONTORRENT_LENGTH)) {
        return false;
      }

    } else if (mt == MachineType.FOURFIVEFOUR_PE || mt == MachineType.FOURFIVEFOUR_SE) {
      if (!cflags.checkRequired(ReadSimCli.MIN_TOTAL_454_LENGTH, ReadSimCli.MAX_TOTAL_454_LENGTH)) {
        return false;
      }
      if (!cflags.checkBanned(ReadSimCli.READLENGTH, ReadSimCli.LEFT_READLENGTH, ReadSimCli.RIGHT_READLENGTH, ReadSimCli.MIN_TOTAL_IONTORRENT_LENGTH, ReadSimCli.MAX_TOTAL_IONTORRENT_LENGTH)) {
        return false;
      }
      if ((Integer) cflags.getValue(ReadSimCli.MAX_TOTAL_454_LENGTH) > (Integer) cflags.getValue(ReadSimCli.MIN_FRAGMENT)) {
        cflags.error("Read length is too large for selected fragment size");
        return false;
      }
    } else if (mt == MachineType.IONTORRENT) {
      if (!cflags.checkRequired(ReadSimCli.MIN_TOTAL_IONTORRENT_LENGTH, ReadSimCli.MAX_TOTAL_IONTORRENT_LENGTH)) {
        return false;
      }
      if (!cflags.checkBanned(ReadSimCli.READLENGTH, ReadSimCli.LEFT_READLENGTH, ReadSimCli.RIGHT_READLENGTH, ReadSimCli.MIN_TOTAL_454_LENGTH, ReadSimCli.MAX_TOTAL_454_LENGTH)) {
        return false;
      }
      if ((Integer) cflags.getValue(ReadSimCli.MAX_TOTAL_IONTORRENT_LENGTH) > (Integer) cflags.getValue(ReadSimCli.MIN_FRAGMENT)) {
        cflags.error("Read length is too large for selected fragment size");
        return false;
      }
    } else {
      throw new IllegalArgumentException("Unhandled machine type: " + mt);
    }
    return true;
  }
}
