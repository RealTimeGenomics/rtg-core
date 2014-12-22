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
package com.rtg.sam;

import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Locale;

import com.rtg.ngs.MapFlags;
import com.rtg.reader.FormatCli;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMReadGroupRecord;

/**
 * Contains methods for helping with command line flags and params related to SAM
 */
public final class SamCommandHelper {

  private SamCommandHelper() {
    super();
  }

  /** flag to capture SAM read group information */
  public static final String SAM_RG = "sam-rg";
  static final String STRING_OR_FILE = "STRING|FILE";

  /**
   * Method adds SAM read group flag to given input
   * @param flags flags to add <code>sam-rg</code> flag.
   */
  public static void initSamRg(CFlags flags) {
    initSamRg(flags, "ILLUMINA", CommonFlagCategories.REPORTING);
  }

  /**
   * Method adds SAM read group flag to given input
   * @param flags flags to add <code>sam-rg</code> flag.
   * @param examplePlatform the example platform to use.
   * @param category the category to register the flag under.
   */
  public static void initSamRg(CFlags flags, String examplePlatform, String category) {
    flags.registerOptional(SAM_RG, String.class, STRING_OR_FILE, "file containing a single valid read group SAM header line or a string in the form \"@RG\\tID:READGROUP1\\tSM:BACT_SAMPLE\\tPL:" + examplePlatform + "\"").setCategory(category);
  }

  /**
   * Validate SAM read group flag
   * @param flags flags to validate from
   * @return true if all OK, false otherwise
   */
  public static boolean validateSamRg(CFlags flags) {
    final String rg = (String) flags.getValue(SAM_RG);
    final File in = new File(rg);
    if (!in.exists() && rg.indexOf('\t') != -1) {
      flags.setParseMessage("given string \"" + rg + "\" for --" + SAM_RG + " contains literal tab characters, please use \\t instead");
      return false;
    } else if (in.isDirectory()) {
      flags.setParseMessage("given input file \"" + in.getPath() + "\" for --" + SAM_RG + " is a directory, must be a file");
      return false;
    }
    return true;
  }

  /**
   * Validate and create SAM read group record.
   * @param value the name of the file or the read group string to be parsed.
   * @param strict if true, throw exceptions when read group information is Invalid, if false return null instead.
   * @return the read group record or null when not in strict mode and input is invalid.
   * @throws InvalidParamsException if there is no read group or more than one read groups in the input.
   * @throws java.io.IOException if there is an IO problem.
   */
  public static SAMReadGroupRecord validateAndCreateSamRG(String value, ReadGroupStrictness strict) throws InvalidParamsException, IOException {
    final File rgFile = new File(value);
    final boolean fileMode = rgFile.exists();
    final BufferedInputStream bis;
    final String errorSubString;
    if (fileMode) {
      errorSubString = "file \"" + rgFile.getPath() + "\", please provide file";
      bis = FileUtils.createInputStream(rgFile, false);
    } else {
      errorSubString = "string \"" + value + "\", please provide string";
      final String convertedValue = value.replaceAll("\\\\t", "\t");
      bis = new BufferedInputStream(new ByteArrayInputStream(convertedValue.getBytes()));
    }
    try {
      final SAMFileReader sfr = new SAMFileReader(bis);
      final List<SAMReadGroupRecord> readGroups = sfr.getFileHeader().getReadGroups();
      final int readGroupCount = readGroups.size();
      if (readGroupCount == 0) {
        if (strict == ReadGroupStrictness.REQUIRED) {
            throw new InvalidParamsException("No read group information present in the input " + errorSubString + " with single read group line");
        } else {
          return null;
        }
      }
      if (readGroupCount > 1) {
        if (strict == ReadGroupStrictness.REQUIRED || strict == ReadGroupStrictness.AT_MOST_ONE) {
          throw new InvalidParamsException("Multiple read group information present in the input " + errorSubString + " with single read group line");
        } else {
          return null;
        }
      }
      final SAMReadGroupRecord samReadGroupRecord = readGroups.get(0);
      if (samReadGroupRecord.getSample() == null) {
        Diagnostic.warning("Sample not specified in read group, it is recommended to set the sample tag.");
      }
      return samReadGroupRecord;
    } finally {
      bis.close();
    }
  }

  /**
   * @param flags flags to check
   * @return true if input is sam format
   */
  public static boolean isSamInput(CFlags flags) {
    final String format;
    if (flags.isSet(FormatCli.FORMAT_FLAG)) {
      format = flags.getValue(FormatCli.FORMAT_FLAG).toString().toLowerCase(Locale.getDefault());
    } else {
      format = MapFlags.SDF_FORMAT;
    }
    return MapFlags.SAM_PE_FORMAT.equals(format) || MapFlags.SAM_SE_FORMAT.equals(format);
  }

  /**
   * This validates an retrieves from a sam file the read group matching the {@code selectReadGroup} parameter
   * @param rgFile sam file containing the read group
   * @param selectReadGroup the read group ID to locate
   * @return a {@code SAM ReadGroupRecord} that corresponds to the requested id
   * @throws java.io.IOException if IO falls over when reading the SAM file
   */
  public static SAMReadGroupRecord validateSelectedSamRG(File rgFile, String selectReadGroup) throws IOException {
    final BufferedInputStream bis;
    bis = FileUtils.createInputStream(rgFile, false);
    try {
      final SAMFileReader sfr = new SAMFileReader(bis);
      final List<SAMReadGroupRecord> readGroups = sfr.getFileHeader().getReadGroups();
      final int readGroupCount = readGroups.size();
      if (readGroupCount == 0) {
        throw new InvalidParamsException("No read group information matching \"" + selectReadGroup + "\" present in the input file \"" + rgFile.getPath() + "\"");
      }
      for (SAMReadGroupRecord r : readGroups) {
        if (selectReadGroup.equals(r.getId())) {
          if (r.getSample() == null) {
            Diagnostic.warning("Sample not specified in read group, it is recommended to set the sample tag.");
          }
          return r;
        }
      }
      throw new InvalidParamsException("No read group information matching \"" + selectReadGroup + "\" present in the input file \"" + rgFile.getPath() + "\"");
    } finally {
      bis.close();
    }

  }

  /**
   * How strict would you like read group validation to be
   */
  public enum ReadGroupStrictness {
    /** A single read group must be present */
    REQUIRED
    /** Either a single read group will be present or null will be returned */
    , OPTIONAL
    /** Null if none, the read group if single, an exception on more than one */
    , AT_MOST_ONE
  }
}
