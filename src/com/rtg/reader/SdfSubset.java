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
package com.rtg.reader;

import static com.rtg.reader.Sdf2Fasta.END_SEQUENCE;
import static com.rtg.reader.Sdf2Fasta.ID_FILE_FLAG;
import static com.rtg.reader.Sdf2Fasta.INPUT;
import static com.rtg.reader.Sdf2Fasta.NAMES_FLAG;
import static com.rtg.reader.Sdf2Fasta.START_SEQUENCE;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Collection;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.LoggedCli;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.io.LogStream;

/**
 * Pulls out a subset of sequences from one SDF into another SDF.
 *
 */
public final class SdfSubset extends LoggedCli {

  static final String OUTPUT = "output";

  private static final Validator VALIDATOR = new Validator() {

    @Override
    public boolean isValid(CFlags flags) {
      if (!CommonFlags.validateOutputDirectory((File) flags.getValue(OUTPUT))) {
        return false;
      }
      if (!flags.isSet(ID_FILE_FLAG) && !flags.getAnonymousFlag(0).isSet() && !(flags.isSet(START_SEQUENCE) || flags.isSet(END_SEQUENCE))) {
        flags.setParseMessage("Sequences to extract must be specified, either explicitly, or using --" + ID_FILE_FLAG + ", or via --" + START_SEQUENCE + "/--" + END_SEQUENCE);
        return false;
      }
      return Sdf2Fasta.validateExtractorFlags(flags);
    }
  };

  @Override
  protected void initFlags() {
    CommonFlagCategories.setCategories(mFlags);
    mFlags.setDescription("Extracts a subset of sequences from one SDF and outputs them to another SDF.");
    Sdf2Fasta.registerExtractorFlags(mFlags);
    mFlags.registerRequired('o', OUTPUT, File.class, "SDF", "output SDF").setCategory(INPUT_OUTPUT);

    mFlags.setValidator(VALIDATOR);
  }

  @Override
  public String moduleName() {
    return "sdfsubset";
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(OUTPUT);
  }


  @Override
  protected int mainExec(OutputStream out, LogStream log) throws IOException {
    try (final SdfReaderWrapper reader = new SdfReaderWrapper((File) mFlags.getValue(INPUT), false, false)) {
      try (final SdfWriterWrapper writer = new SdfWriterWrapper((File) mFlags.getValue(OUTPUT), reader, false)) {
        writer.copySourceTemplatesFile(reader);

        final WrapperFilter filter;
        if (mFlags.isSet(NAMES_FLAG)) {
          filter = new NameWrapperFilter(reader, writer);
        } else {
          filter = new WrapperFilter(reader, writer);
        }

        if (mFlags.getAnonymousFlag(0).isSet()) {
          final Collection<Object> seqs = mFlags.getAnonymousValues(0);
          for (final Object oi : seqs) {
            filter.transfer((String) oi);
          }
        }
        if (mFlags.isSet(ID_FILE_FLAG)) {
          filter.transferFromFile((File) mFlags.getValue(ID_FILE_FLAG));
        }
        if (mFlags.isSet(START_SEQUENCE) || mFlags.isSet(END_SEQUENCE)) {
          final long startId = mFlags.isSet(START_SEQUENCE) ? (Long) mFlags.getValue(START_SEQUENCE) : LongRange.MISSING;
          final long endId = mFlags.isSet(END_SEQUENCE) ? (Long) mFlags.getValue(END_SEQUENCE) : LongRange.MISSING;
          filter.transfer(new LongRange(startId, endId));
        }

        Diagnostic.progress("Extracted " + filter.getWritten() + " sequences");
      }
    }
    return 0;
  }
}
