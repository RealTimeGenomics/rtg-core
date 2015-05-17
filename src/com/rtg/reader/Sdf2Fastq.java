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


import static com.rtg.launcher.CommonFlags.NO_GZIP;
import static com.rtg.reader.Sdf2Fasta.END_SEQUENCE;
import static com.rtg.reader.Sdf2Fasta.ID_FILE_FLAG;
import static com.rtg.reader.Sdf2Fasta.INPUT;
import static com.rtg.reader.Sdf2Fasta.LINE_LENGTH;
import static com.rtg.reader.Sdf2Fasta.NAMES_FLAG;
import static com.rtg.reader.Sdf2Fasta.OUTPUT;
import static com.rtg.reader.Sdf2Fasta.RENAME;
import static com.rtg.reader.Sdf2Fasta.START_SEQUENCE;
import static com.rtg.reader.Sdf2Fasta.registerTextExtractorFlags;
import static com.rtg.reader.Sdf2Fasta.validateTextExtractorFlags;
import static com.rtg.util.cli.CommonFlagCategories.UTILITY;
import static com.rtg.util.cli.CommonFlagCategories.setCategories;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Collection;

import com.rtg.launcher.AbstractCli;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.intervals.LongRange;

/**
 * This class takes format directory and converts to FASTQ format
 */
public final class Sdf2Fastq extends AbstractCli {

  static final String MODULE_NAME = "sdf2fastq";

  static final String DEFAULT_QUALITY = "default-quality";

  /**
   * @return current name of the module
   */
  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  @Override
  protected void initFlags() {
    mFlags.registerExtendedHelp();
    mFlags.setDescription("Converts SDF data into FASTQ file(s).");
    setCategories(mFlags);

    registerTextExtractorFlags(mFlags);

    mFlags.registerOptional('q', DEFAULT_QUALITY, Integer.class, "INT", "default quality value to use if the SDF does not contain quality data (0-63)").setCategory(UTILITY);

    mFlags.setValidator(VALIDATOR);
  }

  private static final Validator VALIDATOR = new Validator() {
    @Override
    public boolean isValid(final CFlags flags) {
      if (flags.isSet(DEFAULT_QUALITY)) {
        final int qual = (Integer) flags.getValue(DEFAULT_QUALITY);
        if (qual < 0) {
          Diagnostic.error(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, DEFAULT_QUALITY, Integer.toString(qual), "0");
          return false;
        }
        if (qual > 63) {
          Diagnostic.error(ErrorType.INVALID_MAX_INTEGER_FLAG_VALUE, DEFAULT_QUALITY, Integer.toString(qual), "63");
          return false;
        }
      }

      return validateTextExtractorFlags(flags);
    }
  };


  @Override
  protected int mainExec(final OutputStream out, final PrintStream err) throws IOException {
    final PrintStream outStream = new PrintStream(out);
    try {
      final int lineLength = (Integer) mFlags.getValue(LINE_LENGTH);
      final boolean gzip = !mFlags.isSet(NO_GZIP);
      final boolean rename = mFlags.isSet(RENAME);
      final int def = mFlags.isSet(DEFAULT_QUALITY) ? (Integer) mFlags.getValue(DEFAULT_QUALITY) + (int) '!' : -1;

      try (SdfReaderWrapper reader = new SdfReaderWrapper((File) mFlags.getValue(INPUT), false, false)) {
        try (WriterWrapper writer = new FastqWriterWrapper((File) mFlags.getValue(OUTPUT), reader, lineLength, rename, gzip, def)) {
          final WrapperFilter filter;
          if (mFlags.isSet(NAMES_FLAG)) {
            filter = new NameWrapperFilter(reader, writer);
          } else {
            filter = new WrapperFilter(reader, writer);
          }

          boolean doAll = true;
          if (mFlags.getAnonymousFlag(0).isSet()) {
            doAll = false;
            final Collection<Object> seqs = mFlags.getAnonymousValues(0);
            for (final Object oi : seqs) {
              filter.transfer((String) oi);
            }
          }
          if (mFlags.isSet(ID_FILE_FLAG)) {
            doAll = false;
            try (BufferedReader br = new BufferedReader(new FileReader((File) mFlags.getValue(ID_FILE_FLAG)))) {
              String line;
              while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.length() > 0) {
                  filter.transfer(line);
                }
              }
            }
          }
          if (doAll || mFlags.isSet(START_SEQUENCE) || mFlags.isSet(END_SEQUENCE)) {
            final long startId = mFlags.isSet(START_SEQUENCE) ? (Long) mFlags.getValue(START_SEQUENCE) : LongRange.MISSING;
            final long endId = mFlags.isSet(END_SEQUENCE) ? (Long) mFlags.getValue(END_SEQUENCE) : LongRange.MISSING;
            final LongRange r = SequencesReaderFactory.resolveRange(new LongRange(startId, endId), reader.numberSequences());
            for (long seq = r.getStart(); seq < r.getEnd(); seq++) {
              filter.transfer(seq);
            }
          }
        }
      } catch (final IOException e) {
        // Ignore broken pipe error so we don't die on | head etc.
        if (!e.getMessage().contains("Broken pipe")) {
          throw e;
        }
      } catch (InvalidParamsException e) {
        e.printErrorNoLog();
        return 1;
      }
      return 0;
    } finally {
      outStream.flush();
    }
  }
}
