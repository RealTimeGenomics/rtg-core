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

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.NoSuchElementException;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.WarningType;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFormatException;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMValidationError;
import net.sf.samtools.util.RuntimeIOException;

/**
 * Connects between SamFileReaderAdaptor and <code>RecordIterator</code>, skipping (and warns for) records
 * that are invalid according to Picard. Use this only if you have a single
 * file to process and do not need additional filtering.
 */
public class SkipInvalidRecordsIterator extends AbstractSamRecordIterator {

  private static final int NUM_RECORDS_TO_WARN = 5;

  private final String mPath;
  private final RecordIterator<SAMRecord> mWrapped;
  private final boolean mSilent;

  private boolean mIsClosed = false;
  protected SAMRecord mRecord;

  /**
   * Constructor
   * @param path filename to use for errors and warnings
   * @param reader supplier of <code>SAMRecords</code>
   * @param silent true to not report warnings
   * @throws IOException if an IO error occurs
   */
  public SkipInvalidRecordsIterator(String path, RecordIterator<SAMRecord> reader, boolean silent) throws IOException {
    super(reader.header());
    mSilent = silent;
    mPath = path;
    mWrapped = reader;

    try {
      nextRecord();
    } catch (final SAMFormatException e) {
      mWrapped.close();
      throw new NoTalkbackSlimException(ErrorType.SAM_BAD_FORMAT, path, e.getMessage());
    } catch (final IllegalArgumentException e) {
      mWrapped.close();
      throw new NoTalkbackSlimException(e, ErrorType.SAM_BAD_FORMAT, path, e.getMessage());
    } catch (final RuntimeException e) {
      if (e.getMessage() != null && e.getMessage().contains("No M operator between pair of IDN operators in CIGAR")) {
        maybeWarn(e);
      } else {
        mWrapped.close();
        throw e;
      }
    }
  }

  /**
   * Constructor
   * @param path filename to use for errors and warnings
   * @param reader supplier of <code>SAMRecords</code>
   * @throws IOException if an IO error occurs
   */
  public SkipInvalidRecordsIterator(String path, RecordIterator<SAMRecord> reader) throws IOException {
    this(path, reader, false);
  }

  /**
   * Constructor
   * @param path filename to use for errors and warnings
   * @param reader to obtain base iterator from
   * @param silent true to not report warnings
   * @throws IOException if an IO error occurs
   */
  public SkipInvalidRecordsIterator(String path, final SAMFileReader reader, boolean silent) throws IOException {
    this(path, new SamFileReaderAdaptor(reader, null), silent);
  }

  /**
   * Constructor
   * @param path filename to use for errors and warnings
   * @param reader to obtain base iterator from
   * @throws IOException if an IO error occurs
   */
  public SkipInvalidRecordsIterator(String path, SAMFileReader reader) throws IOException {
    this(path, new SamFileReaderAdaptor(reader, null));
  }

  /**
   * Constructs an iterator on given file.
   * @param samFile file containing SAM data
   * @throws IOException if an IO Error occurs
   */
  public SkipInvalidRecordsIterator(File samFile) throws IOException {
    this(samFile.getPath(), makeSAMFileReader(samFile));
  }

  private static SAMFileReader makeSAMFileReader(File samFile) throws IOException {
    final String samPath = samFile.getPath();
    try {
      return new SAMFileReader(samFile);
    } catch (final RuntimeIOException e) {
      if (e.getCause() instanceof FileNotFoundException) {
        throw new NoTalkbackSlimException(ErrorType.FILE_NOT_FOUND, samPath);
      } else {
        throw (IOException) e.getCause();
      }
    }
  }

  private void maybeWarn(Throwable t) {
    mNumInvalidRecords++;
    if (!mSilent && mNumInvalidRecords <= NUM_RECORDS_TO_WARN) {
      Diagnostic.warning(WarningType.SAM_BAD_FORMAT_WARNING, mPath, LS + t.getMessage() + " at data line " + mRecordCount);
    }
  }

  private void nextRecord() {
    while (mWrapped.hasNext()) {
      mRecordCount++;
      try {
        final SAMRecord current = mWrapped.next();
        if (current != null) {
          if (current.getIsValid()) {
            mTotalNucleotides += current.getReadLength();
            mRecord = current;
            return;
          } else {
            // invalid records but we continue
            mNumInvalidRecords++;
            if (!mSilent && mNumInvalidRecords <= NUM_RECORDS_TO_WARN) {

              final List<SAMValidationError> valid = current.isValid();
              if (valid == null) {
                Diagnostic.warning(WarningType.SAM_BAD_FORMAT_WARNING, mPath, LS + "At data line " + mRecordCount);
              } else {
                Diagnostic.warning(WarningType.SAM_BAD_FORMAT_WARNING, mPath, ""
                    + LS + current.getSAMString().trim()
                    + LS + valid
                    + LS + "At data line " + mRecordCount
                    );
              }
            }
          }
        }
      } catch (final SAMFormatException e) {
        throw new NoTalkbackSlimException(e, ErrorType.SAM_BAD_FORMAT, mPath, e.getMessage());
      } catch (final IllegalArgumentException e) {
        if (e.getMessage() != null && e.getMessage().contains("Malformed CIGAR string:")) {
          maybeWarn(e);
        } else {
          throw new NoTalkbackSlimException(e, ErrorType.SAM_BAD_FORMAT, mPath, e.getMessage());
        }
      } catch (final RuntimeException e) {
        if (e.getMessage() != null && e.getMessage().contains("No M operator between pair of IDN operators in CIGAR")) {
          maybeWarn(e);
        } else {
          throw e;
        }
      }
    }
    mRecord = null;
  }


  @Override
  public boolean hasNext() {
    return mRecord != null;
  }

  @Override
  public SAMRecord next() {
    if (mRecord == null) {
      throw new NoSuchElementException();
    }
    final SAMRecord ret = mRecord;
    nextRecord();
    return ret;
  }

  @Override
  public void close() throws IOException {
    if (mIsClosed) {
      return;
    }
    if (!mSilent && mNumInvalidRecords > 0) {
      Diagnostic.warning(WarningType.SAM_IGNORED_RECORDS, "" + mNumInvalidRecords, mPath);
    }
    mIsClosed = true;
    mWrapped.close();
  }

}
